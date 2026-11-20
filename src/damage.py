"""Baseline, change, NDVI damage classes, citrus mask and kill zones.

Damage classes (NDVI change vs baseline, vegetated pixels only):
    0 none      change >= -0.10
    1 mild      -0.10 > change >= -0.20
    2 moderate  -0.20 > change >= -0.30
    3 severe    change < -0.30
    255         no data (not vegetated, cloudy, or no baseline)

Kill zone = citrus that is moderate or worse on 21 Feb AND already damaged on
6 Feb (damage that appeared after the freeze and did not recover), in connected
patches of at least 0.1 ha.
"""
import math
import re
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
import rioxarray
import xarray as xr
from odc.geo.xr import xr_reproject
from pyproj import Transformer
from rasterio.features import shapes
from scipy import ndimage
from shapely.geometry import shape

NODATA = 255
CLASS_NAMES = {0: "none", 1: "mild", 2: "moderate", 3: "severe"}


# ---------- baseline, change, damage classes ----------
def baseline(arrays):
    """Per-pixel median of the (cloud-masked) baseline dates, and the number of clear dates."""
    stack = xr.concat(arrays, dim="time")
    return stack.median("time", skipna=True), stack.notnull().sum("time").astype("uint8")


def classify(change, veg, cfg):
    """NDVI change -> damage classes (uint8)."""
    t = cfg["thresholds"]
    c = change.values
    classes = np.select([c < t["severe"], c < t["moderate"], c < t["mild"], c >= t["mild"]],
                        [3, 2, 1, 0], default=NODATA).astype("uint8")   # NaN -> no data
    classes[~veg.values] = NODATA
    return xr.DataArray(classes, coords=change.coords, dims=change.dims, name="damage_class")


# ---------- citrus mask (USDA Cropland Data Layer via CropScape) ----------
def _get(url, c):
    """GET with timeout and retries (the CropScape server is slow)."""
    for attempt in range(1, c["retries"] + 1):
        try:
            r = requests.get(url, timeout=c["timeout_s"])
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt == c["retries"]:
                raise
            time.sleep(c["retry_wait_s"])


def download_cdl(cfg):
    """CDL clip for the AOI (cached GeoTIFF). The returned link expires quickly."""
    c = cfg["citrus"]
    path = cfg["run"]["cache"].parent / f"cdl_{c['cdl_year']}.tif"
    if path.exists():
        return path
    # AOI corners -> CDL projection (all 4 corners: lon/lat boxes are not boxes in Albers)
    w, s, e, n = cfg["aoi"]["bbox"]
    xs, ys = Transformer.from_crs("EPSG:4326", c["cdl_crs"], always_xy=True).transform(
        [w, e, e, w], [s, s, n, n])
    bbox = f"{min(xs):.0f},{min(ys):.0f},{max(xs):.0f},{max(ys):.0f}"
    xml = _get(f"{c['api_url']}?year={c['cdl_year']}&bbox={bbox}", c).text
    link = re.search(r"<returnURL>(.*?)</returnURL>", xml)
    if not link:
        raise RuntimeError(f"CropScape returned no download link:\n{xml[:300]}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(_get(link.group(1), c).content)
    return path


def citrus_mask(template, cfg):
    """Boolean citrus mask on the analysis grid (nearest neighbour keeps class codes intact)."""
    cdl = rioxarray.open_rasterio(download_cdl(cfg)).squeeze("band", drop=True)
    cdl = xr_reproject(cdl, template.odc.geobox, resampling="nearest")
    return (cdl == cfg["citrus"]["cdl_class"]).rename("citrus")


# ---------- kill zones ----------
def kill_zones(dmg_early, dmg_late, citrus, cfg):
    """Zone-ID raster (0 = none): connected kill-zone patches of at least min_area_ha."""
    k = cfg["kill_zones"]
    e, l = dmg_early.values, dmg_late.values
    mask = ((l != NODATA) & (l >= k["late_min_class"]) & (e != NODATA) & (e >= k["early_min_class"])
            & citrus.values)
    labels, _ = ndimage.label(mask, structure=np.ones((3, 3)) if k["connectivity"] == 8 else None)
    min_px = max(1, math.ceil(k["min_area_ha"] / cfg["run"]["pixel_ha"]))
    sizes = np.bincount(labels.ravel())
    keep = sizes >= min_px
    keep[0] = False                                     # 0 = background
    new_id = np.zeros(len(sizes), dtype="uint32")
    new_id[keep] = np.arange(1, keep.sum() + 1)
    return xr.DataArray(new_id[labels], coords=dmg_late.coords, dims=dmg_late.dims, name="kill_zone_id")


def zone_polygons(ids, change_early, change_late, cfg):
    """Kill-zone polygons with area, mean NDVI change and centre (lon/lat), largest first."""
    v = ids.values
    if v.max() == 0:
        return gpd.GeoDataFrame(columns=["zone_id", "geometry"], geometry="geometry", crs=cfg["run"]["crs"])
    geoms = [(int(z), shape(g)) for g, z in shapes(v.astype("int32"), mask=v > 0,
                                                    transform=ids.odc.geobox.affine)]
    gdf = gpd.GeoDataFrame(geoms, columns=["zone_id", "geometry"], crs=cfg["run"]["crs"])
    gdf = gdf.dissolve("zone_id").reset_index()
    index = np.arange(1, v.max() + 1)
    n_px = ndimage.sum(np.ones_like(v), v, index)
    gdf = gdf.merge(pd.DataFrame({
        "zone_id": index,
        "area_ha": (n_px * cfg["run"]["pixel_ha"]).round(2),
        "ndvi_change_6feb": ndimage.mean(np.nan_to_num(change_early.values), v, index).round(3),
        "ndvi_change_21feb": ndimage.mean(np.nan_to_num(change_late.values), v, index).round(3),
    }), on="zone_id")
    centre = gdf.geometry.centroid.to_crs("EPSG:4326")
    gdf["lon"], gdf["lat"] = centre.x.round(5), centre.y.round(5)
    return gdf.sort_values("area_ha", ascending=False)


# ---------- summary tables ----------
def headline(citrus, veg, damages, kz, cfg):
    """Citrus area and NDVI damage figures (one row per item)."""
    ha = cfg["run"]["pixel_ha"]
    c = citrus.values
    rows = [("CDL citrus in AOI (ha)", round(c.sum() * ha, 1)),
            ("Assessed citrus: vegetated, with baseline (ha)", round((c & veg.values).sum() * ha, 1))]
    for date, dmg in damages.items():
        d = dmg.values[c]
        d = d[d != NODATA]
        rows += [(f"{date}: damaged citrus, mild or worse (%)", round(100 * np.mean(d >= 1), 1)),
                 (f"{date}: moderate or worse (%)", round(100 * np.mean(d >= 2), 1)),
                 (f"{date}: severe (ha)", round((d == 3).sum() * ha, 2))]
    rows += [("Kill zones (number)", len(kz)),
             ("Kill zones (ha)", round(float(kz.area_ha.sum()) if len(kz) else 0.0, 2))]
    return pd.DataFrame(rows, columns=["item", "value"])


def index_change_table(stats, cfg):
    """How much each index changed for citrus, also in units of normal date-to-date noise."""
    df = pd.DataFrame(stats)
    late = cfg["dates"]["damage_late"]
    df["change_vs_noise_21feb"] = (df[f"median_change_{late}"] / df["noise_sd"]).round(1)
    return df
