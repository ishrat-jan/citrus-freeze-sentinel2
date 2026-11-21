"""Frost-prone (low-lying) areas from the Copernicus DEM.

On calm, clear freezing nights cold air drains downhill and pools in low ground.
Low-lying index = elevation minus the mean elevation within 300 m
(negative = lower than the surroundings). Classes for the map:
    0 not low-lying    index >= hollow_m (0 m)
    1 hollow           deep_hollow_m (-4 m) <= index < hollow_m
    2 deep hollow      index < deep_hollow_m
    255                lake or no data
Citrus damage on the damage date is then compared between the classes.
"""
import numpy as np
import pandas as pd
import rioxarray
import xarray as xr
from odc.stac import configure_rio, load as stac_load
from pystac_client import Client
from scipy import ndimage

from src.damage import NODATA

FROST_NAMES = {0: "not low-lying", 1: "hollow", 2: "deep hollow"}


def load_dem(template, cfg):
    """Copernicus DEM (m) on the analysis grid, cached as a GeoTIFF."""
    path = cfg["run"]["cache"] / "dem.tif"
    if not path.exists():
        items = list(Client.open(cfg["sentinel2"]["stac_url"]).search(
            collections=[cfg["terrain"]["dem_collection"]], bbox=cfg["aoi"]["bbox"]).items())
        s2 = cfg["sentinel2"]
        configure_rio(cloud_defaults=True, aws={"aws_unsigned": True},   # public bucket, no login
                      GDAL_HTTP_MAX_RETRY=str(s2["gdal_http_retries"]),
                      GDAL_HTTP_RETRY_DELAY=str(s2["gdal_http_retry_delay_s"]))
        print(f"  downloading DEM ({len(items)} tiles)", flush=True)
        dem = stac_load(items, bands=["data"], geobox=template.odc.geobox, resampling="bilinear")
        dem = dem["data"].isel(time=0).astype("float32").compute()
        tmp = path.with_name("dem.part.tif")
        dem.rio.write_crs(cfg["run"]["crs"]).rio.to_raster(tmp, driver="GTiff", compress="deflate")
        tmp.rename(path)
    return rioxarray.open_rasterio(path).squeeze("band", drop=True).load()


def lowlying_index(dem, cfg):
    """Elevation minus the mean elevation in a square window of +/- radius."""
    px = int(round(cfg["terrain"]["lowlying_radius_m"] / cfg["run"]["res"]))
    mean = ndimage.uniform_filter(dem.values, size=2 * px + 1, mode="nearest")
    return dem.copy(data=dem.values - mean).rename("lowlying_index")


def frost_classes(low, water, cfg):
    """Frost-prone class raster (uint8) from the low-lying index; lakes = no data."""
    t, li = cfg["terrain"], low.values
    fp = np.zeros(li.shape, dtype="uint8")
    fp[(li < t["hollow_m"]) & (li >= t["deep_hollow_m"])] = 1
    fp[li < t["deep_hollow_m"]] = 2
    fp[np.isnan(li) | water.values] = NODATA
    return xr.DataArray(fp, coords=low.coords, dims=low.dims, name="frost_prone")


def damage_rates(frost, low, dmg, citrus, cfg):
    """Citrus damage rate per frost-prone class and per low-lying index bin."""
    ha = cfg["run"]["pixel_ha"]
    valid = (dmg.values != NODATA) & citrus.values & np.isfinite(low.values)
    damaged = dmg.values[valid] >= 1

    def table(groups, order):
        df = pd.DataFrame({"group": groups, "damaged": damaged})
        t = df.groupby("group").agg(citrus_ha=("damaged", "size"), damage_rate_pct=("damaged", "mean"))
        t["citrus_ha"] = (t.citrus_ha * ha).round(1)
        t["damage_rate_pct"] = (100 * t.damage_rate_pct).round(1)
        return t.reindex([o for o in order if o in t.index]).reset_index()

    fc = frost.values[valid]
    by_class = table(np.array([FROST_NAMES.get(int(k), "lake / no data") for k in fc]),
                     list(FROST_NAMES.values())).rename(columns={"group": "frost_class"})
    hollow = np.isin(fc, (1, 2))                       # hollow + deep hollow together
    by_class.loc[len(by_class)] = ["all hollows (below surroundings)", round(hollow.sum() * ha, 1),
                                   round(100 * damaged[hollow].mean(), 1) if hollow.any() else np.nan]
    bins = cfg["terrain"]["lowlying_bins"]
    names = [f"< {bins[0]} m"] + [f"{a} to {b} m" for a, b in zip(bins[:-1], bins[1:])] + [f">= {bins[-1]} m"]
    labels = np.array(names, dtype=object)[np.digitize(low.values[valid], bins)]
    by_bin = table(labels, names).rename(columns={"group": "lowlying_index"})
    return by_class, by_bin
