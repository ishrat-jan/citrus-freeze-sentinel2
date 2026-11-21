"""Sentinel-2 data: search, download (cached), reflectance and cloud mask.

Each band of each date is cached as its own GeoTIFF (raw digital numbers):
    data/cache/<res>m/<YYYY-MM-DD>/<band>.tif
so a re-run only downloads what is missing, an interrupted download resumes,
and only one date is held in memory at a time.

BOA offset: newer L2A products can carry a +1000 offset in the digital numbers.
Earth Search's flag for it is not reliable, so it is decided from the data:
the offset is present only if almost no clear-land red value is below 1000
(with an offset, such values would mean negative reflectance).
"""
import ctypes
import time

import numpy as np
import rioxarray
from odc.stac import configure_rio, load as stac_load
from pystac_client import Client


def search(cfg):
    """All Sentinel-2 L2A scenes of the tile that cover the AOI in the date window."""
    s2 = cfg["sentinel2"]
    items = list(Client.open(s2["stac_url"]).search(
        collections=[s2["collection"]],
        bbox=cfg["aoi"]["bbox"],
        datetime=f"{s2['start_date']}/{s2['end_date']}",
        query={"grid:code": {"eq": s2["mgrs_tile"]}},
    ).items())
    if not items:
        raise RuntimeError("No Sentinel-2 scenes found - check config.yaml")
    return items


def day(item):
    return item.datetime.strftime("%Y-%m-%d")


def band_path(cfg, date, band):
    return cfg["run"]["cache"] / date / f"{band}.tif"


def _clear_gdal_cache():
    """After a failed web read GDAL may cache the broken result; forget it before retrying."""
    try:
        lib = [l.split()[-1] for l in open("/proc/self/maps") if "libgdal" in l][0]
        ctypes.CDLL(lib).VSICurlClearCache()
    except Exception:
        pass  # not possible on this system: the retry still runs


def download(cfg, items, date, bands):
    """Download the bands of one date that are not cached yet."""
    missing = [b for b in bands if not band_path(cfg, date, b).exists()]
    if not missing:
        return
    s2 = cfg["sentinel2"]
    configure_rio(cloud_defaults=True,
                  GDAL_HTTP_MAX_RETRY=str(s2["gdal_http_retries"]),
                  GDAL_HTTP_RETRY_DELAY=str(s2["gdal_http_retry_delay_s"]))
    wait = s2["retry_wait_s"]
    for attempt in range(1, s2["retries"] + 1):
        try:
            print(f"  downloading {date}: {', '.join(missing)}", flush=True)
            ds = stac_load(
                [it for it in items if day(it) == date],
                bands=missing,
                bbox=cfg["aoi"]["bbox"],
                crs=cfg["run"]["crs"],
                resolution=cfg["run"]["res"],
                groupby="solar_day",            # merge duplicate scenes of the same day
                chunks={"x": s2["chunk_size"], "y": s2["chunk_size"]},
                resampling={b: "nearest" if b == "scl" else s2["resampling"] for b in missing},
                fail_on_error=True,
            ).isel(time=0).compute()
            break
        except Exception as err:                 # dropped connection: wait and try again
            if attempt == s2["retries"]:
                raise
            print(f"  failed (try {attempt}/{s2['retries']}): {str(err).splitlines()[-1][:120]}"
                  f" - retrying in {wait} s", flush=True)
            time.sleep(wait)
            wait *= 2
            _clear_gdal_cache()

    folder = cfg["run"]["cache"] / date
    folder.mkdir(parents=True, exist_ok=True)
    for b in missing:
        path = band_path(cfg, date, b)
        tmp = path.with_name(path.stem + ".part.tif")     # never leave a half-written file
        ds[b].rio.write_crs(cfg["run"]["crs"]).rio.to_raster(tmp, driver="GTiff", compress="deflate")
        tmp.rename(path)


def read_band(cfg, date, band):
    """One cached band (digital numbers) as a 2-D DataArray (y, x)."""
    return rioxarray.open_rasterio(band_path(cfg, date, band)).squeeze("band", drop=True).load()


def boa_offset(red, scl, cfg):
    """0 or the BOA offset, decided from the red band's dark clear-land pixels."""
    s2 = cfg["sentinel2"]
    clear = scl.isin(s2["clear_scl_classes"]) & (red > 0)
    below = int(((red < -s2["boa_offset"]) & clear).sum()) / max(int(clear.sum()), 1)
    return s2["boa_offset"] if below < s2["offset_check_max_below"] else 0


def load_date(cfg, date, bands, mask=True):
    """Reflectance (float32) of `bands` for one date, NaN where the pixel is not clear
    (mask=False keeps every pixel, e.g. for a true-colour picture).

    Clear = SCL vegetation or bare soil, and valid red + nir.
    Also returns the share of pixels with data that are clear.
    """
    s2 = cfg["sentinel2"]
    dn = {b: read_band(cfg, date, b) for b in set(bands) | set(s2["base_bands"])}
    offset = np.float32(boa_offset(dn["red"], dn["scl"], cfg))
    scale = np.float32(s2["scale"])
    refl = {b: (dn[b].astype("float32").where(dn[b] > 0) + offset) / scale
            for b in dn if b != "scl"}
    clear = dn["scl"].isin(s2["clear_scl_classes"]) & ((refl["nir"] + refl["red"]) > 0)
    clear_pct = 100 * float(clear.sum()) / max(int((dn["scl"] > 0).sum()), 1)
    return {b: refl[b].where(clear) if mask else refl[b] for b in bands}, clear_pct
