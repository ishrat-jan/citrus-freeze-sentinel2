"""Freeze damage assessment for citrus - Lake County (Clermont), Florida, Jan-Feb 2026.

    python run.py           # full 10 m run
    python run.py --test    # quick 60 m run

All settings are in config.yaml.
"""
import argparse
import os
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src import damage, data, indices, plots


def load_config(path, test):
    cfg = yaml.safe_load(open(path))
    r = cfg["run"]
    r["test"] = test
    r["res"] = r["test_resolution_m"] if test else r["resolution_m"]
    r["pixel_ha"] = r["res"] ** 2 / 10_000
    r["cache"] = Path(r["cache_dir"]) / f"{r['res']}m"
    r["out"] = Path(r["output_dir"]) / f"{r['res']}m"
    r["figures"] = r["out"] / "figures" if test else Path(r["figures_dir"])
    for p in (r["cache"], r["out"] / "tables", r["out"] / "rasters", r["figures"]):
        p.mkdir(parents=True, exist_ok=True)
    return cfg


def step(msg):
    print(f"\n[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def save_tif(da, path, nodata=None):
    """Write a 2-D DataArray as a compressed GeoTIFF."""
    out = da.copy()
    out.attrs, out.encoding = {}, {}
    if nodata is not None:
        out = out.rio.write_nodata(nodata)
    out.rio.to_raster(path, compress="deflate")


def main(config="config.yaml", test=False):
    cfg = load_config(config, test)
    d, r, tp = cfg["dates"], cfg["run"], cfg["timeseries"]
    early, late = d["damage_early"], d["damage_late"]
    names = cfg["indices"]
    print(f"Resolution {r['res']} m -> {r['out']}, figures -> {r['figures']}")

    step("1. Search Sentinel-2 scenes")
    items = data.search(cfg)
    dates = sorted({data.day(it) for it in items})
    print(f"  {len(dates)} dates: {', '.join(dates)}")

    step("2. Download bands (cached) and find the clear dates")
    clear = {}
    for date in dates:
        data.download(cfg, items, date, cfg["sentinel2"]["base_bands"])
        clear[date] = data.load_date(cfg, date, ["red"])[1]
    must = set(d["baseline"]) | {early, late} | set(d["noise_pair"])
    index_dates = [x for x in dates if clear[x] >= tp["min_clear_pct"] or x in must]
    for date in index_dates:
        data.download(cfg, items, date, indices.bands_for(names))
    clear_tab = pd.DataFrame({"date": dates, "clear_pct": [round(clear[x], 1) for x in dates],
                              "used_for_indices": [x in index_dates for x in dates]})
    clear_tab.to_csv(r["out"] / "tables" / "clear_fraction_by_date.csv", index=False)
    print(clear_tab.to_string(index=False))

    def index_on(date, name):
        refl, _ = data.load_date(cfg, date, indices.BANDS[name])
        return indices.compute(name, refl)

    step("3. NDVI baseline, vegetation and citrus masks")
    base_ndvi, base_count = damage.baseline([index_on(x, "NDVI") for x in d["baseline"]])
    veg = base_ndvi >= cfg["thresholds"]["vegetation_ndvi"]
    land = base_ndvi >= 0                      # lakes have negative NDVI (display only)
    citrus = damage.citrus_mask(base_ndvi, cfg)
    groups = {"citrus": (citrus & veg).values, "other vegetation": (~citrus & veg).values}

    step("4. Each index: baseline, change on the damage dates, maps")
    stats, ndvi_change = [], {}
    for name in names:
        base = base_ndvi if name == "NDVI" else damage.baseline([index_on(x, name) for x in d["baseline"]])[0]
        a, b = (index_on(x, name) for x in d["noise_pair"])
        noise = (b - a).values[veg.values]
        noise = noise[np.isfinite(noise)].astype("float64")
        row = {"index": name,
               "baseline_median": round(float(np.nanmedian(base.values[groups["citrus"]])), 3),
               # robust SD (1.4826 x median absolute deviation): not inflated by the few extreme
               # pixels that ratio indices such as EVI and MCARI produce
               "noise_sd": round(1.4826 * float(np.median(np.abs(noise - np.median(noise)))), 4)}
        maps, changes = {"Baseline (median, 4-27 Jan)": plots.thumb(base, cfg)}, {}
        for date in (early, late):
            val = index_on(date, name)
            ch = (val - base).where(veg)
            v = ch.values[groups["citrus"]]
            row[f"median_change_{date}"] = round(float(np.nanmedian(v)), 3)
            maps[date] = plots.thumb(val, cfg)
            changes[date] = plots.thumb(ch, cfg)
            if name == "NDVI":
                ndvi_change[date] = ch
        stats.append(row)
        lo, hi = np.nanpercentile(base.values[veg.values], cfg["plots"]["map_percentiles"])
        lim = np.nanpercentile(np.abs(changes[late][0]), cfg["plots"]["change_percentile"])
        plots.index_maps(name, maps, lo, hi, cfg, r["figures"])
        plots.change_maps(name, changes, float(f"{lim:.2g}"), cfg, r["figures"])
        print(f"  {name}: citrus median change {row[f'median_change_{early}']:+.3f} ({early}), "
              f"{row[f'median_change_{late}']:+.3f} ({late})")
    comparison = damage.index_change_table(stats, cfg)
    comparison.to_csv(r["out"] / "tables" / "index_change_citrus.csv", index=False)

    step("5. Time series of all indices (citrus vs other vegetation)")
    rows = []
    for date in index_dates:
        refl, _ = data.load_date(cfg, date, indices.bands_for(names))
        for name in names:
            v = indices.compute(name, refl).values
            for group, m in groups.items():
                ok = m & np.isfinite(v)
                pct = 100 * ok.sum() / max(m.sum(), 1)
                if pct >= tp["min_group_clear_pct"]:
                    rows.append({"date": date, "index": name, "group": group,
                                 "median": round(float(np.median(v[ok])), 4), "clear_pct": round(pct, 1)})
        del refl
    ts = pd.DataFrame(rows)
    ts.to_csv(r["out"] / "tables" / "timeseries_indices.csv", index=False)
    plots.timeseries(ts, cfg, r["figures"])

    step("6. NDVI damage classes and kill zones")
    damages = {x: damage.classify(ndvi_change[x], veg, cfg) for x in (early, late)}
    ids = damage.kill_zones(damages[early], damages[late], citrus, cfg)
    kz = damage.zone_polygons(ids, ndvi_change[early], ndvi_change[late], cfg)
    head = damage.headline(citrus, veg, damages, kz, cfg)
    head.to_csv(r["out"] / "tables" / "headline_results.csv", index=False)
    kz.drop(columns="geometry").to_csv(r["out"] / "tables" / "kill_zones.csv", index=False)
    if len(kz):
        kz.to_file(r["out"] / "kill_zones.gpkg", driver="GPKG")
    plots.damage_map(damages, citrus, land, ids, kz, cfg, r["figures"])
    print(head.to_string(index=False))

    step("7. GeoTIFFs")
    save_tif(base_ndvi, r["out"] / "rasters" / "ndvi_baseline.tif")
    save_tif(base_count, r["out"] / "rasters" / "baseline_clear_dates.tif")
    for x in (early, late):
        save_tif(ndvi_change[x], r["out"] / "rasters" / f"ndvi_change_{x}.tif")
        save_tif(damages[x], r["out"] / "rasters" / f"damage_class_{x}.tif", nodata=damage.NODATA)
    save_tif(ids, r["out"] / "rasters" / "kill_zone_id.tif", nodata=0)
    save_tif(citrus.astype("uint8"), r["out"] / "rasters" / "citrus_mask.tif")

    print(f"\nDone. Tables and rasters in {r['out']}, figures in {r['figures']}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)      # e.g. median of all-NaN pixels
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--test", action="store_true", help="quick 60 m run")
    a = p.parse_args()
    try:
        main(a.config, a.test)
    except Exception:
        traceback.print_exc()
        sys.stdout.flush()
        os._exit(1)        # exit at once: background download threads could otherwise keep it hanging
