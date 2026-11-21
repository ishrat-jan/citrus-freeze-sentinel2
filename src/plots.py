"""All figures (PNG).

Maps are drawn with matplotlib imshow on a downsampled array (every n-th pixel),
which is fast even at 10 m. Axes are UTM 17N coordinates in km.
"""
import math

import matplotlib

matplotlib.use("Agg")                       # no screen needed
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import geopandas as gpd  # noqa: E402
from rasterio.features import shapes  # noqa: E402
from scipy import ndimage  # noqa: E402
from shapely.geometry import shape  # noqa: E402

from src.damage import NODATA  # noqa: E402


# ---------- helpers ----------
def thumb(da, cfg):
    """Downsampled array and its imshow extent in km (left, right, bottom, top)."""
    s = max(1, math.ceil(max(da.shape) / cfg["plots"]["max_pixels"]))
    x, y = da.x.values, da.y.values
    dx, dy = abs(x[1] - x[0]) / 2, abs(y[1] - y[0]) / 2
    ext = [(x.min() - dx) / 1e3, (x.max() + dx) / 1e3, (y.min() - dy) / 1e3, (y.max() + dy) / 1e3]
    return da.values[::s, ::s], ext


def _style(ax, cfg, title):
    c = cfg["plots"]["colors"]
    ax.set_facecolor(c["surface"])
    ax.set_title(title, fontsize=10, color=c["text"], loc="left")
    ax.tick_params(labelsize=7, colors=c["muted"])
    for s in ax.spines.values():
        s.set_visible(False)


def _save(fig, folder, name, cfg):
    path = folder / f"{name}.png"
    fig.savefig(path, dpi=cfg["plots"]["dpi"], bbox_inches="tight",
                facecolor=cfg["plots"]["colors"]["surface"])
    plt.close(fig)
    return path


def _map_row(n, cfg):
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 5.6), sharey=True, layout="constrained",
                             facecolor=cfg["plots"]["colors"]["surface"])
    axes[0].set_ylabel("Northing (km, UTM 17N)", fontsize=8, color=cfg["plots"]["colors"]["muted"])
    for ax in axes:
        ax.set_xlabel("Easting (km)", fontsize=8, color=cfg["plots"]["colors"]["muted"])
    return fig, axes


# ---------- time series ----------
def timeseries(ts, cfg, folder):
    """Median of each index over clear pixels: citrus vs other vegetation, freeze shaded."""
    c, d = cfg["plots"]["colors"], cfg["dates"]
    names = cfg["indices"]
    ncol = 2
    nrow = math.ceil((len(names) + 1) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 3 * nrow), sharex=True, facecolor=c["surface"])
    axes = axes.ravel()
    f0, f1 = pd.to_datetime(d["freeze_period"])
    for ax, name in zip(axes, names):
        _style(ax, cfg, name)
        ax.axvspan(f0, f1 + pd.Timedelta(days=1), color=c["freeze"], lw=0, zorder=0)
        for day in (d["damage_early"], d["damage_late"]):
            ax.axvline(pd.Timestamp(day), color=c["muted"], lw=0.8, ls=":", zorder=1)
        for group, col in (("citrus", c["citrus"]), ("other vegetation", c["other_veg"])):
            g = ts[(ts["index"] == name) & (ts.group == group)]
            ax.plot(pd.to_datetime(g.date), g["median"], color=col, lw=2, marker="o", ms=5, zorder=3)
        ax.grid(axis="y", color=c["freeze"], lw=0.5)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        ax.tick_params(labelbottom=True)          # dates under every panel
    # last cell: legend and notes
    leg = axes[len(names)]
    leg.axis("off")
    handles = [Line2D([], [], color=c["citrus"], lw=2, marker="o", label="citrus"),
               Line2D([], [], color=c["other_veg"], lw=2, marker="o", label="other vegetation"),
               Patch(color=c["freeze"], label="main freeze (26 Jan - 2 Feb)"),
               Line2D([], [], color=c["muted"], lw=0.8, ls=":", label="damage dates (6 and 21 Feb)")]
    leg.legend(handles=handles, loc="center left", frameon=False, fontsize=9)
    for ax in axes[len(names) + 1:]:
        ax.axis("off")
    fig.suptitle("Median index value per clear date - citrus vs all other vegetation", x=0.01, ha="left",
                 fontsize=12, color=c["text"])
    fig.tight_layout()
    return _save(fig, folder, "timeseries_indices", cfg)


# ---------- index and change maps ----------
def index_maps(name, panels, vmin, vmax, cfg, folder):
    """One index on the baseline, early and late date (shared colour scale)."""
    fig, axes = _map_row(len(panels), cfg)
    for ax, (title, (arr, ext)) in zip(axes, panels.items()):
        _style(ax, cfg, title)
        im = ax.imshow(arr, extent=ext, cmap=cfg["plots"]["index_cmaps"][name], vmin=vmin, vmax=vmax,
                       interpolation="nearest")
    cb = fig.colorbar(im, ax=axes, shrink=0.7)
    cb.set_label(name, fontsize=9)
    fig.suptitle(f"{name} - all vegetation (baseline NDVI >= 0.3), cloud-masked", x=0.01, ha="left",
                 fontsize=12)
    return _save(fig, folder, f"map_{name.lower()}", cfg)


def change_maps(name, panels, lim, cfg, folder):
    """Index change vs baseline on both damage dates (red = decrease, blue = increase)."""
    fig, axes = _map_row(len(panels), cfg)
    for ax, (title, (arr, ext)) in zip(axes, panels.items()):
        _style(ax, cfg, title)
        im = ax.imshow(arr, extent=ext, cmap="RdBu", vmin=-lim, vmax=lim, interpolation="nearest")
    cb = fig.colorbar(im, ax=axes, shrink=0.7)
    cb.set_label(f"{name} change", fontsize=9)
    fig.suptitle(f"{name} change vs baseline - all vegetation (red = decrease)", x=0.01, ha="left",
                 fontsize=12)
    return _save(fig, folder, f"change_{name.lower()}", cfg)


# ---------- NDVI damage classes + kill zones ----------
def _damage_rgb(dmg, citrus, land, cfg):
    """Class codes for display: 0 not shown, 1 land, 2 citrus without damage, 3-5 damage classes."""
    c = cfg["plots"]["colors"]
    out = np.where(land, 1, 0)
    cit = citrus & (dmg != NODATA)
    out = np.where(cit & (dmg == 0), 2, out)
    for k in (1, 2, 3):
        out = np.where(cit & (dmg == k), 2 + k, out)
    cmap = ListedColormap([c["surface"], c["land"], c["citrus_ok"], *c["damage"]])
    return out, cmap


def damage_map(damages, citrus, land, ids, kz, cfg, folder):
    """Citrus damage classes on both dates, kill zones outlined, plus a close-up."""
    c = cfg["plots"]["colors"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 6.2), facecolor=c["surface"],
                             gridspec_kw={"width_ratios": [1, 1, 1.05]})
    for ax, (date, dmg) in zip(axes[:2], damages.items()):
        codes, cmap = _damage_rgb(dmg.values, citrus.values, land.values, cfg)
        small, ext = thumb(dmg.copy(data=codes), cfg)
        _style(ax, cfg, f"Citrus damage (NDVI), {date}")
        ax.imshow(small, extent=ext, cmap=cmap, vmin=-0.5, vmax=5.5, interpolation="nearest")
        ax.set_xlabel("Easting (km)", fontsize=8, color=c["muted"])
    axes[0].set_ylabel("Northing (km, UTM 17N)", fontsize=8, color=c["muted"])
    if len(kz):
        g = kz.copy()
        g["geometry"] = g.geometry.scale(1e-3, 1e-3, origin=(0, 0))     # m -> km
        g.boundary.plot(ax=axes[1], color=c["kill_zone"], lw=0.6)
        pts = g.geometry.centroid
        axes[1].scatter(pts.x, pts.y, s=70, facecolors="none", edgecolors=c["kill_zone"], lw=1.3)
        # close-up around the largest zone, at full resolution
        half = cfg["plots"]["zoom_km"] * 500
        cx, cy = kz.geometry.iloc[0].centroid.coords[0]
        win = dict(x=slice(cx - half, cx + half), y=slice(cy + half, cy - half))
        late = list(damages)[-1]
        codes, cmap = _damage_rgb(damages[late].sel(**win).values, citrus.sel(**win).values,
                                  land.sel(**win).values, cfg)
        arr, ext = thumb(damages[late].sel(**win).copy(data=codes), cfg)
        _style(axes[2], cfg, f"Close-up of the largest kill zone ({kz.area_ha.iloc[0]:.2f} ha)")
        axes[2].imshow(arr, extent=ext, cmap=cmap, vmin=-0.5, vmax=5.5, interpolation="nearest")
        g.boundary.plot(ax=axes[2], color=c["kill_zone"], lw=0.7)
        axes[2].set_xlim(ext[0], ext[1])
        axes[2].set_ylim(ext[2], ext[3])
        axes[2].set_xlabel("Easting (km)", fontsize=8, color=c["muted"])
    else:
        axes[2].axis("off")
    axes[1].set_ylabel("")                       # geopandas adds its own axis labels
    axes[2].set_ylabel("Northing (km)", fontsize=8, color=c["muted"])
    handles = [Patch(color=c["citrus_ok"], label="citrus, no damage")] + \
              [Patch(color=col, label=lab) for col, lab in zip(c["damage"], ["mild", "moderate", "severe"])] + \
              [Patch(facecolor="none", edgecolor=c["kill_zone"], lw=1.5, label="kill zone"),
               Patch(color=c["land"], label="other land")]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=9)
    fig.suptitle("NDVI damage classes - citrus only (whole area, and close-up of the largest kill zone)",
                 x=0.01, ha="left", fontsize=12, color=c["text"])
    fig.subplots_adjust(bottom=0.12, wspace=0.15)
    return _save(fig, folder, "damage_ndvi_killzones", cfg)


# ---------- where the citrus is, and citrus damage in the main citrus area ----------
def main_citrus_window(citrus, cfg):
    """Square window (zoom_citrus_km wide) centred where citrus is densest. Returns .sel() slices."""
    res = cfg["run"]["res"]
    size = max(1, int(cfg["plots"]["zoom_citrus_km"] * 1000 / res))
    density = ndimage.uniform_filter(citrus.values.astype("float32"), size=size, mode="constant")
    iy, ix = np.unravel_index(np.argmax(density), density.shape)
    cx, cy, half = float(citrus.x[ix]), float(citrus.y[iy]), size * res / 2
    return dict(x=slice(cx - half, cx + half), y=slice(cy + half, cy - half))


def _outlines(mask, cfg):
    """Outlines of a boolean raster (e.g. the citrus mask) as polygons in km."""
    v = mask.values.astype("uint8")
    geoms = [shape(g) for g, _ in shapes(v, mask=v > 0, transform=mask.odc.geobox.affine)]
    g = gpd.GeoSeries(geoms, crs=cfg["run"]["crs"])
    return g.scale(1e-3, 1e-3, origin=(0, 0))                    # m -> km


def _rgb(bands, cfg):
    """True-colour RGB array (contrast-stretched), NaN shown as white."""
    rgb = np.dstack([bands[b] for b in ("red", "green", "blue")])
    lo, hi = np.nanpercentile(rgb, cfg["plots"]["rgb_stretch"])
    return np.nan_to_num(np.clip((rgb - lo) / (hi - lo), 0, 1), nan=1.0)


def citrus_location(rgb_bands, citrus, win, cfg, folder):
    """True-colour image with the CDL citrus mask outlined in black: whole area and main citrus area."""
    c = cfg["plots"]["colors"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.6), facecolor=c["surface"],
                             gridspec_kw={"width_ratios": [1, 1.25]})
    # whole area (downsampled)
    small = {b: thumb(da, cfg)[0] for b, da in rgb_bands.items()}
    ext = thumb(rgb_bands["red"], cfg)[1]
    _style(axes[0], cfg, "Whole study area")
    axes[0].imshow(_rgb(small, cfg), extent=ext, interpolation="nearest")
    _outlines(citrus, cfg).boundary.plot(ax=axes[0], color=c["citrus_outline"], lw=0.4)
    x0, x1, y0, y1 = (win["x"].start / 1e3, win["x"].stop / 1e3, win["y"].stop / 1e3, win["y"].start / 1e3)
    axes[0].plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color=c["zoom_box"], lw=1.5)
    axes[0].set_xlim(ext[0], ext[1])
    axes[0].set_ylim(ext[2], ext[3])
    # main citrus area (full resolution)
    part = {b: da.sel(**win) for b, da in rgb_bands.items()}
    ext2 = thumb(part["red"], cfg)[1]
    _style(axes[1], cfg, "Main citrus area (box on the left)")
    axes[1].imshow(_rgb({b: thumb(da, cfg)[0] for b, da in part.items()}, cfg), extent=ext2,
                   interpolation="nearest")
    _outlines(citrus.sel(**win), cfg).boundary.plot(ax=axes[1], color=c["citrus_outline"], lw=0.8)
    axes[1].set_xlim(ext2[0], ext2[1])
    axes[1].set_ylim(ext2[2], ext2[3])
    for ax in axes:
        ax.set_xlabel("Easting (km, UTM 17N)", fontsize=8, color=c["muted"])
        ax.set_ylabel("Northing (km)", fontsize=8, color=c["muted"])
    handles = [Patch(facecolor="none", edgecolor=c["citrus_outline"], lw=1.2, label="citrus (USDA CDL 2025)"),
               Patch(facecolor="none", edgecolor=c["zoom_box"], lw=1.5, label="main citrus area")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"Where the citrus is - true colour, {cfg['plots']['rgb_date']}", x=0.01, ha="left",
                 fontsize=12, color=c["text"])
    fig.subplots_adjust(bottom=0.1, wspace=0.18)
    return _save(fig, folder, "citrus_location", cfg)


def damage_citrus_zoom(damages, citrus, land, kz, win, cfg, folder):
    """Citrus-only damage classes in the main citrus area, both dates, kill zones outlined."""
    c = cfg["plots"]["colors"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), sharey=True, facecolor=c["surface"])
    for ax, (date, dmg) in zip(axes, damages.items()):
        d = dmg.sel(**win)
        codes, cmap = _damage_rgb(d.values, citrus.sel(**win).values, land.sel(**win).values, cfg)
        arr, ext = thumb(d.copy(data=codes), cfg)
        _style(ax, cfg, date)
        ax.imshow(arr, extent=ext, cmap=cmap, vmin=-0.5, vmax=5.5, interpolation="nearest")
        if len(kz):
            g = kz.copy()
            g["geometry"] = g.geometry.scale(1e-3, 1e-3, origin=(0, 0))
            g.boundary.plot(ax=ax, color=c["kill_zone"], lw=0.8)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_xlabel("Easting (km, UTM 17N)", fontsize=8, color=c["muted"])
        ax.set_ylabel("")
    axes[0].set_ylabel("Northing (km)", fontsize=8, color=c["muted"])
    handles = [Patch(color=c["citrus_ok"], label="citrus, no damage")] + \
              [Patch(color=col, label=lab) for col, lab in zip(c["damage"], ["mild", "moderate", "severe"])] + \
              [Patch(facecolor="none", edgecolor=c["kill_zone"], lw=1.5, label="kill zone"),
               Patch(color=c["land"], label="not citrus / not assessed")]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=9)
    fig.suptitle("NDVI damage classes - citrus only, main citrus area", x=0.01, ha="left",
                 fontsize=12, color=c["text"])
    fig.subplots_adjust(bottom=0.12, wspace=0.08)
    return _save(fig, folder, "damage_citrus_main_area", cfg)
