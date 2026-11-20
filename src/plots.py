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
    fig.suptitle("Median index value per clear date (vegetated pixels)", x=0.01, ha="left",
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
    fig.suptitle(f"{name}: baseline and after the freeze (cloud-masked)", x=0.01, ha="left", fontsize=12)
    return _save(fig, folder, f"map_{name.lower()}", cfg)


def change_maps(name, panels, lim, cfg, folder):
    """Index change vs baseline on both damage dates (red = decrease, blue = increase)."""
    fig, axes = _map_row(len(panels), cfg)
    for ax, (title, (arr, ext)) in zip(axes, panels.items()):
        _style(ax, cfg, title)
        im = ax.imshow(arr, extent=ext, cmap="RdBu", vmin=-lim, vmax=lim, interpolation="nearest")
    cb = fig.colorbar(im, ax=axes, shrink=0.7)
    cb.set_label(f"{name} change", fontsize=9)
    fig.suptitle(f"{name} change vs baseline (vegetated pixels)", x=0.01, ha="left", fontsize=12)
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
    fig.subplots_adjust(bottom=0.12, wspace=0.15)
    return _save(fig, folder, "damage_ndvi_killzones", cfg)
