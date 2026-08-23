"""
Plot the FIR fraction for both the PREFIRE SAT2 PRIME mission (CESM2 and obs) and the 21st century generally.
"""
# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from pathlib import Path

# Ensure project-root modules are importable when this file is run directly.
try:
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
except NameError:
    PROJECT_ROOT = Path.cwd().resolve()  # or hardcode/adjust as needed
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from data_config import DATA_ROOT, OUTPUT_ROOT, get_data_file

# %%
def to_png(file, filename, loc=None, dpi=200, ext='png', **kwargs):
    '''
    Simple function for one-line saving.
    Saves to OUTPUT_ROOT by default
    '''
    if loc is None:
        loc = OUTPUT_ROOT
    output_dir = loc
    full_path = '%s/%s.%s' % (output_dir, filename, ext)

    if not os.path.exists(full_path):
        file.savefig(
            full_path,
            format=ext,
            dpi=dpi,
            **kwargs,
        )


# %%

if __name__ == "__main__":
    data_loaddir = DATA_ROOT
    prefireyear_firfrac_ds = xr.open_dataset(str(data_loaddir / "FIRfraction_PREFIREprime.nc"))
    century_firfrac_ds = xr.open_dataset(str(data_loaddir / "FIRfraction_PREFIREssp585.nc")).sel(time=slice("2015", "2089"))
    obs_firfrac_ds = xr.open_dataset(str(data_loaddir / "PREFIRE_Polar_Fraction.nc"))

    # %%
    # Create first panel
    fontsize = 14
    cfg1 = {
        "LU_TOA": {
            "label": "CESM2 All-sky",
            "color": "blue",
            "linestyle": "solid",
        },
        "LUC_TOA": {
            "label": "CESM2 Clear-sky",
            "color": "blue",
            "linestyle": "dashed",
        },
    }

    fig, axs = plt.subplots(1, 3, figsize=(10, 4), sharey=True)
    # First panel first year PREFIRE
    ax = axs[0]
    # Plot obs
    ax.plot(
        obs_firfrac_ds["FRAC_Arctic_Sat1"].lat,
        obs_firfrac_ds["FRAC_Arctic_Sat1"],
        label="PREFIRE All-sky",
        color="red",
        linestyle="solid",
    )
    ax.plot(
        obs_firfrac_ds["FRAC_Anta_Sat2"].lat,
        obs_firfrac_ds["FRAC_Anta_Sat2"],
        color="red",
        linestyle="solid",
    )
    for _var in prefireyear_firfrac_ds:
        _data = prefireyear_firfrac_ds[_var]
        _fir_ratio = _data.sel(LW_band="FIR2") / _data.sel(LW_band="OLR")

        ax.plot(
            _fir_ratio.lat,
            _fir_ratio,
            label=cfg1[_var]["label"],
            color=cfg1[_var]["color"],
            linestyle=cfg1[_var]["linestyle"],
        )
    ax.set_ylim(0.4, 0.7)
    ax.set_ylabel("Far-Infrared Fraction", fontsize=fontsize)

    # Second panel seasonal patterns
    ax = axs[1]
    for _var in century_firfrac_ds:
        _data = century_firfrac_ds[_var].sel(time=slice("2015", "2035"))
        _data_monthly = _data.groupby("time.season").mean(dim="time").sel(season=["DJF", "JJA"])
        seasons = _data_monthly["season"].values
        cmap = sns.color_palette("colorblind", n_colors=4)
        for i,_season in enumerate(seasons):
            _fir_ratio = _data_monthly.sel(season=_season)
            color = cmap[i]

            label = None
            if _var == "LU_TOA":
                label = _season

            ax.plot(
            _fir_ratio.lat,
            _fir_ratio,
            color=color,
            linestyle=cfg1[_var]["linestyle"],
            label=label,
            alpha=1,
            linewidth=1,
            )

    # Third panel 21 century change
    ax = axs[2]
    for _var in century_firfrac_ds:
        _data = century_firfrac_ds[_var]
        _data_annual = _data.groupby("time.year").mean(dim="time")
        years = _data_annual["year"].values
        cmap = plt.get_cmap("viridis")
        norm = plt.Normalize(vmin=years.min(), vmax=years.max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        for _year in years[::5]:
            _fir_ratio = _data_annual.sel(year=_year)
            color = cmap(norm(_year))

            ax.plot(
            _fir_ratio.lat,
            _fir_ratio,
            color=color,
            linestyle=cfg1[_var]["linestyle"],
            alpha=0.85,
            linewidth=0.7,
            )
        ax.plot(
            _fir_ratio.lat,
            np.nan * _fir_ratio.lat,
            color=cmap(norm(years[0])),
            linestyle=cfg1[_var]["linestyle"],
            alpha=0.85,
            linewidth=0.7,
            label=cfg1[_var]["label"].split(" ")[-1],
        )
        # add a single colorbar for the year colormap (attach once to the figure)
        if not getattr(fig, "_fir_cb_added", False):
            cbar = fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.046, pad=0.04)
            cbar.ax.invert_yaxis()
            cbar.set_label("Year", fontsize=fontsize)
            cbar.ax.set_yticks([2015, 2030, 2045, 2060, 2075, 2090])
            fig._fir_cb_added = True

    for ax in axs:
        ax.legend(loc="upper right")
        ax.set_xlim(-90, 90)
        ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
        ax.set_xticklabels(["90S", "60S", "30S", "Eq.", "30N", "60N", "90N"])
        ax.set_xlabel("Latitude", fontsize=fontsize)
        ax.grid(True, alpha=0.5, color='black')
        ax.set_facecolor("#CCCCCC")

    panel_letters = ["a", "b", "c"]
    for ax, letter in zip(axs.flat[:6], panel_letters):
        ax.text(0.03, 0.98, f"({letter})", transform=ax.transAxes, ha="left", va="top", fontsize=fontsize, weight="bold")

    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig07_FIRfraction_2015start",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    # %%
    fir_fraction_global = century_firfrac_ds.weighted(np.cos(np.deg2rad(century_firfrac_ds["lat"]))).mean(dim="lat")
    fir_fraction_global_annual = fir_fraction_global.groupby("time.year").mean(dim="time")

    # Regionally    
    fir_fraction_arctic = century_firfrac_ds.sel(lat=slice(60, None)).weighted(np.cos(np.deg2rad(century_firfrac_ds["lat"]))).mean(dim="lat")
    fir_fraction_arctic_annual = fir_fraction_arctic.groupby("time.year").mean(dim="time")
    fir_fraction_antarctic = century_firfrac_ds.sel(lat=slice(None, -60)).weighted(np.cos(np.deg2rad(century_firfrac_ds["lat"]))).mean(dim="lat")
    fir_fraction_antarctic_annual = fir_fraction_antarctic.groupby("time.year").mean(dim="time")
    fir_fraction_tropic = century_firfrac_ds.sel(lat=slice(-30, 30)).weighted(np.cos(np.deg2rad(century_firfrac_ds["lat"]))).mean(dim="lat")
    fir_fraction_tropic_annual = fir_fraction_tropic.groupby("time.year").mean(dim="time")

    midlat_mask = (np.abs(century_firfrac_ds["lat"]) < 60) & (np.abs(century_firfrac_ds["lat"]) > 30)
    fir_fraction_midlat = century_firfrac_ds.sel(lat=midlat_mask).weighted(np.cos(np.deg2rad(century_firfrac_ds["lat"]))).mean(dim="lat")
    fir_fraction_midlat_annual = fir_fraction_midlat.groupby("time.year").mean(dim="time")

    # %%
    for region_label, region_data in zip(
        ["Global", "Arctic", "Antarctic", "Tropic", "Midlatitudes"],
        [fir_fraction_global_annual, fir_fraction_arctic_annual, fir_fraction_antarctic_annual, fir_fraction_tropic_annual, fir_fraction_midlat_annual],
    ):
        region_fir_fraction_2015_2034 = region_data.sel(year=slice("2015", "2034")).mean(dim="year")
        region_fir_fraction_2070_2089 = region_data.sel(year=slice("2070", "2089")).mean(dim="year")
        percent_change = 100 * (region_fir_fraction_2070_2089 - region_fir_fraction_2015_2034) / region_fir_fraction_2015_2034

        luc_2015_2034 = 100*region_fir_fraction_2015_2034["LUC_TOA"].values
        luc_2070_2089 = 100*region_fir_fraction_2070_2089["LUC_TOA"].values
        luc_pct_change = percent_change["LUC_TOA"]
        print(f"{region_label} Clear-sky FIR fraction went from {luc_2015_2034:.1f} to {luc_2070_2089:.1f} ({luc_pct_change:.1f} % change)")
        
        lu_2015_2034 = 100*region_fir_fraction_2015_2034["LU_TOA"].values
        lu_2070_2089 = 100*region_fir_fraction_2070_2089["LU_TOA"].values
        lu_pct_change = percent_change["LU_TOA"]
        print(f"{region_label} All-sky FIR fraction went from {lu_2015_2034:.1f} to {lu_2070_2089:.1f} ({lu_pct_change:.1f} % change)")
    print()

    # %%