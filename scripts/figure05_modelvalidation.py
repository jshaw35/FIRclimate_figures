"""
Figure 2: Demonstrate that CESM is adequately representative of year 1 PREFIRE observations
"""
# %%
import os
import glob
import xarray as xr
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from data_config import DATA_ROOT, OUTPUT_ROOT, get_data_file

def label_PREFIRE(
    ds,
):
    # Label with the wavelength
    channel_wavelengths = [5.9, 11.8, 23.7, 24.5]
    ds["channel"] = channel_wavelengths
    ds = ds.assign_coords(wnum=("channel", 1e4 / np.array(channel_wavelengths)))

    return ds


def to_png(file, filename, loc=None, dpi=200, ext='png', **kwargs):
    '''
    Simple function for one-line saving.
    Saves to OUTPUT_ROOT by default
    '''
    if loc is None:
        loc = OUTPUT_ROOT
    output_dir = loc
    full_path = '%s%s.%s' % (output_dir, filename, ext)

    if not os.path.exists(output_dir + filename):
        file.savefig(
            full_path,
            format=ext,
            dpi=dpi,
            **kwargs,
        )


# %%
if __name__ == "__main__":
    # Load PREFIRE observations
    load_dir = DATA_ROOT
    prefire_filepath = load_dir / "PREFIRE" / "PREFIRE_zonal_stats.zarr"
    prefire_ds = xr.open_zarr(str(prefire_filepath))
    cesm_rttov_vars = [
        "rttov_rad_clear_inst001",
        "rttov_rad_total_inst001",
    ]
    cesm_cases = [
        "20241014_154735.FHIST.f09_f09_mg17.cesm2.1.5_port",
        "20250616_103133.FHIST.f09_f09_mg17.cesm2.1.5_port_SSP585branch_PREFIRE",
    ]
    file_dict = {}
    for _var in cesm_rttov_vars:
        cesm_files = []
        for case in cesm_cases:
            cesm_files += glob.glob(str(load_dir / "CESM2" / f"{case}*{_var}*.zarr"))
        file_dict[_var] = cesm_files
    cesm_total_ds = xr.open_mfdataset(file_dict["rttov_rad_total_inst001"], combine="nested")
    cesm_clear_ds = xr.open_mfdataset(file_dict["rttov_rad_clear_inst001"], combine="nested")

    # %%
    # Load zonal mean PREFIRE' data
    # RRTMG fluxes
    zonal_load_dir = DATA_ROOT
    prefireprime_zonal_rrtmg_file = zonal_load_dir / "FIRfraction_PREFIREprime.nc"
    prefireprime_zonal_rttov_file = zonal_load_dir / "zonalmeans_PREFIREprime.zarr"
    prefireprime_monthly_allsky_file = zonal_load_dir / "monthlymeans_rttov_bt_total_inst001_PREFIREprime.zarr"

    zonal_rrtmg_prefireprime_ds = xr.open_dataset(str(prefireprime_zonal_rrtmg_file))
    zonal_rttov_prefireprime_ds = xr.open_zarr(str(prefireprime_zonal_rttov_file))
    prefireprime_monthly_allsky_rttov_ds = xr.open_zarr(str(prefireprime_monthly_allsky_file))
    # Fix naming for consistency...
    zonal_rrtmg_prefireprime_ds["LW_band"] = ['AW', 'FIR', 'MIR', 'OLR']
    
    # %%
    # Compute fields for PREFIRE' simulation
    cesm_ds_monthly = prefireprime_monthly_allsky_rttov_ds["rttov_bt_total_inst001"].sel(lat=slice(-84, 84))
    cesm_ds_year = cesm_ds_monthly.mean(dim=["month"])
    prefireprime_seasonal_cycle_ds = cesm_ds_monthly.mean(dim="lon").sel(month=cesm_ds_monthly["month"].isin([6, 7, 8])).mean(dim="month") - cesm_ds_monthly.mean(dim="lon").sel(month=cesm_ds_monthly["month"].isin([12, 1, 2])).mean(dim="month")
    prefireprime_zonalvariance_ds = cesm_ds_year.std(dim="lon")
    prefireprime_anomaly_ds = cesm_ds_year.mean("lon") - cesm_ds_year.weighted(np.cos(np.deg2rad(cesm_ds_year.lat))).mean(dim=["lat", "lon"])
    
    # %%
    # Compute zonal-mean BT anomalies relative to the area-weighted global mean
    prefire_globalmean = prefire_ds["brightness_temperature"].sel(operation="mean").weighted(np.cos(np.deg2rad(prefire_ds.lat))).mean("lat")
    prefire_anomaly = prefire_ds["brightness_temperature"].sel(operation="mean") - prefire_globalmean
    prefire_zonal_stddev = prefire_ds["brightness_temperature"].sel(operation="stddev")

    cesm_globalmean = cesm_total_ds["mean"].weighted(np.cos(np.deg2rad(cesm_total_ds.lat))).mean("lat")
    cesm_anomaly = cesm_total_ds["mean"] - cesm_globalmean
    cesm_zonal_stddev = cesm_total_ds["stddev"]
    # %%
    representative_channels = [
        6,  # Channel 6, 5.89um (WV)
        13, # Channel 13, 11.77um
        28, # Channel 28, 24.72
    ]
    representative_wavelengths = [
        5.9,  # Channel 6, 5.89um (WV)
        11.8, # Channel 13, 11.77um
        24.5, # Channel 28, 24.72
    ]
    # %%
    # Load data conditioned on CWV
    conditional_load_dir = f"{load_dir}/PREFIRE_conditional_correlations/"
    data_paths = glob.glob(f"{conditional_load_dir}/*.binned.zarr/")
    data_paths.sort()

    data_dict = {}
    for path in data_paths:
        _var = Path(path).stem.split('.')[-3]
        # PREFIRE
        if _var in ["rttov_rad_clear_inst001", "rttov_rad_total_inst001", "rttov_rad_cre_inst001"]:
            data_dict[_var] = label_PREFIRE(xr.open_zarr(path))
        # FORUM
        elif _var in ["rttov_rad_clear_inst002", "rttov_rad_total_inst002", "rttov_rad_cre_inst002"]:
            continue
        # RRTMG-LW
        else:
            data_dict[_var] = xr.open_zarr(path)
    data_vars = list(data_dict.keys())

    #%%
    # Figure 1, seasonal cycle and zonal mean values
    # Clear-sky values here on CURC: /pl/active/kaygroup/PREFIRE/PREFIRE_SAT2_1B-RAD_R01_Gridded_and_Masked/
    year_slice = slice(2015, 2035) # bracket the CESM year with 10 years on either side
    colors = sns.color_palette("colorblind")
    fontsize = 12
    channel_labels = ["Mid-Infrared WV", "Atmospheric Window", "Far-Infrared"]
    ylims_mean = [(-17, 10), (-42, 18), (-22, 11)]
    ylims_seasonal = [(-15, 15), (-30, 30), (-20, 20)]
    ylims_stddev = [(0, 7), (0, 12), (0, 9)]
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 6.5))
    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    # Plot zonal-mean BT anomalies in the first row
    axs = axes[0, :]
    for ax, channel, color, ylims in zip(axs, representative_channels, colors, ylims_mean):
        # Plot PREFIRE
        prefire_data = prefire_anomaly.sel(channel=channel)
        ax.plot(
            prefire_data.lat,
            prefire_data,
            label=f"PREFIRE",
            color=color,
            alpha=0.75,
            zorder=10,
        )
        # Plot CESM
        cesm_data = cesm_anomaly.sel(channel=channel, year=year_slice)
        for year in cesm_data.year:
            ax.plot(
                cesm_data.lat,
                cesm_data.sel(year=year),
                color="black",
                alpha=0.5,
                linewidth=0.25,
                zorder=2,
            )
        # Plot PREFIRE'
        prefireprime_data = prefireprime_anomaly_ds.sel(channel=channel)
        ax.plot(
            prefireprime_data.lat,
            prefireprime_data,
            color="black",
            alpha=0.75,
            linewidth=1.0,
            linestyle="dashed",
            zorder=2,
        )
        ax.set_ylabel("BT Anomaly (K)", fontsize=fontsize+2)
        ax.set_ylim(ylims)
        # Create a filler plot for the legend

    # Second row is the zonal seasonal cycle (JJA - DJF)
    prefire_seasonal_cycle_ds = prefire_ds["JJA-DJF"]
    cesm_seasonal_cycle_ds = cesm_total_ds["JJA-DJF"]
    axs = axes[1, :]
    for ax, channel, color, ylims in zip(axs, representative_channels, colors, ylims_seasonal):
        # Plot PREFIRE
        prefire_data = prefire_seasonal_cycle_ds.sel(channel=channel)
        ax.plot(
            prefire_data.lat,
            prefire_data,
            label=f"PREFIRE",
            color=color,
            alpha=0.75,
            zorder=10,
        )
        # Plot CESM
        cesm_data = cesm_seasonal_cycle_ds.sel(channel=channel, year=year_slice)
        for year in cesm_data.year:
            ax.plot(
                cesm_data.lat,
                cesm_data.sel(year=year),
                color="black",
                alpha=0.5,
                linewidth=0.25,
                zorder=2,
            )
        # Plot PREFIRE'
        prefireprime_data = prefireprime_seasonal_cycle_ds.sel(channel=channel)
        ax.plot(
            prefireprime_data.lat,
            prefireprime_data,
            color="black",
            alpha=0.75,
            linewidth=1.0,
            linestyle="dashed",
            zorder=2,
        )
        ax.set_ylabel("Seasonal Cycle \n Amplitude (K)", fontsize=fontsize+2)
        ax.set_ylim(ylims)
        ax.set_xlabel("Latitude", fontsize=fontsize)

    for ax in axes.flat:
        ax.grid(True)
        ax.set_xlim(-90, 90)
        ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
        ax.set_xticklabels(["90S", "60S", "30S", "Eq.", "30N", "60N", "90N"])

    letters = list("abcdefghijkl")
    # repeat each wavelength twice so labels map to axes in order: [w1, w1, w2, w2, w3, w3]
    rrtmg_labels = [""]
    panel_letters = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
    for ax, letter in zip(
        axes.flatten(),
        letters,
    ):
        ax.text(
            0.02,
            0.98,
            f"({letter})",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )
    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")
    for ax, label in zip(axes[0, :].flat, channel_labels):
        ax.set_title(f"{label}", fontsize=fontsize+2, fontweight="bold")

    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig2_modelvalidation",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
    # Figure 2: Zonal variance for supplement
    year_slice = slice(2015, 2035) # bracket the CESM year with 10 years on either side
    colors = sns.color_palette("colorblind")
    fontsize = 12
    channel_labels = ["Mid-Infrared WV", "Atmospheric Window", "Far-Infrared"]
    ylims_stddev = [(0, 7), (0, 12), (0, 9)]
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 3))
    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    # Plot zonal-mean BT anomalies in the first row
    axs = axes#[0, :]
    for ax, channel, color, ylims in zip(axs, representative_channels, colors, ylims_stddev):
        # Plot PREFIRE
        prefire_data = prefire_zonal_stddev.sel(channel=channel)
        ax.plot(
            prefire_data.lat,
            prefire_data,
            label=f"PREFIRE",
            color=color,
            alpha=0.75,
            zorder=10,
        )
        # Plot CESM
        cesm_data = cesm_zonal_stddev.sel(channel=channel, year=year_slice)
        for year in cesm_data.year:
            ax.plot(
                cesm_data.lat,
                cesm_data.sel(year=year),
                color="black",
                alpha=0.5,
                linewidth=0.25,
                zorder=2,
            )
        # Plot PREFIRE'
        prefireprime_data = prefireprime_zonalvariance_ds.sel(channel=channel)
        ax.plot(
            prefireprime_data.lat,
            prefireprime_data,
            color="black",
            alpha=0.75,
            linewidth=1.0,
            linestyle="dashed",
            zorder=2,
        )
        ax.set_ylabel("BT Std. Dev. (K)", fontsize=fontsize+2)
        ax.set_xlabel("Latitude", fontsize=fontsize)
        ax.set_ylim(ylims)
        # Create a filler plot for the legend

    letters = list("abcdefghijkl")
    # repeat each wavelength twice so labels map to axes in order: [w1, w1, w2, w2, w3, w3]
    row_labels = ["MIR", "AW", "FIR"]
    rrtmg_labels = [""]
    panel_letters = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
    for ax, letter in zip(
        axes.flatten(),
        letters,
    ):
        ax.text(
            0.02,
            0.98,
            f"({letter})",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )
    for ax in axes[1:]:
        ax.set_ylabel("")
    for ax, label in zip(axes.flat, channel_labels):
        ax.set_title(f"{label}", fontsize=fontsize+2, fontweight="bold")

    for ax in axes.flat:
        ax.grid(True)
        ax.set_xlim(-90, 90)
        ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
        ax.set_xticklabels(["90S", "60S", "30S", "Eq.", "30N", "60N", "90N"])

    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig_supp_modelvalidation_zonalvariance",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    # %%
