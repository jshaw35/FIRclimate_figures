"""
Figure 3: Show the physics controls on the spectral cloud radiative effect
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

    zonal_rrtmg_prefireprime_ds = xr.open_dataset(str(prefireprime_zonal_rrtmg_file))
    zonal_rttov_prefireprime_ds = xr.open_zarr(str(prefireprime_zonal_rttov_file))
    # Fix naming for consistency...
    zonal_rrtmg_prefireprime_ds["LW_band"] = ['AW', 'FIR', 'MWV', 'OLR']
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
    conditional_load_dir = load_dir / "PREFIRE_conditional_correlations"
    data_paths = glob.glob(str(conditional_load_dir / "*.binned.zarr/"))
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

    # %%
    # Visualize the FIR cloud radiative effect
    # Panel 1: zonal-mean CRE for different spectral regions (units flux or BT?)
    # Panel 2: CRE vs. CWV for the different spectral channels (I think thermal contrast (BT difference) is the most intuitive here)
    # Panel 3: Plot the emission temperature and the cloud temperature vs. CWV
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
    axs = axs.flat
    fig.subplots_adjust(wspace=0.47)

    plot_var = "global_mean"
    year = 2025
    colors = sns.color_palette("colorblind")
    fontsize = 12
    prefire_rttov_vars = ["rttov_rad_clear_inst001", "rttov_rad_total_inst001", "rttov_rad_cre_inst001", "channel"]
    prefire_rttov_channels = [6, 13, 28]
    prefire_rttov_labels = ["MWV", "AW", "FIR"]

    # Panel 1: zonal-mean flux CRE for different spectral regions
    ax = axs[0]
    twin_ax = False

    bt_cre_data = zonal_rrtmg_prefireprime_ds["LUC_TOA"] - zonal_rrtmg_prefireprime_ds["LU_TOA"]
    clt_data = zonal_rttov_prefireprime_ds["CLTMODIS"]

    for _chan, _color, _label in zip(prefire_rttov_labels, colors, prefire_rttov_labels):
        chan_cre_data = bt_cre_data.sel({"LW_band":_chan})
        ax.plot(
            chan_cre_data.lat,
            chan_cre_data,
            c=_color,
            label=_label,
        )

    ax.set_xlim(-90, 90)
    ax.set_ylim(0, 12)
    ax.set_xlabel("Latitude", fontsize=fontsize)
    ax.set_xticks([-90, -60, -30, 0 , 30, 60, 90])
    ax.set_ylabel("Cloud Radiative Effect (Wm$^{-2}$)", fontsize=fontsize)
    ax.legend()

    if twin_ax:
        axb = ax.twinx()
        axb.plot(
            clt_data.lat,
            clt_data,
            linestyle="dashed",
            color="grey",
        )
        # make ticks and tick labels grey to match the label
        axb.set_ylim(0, 100)
        axb.set_ylabel("Cloud Fraction (%)", fontsize=fontsize, color='grey')
        axb.tick_params(axis='y', colors='grey')

    # Panel 2: zonal-mean CRE for different spectral regions (units flux or BT?)
    ax = axs[1]
    axb = ax.twinx()

    bt_cre_data = zonal_rttov_prefireprime_ds["rttov_bt_clear_inst001"] - zonal_rttov_prefireprime_ds["rttov_bt_total_inst001"]
    clt_data = zonal_rttov_prefireprime_ds["CLTMODIS"]

    axb.plot(
        clt_data.lat,
        clt_data,
        linestyle="dashed",
        color="grey",
        label="Cloud Fraction"
    )

    for _chan, _color, _label in zip(prefire_rttov_channels, colors, prefire_rttov_labels):
        chan_bt_cre_data = bt_cre_data.sel({prefire_rttov_vars[-1]:_chan})
        ax.plot(
            chan_bt_cre_data.lat,
            chan_bt_cre_data,
            c=_color,
            label=_label,
        )
    ax.set_xlim(-90, 90)
    ax.set_ylim(0, 20)
    ax.set_yticks(np.arange(0, 21 ,4))
    ax.set_xlabel("Latitude", fontsize=fontsize)
    ax.set_xticks([-90, -60, -30, 0 , 30, 60, 90])
    ax.set_ylabel("Cloud Radiative Effect (K)", fontsize=fontsize)

    axb.set_ylim(0, 100)
    axb.set_ylabel("Cloud Fraction (%)", fontsize=fontsize, color='grey')
    axb.tick_params(axis='y', colors='grey')

    # Panel 3: CRE vs. CWV for the different spectral channels (I think thermal contrast (BT difference) is the most intuitive here)
    ax = axs[2]
    axb = ax.twinx()
    data = data_dict["rttov_bt_clear_inst001"] - data_dict["rttov_bt_total_inst001"]
    clt_data = data_dict["CLTMODIS"]

    axb.plot(
        clt_data.bin_center,
        clt_data[plot_var],
        linestyle="dashed",
        color="grey",
        label="Cloudtop Temperature"
    )
    for _chan, _color in zip(prefire_rttov_channels, colors):
        cre_data = data[plot_var].sel({prefire_rttov_vars[-1]:_chan})
        ax.plot(
            cre_data.bin_center,
            cre_data,
            linestyle="solid",
            color=_color
        )

    ax.set_xlim(0, 70)
    ax.set_ylim(0, 52)
    ax.set_xlabel("PWV (mm)", fontsize=fontsize)
    ax.set_ylabel("Cloud Radiative Effect (K)", fontsize=fontsize)

    axb.set_ylim(0, 104)
    axb.set_ylabel("Cloud Fraction (%)", fontsize=fontsize, color='grey')
    axb.tick_params(axis='y', colors='grey')

    # Panel 4: Plot the emission temperature and the cloud temperature vs. CWV
    ax = axs[3]
    axb = ax.twinx()
    tct_data = data_dict["TCTMODIS_norm"]
    bt_data = data_dict["rttov_bt_clear_inst001"]

    axb.plot(
        tct_data.bin_center,
        tct_data[plot_var],
        linestyle="dashed",
        color="grey",
    )
    for _chan, _color, _label in zip(prefire_rttov_channels, colors, prefire_rttov_labels):
        data = bt_data[plot_var].sel({prefire_rttov_vars[-1]:_chan})
        ax.plot(
            data.bin_center,
            data,
            linestyle="solid",
            color=_color,
            label=_label,
        )

    ax.set_xlim(0, 70)
    ax.set_ylim(210, 300)
    ax.set_xlabel("PWV (mm)", fontsize=fontsize)
    ax.set_ylabel("Clear-sky Brightness Temperature (K)", fontsize=fontsize)

    axb.set_ylim(210, 300)
    axb.set_ylabel("Cloudtop Temperature (K)", fontsize=fontsize, color='grey')
    # make ticks and tick labels grey to match the label
    axb.tick_params(axis='y', colors='grey')

    panel_letters = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"]
    for ax, letter in zip(
        axs,
        panel_letters,
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
        ax.grid(True)

    #%%

    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig3_cloudintuition",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    # %%