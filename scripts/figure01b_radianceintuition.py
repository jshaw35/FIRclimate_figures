"""
Figure 1 (lower panel): Show how the different representative channels
"""
# %%
import os
import glob
import xarray as xr
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Ensure project-root modules are importable when this file is run directly.
try:
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
except NameError:
    PROJECT_ROOT = Path.cwd().resolve()  # or hardcode/adjust as needed
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from data_config import DATA_ROOT, OUTPUT_ROOT

# %%
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
    full_path = '%s/%s.%s' % (output_dir, filename, ext)

    if not os.path.exists(full_path):
        file.savefig(
            full_path,
            format=ext,
            dpi=dpi,
            **kwargs,
        )


def transform_cwv_axis(cwv_values, threshold=5, low_scale=2.0, high_scale=0.4):
    """
    Apply piecewise linear transformation to CWV values for better visualization.
    
    Parameters
    ----------
    cwv_values : array-like
        CWV values in mm
    threshold : float
        Breakpoint between low and high scale regions (default: 5 mm)
    low_scale : float
        Expansion factor for values 0 to threshold (default: 2.0)
    high_scale : float
        Compression factor for values above threshold (default: 0.4)
    
    Returns
    -------
    transformed : array-like
        Transformed CWV values
    """
    cwv_array = np.asarray(cwv_values)
    transformed = np.zeros_like(cwv_array, dtype=float)
    
    # Low CWV region (0 to threshold): stretched
    low_mask = cwv_array <= threshold
    transformed[low_mask] = cwv_array[low_mask] * low_scale
    
    # High CWV region (threshold to 70): compressed
    high_mask = cwv_array > threshold
    max_cwv = np.max(cwv_array)
    high_range = max_cwv - threshold
    transformed[high_mask] = (threshold * low_scale) + (cwv_array[high_mask] - threshold) * high_scale
    
    return transformed


def rad2bt(fr, rad):
    """
    Convert radiance to brightness temperature.
    Based on matlab function rad2bt from Sergio Desouza-Machado (UMBC)
    The scalar calculation, for reference:
    bt = c2 * fr / log(1 + c1 * fr^3 / rad)

    Parameters
    ----------
    fr : array-like
        Wavenumbers in cm^-1 (1/cm).
    rad : array-like
        Radiances in mW/m^2/sr.

    Returns
    -------
    bt : array-like
        Brightness temperatures in K.
    """

    # Constants; values from NIST (CODATA98)
    c = 2.99792458e+08  # speed of light      299 792 458 m s-1
    h = 6.62606876e-34  # Planck constant     6.626 068 76 x 10-34 J s
    k = 1.3806503e-23   # Boltzmann constant  1.380 6503 x 10-23 J K-1

    # Compute radiation constants c1 and c2
    c1 = 2 * h * c**2 * 1e+11  # in mW m^2 / sr
    c2 = (h * c / k) * 100      # in K cm

    # Calculate brightness temperature
    with np.errstate(divide='ignore', invalid='ignore'):
        bt = c2 * fr / np.log(1 + (c1 * fr**3) / rad)

    return bt


# Use computed PREFIRE channel central wavelengths to get wavenumbers
PREFIRE_CHANNEL_WAVELENGTHS_UM = xr.open_dataarray(
    str(DATA_ROOT / "PREFIRE_TIRS2_SRF_central_wavelengths.nc")
)

def _prefire_channel_wavenumbers(channel_coord: xr.DataArray) -> xr.DataArray:
    """Return central wavenumbers (cm^-1) for PREFIRE channels."""

    channels = [int(ch) for ch in channel_coord.values]
    try:
        wavelengths = np.array([PREFIRE_CHANNEL_WAVELENGTHS_UM.sel(channel=ch) for ch in channels])
    except KeyError as exc:
        raise ValueError(f"Unsupported PREFIRE channel {exc.args[0]} for BT conversion") from exc
    wavenumbers = 10000.0 / wavelengths
    return xr.DataArray(
        wavenumbers,
        dims=["channel"],
        coords={"channel": channel_coord},
        name="wavenumber",
    )


def get_channel_wavenumbers(data_arr: xr.DataArray, instrument: str) -> xr.DataArray:
    """Infer per-channel wavenumbers for a radiance field."""

    instrument = instrument.lower()
    if "channel" not in data_arr.coords:
        raise ValueError("Radiance field lacks 'channel' coordinate required for BT conversion")
    if instrument == "prefire":
        return _prefire_channel_wavenumbers(data_arr["channel"])
    if instrument == "forum":
        return xr.DataArray(
            data_arr["channel"].astype(np.float64),
            dims=["channel"],
            coords={"channel": data_arr["channel"]},
            name="wavenumber",
        )
    raise ValueError(f"Unsupported instrument '{instrument}' for BT conversion")


# %%
if __name__ == "__main__":
    # Load PREFIRE observations
    load_dir = DATA_ROOT
    prefire_filepath = load_dir / "PREFIRE_zonal_stats.zarr"
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
            cesm_files += glob.glob(str(load_dir / f"{case}*{_var}*.zarr"))
        file_dict[_var] = cesm_files
    cesm_total_ds = xr.open_mfdataset(file_dict["rttov_rad_total_inst001"], combine="nested")
    cesm_clear_ds = xr.open_mfdataset(file_dict["rttov_rad_clear_inst001"], combine="nested")

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
    conditional_load_dir = load_dir
    mask_path = glob.glob(str(load_dir / "*binmasks.zarr"))
    data_paths = glob.glob(str(load_dir / "*.binned.zarr"))
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
    # Figure 1 lower panels
    # Plot with brightness temperatures
    colors = sns.color_palette("colorblind")
    fontsize = 12
    prefire_rttov_channels = [5.9, 11.8, 24.5]
    prefire_rttov_channel_indices = [6, 13, 28]
    prefire_rttov_channel_lambdas = PREFIRE_CHANNEL_WAVELENGTHS_UM.sel(channel=prefire_rttov_channel_indices).values
    prefire_rttov_channel_wnums = 1e4 / prefire_rttov_channel_lambdas
    rrtmg_channels = ["NIR", "AW", "FIR2"]
    panel_labels = ["Clear-sky", "All-sky", "CRE"]
    channel_labels = ["MWV", "AW", "FIR"]

    xlims = [0.125, 70]
    xticks = [0.125, 0.5, 1, 2, 5, 10, 20, 40, 70]

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(13, 4))
    fig.subplots_adjust(hspace=0.25, wspace=0.35)

    clearsky_prefire = data_dict["rttov_rad_clear_inst001"]
    allsky_prefire = data_dict["rttov_rad_total_inst001"]
    # Pre-process to BTs:
    clearsky_bt_prefire = rad2bt(prefire_rttov_channel_wnums, clearsky_prefire["global_mean"].isel(channel=[0, 1, 3]))
    allsky_bt_prefire = rad2bt(prefire_rttov_channel_wnums, allsky_prefire["global_mean"].isel(channel=[0, 1, 3]))

    for ax, data in zip(axes, [clearsky_bt_prefire, allsky_bt_prefire]):
        for prefire_chan, rrtmg_chan, wnum, label, color in zip(prefire_rttov_channels, rrtmg_channels, prefire_rttov_channel_wnums, channel_labels, colors):
            # Plot PREFIRE clear-sky radiance on primary y-axis
            chan_data = data.sel({"channel": prefire_chan})

            ax.plot(
                chan_data.bin_center,
                chan_data,
                label=label,
                linestyle="solid",
                color=color,
                linewidth=2
            )
        ax.set_xlim(*xlims)
        ax.set_xlabel("PWV (mm)", fontsize=fontsize)
        ax.set_ylabel("Brightness Temperature (K)", fontsize=fontsize+2)
        ax.set_xscale('log', base=2)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(i) for i in xticks])
        ax.grid(alpha=0.7, axis="both")

    axes[0].set_ylim(215, 295)
    axes[1].set_ylim(215, 295)
    axes[0].legend()


    # Add panel labels
    letters = ["c", "d"]
    for ax, letter, chan in zip(axes, letters, panel_labels):
        ax.text(
            0.0,
            1.06,
            f"({letter}) {chan}",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )

    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig1_lowerpanels_BTalt2",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
    # Figure 4: Clear-sky radiance and flux on same plot
    colors = sns.color_palette("colorblind")
    fontsize = 12
    prefire_rttov_channels = [5.9, 11.8, 24.5]
    rrtmg_channels = ["NIR", "AW", "FIR2"]
    channel_labels = ["Mid-Infrared WV", "Atmospheric Window", "Far-Infrared"]

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(13, 4))
    fig.subplots_adjust(hspace=0.25, wspace=0.45)

    clearsky_prefire = data_dict["rttov_rad_clear_inst001"]
    clearsky_rrtmg = data_dict["LUC_TOA"]

    for ax, prefire_chan, rrtmg_chan, color in zip(axes, prefire_rttov_channels, rrtmg_channels, colors):
        # Plot PREFIRE clear-sky radiance on primary y-axis
        chan_clear_data = clearsky_prefire["global_mean"].sel({"channel": prefire_chan})
        ax.plot(
            chan_clear_data.bin_center,
            chan_clear_data,
            label="Radiance",
            linestyle="solid",
            color=color,
            linewidth=2
        )
        ax.set_xlim(0, 70)
        ax.set_xlabel("PWV (mm)", fontsize=fontsize)
        ax.set_ylabel("Radiance (mWm$^{-2}$sr$^{-1}$ cm$^{-1}$)", fontsize=fontsize+2, color=color)
        ax.tick_params(axis='y', labelcolor=color)
        ax.grid(alpha=0.7, axis="both")

        # Create secondary y-axis for RRTMG flux
        ax2 = ax.twinx()
        chan_flux_data = clearsky_rrtmg["global_mean"].sel({"LW_band": rrtmg_chan})
        ax2.plot(
            chan_flux_data.bin_center,
            chan_flux_data,
            label="Flux",
            linestyle="dashed",
            color=color,
            linewidth=2
        )
        ax2.set_ylabel("Spectral Flux (Wm$^{-2}$)", fontsize=fontsize+2, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
        # Add legend
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc=[0.45, 0.08], fontsize=fontsize)

    # Add panel labels
    letters = ["a", "b", "c"]
    for ax, letter, chan in zip(axes, letters, channel_labels):
        ax.text(
            0.0,
            1.06,
            f"({letter}) {chan}",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )

    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig_supp_radflux_comparison",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
    # Unused alternative to figure 1 lower panels
    # Plot with brightness temperatures
    colors = sns.color_palette("colorblind")
    fontsize = 12
    prefire_rttov_channels = [5.9, 11.8, 24.5]
    prefire_rttov_channel_indices = [6, 13, 28]
    prefire_rttov_channel_lambdas = PREFIRE_CHANNEL_WAVELENGTHS_UM.sel(channel=prefire_rttov_channel_indices).values
    prefire_rttov_channel_wnums = 1e4 / prefire_rttov_channel_lambdas
    rrtmg_channels = ["NIR", "AW", "FIR2"]
    panel_labels = ["Clear-sky", "All-sky", "CRE"]
    channel_labels = ["MWV", "AW", "FIR"]

    xlims = [0.125, 70]
    xticks = [0.125, 0.5, 1, 2, 5, 10, 20, 40, 70]

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(13, 4))
    fig.subplots_adjust(hspace=0.25, wspace=0.45)

    clearsky_prefire = data_dict["rttov_rad_clear_inst001"]
    allsky_prefire = data_dict["rttov_rad_total_inst001"]
    # Pre-process to BTs:
    clearsky_bt_prefire = rad2bt(prefire_rttov_channel_wnums, clearsky_prefire["global_mean"].isel(channel=[0, 1, 3]))
    allsky_bt_prefire = rad2bt(prefire_rttov_channel_wnums, allsky_prefire["global_mean"].isel(channel=[0, 1, 3]))
    cre_bt_prefire = clearsky_bt_prefire - allsky_bt_prefire

    for ax, data in zip(axes, [clearsky_bt_prefire, allsky_bt_prefire, cre_bt_prefire]):
        for prefire_chan, rrtmg_chan, wnum, label, color in zip(prefire_rttov_channels, rrtmg_channels, prefire_rttov_channel_wnums, channel_labels, colors):
            # Plot PREFIRE clear-sky radiance on primary y-axis
            chan_data = data.sel({"channel": prefire_chan})

            ax.plot(
                chan_data.bin_center,
                chan_data,
                label=label,
                linestyle="solid",
                color=color,
                linewidth=2
            )
        ax.set_xlim(*xlims)
        ax.set_xlabel("PWV (mm)", fontsize=fontsize)
        ax.set_ylabel("Brightness Temperature (K)", fontsize=fontsize+2)
        ax.set_xscale('log', base=2)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(i) for i in xticks])
        ax.grid(alpha=0.7, axis="both")

    axes[0].set_ylim(215, 295)
    axes[1].set_ylim(215, 295)
    axes[0].legend()


    # Add panel labels
    letters = ["c", "d", "e"]
    for ax, letter, chan in zip(axes, letters, panel_labels):
        ax.text(
            0.0,
            1.06,
            f"({letter}) {chan}",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )

    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig,
        "fig1_lowerpanels_BTalt",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%