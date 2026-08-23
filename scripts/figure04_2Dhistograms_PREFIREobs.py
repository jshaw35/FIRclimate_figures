"""
Create plots showing 2D histograms conditioned on PWV and TS.

"""
# %%
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import xarray as xr
import glob
from pathlib import Path
import seaborn as sns
import os
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


def transform_cwv_axis(cwv_values, threshold=5, low_scale=2.0, high_scale=0.4):
    """
    Apply piecewise linear transformation to PWV values for better visualization.
    
    Parameters
    ----------
    cwv_values : array-like
        PWV values in mm
    threshold : float
        Breakpoint between low and high scale regions (default: 5 mm)
    low_scale : float
        Expansion factor for values 0 to threshold (default: 2.0)
    high_scale : float
        Compression factor for values above threshold (default: 0.4)
    
    Returns
    -------
    transformed : array-like
        Transformed PWV values
    """
    cwv_array = np.asarray(cwv_values)
    transformed = np.zeros_like(cwv_array, dtype=float)
    
    # Low PWV region (0 to threshold): stretched
    low_mask = cwv_array <= threshold
    transformed[low_mask] = cwv_array[low_mask] * low_scale
    
    # High PWV region (threshold to 70): compressed
    high_mask = cwv_array > threshold
    max_cwv = np.max(cwv_array)
    high_range = max_cwv - threshold
    transformed[high_mask] = (threshold * low_scale) + (cwv_array[high_mask] - threshold) * high_scale
    
    return transformed


def plot_2d_histogram_with_contours(
    data_dict,
    counts_2d_ds,
    cwv_bin_bounds,
    ts_bin_bounds,
    axs,
    cax,
    var="rttov_bt_clear_inst001",
    subvar="global_mean",
    sel_dict=None,
    fontsize=14,
    counts_threshold=300,
    contour_levels_bt=None,
    mean_ts_for_cwv_bins=None,
    ylims=(210, 305),
    xlims=(0, 70),
    xticks=None,
    threshold=5,
    low_scale=2.0,
    high_scale=0.4,
    norm=None,
    cmap=None,
):
    """
    Create 2D histogram plots with contours and annotations.
    
    Parameters
    ----------
    data_dict : dict
        Dictionary containing the data arrays to plot
    counts_2d_ds : xr.Dataset
        Dataset containing bin counts and metadata
    cwv_bin_bounds : array-like
        Bin boundaries for column water vapor
    ts_bin_bounds : array-like
        Bin boundaries for surface temperature
    mean_ts_global : array-like
        Mean surface temperature per PWV bin
    axs : np.ndarray or matplotlib.axes.Axes
        Array of axes objects to plot on
    cax: matplotlib.axes.Axes
        Axes object for the shared colorbar
    var : str, optional
        Variable name to plot (default: "rttov_bt_clear_inst001")
    subvar : str, optional
        Sub-variable name within the data (default: "global_mean")
    sel_dict : dict, optional
        Dictionary with "channel" key containing list of channels to plot
        (default: {"channel": [6, 13, 27]})
    fontsize : int, optional
        Font size for titles and labels (default: 14)
    counts_threshold : int, optional
        Minimum count threshold for masking (default: 300)
    contour_levels_bt : array-like, optional
        Contour levels for brightness temperature
        (default: np.arange(210, 310, 5))
    ylims : tuple, optional
        Y-axis limits (default: (210, 305))
    xlims : tuple, optional
        X-axis limits in original PWV units (default: (0, 70))
    xticks : array-like, optional
        X-tick positions in original PWV units
        (default: [0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
    threshold : float, optional
        Breakpoint for piecewise linear PWV transformation (default: 5)
    low_scale : float, optional
        Expansion factor for low PWV values (default: 2.0)
    high_scale : float, optional
        Compression factor for high PWV values (default: 0.4)
    
    Returns
    -------
    axs : np.ndarray
        Array of axes objects
    """
    if sel_dict is None:
        sel_dict = {"channel": [6, 13, 27]}
    if contour_levels_bt is None:
        contour_levels_bt = np.arange(210, 310, 5)
    if xticks is None:
        xticks = np.array([0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
    
    # Create mask based on counts threshold
    if isinstance(counts_2d_ds, xr.Dataset):
        counts_mask = counts_2d_ds["global_counts"] > counts_threshold
    elif isinstance(counts_2d_ds, xr.DataArray):
        counts_mask = counts_2d_ds > counts_threshold
    else:
        raise ValueError("counts_2d_ds must contain 'global_counts' variable or be an xarray.DataArray with counts data")

    # Transform PWV bin bounds to use piecewise linear scaling
    cwv_bin_bounds_transformed = transform_cwv_axis(
        cwv_bin_bounds, threshold=threshold, low_scale=low_scale, high_scale=high_scale
    )
    xlims_transformed = transform_cwv_axis(
        xlims, threshold=threshold, low_scale=low_scale, high_scale=high_scale
    )

    # Create discrete colorbar for brightness temperature
    if norm is None:
        norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    if cmap is None:
        cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)
    
    # Plot each channel
    sel_dim = list(sel_dict.keys())[0]
    for i, (ax, sel) in enumerate(zip(axs, sel_dict[sel_dim])):
        # Select and mask data
        if subvar is not None:
            hist_data = data_dict[var][subvar].sel({sel_dim:sel})
        else:
            hist_data = data_dict[var].sel({sel_dim:sel})
        
        # Create pcolormesh
        im = ax.pcolormesh(
            cwv_bin_bounds_transformed,
            ts_bin_bounds,
            hist_data.values.T,
            norm=norm,
            cmap=cmap,
            shading="auto",
        )
        if isinstance(cax, np.ndarray):
            cbar = plt.colorbar(im, cax=cax[i], orientation="horizontal")
            
        # Set limits and labels
        ax.set_ylim(*ylims)
        ax.set_xlim(*xlims_transformed)
        ax.set_xlabel("PWV (mm)")
        ax.set_ylabel("Surface Temperature (K)")
        
        # Calculate bin centers for contours
        cwv_centers = (cwv_bin_bounds[:-1] + cwv_bin_bounds[1:]) / 2
        ts_centers = (ts_bin_bounds[:-1] + ts_bin_bounds[1:]) / 2
        cwv_centers_transformed = transform_cwv_axis(
            cwv_centers, threshold=threshold, low_scale=low_scale, high_scale=high_scale
        )
        
        # Create meshgrid
        cwv_mesh, ts_mesh = np.meshgrid(cwv_centers_transformed, ts_centers)
        
        # Add contours
        ax.contour(
            cwv_mesh,
            ts_mesh,
            hist_data.values.T,
            levels=contour_levels_bt,
            colors="black",
            linewidths=0.5,
            alpha=0.5,
            linestyles="-",
        )
        
        # Set x-ticks with original PWV labels
        transformed_ticks = transform_cwv_axis(
            xticks, threshold=threshold, low_scale=low_scale, high_scale=high_scale
        )
        valid_ticks = transformed_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]
        valid_labels = xticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]
        
        ax.set_xticks(valid_ticks)
        ax.set_xticklabels([f"{int(val)}" for val in valid_labels], rotation=0)
        
        # Add hatching for NaN regions
        if counts_mask is not None:
            ax.contourf(
                cwv_mesh,
                ts_mesh,
                (~counts_mask).T,
                levels=[0.5, 1.5],
                colors="none",
                hatches=["///"],
            )
        
        # Plot mean TS line
        if mean_ts_for_cwv_bins is not None:
            ax.plot(
                cwv_bin_bounds_transformed[:-1] + np.diff(cwv_bin_bounds_transformed) / 2,
                mean_ts_for_cwv_bins,
                color="black",
                linestyle="--",
                label="Mean Surface Temperature per PWV bin",
            )
    
    # Add shared colorbar
    if isinstance(cax, plt.Axes):
        cbar = plt.colorbar(
            im,
            cax=cax,
            orientation="horizontal",
            boundaries=contour_levels_bt,
        )

    return axs, cbar


# %%
# Load preprocessed data for 2D histograms.
if __name__ == "__main__":
    load_dir = DATA_ROOT
    fig_save_dir = OUTPUT_ROOT

    obs_ds = xr.open_dataset(str(load_dir / "PREFIRE_BT_grids.nc"))
    obs_annualmean_ds = obs_ds[["MWV_mean", "AW_mean", "FIR_mean"]].mean("time")
    obs_counts_ds = obs_ds["FIR_count"].sum("time")

    # Reshape to work with the plotting function
    obs_annualmean_channels = xr.concat(
        [obs_annualmean_ds["MWV_mean"], obs_annualmean_ds["AW_mean"], obs_annualmean_ds["FIR_mean"]],
        dim=xr.DataArray([6, 13, 27], dims="channel", name="channel"),
    )

    # %%
    # Reconstruct the bin centers and bounds for PWV and TS.
    low_bounds = np.linspace(0, 0.75, 4)
    mid_bounds = np.linspace(1, 29, 29)
    high_bounds = np.linspace(30, 70, 21)
    top_bound = np.array([120])
    cwv_bin_bounds = np.concat([low_bounds, mid_bounds, high_bounds, top_bound])
    cwv_condition_bins = np.array([cwv_bin_bounds[:-1], cwv_bin_bounds[1:]]).T

    bottom_bound = np.array([0])
    low_bounds = np.linspace(210, 285, 16)
    mid_bounds = np.linspace(290, 305, 7)
    high_bounds = np.array([315, 350])
    ts_bin_bounds = np.concat([bottom_bound, low_bounds, mid_bounds, high_bounds])
    ts_condition_bins = np.array([ts_bin_bounds[:-1], ts_bin_bounds[1:]]).T

    # Convert bin bounds to bin centers
    cwv_centers = (cwv_bin_bounds[:-1] + cwv_bin_bounds[1:]) / 2
    ts_centers = (ts_bin_bounds[:-1] + ts_bin_bounds[1:]) / 2

    # Compute mean TS for each PWV bin to plot as a reference line
    # Reshape the counts and mean datasets to have dimensions bin_center_TS, bin_center_PWV
    obs_counts_ds = obs_counts_ds.rename({"y_idx": "bin_center_TS", "x_idx": "bin_center_TMQ"})
    obs_counts_ds.coords["bin_center_TS"] = ts_centers
    obs_counts_ds.coords["bin_center_TMQ"] = cwv_centers
    obs_annualmean_channels = obs_annualmean_channels.rename({"x_idx": "bin_center_TMQ", "y_idx": "bin_center_TS"})
    obs_annualmean_channels.coords["bin_center_TS"] = ts_centers
    obs_annualmean_channels.coords["bin_center_TMQ"] = cwv_centers

    # Determine the mean PWV and TS values for each TS and PWV bin, respectively.
    prefire_mean_ts_global = obs_counts_ds["bin_center_TS"].weighted(obs_counts_ds).mean("bin_center_TS")
    prefire_mean_tmq_global = obs_counts_ds["bin_center_TMQ"].weighted(obs_counts_ds).mean("bin_center_TMQ")

    # %%
    load_dir = DATA_ROOT
    cesm_case = "20251121_081336.FHIST.f09_f09_mg17.PREFIREPRIME"
    fig_save_dir = OUTPUT_ROOT

    # Load histogram counts for pdf
    counts_2d_file = glob.glob(str(load_dir / f"{cesm_case}*.TMQ.TS.bincounts.zarr"))
    assert len(counts_2d_file) == 1, f"Expected exactly one counts file, found {len(counts_2d_file)}"
    counts_2d_file = counts_2d_file[0]
    counts_2d_ds = xr.open_zarr(counts_2d_file).compute()

    # # Modify dimensions so it will broadcast correctly with the other datasets
    counts_2d_ds = counts_2d_ds.drop_vars(["bin_center_TMQ", "bin_center_TS"]).rename({"bin_TMQ": "bin_center_TMQ", "bin_TS": "bin_center_TS"})
    counts_2d_ds["bin_center_TMQ"] = cwv_centers
    counts_2d_ds["bin_center_TS"] = ts_centers

    # Load binned 2D histogram data for specified variables, if files exist
    load_vars = ["rttov_rad_total_inst001"]
    cesm_data_dict = {}
    for var in load_vars:
        file = glob.glob(str(Path(load_dir) / f"{cesm_case}*.{var}.*binned2D.zarr"))
        if file:
            cesm_data_dict[var] = xr.open_zarr(*file).compute()
        else:
            print(f"Warning: No file found for variable {var} in {load_dir}. Skipping.")

    # %%
    # Convert radiances to brightness temperatures for the radiance variables.
    for var in load_vars:
        if var in cesm_data_dict:
            wavenumbers = get_channel_wavenumbers(cesm_data_dict[var], instrument="PREFIRE")
            bt_var = var.replace("rad", "bt")
            bt_data = rad2bt(wavenumbers.values, cesm_data_dict[var])
            cesm_data_dict[bt_var] = bt_data
    
    # Compute the mean BT for each PWV bin to plot as a reference line
    mean_ts_for_cwv_bins_dict = {}
    for var in cesm_data_dict:
        bt_data = cesm_data_dict[var]
        mean_ts_for_cwv_bins = bt_data.weighted(counts_2d_ds["global_counts"]).mean("bin_center_TS")
        mean_ts_for_cwv_bins_dict[f"{var}_mean_ts_for_cwv_bins"] = mean_ts_for_cwv_bins

    # %%
    # Create a 3-row figure comparing with the previous model results
    colors = sns.color_palette("colorblind")
    fontsize = 16
    prefire_rttov_channels = [6, 13, 27]
    channel_labels = ["Mid-Infrared WV", "Atmospheric Window", "Far-Infrared"]

    fig = plt.figure(figsize=(14, 12))
    
    # Create separate GridSpecs for different row groups with different spacing
    gs_top = GridSpec(2, 3, figure=fig, top=0.95, bottom=0.41, hspace=0.25, wspace=0.2)
    gs_bottom = GridSpec(1, 3, figure=fig, top=0.28, bottom=0.04, hspace=0.25, wspace=0.2)
    
    # Create axes array
    axes = np.empty((3, 3), dtype=object)
    
    # Fill top 2 rows
    for i in range(2):
        for j in range(3):
            axes[i, j] = fig.add_subplot(gs_top[i, j])
    
    # Fill bottom row
    for j in range(3):
        axes[2, j] = fig.add_subplot(gs_bottom[0, j])

    cax1 = fig.add_axes([0.24, 0.34, 0.55, 0.015])
    cax2 = fig.add_axes([0.24, -0.03, 0.55, 0.015])

    obs_allsky_bt_data = obs_annualmean_channels
    cesm_allsky_bt_data = cesm_data_dict["rttov_bt_total_inst001"]

    counts_threshold=300
    ylims=(210, 305)
    xlims=(0, 70)
    xticks=np.array([0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
    threshold=5
    low_scale=2.0
    high_scale=0.4

    for ax, prefire_chan, color in zip(axes[0,:], prefire_rttov_channels, colors):
        # Plot PREFIRE clear-sky radiance on primary y-axis
        obs_data = obs_allsky_bt_data.sel({"channel": prefire_chan})
        cesm_data = cesm_allsky_bt_data["global_mean"].sel({"channel": prefire_chan})

        # Row 1: PREFIRE vs CESM BT conditioned on PWV
        # Determine the mean BT for each PWV bin
        obs_marginal_mean_cwv = obs_data.weighted(obs_counts_ds).mean("bin_center_TS")
        cesm_marginal_mean_cwv = cesm_data.weighted(counts_2d_ds["global_counts"]).mean("bin_center_TS")
        cwv_bin_centers_transformed = transform_cwv_axis(obs_marginal_mean_cwv.bin_center_TMQ, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        cwv_ticks_transformed = transform_cwv_axis(xticks, threshold=threshold, low_scale=low_scale, high_scale=high_scale)

        ax.plot(
            cwv_bin_centers_transformed,
            obs_marginal_mean_cwv,
            label="PREFIRE",
            linestyle="solid",
            color=color,
            linewidth=2
        )
        ax.plot(
            cwv_bin_centers_transformed,
            cesm_marginal_mean_cwv,
            label="CESM2",
            linestyle="dashed",
            color=color,
            linewidth=2
        )
        ax.set_xticks(cwv_ticks_transformed)
        ax.set_xticklabels(xticks)
        ax.set_xlim(*transform_cwv_axis(xlims, threshold=threshold, low_scale=low_scale, high_scale=high_scale))
        ax.set_xlabel("PWV (mm)", fontsize=fontsize)
        ax.set_ylabel("Brightness Temperature (K)", fontsize=fontsize, color=color)
        ax.tick_params(axis='y', labelcolor=color)
        ax.grid(alpha=0.3, axis="both")
        ax.legend(fontsize=fontsize-2, loc=[0.25, 0.02])

    # Row 2: PREFIRE raw histograms
    obs_data = obs_allsky_bt_data

    # counts_threshold=300

    contour_levels_bt=np.arange(210, 310, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)

    _axs, cbar = plot_2d_histogram_with_contours(
        data_dict={"allsky": obs_data},
        counts_2d_ds=obs_counts_ds,
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        mean_ts_for_cwv_bins=prefire_mean_ts_global, # None
        axs=axes[1,:],
        cax=cax1,
        var="allsky",
        subvar=None,
        sel_dict={"channel": [6, 13, 27]},
        fontsize=fontsize,
        counts_threshold=counts_threshold,
        contour_levels_bt=contour_levels_bt,
        ylims=ylims,
        xlims=xlims,
        xticks=xticks,
        threshold=threshold,
        low_scale=low_scale,
        high_scale=high_scale,
        norm=norm,
        cmap=cmap,
    )
    cax1.set_xlabel("PREFIRE Brightness Temperature (K)", fontsize=fontsize-3)

    # Row 3: CESM2 histogram error
    obs_data = cesm_allsky_bt_data["global_mean"] - obs_allsky_bt_data

    contour_levels_bt=np.arange(-10, 10+0.1, 2)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("seismic", len(contour_levels_bt) - 1)

    _axs, cbar = plot_2d_histogram_with_contours(
        data_dict={"cesm_error": obs_data},
        counts_2d_ds=obs_counts_ds,
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        mean_ts_for_cwv_bins=None,
        axs=axes[2,:],
        cax=cax2,
        var="cesm_error",
        subvar=None,
        sel_dict={"channel": [6, 13, 27]},
        fontsize=fontsize,
        counts_threshold=counts_threshold,
        contour_levels_bt=contour_levels_bt,
        ylims=ylims,
        xlims=xlims,
        xticks=xticks,
        threshold=threshold,
        low_scale=low_scale,
        high_scale=high_scale,
        norm=norm,
        cmap=cmap,
    )
    cax2.set_xlabel("CESM2 minus PREFIRE Brightness Temperature (K)", fontsize=fontsize - 3)

    # Set tick label fontsize for all axes
    for ax in [*axes.flat] + [cax1, cax2]:
        ax.tick_params(axis='both', which='major', labelsize=fontsize-2)
        ax.xaxis.label.set_fontsize(fontsize-2)
        ax.yaxis.label.set_fontsize(fontsize-2)

    for ax, label in zip(axes[0, :], channel_labels):
        ax.set_title(label, fontsize=fontsize, fontweight="bold", pad=5)

    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")

    for ax, label in zip(axes.flat, ["a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i."]):
        ax.text(
            0.97, 0.03, label,
            ha="right",
            va="bottom",
            transform=ax.transAxes,
            fontsize=fontsize-2,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )

    # %%
    to_png(
        fig,
        "fig04_PREFIREobs_CESM_comparison_BT",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
    # Individual plots not used in manuscript
    # All-sky plot
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    cax = fig.add_axes([0.24, -0.07, 0.55, 0.04])

    labels = ["MWV", "AW", "FIR"]
    fontsize = 14

    data = obs_annualmean_channels
    data_dict = {"allsky": data}

    counts_threshold=300
    ylims=(210, 305)
    xlims=(0, 70)
    xticks=np.array([0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
    threshold=5
    low_scale=2.0
    high_scale=0.4

    contour_levels_bt=np.arange(210, 310, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)

    _axs, cbar = plot_2d_histogram_with_contours(
        data_dict=data_dict,
        counts_2d_ds=obs_counts_ds,
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        mean_ts_for_cwv_bins=prefire_mean_ts_global, # None
        axs=axs,
        cax=cax,
        var="allsky",
        subvar=None,
        sel_dict={"channel": [6, 13, 27]},
        fontsize=14,
        counts_threshold=counts_threshold,
        contour_levels_bt=contour_levels_bt,
        ylims=ylims,
        xlims=xlims,
        xticks=xticks,
        threshold=threshold,
        low_scale=low_scale,
        high_scale=high_scale,
        norm=norm,
        cmap=cmap,
    )

    cax.set_xticks(contour_levels_bt[::4])
    cax.set_xlabel("Brightness Temperature (K)")

    # Remove y-labels from subplots after the first
    for ax in axs[1:]:
        ax.set_ylabel("")
    for ax, label in zip(axs, labels):
        ax.set_title(label, fontsize=fontsize, fontweight='bold', pad=15)

    axs[0].text(
        0.0,
        1.06,
        "All-sky",
        transform=axs[0].transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )

    # %%
    to_png(
        fig,
        "fig_2D_histograms_PREFIREobs_BT",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    # %%