"""
Create plots showing 2D histograms conditioned on CWV and TS.

"""
# %%
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob
from pathlib import Path
import seaborn as sns
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
    full_path = '%s%s.%s' % (output_dir, filename, ext)

    if not os.path.exists(output_dir + filename):
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
        Mean surface temperature per CWV bin
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
        X-axis limits in original CWV units (default: (0, 70))
    xticks : array-like, optional
        X-tick positions in original CWV units
        (default: [0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
    threshold : float, optional
        Breakpoint for piecewise linear CWV transformation (default: 5)
    low_scale : float, optional
        Expansion factor for low CWV values (default: 2.0)
    high_scale : float, optional
        Compression factor for high CWV values (default: 0.4)
    
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

    # Transform CWV bin bounds to use piecewise linear scaling
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
        hist_data = data_dict[var][subvar].sel({sel_dim:sel})
        
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
        
        # Set x-ticks with original CWV labels
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
    load_dir = DATA_ROOT / "PREFIRE_conditional_correlations/"
    cesm_case = "20251121_081336.FHIST.f09_f09_mg17.PREFIREPRIME"
    fig_save_dir = OUTPUT_ROOT

    # Load histogram counts for pdf
    counts_2d_file = glob.glob(str(load_dir / f"{cesm_case}*.TMQ.TS.bincounts.zarr"))
    assert len(counts_2d_file) == 1, f"Expected exactly one counts file, found {len(counts_2d_file)}"
    counts_2d_file = counts_2d_file[0]
    counts_2d_ds = xr.open_zarr(counts_2d_file).compute()

    # Determine the mean CWV and TS values for each TS and CWV bin, respectively.
    mean_ts_global = counts_2d_ds["bin_center_TS"].weighted(counts_2d_ds["global_counts"]).mean("bin_TS")
    mean_tmq_global = counts_2d_ds["bin_center_TMQ"].weighted(counts_2d_ds["global_counts"]).mean("bin_TMQ")
    mean_ts_polar = counts_2d_ds["bin_center_TS"].weighted(counts_2d_ds["polar_counts"]).mean("bin_TS")
    mean_tmq_polar = counts_2d_ds["bin_center_TMQ"].weighted(counts_2d_ds["polar_counts"]).mean("bin_TMQ")

    # Modify dimensions so it will broadcast correctly with the other datasets
    counts_2d_ds = counts_2d_ds.drop_vars(["bin_center_TMQ", "bin_center_TS"]).rename({"bin_TMQ": "bin_center_TMQ", "bin_TS": "bin_center_TS"})
    counts_2d_ds["bin_center_TMQ"] = counts_2d_ds["bin_center_TMQ"] - 1
    counts_2d_ds["bin_center_TS"] = counts_2d_ds["bin_center_TS"] - 1

    # Load binned 2D histogram data for specified variables, if files exist
    load_vars = ["rttov_rad_total_inst001", "rttov_rad_clear_inst001", "LU_TOA", "LUC_TOA", "CLTMODIS"]
    data_dict = {}
    for var in load_vars:
        file = glob.glob(str(Path(load_dir) / f"{cesm_case}*.{var}.*binned2D.zarr"))
        if file:
            data_dict[var] = xr.open_zarr(*file).compute()
        else:
            print(f"Warning: No file found for variable {var} in {load_dir}. Skipping.")

    # %%
    # Convert radiances to brightness temperatures for the radiance variables.
    for var in ["rttov_rad_total_inst001", "rttov_rad_clear_inst001"]:
        if var in data_dict:
            wavenumbers = get_channel_wavenumbers(data_dict[var], instrument="PREFIRE")
            bt_var = var.replace("rad", "bt")
            bt_data = rad2bt(wavenumbers.values, data_dict[var])
            data_dict[bt_var] = bt_data

    # %%
    # Add the cloud radiative effect (CRE) variables by taking the difference between total and clear-sky brightness temperatures.
    data_dict["rttov_rad_CRE_inst001"] = data_dict["rttov_rad_clear_inst001"] - data_dict["rttov_rad_total_inst001"]
    data_dict["rttov_bt_CRE_inst001"] = data_dict["rttov_bt_clear_inst001"] - data_dict["rttov_bt_total_inst001"]
    data_dict["LU_CRE"] = data_dict["LUC_TOA"] - data_dict["LU_TOA"]

    # %%
    # Reconstruct the bin centers and bounds for CWV and TS.
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

    # %%
    # Combine clear-sky, all-sky, and CRE plots
    fig, axes = plt.subplots(3, 3, figsize=(12, 14))
    plt.subplots_adjust(hspace=0.45)
    cax1 = fig.add_axes([0.24, 0.635, 0.55, 0.01])
    cax2 = fig.add_axes([0.24, 0.35, 0.55, 0.01])
    cax3 = fig.add_axes([0.24, 0.06, 0.55, 0.01])

    labels = ["Mid-Infrared WV", "Atmospheric Window", "Far-Infrared"]
    row_labels = ["Clear-sky", "All-sky", "Cloud Radiative Effect"]
    fontsize = 14

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

    axs, cbar = plot_2d_histogram_with_contours(
        data_dict=data_dict,
        counts_2d_ds=counts_2d_ds,
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        mean_ts_for_cwv_bins=mean_ts_global,
        axs=axes[0, :],
        cax=cax1,
        var="rttov_bt_clear_inst001",
        subvar="global_mean",
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

    cax1.set_xticks(contour_levels_bt[::4])
    cax1.set_xlabel("Brightness Temperature (K)", fontsize=fontsize-2)

    contour_levels_bt=np.arange(210, 310, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)

    axs, cbar = plot_2d_histogram_with_contours(
        data_dict=data_dict,
        counts_2d_ds=counts_2d_ds,
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        mean_ts_for_cwv_bins=mean_ts_global,
        axs=axes[1, :],
        cax=cax2,
        var="rttov_bt_total_inst001",
        subvar="global_mean",
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

    cax2.set_xticks(contour_levels_bt[::4])
    cax2.set_xlabel("Brightness Temperature (K)", fontsize=fontsize-2)

    contour_levels_bt=np.arange(-5, 55, 5)
    contour_levels_bt_fake=np.arange(-50, 55, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap1 = sns.color_palette("coolwarm", n_colors=len(contour_levels_bt_fake) - 1)[9:]
    cmap = ListedColormap(cmap1)

    axs, cbar = plot_2d_histogram_with_contours(
        data_dict=data_dict,
        counts_2d_ds=counts_2d_ds,
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        mean_ts_for_cwv_bins=mean_ts_global,
        axs=axes[2, :],
        cax=cax3,
        var="rttov_bt_CRE_inst001",
        subvar="global_mean",
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

    cax3.set_xticks(contour_levels_bt[::2])
    cax3.set_xlabel("Cloud Radiative Effect (K)", fontsize=fontsize-2)

    # Remove y-labels from subplots after the first
    for ax in axes[:, 1:].flatten():
        ax.set_ylabel("")

    for ax, label in zip(axes[0, :], labels):
        ax.set_title(label, fontsize=fontsize, fontweight="bold", pad=20)

    # Standardize font sizes for all axes
    for ax in [*axes.flat] + [cax1, cax2]:
        ax.tick_params(axis='both', which='major', labelsize=fontsize-2)
        ax.xaxis.label.set_fontsize(fontsize-2)
        ax.yaxis.label.set_fontsize(fontsize-2)

    for ax, label in zip(axes[:, 0], row_labels):
        ax.text(
            0.0,
            1.06,
            label,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )

    for ax, label in zip(axes.flat, ["a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i."]):
        ax.text(
            0.97, 0.03, label,
            ha="right",
            va="bottom",
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )

    # %%
    to_png(
        fig,
        "fig_2D_histograms_PREFIRE_CESM",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
