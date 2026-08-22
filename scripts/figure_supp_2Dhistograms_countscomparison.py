"""
Create plots showing 2D histograms conditioned on PWV and TS.

"""
# %%
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob
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
    obs_load_dir = DATA_ROOT
    cesm_load_dir = DATA_ROOT / "PREFIRE_conditional_correlations"
    cesm_case = "20251121_081336.FHIST.f09_f09_mg17.PREFIREPRIME"
    fig_save_dir = OUTPUT_ROOT

    obs_ds = xr.open_dataset(str(obs_load_dir / "PREFIRE_BT_grids.nc"))
    obs_counts_ds = obs_ds["FIR_count"].sum("time")

    # Load histogram counts for pdf
    cesm_counts_file = glob.glob(f"{cesm_load_dir}{cesm_case}*.TMQ.TS.bincounts.zarr")
    assert len(cesm_counts_file) == 1, f"Expected exactly one counts file, found {len(cesm_counts_file)}"
    cesm_counts_file = cesm_counts_file[0]
    cesm_counts_ds = xr.open_zarr(cesm_counts_file)["global_counts"].compute()

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

    # %%
    # Compute mean TS for each PWV bin to plot as a reference line
    # Reshape the counts and mean datasets to have dimensions bin_center_TS, bin_center_PWV
    obs_counts_ds = obs_counts_ds.rename({"y_idx": "bin_center_TS", "x_idx": "bin_center_TMQ"})
    obs_counts_ds.coords["bin_center_TS"] = ts_centers
    obs_counts_ds.coords["bin_center_TMQ"] = cwv_centers

    # Determine the mean PWV and TS values for each TS and PWV bin, respectively.
    prefire_mean_ts_global = obs_counts_ds["bin_center_TS"].weighted(obs_counts_ds).mean("bin_center_TS")
    prefire_mean_tmq_global = obs_counts_ds["bin_center_TMQ"].weighted(obs_counts_ds).mean("bin_center_TMQ")

    # Modify dimensions so it will broadcast correctly with the other datasets
    cesm_counts_ds = cesm_counts_ds.drop_vars(["bin_center_TMQ", "bin_center_TS"]).rename({"bin_TMQ": "bin_center_TMQ", "bin_TS": "bin_center_TS"})
    cesm_counts_ds.coords["bin_center_TS"] = ts_centers
    cesm_counts_ds.coords["bin_center_TMQ"] = cwv_centers

    obs_counts_norm = obs_counts_ds / obs_counts_ds.sum()
    cesm_counts_norm = cesm_counts_ds / cesm_counts_ds.sum()
    norm_counts_diff = cesm_counts_norm - obs_counts_norm

    # %%
    # Plot the counts in each histogram bin
    # Transform PWV bin bounds to use piecewise linear scaling
    threshold = 5
    low_scale = 2.0
    high_scale = 0.4
    cwv_bin_bounds_transformed = transform_cwv_axis(cwv_bin_bounds, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
    upper_limit_cwv = 70
    upper_limit_transformed = transform_cwv_axis(upper_limit_cwv, threshold=threshold, low_scale=low_scale, high_scale=high_scale)

    # Only plot data with counts > threshold
    counts_threshold = 300
    fontsize = 12

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    for data, label, cbarlabel, ax in zip(
        [obs_counts_norm, cesm_counts_norm, norm_counts_diff],
        ["PREFIRE Observations", "CESM2", "CESM2 minus PREFIRE"],
        ["Frequency", "Frequency", "Frequency Difference"],
        axs,
    ):

        # data = obs_counts_ds
        counts_mask = obs_counts_ds > counts_threshold

        # Use a logarithmic color scale for the counts to better visualize the range of values
        if data.min() < 0:
            max_abs = np.max(np.abs(data))
            norm_counts = cm.colors.SymLogNorm(linthresh=0.0001, vmin=-max_abs, vmax=max_abs)
            cmap = plt.get_cmap("seismic")
            counts_mask = None
        else:
            norm_counts = cm.colors.LogNorm(vmin=1e-6, vmax=data.max().compute())
            cmap = plt.get_cmap("Reds")

        cax = ax.inset_axes([0.1, -0.2, 0.8, 0.04])

        # First plot the PREFIRE observations
        im = ax.pcolormesh(
            cwv_bin_bounds_transformed, ts_bin_bounds, data.values.T,
            cmap=cmap, norm=norm_counts, shading="auto",
        )
        # Add hatching for matched
        cwv_centers_transformed = transform_cwv_axis(cwv_centers, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        cwv_mesh, ts_mesh = np.meshgrid(cwv_centers_transformed, ts_centers)
        if counts_mask is not None:
            ax.contourf(
                cwv_mesh, ts_mesh, (~counts_mask).T.astype(float),
                levels=[0.5, 1.5], colors='none', hatches=['///'],
            )
        plt.colorbar(im, cax=cax, label=cbarlabel, orientation="horizontal")

        # Create custom x-ticks at important PWV values
        original_ticks = np.array([0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
        transformed_ticks = transform_cwv_axis(original_ticks, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        valid_ticks = transformed_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]
        valid_labels = original_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]

        ax.set_xticks(valid_ticks)
        ax.set_xticklabels([f"{int(val)}" for val in valid_labels], rotation=0)

        ax.set_xlabel("PWV (mm)")
        ax.set_ylabel("Surface Temperature (K)")
        ax.set_ylim(210, 305)
        ax.set_xlim(cwv_bin_bounds_transformed[0], upper_limit_transformed)
        ax.set_title(label, fontsize=fontsize+3, fontweight="bold")

    for ax, label in zip(axs.flat, ["a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i."]):
        ax.text(
            0.97, 0.03, label,
            ha="right",
            va="bottom",
            transform=ax.transAxes,
            fontsize=fontsize+2,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )
    # %%
    to_png(
        fig,
        "fig_2Dhistograms_countscomparison",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%