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
from pathlib import Path
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


def draw_line_with_endpoints_custom(ax, point1, point2, label=None, line_kwargs=None, 
                                    marker1_kwargs=None, marker2_kwargs=None):
    """
    Draw a line between two points with different icons at each endpoint.
    Legend is applied to the marker, not the line.
    
    Parameters
    ----------
    marker1_kwargs : dict
        Marker settings for point1
    marker2_kwargs : dict
        Marker settings for point2 (can differ from marker1_kwargs)
    """
    
    if line_kwargs is None:
        line_kwargs = {'color': 'black', 'linestyle': '--', 'linewidth': 2}
    if marker1_kwargs is None:
        marker1_kwargs = {'marker': 'o', 'markersize': 8, 'color': 'blue'}
    if marker2_kwargs is None:
        marker2_kwargs = {'marker': 's', 'markersize': 8, 'color': 'red'}
    
    # Draw line without label
    x_points = [point1[0], point2[0]]
    y_points = [point1[1], point2[1]]
    line = ax.plot(x_points, y_points, **line_kwargs)[0]
    
    # Plot markers with different styles - apply label to first marker
    marker1 = ax.plot(point1[0], point1[1], label=label, **marker1_kwargs)[0]
    marker2 = ax.plot(point2[0], point2[1], **marker2_kwargs)[0]
    
    return line, [marker1, marker2]


# %%
# Load preprocessed data for 2D histograms.
if __name__ == "__main__":

    load_dir = DATA_ROOT
    cesm_case = "20250616_103133.FHIST.f09_f09_mg17.cesm2.1.5_port_SSP585branch_PREFIRE"
    fig_save_dir = OUTPUT_ROOT

    # Load histogram counts for pdf
    counts_2d_file = f"{cesm_case}.h0.20150101-20891231.TMQ.TS.bincounts.zarr"
    counts_2d_ds = xr.open_zarr(str(load_dir / counts_2d_file))

    # Determine the mean PWV and TS values for each TS and PWV bin, respectively.
    mean_ts_global = counts_2d_ds["bin_center_TS"].weighted(counts_2d_ds["global_counts"]).mean("bin_TS")
    mean_tmq_global = counts_2d_ds["bin_center_TMQ"].weighted(counts_2d_ds["global_counts"]).mean("bin_TMQ")
    mean_ts_polar = counts_2d_ds["bin_center_TS"].weighted(counts_2d_ds["polar_counts"]).mean("bin_TS")
    mean_tmq_polar = counts_2d_ds["bin_center_TMQ"].weighted(counts_2d_ds["polar_counts"]).mean("bin_TMQ")

    # Modify dimensions so it will broadcast correctly with the other datasets
    counts_2d_ds = counts_2d_ds.drop_vars(["bin_center_TMQ", "bin_center_TS"]).rename({"bin_TMQ": "bin_center_TMQ", "bin_TS": "bin_center_TS"})
    counts_2d_ds["bin_center_TMQ"] = counts_2d_ds["bin_center_TMQ"] - 1
    counts_2d_ds["bin_center_TS"] = counts_2d_ds["bin_center_TS"] - 1

    # Load binned 2D histogram data for specified variables, if files exist
    load_vars = ["rttov_rad_total_inst001", "rttov_rad_clear_inst001", "LUC_TOA"]#, "LU_TOA", "LUC_TOA", "CLTMODIS"]
    data_dict = {}
    for var in load_vars:
        file = glob.glob(str(Path(load_dir) / f"{cesm_case}*.{var}.*binned2D.lon.zarr"))
        if file:
            data_dict[var] = xr.open_zarr(*file)
        else:
            print(f"Warning: No file found for variable {var} in {load_dir}. Skipping.")

    # %%
    # Convert radiances to brightness temperatures for the radiance variables.
    for var in ["rttov_rad_total_inst001", "rttov_rad_clear_inst001"]:
        if var in data_dict:
            wavenumbers = get_channel_wavenumbers(data_dict[var], instrument="PREFIRE")
            bt_var = var.replace("rad", "bt")
            bt_data = rad2bt(wavenumbers, data_dict[var].drop_vars(["all_climo_count", "early_climo_count", "late_climo_count"]))
            data_dict[bt_var] = bt_data

    # Compute Arctic and Antarctic specific values
    arctic_data_dict = {}
    antarctic_data_dict = {}
    for var in data_dict:
        _data = data_dict[var]
        bt_data_arctic = _data.sel(lat=slice(70, 90)).weighted(np.cos(np.deg2rad(_data["lat"]))).mean(dim="lat")
        bt_data_antarctic = _data.sel(lat=slice(-90, -70)).weighted(np.cos(np.deg2rad(_data["lat"]))).mean(dim="lat")
        arctic_data_dict[var] = bt_data_arctic
        antarctic_data_dict[var] = bt_data_antarctic

    # %%
    # Open TS and TMQ fields in order to compute the Arctic and Antarctic means
    case = "20250616_103133.FHIST.f09_f09_mg17.cesm2.1.5_port_SSP585branch_PREFIRE"
    ts_files = glob.glob(str(DATA_ROOT / f"{case}.h0.TS.?????????????.zarr"))
    tmq_files = glob.glob(str(DATA_ROOT / f"{case}.h0.TMQ*.nc"))
    ts_files.sort()
    tmq_files.sort()
    early_timesel = slice("2015", "2034")
    late_timesel = slice("2070", "2089")

    ts_ds = xr.concat([xr.open_zarr(i) for i in ts_files], dim="time")["TS"]
    tmq_ds = xr.open_mfdataset(tmq_files, chunks={"time": 120})["TMQ"]

    ts_ds["time"] = xr.CFTimeIndex(ts_ds["time"].values).shift(-1, "MS").shift(14, "1D")
    tmq_ds["time"] = xr.CFTimeIndex(tmq_ds["time"].values).shift(-1, "MS").shift(14, "1D")

    ts_arctic_mean = ts_ds.sel(lat=slice(70, 90)).weighted(np.cos(np.deg2rad(ts_ds["lat"]))).mean("lat").mean("lon")
    ts_antarctic_mean = ts_ds.sel(lat=slice(-90, -70)).weighted(np.cos(np.deg2rad(ts_ds["lat"]))).mean(dim=["lat", "lon"])

    tmq_arctic_mean = tmq_ds.sel(lat=slice(70, 90)).weighted(np.cos(np.deg2rad(tmq_ds["lat"]))).mean(dim=["lat", "lon"])
    tmq_antarctic_mean = tmq_ds.sel(lat=slice(-90, -70)).weighted(np.cos(np.deg2rad(tmq_ds["lat"]))).mean(dim=["lat", "lon"])

    ts_arctic_start = ts_arctic_mean.sel(time=early_timesel).groupby("time.season").mean("time")
    ts_arctic_end = ts_arctic_mean.sel(time=late_timesel).groupby("time.season").mean("time")
    tmq_arctic_start = tmq_arctic_mean.sel(time=early_timesel).groupby("time.season").mean("time")
    tmq_arctic_end = tmq_arctic_mean.sel(time=late_timesel).groupby("time.season").mean("time")

    print("Arctic start (TS, PWV):", ts_arctic_start.values, tmq_arctic_start.values)
    print("Arctic end (TS, PWV):", ts_arctic_end.values, tmq_arctic_end.values)

    ts_antarctic_start = ts_antarctic_mean.sel(time=early_timesel).groupby("time.season").mean("time")
    ts_antarctic_end = ts_antarctic_mean.sel(time=late_timesel).groupby("time.season").mean("time")
    tmq_antarctic_start = tmq_antarctic_mean.sel(time=early_timesel).groupby("time.season").mean("time")
    tmq_antarctic_end = tmq_antarctic_mean.sel(time=late_timesel).groupby("time.season").mean("time")

    print("Antarctic start (TS, PWV):", ts_antarctic_start.values, tmq_antarctic_start.values)
    print("Antarctic end (TS, PWV):", ts_antarctic_end.values, tmq_antarctic_end.values)

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


    # %%
    # Only plot data with counts > threshold (actually hash out the other points)
    counts_threshold = 300
    counts_max = 1e5
    weights = np.cos(np.deg2rad(data_dict["rttov_rad_total_inst001"]["lat"]))
    counts_arctic = data_dict["rttov_rad_total_inst001"][["all_climo_count", "early_climo_count", "late_climo_count"]].sel(lat=slice(70, 90)).weighted(weights).sum("lat")
    counts_antarctic = data_dict["rttov_rad_total_inst001"][["all_climo_count", "early_climo_count", "late_climo_count"]].sel(lat=slice(-90, -70)).weighted(weights).sum("lat")
    counts_mask_arctic = counts_arctic > counts_threshold
    counts_mask_antarctic = counts_antarctic > counts_threshold

    # Transform PWV bin bounds to use piecewise linear scaling
    threshold = 5
    low_scale = 2.0
    high_scale = 0.4
    cwv_bin_bounds_transformed = transform_cwv_axis(cwv_bin_bounds, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
    upper_limit_cwv = 70
    upper_limit_transformed = transform_cwv_axis(upper_limit_cwv, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
    
    # %%
    # Plot distribution change and 2D histograms for Arctic and Antarctic regions
    fig, axes = plt.subplots(2, 4, figsize=(16, 10))
    plt.subplots_adjust(hspace=0.28)

    caxa = fig.add_axes([0.13, 0.48, 0.16, 0.01])
    caxb = fig.add_axes([0.13, 0.04, 0.16, 0.01])
    cax1 = fig.add_axes([0.335, 0.04, 0.16, 0.01])
    cax2 = fig.add_axes([0.535, 0.04, 0.16, 0.01])
    cax3 = fig.add_axes([0.737, 0.04, 0.16, 0.01])

    labels = ["Distribution Change", "Mid-Infrared WV", "Atmospheric Window", "Far-Infrared", ]
    row_labels = ["Arctic", "Antarctic"]
    counts_climo_var = "all_climo_count"
    fontsize = 14

    counts_threshold=300
    ylims=(210, 295)
    xlims=(0, 35)
    xticks=np.array([0, 1, 2, 3, 4, 5, 10, 20, 30])
    threshold=5
    low_scale=2.0
    high_scale=0.4

    contour_levels_bt=np.arange(210, 310, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)

    plot_params = {
        "MWV": {
            "band": "NIR",
            "contour_levels": np.arange(2.0, 7.0+0.1, 0.5),
            "cmap": "RdYlBu_r",
        },
        "AW": {
            "band": "AW",
            "contour_levels": np.arange(5, 50+0.1, 5),
            "cmap": "RdYlBu_r",
        },
        "FIR": {
            "band": "FIR2",
            "contour_levels": np.arange(75, 115+0.1, 5),
            "cmap": "RdYlBu_r",
        },
    }
    
    for key, axs, cax in zip(plot_params, axes[:, 1:].T, [cax1, cax2, cax3]):
        band = plot_params[key]["band"]
        contour_levels = plot_params[key]["contour_levels"]
        norm = BoundaryNorm(contour_levels, ncolors=len(contour_levels) - 1)
        cmap = plt.get_cmap(plot_params[key]["cmap"], len(contour_levels) - 1)

        _axs, cbar = plot_2d_histogram_with_contours(
            data_dict=arctic_data_dict,
            counts_2d_ds=counts_arctic[counts_climo_var],
            cwv_bin_bounds=cwv_bin_bounds,
            ts_bin_bounds=ts_bin_bounds,
            axs=[axs[0]],
            cax=cax,
            var="LUC_TOA",
            subvar="all_climo_mean",
            sel_dict={"LW_band": [band]},
            fontsize=14,
            counts_threshold=counts_threshold,
            contour_levels_bt=contour_levels,
            ylims=ylims,
            xlims=xlims,
            xticks=xticks,
            threshold=threshold,
            low_scale=low_scale,
            high_scale=high_scale,
            norm=norm,
            cmap=cmap,
        )
        _axs, cbar = plot_2d_histogram_with_contours(
            data_dict=antarctic_data_dict,
            counts_2d_ds=counts_antarctic[counts_climo_var],
            cwv_bin_bounds=cwv_bin_bounds,
            ts_bin_bounds=ts_bin_bounds,
            axs=[axs[1]],
            cax=cax,
            var="LUC_TOA",
            subvar="all_climo_mean",
            sel_dict={"LW_band": [band]},
            fontsize=14,
            counts_threshold=counts_threshold,
            contour_levels_bt=contour_levels,
            ylims=ylims,
            xlims=xlims,
            xticks=xticks,
            threshold=threshold,
            low_scale=low_scale,
            high_scale=high_scale,
            norm=norm,
            cmap=cmap,
        )

        cax.set_xticks(contour_levels[::2])
        cax.set_xlabel("Spectral Flux (Wm$^{-2}$)")

    # Add points showing how the seasons change
    markers = ["o", "s", "v", "X"]  # Circle, square, triangle, cross
    start_color = 'white'
    end_color = 'black'

    tmq_arctic_start_scaled = xr.DataArray(transform_cwv_axis(tmq_arctic_start, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_arctic_start.coords, dims=tmq_arctic_start.dims)
    tmq_arctic_end_scaled = xr.DataArray(transform_cwv_axis(tmq_arctic_end, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_arctic_end.coords, dims=tmq_arctic_end.dims)
    tmq_antarctic_start_scaled = xr.DataArray(transform_cwv_axis(tmq_antarctic_start, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_antarctic_start.coords, dims=tmq_antarctic_start.dims)
    tmq_antarctic_end_scaled = xr.DataArray(transform_cwv_axis(tmq_antarctic_end, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_antarctic_end.coords, dims=tmq_antarctic_end.dims)
    for season, marker in zip(tmq_arctic_start_scaled["season"], markers):
        arctic_point_start = (tmq_arctic_start_scaled.sel(season=season).values, ts_arctic_start.sel(season=season).values)
        arctic_point_end = (tmq_arctic_end_scaled.sel(season=season).values, ts_arctic_end.sel(season=season).values)
        antarctic_point_start = (tmq_antarctic_start_scaled.sel(season=season).values, ts_antarctic_start.sel(season=season).values)
        antarctic_point_end = (tmq_antarctic_end_scaled.sel(season=season).values, ts_antarctic_end.sel(season=season).values)
        for (arctic_ax, antarctic_ax) in zip(axes[0, 1:], axes[1, 1:]):  # Add to Arctic row
            draw_line_with_endpoints_custom(
                arctic_ax, 
                arctic_point_start, 
                arctic_point_end,
                label=str(season.values),
                line_kwargs={'color': 'black', 'linestyle': '-', 'linewidth': 2},
                marker1_kwargs={'marker': marker, 'markersize': 7, 'color': start_color, 'markeredgecolor': 'black', 'markeredgewidth': 1},
                marker2_kwargs={'marker': marker, 'markersize': 7, 'color': end_color, 'markeredgecolor': 'black', 'markeredgewidth': 1}
            )
            draw_line_with_endpoints_custom(
                antarctic_ax, 
                antarctic_point_start, 
                antarctic_point_end,
                label="Seasonal change",
                line_kwargs={'color': 'black', 'linestyle': '-', 'linewidth': 2},
                marker1_kwargs={'marker': marker, 'markersize': 7, 'color': start_color, 'markeredgecolor': 'black', 'markeredgewidth': 1},
                marker2_kwargs={'marker': marker, 'markersize': 7, 'color': end_color, 'markeredgecolor': 'black', 'markeredgewidth': 1}
            )

    # Show the distribution change explicitly in the rightmost panel
    counts_cmap = plt.get_cmap("seismic")
    for ax, cax, data, mask in zip(axes[:, 0], [caxa, caxb], [counts_arctic, counts_antarctic], [counts_mask_arctic, counts_mask_antarctic]):

        data = data["late_climo_count"] - data["early_climo_count"]
        mask = mask["all_climo_count"]
        # Symmetric diverging norm with logarithmic scaling
        max_abs = np.max(np.abs(data.values))
        norm_counts = cm.colors.SymLogNorm(linthresh=100, vmin=-max_abs, vmax=max_abs)

        im = ax.pcolormesh(
            cwv_bin_bounds_transformed, ts_bin_bounds, data.values.T,
            cmap=counts_cmap, norm=norm_counts, shading="auto",
        )
        # Add hatching for matched
        # Convert bin bounds to bin centers for contour plot
        cwv_centers = (cwv_bin_bounds[:-1] + cwv_bin_bounds[1:]) / 2
        ts_centers = (ts_bin_bounds[:-1] + ts_bin_bounds[1:]) / 2
        cwv_centers_transformed = transform_cwv_axis(cwv_centers, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        cwv_mesh, ts_mesh = np.meshgrid(cwv_centers_transformed, ts_centers)
        ax.contourf(
            cwv_mesh, ts_mesh, (data==0).T.astype(float),
            levels=[0.5, 1.5], colors='grey',
        )

        cbar = plt.colorbar(im, cax=cax, label="Area-weighted Counts Change", orientation="horizontal")

        # Create custom x-ticks at important PWV values
        transformed_ticks = transform_cwv_axis(xticks, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        valid_ticks = transformed_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]
        valid_labels = xticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]

        ax.set_xticks(valid_ticks)
        ax.set_xticklabels([f"{int(val)}" for val in valid_labels], rotation=0)

        ax.set_xlabel("PWV (mm)")
        ax.set_ylabel("Surface Temperature (K)")
        ax.set_ylim(*ylims)
        ax.set_xlim(*transform_cwv_axis(xlims, threshold=threshold, low_scale=low_scale, high_scale=high_scale))
    caxa.set_xlabel("")

    # Remove y-labels from subplots after the first
    for ax in axes[:, 1:].flatten():
        ax.set_ylabel("")

    for ax, label in zip(axes[0, :], labels):
        ax.set_title(label, fontsize=fontsize, fontweight="bold", pad=15)

    for ax, label in zip(axes.flat, ["a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."]):
        ax.text(
            0.09, 0.98, label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )
    for ax, label in zip(axes[:, 0], row_labels):
        ax.text(
            -0.28,
            0.5,
            label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
        )

    axes[0, 1].legend(loc="lower right", fontsize=10, frameon=True, framealpha=1)

    # %%
    to_png(
        fig,
        "fig_2D_histograms_AAexplanation",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    
    # %%
    # Show histograms of the BT values for the Arctic and Antarctic, with contours and points showing how the seasons change.
    fig, axes = plt.subplots(2, 4, figsize=(16, 10))
    plt.subplots_adjust(hspace=0.40)

    cax1 = fig.add_axes([0.4, 0.49, 0.41, 0.01])
    cax2 = fig.add_axes([0.4, 0.04, 0.41, 0.01])

    cax1b = fig.add_axes([0.13, 0.49, 0.15, 0.01])
    cax2b = fig.add_axes([0.13, 0.04, 0.15, 0.01])

    labels = ["Distribution Change", "Mid-Infrared WV", "Atmospheric Window", "Far-Infrared"]
    row_labels = ["Arctic", "Antarctic"]
    counts_climo_var = "all_climo_count"
    fontsize = 14

    counts_threshold=300
    ylims=(210, 295)
    xlims=(0, 35)
    xticks=np.array([0, 1, 2, 3, 4, 5, 10, 20, 30])
    threshold=5
    low_scale=2.0
    high_scale=0.4

    contour_levels_bt=np.arange(210, 310, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)

    axs, cbar = plot_2d_histogram_with_contours(
        data_dict=arctic_data_dict,
        counts_2d_ds=counts_arctic[counts_climo_var],
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        axs=axes[0, 1:],
        cax=cax1,
        var="rttov_bt_clear_inst001",
        subvar="all_climo_mean",
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

    # Antarctic
    contour_levels_bt=np.arange(210, 310, 5)
    norm = BoundaryNorm(contour_levels_bt, ncolors=len(contour_levels_bt) - 1)
    cmap = plt.get_cmap("RdYlBu_r", len(contour_levels_bt) - 1)

    axs, cbar = plot_2d_histogram_with_contours(
        data_dict=antarctic_data_dict,
        counts_2d_ds=counts_antarctic[counts_climo_var],
        cwv_bin_bounds=cwv_bin_bounds,
        ts_bin_bounds=ts_bin_bounds,
        axs=axes[1, 1:],
        cax=cax2,
        var="rttov_bt_clear_inst001",
        subvar="all_climo_mean",
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

    for cax in [cax1, cax2]:
        cax.set_xticks(contour_levels_bt[::4])
        cax.set_xlabel("Brightness Temperature (K)")

    # Add points showing how the seasons change
    markers = ["o", "s", "v", "X"]  # Circle, square, triangle, cross
    start_color = 'white'
    end_color = 'black'

    tmq_arctic_start_scaled = xr.DataArray(transform_cwv_axis(tmq_arctic_start, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_arctic_start.coords, dims=tmq_arctic_start.dims)
    tmq_arctic_end_scaled = xr.DataArray(transform_cwv_axis(tmq_arctic_end, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_arctic_end.coords, dims=tmq_arctic_end.dims)
    tmq_antarctic_start_scaled = xr.DataArray(transform_cwv_axis(tmq_antarctic_start, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_antarctic_start.coords, dims=tmq_antarctic_start.dims)
    tmq_antarctic_end_scaled = xr.DataArray(transform_cwv_axis(tmq_antarctic_end, threshold=threshold, low_scale=low_scale, high_scale=high_scale), coords=tmq_antarctic_end.coords, dims=tmq_antarctic_end.dims)
    for season, marker in zip(tmq_arctic_start_scaled["season"], markers):
        arctic_point_start = (tmq_arctic_start_scaled.sel(season=season).values, ts_arctic_start.sel(season=season).values)
        arctic_point_end = (tmq_arctic_end_scaled.sel(season=season).values, ts_arctic_end.sel(season=season).values)
        antarctic_point_start = (tmq_antarctic_start_scaled.sel(season=season).values, ts_antarctic_start.sel(season=season).values)
        antarctic_point_end = (tmq_antarctic_end_scaled.sel(season=season).values, ts_antarctic_end.sel(season=season).values)
        for (arctic_ax, antarctic_ax) in zip(axes[0, 1:], axes[1, 1:]):  # Add to Arctic row
            draw_line_with_endpoints_custom(
                arctic_ax, 
                arctic_point_start, 
                arctic_point_end,
                label=str(season.values),
                line_kwargs={'color': 'black', 'linestyle': '-', 'linewidth': 2},
                marker1_kwargs={'marker': marker, 'markersize': 7, 'color': start_color, 'markeredgecolor': 'black', 'markeredgewidth': 1},
                marker2_kwargs={'marker': marker, 'markersize': 7, 'color': end_color, 'markeredgecolor': 'black', 'markeredgewidth': 1}
            )
            draw_line_with_endpoints_custom(
                antarctic_ax, 
                antarctic_point_start, 
                antarctic_point_end,
                label="Seasonal change",
                line_kwargs={'color': 'black', 'linestyle': '-', 'linewidth': 2},
                marker1_kwargs={'marker': marker, 'markersize': 7, 'color': start_color, 'markeredgecolor': 'black', 'markeredgewidth': 1},
                marker2_kwargs={'marker': marker, 'markersize': 7, 'color': end_color, 'markeredgecolor': 'black', 'markeredgewidth': 1}
            )

    # Show the distribution change explicitly in the leftmost panel
    counts_cmap = plt.get_cmap("seismic")
    for ax, cax, data, mask in zip(axes[:, 0], [cax1b, cax2b], [counts_arctic, counts_antarctic], [counts_mask_arctic, counts_mask_antarctic]):

        data = data["late_climo_count"] - data["early_climo_count"]
        mask = mask["all_climo_count"]
        # Symmetric diverging norm with logarithmic scaling
        max_abs = np.max(np.abs(data.values))
        norm_counts = cm.colors.SymLogNorm(linthresh=100, vmin=-max_abs, vmax=max_abs)

        im = ax.pcolormesh(
            cwv_bin_bounds_transformed, ts_bin_bounds, data.values.T,
            cmap=counts_cmap, norm=norm_counts, shading="auto",
        )
        # Add hatching for matched
        # Convert bin bounds to bin centers for contour plot
        cwv_centers = (cwv_bin_bounds[:-1] + cwv_bin_bounds[1:]) / 2
        ts_centers = (ts_bin_bounds[:-1] + ts_bin_bounds[1:]) / 2
        cwv_centers_transformed = transform_cwv_axis(cwv_centers, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        cwv_mesh, ts_mesh = np.meshgrid(cwv_centers_transformed, ts_centers)
        ax.contourf(
            cwv_mesh, ts_mesh, (data==0).T.astype(float),
            levels=[0.5, 1.5], colors='grey',
        )

        cbar = plt.colorbar(im, cax=cax, label="Area-weighted Counts Change", orientation="horizontal")

        # Create custom x-ticks at important PWV values
        transformed_ticks = transform_cwv_axis(xticks, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        valid_ticks = transformed_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]
        valid_labels = xticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]

        ax.set_xticks(valid_ticks)
        ax.set_xticklabels([f"{int(val)}" for val in valid_labels], rotation=0)

        ax.set_xlabel("PWV (mm)")
        ax.set_ylabel("Surface Temperature (K)")
        ax.set_ylim(*ylims)
        ax.set_xlim(*transform_cwv_axis(xlims, threshold=threshold, low_scale=low_scale, high_scale=high_scale))

    # Remove y-labels from subplots after the first
    for ax in axes[:, 1:].flatten():
        ax.set_ylabel("")

    for ax, label in zip(axes[0, :], labels):
        ax.set_title(label, fontsize=fontsize, fontweight="bold", pad=15)

    for ax, label in zip(axes[:, 0], row_labels):
        ax.text(
            -0.28,
            0.5,
            label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
        )

    for ax, label in zip(axes.flat, ["a.", "b.", "c.", "d.", "e.", "f.", "g.", "h."]):
        ax.text(
            0.09, 0.98, label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )

    axes[0, 1].legend(loc="lower right", fontsize=10, frameon=True, framealpha=1)

    # %%

    to_png(
        fig,
        "fig_2D_histograms_AAexplanation_BT",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
    # Plot the counts in each histogram bin
    # Transform PWV bin bounds to use piecewise linear scaling
    threshold = 5
    low_scale = 2.0
    high_scale = 0.4
    cwv_bin_bounds_transformed = transform_cwv_axis(cwv_bin_bounds, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
    upper_limit_cwv = 70
    upper_limit_transformed = transform_cwv_axis(upper_limit_cwv, threshold=threshold, low_scale=low_scale, high_scale=high_scale)

    # Use a logarithmic color scale for the counts to better visualize the range of values
    cmap = plt.get_cmap("Reds")
    norm_counts = cm.colors.LogNorm(vmin=100, vmax=counts_max)
    panel_labels = ["a.", "b."]
    panel_titles = ["Arctic", "Antarctic"]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    cax = fig.add_axes([0.25, -0.05, 0.5, 0.04])

    counts_climo_var = "all_climo_count"
    for ax, data, mask in zip(axs, [counts_arctic, counts_antarctic], [counts_mask_arctic, counts_mask_antarctic]):
        data = data[counts_climo_var]
        mask = mask[counts_climo_var]
        im = ax.pcolormesh(
            cwv_bin_bounds_transformed, ts_bin_bounds, data.values.T,
            cmap=cmap, norm=norm_counts, shading="auto",
        )
        # Add hatching for matched
        # Convert bin bounds to bin centers for contour plot
        cwv_centers = (cwv_bin_bounds[:-1] + cwv_bin_bounds[1:]) / 2
        ts_centers = (ts_bin_bounds[:-1] + ts_bin_bounds[1:]) / 2
        cwv_centers_transformed = transform_cwv_axis(cwv_centers, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
        cwv_mesh, ts_mesh = np.meshgrid(cwv_centers_transformed, ts_centers)
        if np.any(~mask):
            ax.contourf(
                cwv_mesh, ts_mesh, (~mask).T.astype(float),
                levels=[0.5, 1.5], colors='none', hatches=['///'],
            )
    plt.colorbar(im, cax=cax, label="Area-weighted Counts", orientation="horizontal")

    # Create custom x-ticks at important PWV values
    original_ticks = np.array([0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70])
    transformed_ticks = transform_cwv_axis(original_ticks, threshold=threshold, low_scale=low_scale, high_scale=high_scale)
    valid_ticks = transformed_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]
    valid_labels = original_ticks[transformed_ticks <= cwv_bin_bounds_transformed[-1]]

    for ax, title, label in zip(axs, panel_titles, panel_labels):
        ax.set_xticks(valid_ticks)
        ax.set_xticklabels([f"{int(val)}" for val in valid_labels], rotation=0)

        ax.set_xlabel("PWV (mm)")
        ax.set_ylabel("Surface Temperature (K)")
        ax.set_ylim(210, 305)
        ax.set_xlim(cwv_bin_bounds_transformed[0], upper_limit_transformed)
        ax.set_title(title, fontsize=14)
        ax.text(
            0.06, 0.98, label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )

    # %%
    to_png(
        fig,
        "fig_2D_histograms_weightedcounts_CESM_polar",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    # %%