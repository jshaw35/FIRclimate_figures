"""
Use clear-sky PREFIRE observations (limited to the Arctic) for additional CESM2 validation.
- Limited spatial (polar) and temporal (requires clear-skies) sampling.
- Might mean that we can only get a single monthly average, in which case we'd just validate the monthly cycle for the Arctic and Antarctic
"""
# %%
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys
import os

# Ensure project-root modules are importable when this file is run directly.
try:
    PROJECT_ROOT = Path(__file__).resolve().parents[1]
except NameError:
    PROJECT_ROOT = Path.cwd().resolve()  # or hardcode/adjust as needed
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from data_config import DATA_ROOT, OUTPUT_ROOT

# %%
def plot_channel_comparison(
    obs_clearsky_subset_ds,
    cesm_clearsky_subset_ds,
    channel_idx: int,
    channel_label: str = "",
    vlims: list = [200, 290],
    prefire_channels: list = [6, 13, 28],
    prefire_wavelengths: list = [5.9, 11.8, 24.5],
    cesm_clearsky_var: str = "rttov_bt_clear_inst001",
    obs_clearsky_var: str = "brightness_temperature"
):
    """
    Create comparison figure for observations vs model for a single channel.
    
    Parameters
    ----------
    obs_clearsky_subset_ds : xr.Dataset
        PREFIRE observations dataset
    cesm_clearsky_subset_ds : xr.Dataset
        CESM2 model dataset
    prefire_channels : list
        List of channel numbers
    prefire_wavelengths : list
        List of wavelengths corresponding to channels
    cesm_clearsky_var : str
        Name of CESM2 brightness temperature variable
    vlims : list
        [min, max] for colorbar limits
    channel_idx : int
        Index of channel to plot (default: 1)
    """
    channel = prefire_channels[channel_idx]
    wavelength = prefire_wavelengths[channel_idx]
    
    fig = plt.figure(figsize=(16, 8))
    axes_flat = []
    for i in range(8):
        proj = ccrs.NorthPolarStereo() if i < 4 else ccrs.SouthPolarStereo()
        ax = fig.add_subplot(2, 4, i + 1, projection=proj)
        axes_flat.append(ax)
    
    months_data = [(7, "July"), (2, "February")]
    
    for month_idx, (month_num, month_name) in enumerate(months_data):
        obs_time = obs_clearsky_subset_ds.sel(
            time=(obs_clearsky_subset_ds["time.month"] == month_num)
        ).time.values[0]
        obs_data = obs_clearsky_subset_ds[obs_clearsky_var].sel(
            time=obs_time, channel=channel
        )
        model_data = cesm_clearsky_subset_ds[cesm_clearsky_var].sel(
            month=month_num, channel=channel
        )
        
        for pole_idx, (extent, pole_name) in enumerate(
            [([-180, 180, 60, 90], "North"), ([-180, 180, -90, -60], "South")]
        ):
            for data_idx, (data, label) in enumerate(
                [(obs_data, "PREFIRE"), (model_data, "CESM2")]
            ):
                ax = axes_flat[pole_idx * 4 + month_idx * 2 + data_idx]
                im = ax.pcolormesh(data.lon, data.lat, data.values,
                                   cmap="RdYlBu_r", vmin=vlims[0], vmax=vlims[1],
                                   transform=ccrs.PlateCarree())
                ax.coastlines()
                ax.set_extent(extent, crs=ccrs.PlateCarree())
                ax.set_title(f"{label}")
                plt.colorbar(im, ax=ax, label="BT (K)")

    for ax, label in zip(axes_flat, ["a.", "b.", "c.", "d.", "e.", "f.", "g.", "h.", "i."]):
        ax.text(
            0.97, 0.03, label,
            ha="right",
            va="bottom",
            transform=ax.transAxes,
            fontsize=14,
            fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="square,pad=0.2"),
        )

    fig.suptitle(f"{channel_label}: PREFIRE Channel {channel} - {wavelength} µm", fontsize=20, fontweight="bold", y=0.99)
    fig.text(0.25, 0.92, "July", fontsize=18, ha='right')
    fig.text(0.75, 0.92, "February", fontsize=18, ha='right')
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    plt.show()
    return fig


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
    load_dir = DATA_ROOT
    complete_clearsky_ds = xr.open_dataset(str(load_dir / "PREFIRE_SAT2_clearsky_202407_202505.nc"))
    prefire_channels = [6, 13, 28] # 5.9, 11.77, 24.5 um
    prefire_wavelengths = [5.9, 11.8, 24.5]

    # July and February have the best clear-sky coverage, so just plot those
    good_months = [7, 2]
    obs_clearsky_subset_ds = complete_clearsky_ds.sel(time=(complete_clearsky_ds["time.month"].isin(good_months)))
    obs_clearsky_subset_ds = obs_clearsky_subset_ds.sortby("time.month")
    vlims = [200, 290]

    # Model data
    data_savedir = DATA_ROOT
    cesm_clearsky_var = "rttov_bt_clear_inst001"
    model_datapath = data_savedir / f"monthlymeans_{cesm_clearsky_var}_PREFIREprime.zarr"

    model_clearsky_ds = xr.open_zarr(str(model_datapath))
    cesm_clearsky_subset_ds = model_clearsky_ds.sel(month=(model_clearsky_ds["month"].isin(good_months)))

    # %%
    fig_ch6 = plot_channel_comparison(
        obs_clearsky_subset_ds,
        cesm_clearsky_subset_ds,
        channel_idx=0,
        channel_label="Mid-Infrared WV",
        vlims=[200,290],
    )
    # %%
    fig_ch13 = plot_channel_comparison(
        obs_clearsky_subset_ds,
        cesm_clearsky_subset_ds,
        channel_idx=1,
        channel_label="Atmospheric Window",
        vlims=[200,290],
    )
    # %%
    fig_ch28 = plot_channel_comparison(
        obs_clearsky_subset_ds,
        cesm_clearsky_subset_ds,
        channel_idx=2,
        channel_label="Far-Infrared",
        vlims=[200,290],
    )
    # %%
    fig_save_dir = OUTPUT_ROOT
    to_png(
        fig_ch6,
        "fig_supp_ch6_clearskyvalidation",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    to_png(
        fig_ch13,
        "fig_supp_ch13_clearskyvalidation",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )
    to_png(
        fig_ch28,
        "fig_supp_ch28_clearskyvalidation",
        loc=fig_save_dir,
        dpi=200,
        ext="png",
        bbox_inches="tight",
    )

    # %%
