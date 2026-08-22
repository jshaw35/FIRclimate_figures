"""
Figure 4: Spectral band zonal trend and variability - 4-column comparison
Columns: Clear-sky, All-sky, CRE, SNR
Rows: MWV, AW, FIR (350-500), Broadband OLR
"""
# %%
import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from data_config import DATA_ROOT, OUTPUT_ROOT, get_data_file

# %%

def load_variable_multicase(data_dir, cases, variable_name):
    """Load a variable across multiple cases and time periods"""
    if isinstance(cases, str):
        cases = [cases]
    
    all_files = []
    for case in cases:
        pattern = f"{data_dir}/{case}*.{variable_name}.*zonalfields.zarr"
        files = sorted(glob.glob(pattern))
        all_files.extend(files)
    
    print(f"\nLoading {variable_name}:")
    print(f"  Found {len(all_files)} files")
    
    if len(all_files) == 0:
        raise FileNotFoundError(f"No files found for {variable_name}")
    
    datasets = []
    for f in all_files:
        ds = xr.open_zarr(f)
        datasets.append(ds)
    
    ds = xr.concat(datasets, dim='year')
    print(f"  Years: {ds.year.min().values} to {ds.year.max().values}")
    
    return ds


def calculate_trend_and_variability2(
    data,
    time_coord='year',
    durations=np.array([20]),
):
    """Calculate linear trend and detrended standard deviation for each latitude."""

    data_polyfit = data.polyfit(dim='year', deg=1)
    data_fit = xr.polyval(data.year, data_polyfit)['polyfit_coefficients']
    data_residual = data - data_fit
    data_trend = data_polyfit['polyfit_coefficients'].sel(degree=1)

    # Compute trend uncertainty follow Weatherhead et al.
    # Compute the stddev and lag-1 autocorrelation of the residuals
    residuals_stddev = data_residual.std(dim=time_coord)
    try:
        residuals_lag1r = xr.corr(data_residual, data_residual.shift({time_coord:1}), dim=time_coord)
    except:
        return data_residual

    # Add the net variability term
    tau = (1 + residuals_lag1r) / (1 - residuals_lag1r)
    residuals_netvar = residuals_stddev * np.sqrt(tau)

    durations_da = xr.DataArray(durations, dims='duration', coords={'duration': durations})
    data_trend_se = (12 * durations_da**(-3.0)) ** (1/2) * residuals_netvar

    return data_trend, data_trend_se


def apply_axis_styling(ax, ylim, vline_refs=None, hline_refs=None, xlabel_text='Latitude', fontsize=12):
    """
    Apply consistent styling to a plot axis.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to style
    ylim : tuple
        (ymin, ymax) limits for spectral flux/SNR values
    vline_refs : list of dict, optional
        Vertical line references. Each dict has keys: x, color, linestyle, linewidth, alpha
    hline_refs : list of dict, optional
        Horizontal line references. Each dict has keys: y, color, linestyle, linewidth, alpha
    xlabel_text : str
        X-axis label text (now latitude by default)
    """
    # Default reference lines
    if vline_refs is None:
        vline_refs = [{'x': 0, 'color': 'black', 'linestyle': '--', 'linewidth': 1, 'alpha': 0.8}]
    if hline_refs is None:
        hline_refs = [{'y': 0, 'color': 'black', 'linestyle': '-', 'linewidth': 1, 'alpha': 0.8}]
    
    # Draw reference lines
    for vref in vline_refs:
        ax.axvline(**vref)
    for href in hline_refs:
        ax.axhline(**href)
    
    # Set limits and ticks
    ax.set_xlim(-90, 90)
    ax.set_ylim(ylim)
    ax.set_xticks(np.arange(-90, 91, 30))
    
    # Grid and labels
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax.tick_params(axis='both', labelsize=fontsize-2)
    ax.set_xlabel(xlabel_text, fontsize=fontsize)

    # Set background color for better visibility
    ax.set_facecolor('whitesmoke') # "grey"


def create_figure_panels(fig, gs, panel_data, lats, fontsize=12):
    """
    Create and populate 4-column comparison figure panels.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to plot on
    gs : matplotlib.gridspec.GridSpec
        GridSpec for layout
    panel_data : list of tuples
        Data for each row: (label, trend_clear, var_clear, trend_all, var_all,
                           trend_cre, var_cre, snr_all, snr_clear, color)
    lats : array-like
        Latitude values
    total_snr_all, total_snr_clear, fir_snr_all, fir_snr_clear,
    aw_snr_all, aw_snr_clear, mwv_snr_all, mwv_snr_clear : array-like
        SNR arrays for calculating global limits
    """
    # Find global SNR limits for consistent scaling in column 4
    snr_data = []
    for i in panel_data:
        snr_data.extend(i[-4:-1])
    all_snr_values = np.concatenate(snr_data)

    snr_min = np.nanmin(all_snr_values)
    snr_max = np.nanmax(all_snr_values)
    snr_ylim = [snr_min * 1.1 if snr_min < 0 else -0.5, snr_max * 1.1]

    trend_styles = {'linestyle': '-', 'linewidth': 2.5, 'alpha': 0.9}
    var_styles = {'linestyle': '--', 'linewidth': 2, 'alpha': 0.7}
    varlabel_height = 1.05
    panellabel_fontsize = fontsize
    legend_fontsize = fontsize - 3.5
    panellabel_height = 0.9

    # Loop through rows
    for row, (label, trend_clear, var_clear, trend_all, var_all,
            trend_cre, var_cre, snr_all, snr_clear, snr_cre, color) in enumerate(panel_data):
        
        # ========== COLUMN 1: CLEAR-SKY ==========
        ax1 = fig.add_subplot(gs[row, 0])
        # Convert to per decade units for better interpretability
        trend_clear = trend_clear * 10
        var_clear = var_clear * 10
        trend_all = trend_all * 10
        var_all = var_all * 10
        trend_cre = trend_cre * 10
        var_cre = var_cre * 10

        # Determine limits
        all_vals = np.concatenate([trend_clear, var_clear])
        val_min = all_vals.min()
        val_max = all_vals.max()
        ylim_min = val_min * 1.2 if val_min < 0 else -0.15 * val_max
        ylim_max = val_max * 1.1

        # Plot and style (swapped: now plot(lats, values))
        ax1.plot(lats, trend_clear, color=color, label='Trend', **trend_styles)
        ax1.plot(lats, var_clear, color=color, label='Variability', **var_styles)
        apply_axis_styling(ax1, (ylim_min, ylim_max), fontsize=fontsize, xlabel_text='')
        
        # Column title for first row
        if row == 0:
            ax1.text(0.5, varlabel_height, 'Clear-sky', transform=ax1.transAxes,
                    fontsize=fontsize+2, fontweight='bold', ha='center', va='bottom')
            ax1.legend(loc=[0.59, 0.01], fontsize=legend_fontsize, frameon=True)

        # Band label
        pos = ax1.get_position()
        y_center = pos.y0 + 0.5 * pos.height
        fig.text(0.02, y_center, label, fontsize=fontsize, fontweight='bold', # JKS
        # fig.text(0.00, y_center, label, fontsize=fontsize, fontweight='bold', # JKS
            ha='left', va='center', rotation=90)

        # Panel label
        ax1.text(0.02, panellabel_height, f'({chr(97+row*4)})', transform=ax1.transAxes, 
                fontsize=panellabel_fontsize, fontweight='bold', va='bottom', ha='left')
        
        # ========== COLUMN 2: ALL-SKY ==========
        ax2 = fig.add_subplot(gs[row, 1])
        
        # Determine limits
        all_vals = np.concatenate([trend_all, var_all])
        val_min = all_vals.min()
        val_max = all_vals.max()
        ylim_min = val_min * 1.2 if val_min < 0 else -0.15 * val_max
        ylim_max = val_max * 1.1
        
        # Plot and style (swapped: now plot(lats, values))
        ax2.plot(lats, trend_all, color=color, label='Trend', **trend_styles)
        ax2.plot(lats, var_all, color=color, label='Variability', **var_styles)
        apply_axis_styling(ax2, (ylim_min, ylim_max), xlabel_text='', fontsize=fontsize)
        ax2.set_xlabel('')
        
        # Column title for first row
        if row == 0:
            ax2.text(0.5, varlabel_height, 'All-sky', transform=ax2.transAxes,
                    fontsize=fontsize+2, fontweight='bold', ha='center', va='bottom')
            ax2.legend(loc=[0.59, 0.01], fontsize=legend_fontsize, frameon=True)
        
        # Panel label
        ax2.text(0.02, panellabel_height, f'({chr(98+row*4)})', transform=ax2.transAxes, 
                fontsize=panellabel_fontsize, fontweight='bold', va='bottom', ha='left')
        
        # ========== COLUMN 3: CRE ==========
        ax3 = fig.add_subplot(gs[row, 2])
        
        # Determine limits
        all_vals = np.concatenate([trend_cre, var_cre])
        val_min = all_vals.min()
        val_max = all_vals.max()
        ylim_min = val_min * 1.2 if val_min < 0 else -0.15 * max(abs(val_min), val_max)
        ylim_max = val_max * 1.2 if val_max > 0 else 0.15 * max(abs(val_min), val_max)
        
        # Plot and style (swapped: now plot(lats, values))
        ax3.plot(lats, trend_cre, color=color, label='Trend', **trend_styles)
        ax3.plot(lats, var_cre, color=color, label='Variability', **var_styles)
        apply_axis_styling(ax3, (ylim_min, ylim_max), xlabel_text='', fontsize=fontsize)
        ax3.set_xlabel('')

        # Column title for first row
        if row == 0:
            ax3.text(0.5, varlabel_height, 'Cloud Radiative Effect', transform=ax3.transAxes,
                    fontsize=fontsize+2, fontweight='bold', ha='center', va='bottom')
            ax3.legend(loc=[0.59, 0.01], fontsize=legend_fontsize, frameon=True)
        
        # Panel label
        ax3.text(0.02, panellabel_height, f'({chr(99+row*4)})', transform=ax3.transAxes, 
                fontsize=panellabel_fontsize, fontweight='bold', va='bottom', ha='left')
        
        # ========== COLUMN 4: SNR (SHARED X-AXIS) ==========
        ax4 = fig.add_subplot(gs[row, 3])
        
        # Plot (swapped: now plot(lats, snr_values))
        ax4.plot(lats, snr_clear, color=color, label='Clear-sky', **trend_styles)
        ax4.plot(lats, snr_all, color=color, label='All-sky', **var_styles)
        ax4.plot(lats, snr_cre, color=color, label='CRE', linestyle=':', linewidth=2, alpha=0.9)
        
        # Special SNR styling with reference lines (now horizontal)
        hline_refs = [
            {'y': -2, 'color': 'black', 'linestyle': '-', 'linewidth': 1, 'alpha': 0.8},
            {'y': 0, 'color': 'black', 'linestyle': '--', 'linewidth': 1, 'alpha': 0.5},
            {'y': 2, 'color': 'black', 'linestyle': '-', 'linewidth': 1, 'alpha': 0.5},
        ]
        apply_axis_styling(ax4, snr_ylim, hline_refs=hline_refs, xlabel_text='', fontsize=fontsize)

        # Fill region where -2 < SNR < 2 (now fill_between instead of fill_betweenx)
        ax4.fill_between(lats, -2, 2, color='lightcoral', alpha=0.3)

        # Column title for first row
        if row == 0:
            ax4.text(0.5, varlabel_height, 'Signal-to-Noise Ratio', transform=ax4.transAxes,
                    fontsize=fontsize+2, fontweight='bold', ha='center', va='bottom')
            ax4.legend(loc=[0.61, 0.72], fontsize=legend_fontsize, frameon=True)

        # Panel label
        ax4.text(0.02, panellabel_height, f'({chr(100+row*4)})', transform=ax4.transAxes, 
                fontsize=panellabel_fontsize, fontweight='bold', va='bottom', ha='left')

        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xticks(np.arange(-90, 91, 30))
            ax.set_xlim(-90, 90)
            ax.set_xticklabels(["90S", "60S", "30S", "0", "30N", "60N", "90N"])
        # Add figure-wide x-axis and y-axis labels
        fig.text(0.41, 0.22, 'Latitude', ha='center', va='center', fontsize=fontsize+4)
        fig.text(0.885, 0.22, 'Latitude', ha='center', va='center', fontsize=fontsize+4)
        fig.text(0.41, 0.935, 'Spectral Flux Trend (W/m²/decade)', ha='center', va='center', fontsize=fontsize+4, rotation='horizontal')

        # Add vertical divider between columns 2 and 3 (CRE and SNR)
        ax_col2 = fig.get_axes()[2]  # Get an axis from column 2
        pos_col2 = ax_col2.get_position()
        divider_x = pos_col2.x1  # Right edge of column 2
        # Draw vertical line from bottom to top of plot area
        line = plt.Line2D([divider_x + 0.014, divider_x + 0.014], [0.215, 0.94], 
                        transform=fig.transFigure, 
                        color='black', linewidth=1.5, linestyle='-', alpha=0.7)
        fig.add_artist(line)


# %%

if __name__ == "__main__":

    data_dir = DATA_ROOT / "CESM2"
    output_dir = OUTPUT_ROOT

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    case_2015_2089 = "20250616_103133.FHIST.f09_f09_mg17.cesm2.1.5_port_SSP585branch_PREFIRE"

    #Load data
    print("="*70)
    print("LOADING DATA")
    print("="*70)

    ds_total_olr_allsky = load_variable_multicase(data_dir, [case_2015_2089], "FLUT")
    ds_spectral_allsky = load_variable_multicase(data_dir, [case_2015_2089], "LU_TOA")
    ds_total_olr_clearsky = load_variable_multicase(data_dir, [case_2015_2089], "FLUTC")
    ds_spectral_clearsky = load_variable_multicase(data_dir, [case_2015_2089], "LUC_TOA")

    #Extract and process data
    print("\n" + "="*70)
    print("PROCESSING DATA (2000-2089)")
    print("="*70)

    year_slice = slice(2000, 2089)
    mwv_slice = slice(1390, 2080)
    aw_slice = slice(820, 980)
    fir_narrow_slice = slice(350, 500)
    fir_full_slice = slice(10, 630)
    fir_med_slice = slice(10, 500) # "FIR3": Exclude edges of CO2 absorption

    total_allsky = ds_total_olr_allsky['mean'].sel(year=year_slice)
    fir_allsky = ds_spectral_allsky['mean'].sel(lw_band=fir_narrow_slice, year=year_slice).sum(dim='lw_band')
    fir_full_allsky = ds_spectral_allsky['mean'].sel(lw_band=fir_full_slice, year=year_slice).sum(dim='lw_band')
    fir_med_allsky = ds_spectral_allsky['mean'].sel(lw_band=fir_med_slice, year=year_slice).sum(dim='lw_band')
    aw_allsky = ds_spectral_allsky['mean'].sel(lw_band=aw_slice, year=year_slice).sum(dim='lw_band')
    mwv_allsky = ds_spectral_allsky['mean'].sel(lw_band=mwv_slice, year=year_slice).sum(dim='lw_band')

    total_clearsky = ds_total_olr_clearsky['mean'].sel(year=year_slice)
    fir_clearsky = ds_spectral_clearsky['mean'].sel(lw_band=fir_narrow_slice, year=year_slice).sum(dim='lw_band')
    fir_full_clearsky = ds_spectral_clearsky['mean'].sel(lw_band=fir_full_slice, year=year_slice).sum(dim='lw_band')
    fir_med_clearsky = ds_spectral_clearsky['mean'].sel(lw_band=fir_med_slice, year=year_slice).sum(dim='lw_band')
    aw_clearsky = ds_spectral_clearsky['mean'].sel(lw_band=aw_slice, year=year_slice).sum(dim='lw_band')
    mwv_clearsky = ds_spectral_clearsky['mean'].sel(lw_band=mwv_slice, year=year_slice).sum(dim='lw_band')

    # CRE (Cloud Radiative Effect)
    # Standard definition: CRE = Clear-sky - All-sky
    total_cre = total_clearsky - total_allsky
    fir_cre = fir_clearsky - fir_allsky
    fir_full_cre = fir_full_clearsky - fir_full_allsky
    fir_med_cre = fir_med_clearsky - fir_med_allsky
    aw_cre = aw_clearsky - aw_allsky
    mwv_cre = mwv_clearsky - mwv_allsky

    lats = total_allsky.lat.values

    # Calculate trends and variability
    print("\n" + "="*70)
    print("CALCULATING TRENDS AND VARIABILITY")
    print("="*70)

    # All-sky
    # total_trend_all2, total_var_all2 = calculate_trend_and_variability(total_allsky)
    total_trend_all, total_var_all = calculate_trend_and_variability2(total_allsky)
    fir_trend_all, fir_var_all = calculate_trend_and_variability2(fir_allsky)
    fir_full_trend_all, fir_full_var_all = calculate_trend_and_variability2(fir_full_allsky)
    fir_med_trend_all, fir_med_var_all = calculate_trend_and_variability2(fir_med_allsky)
    aw_trend_all, aw_var_all = calculate_trend_and_variability2(aw_allsky)
    mwv_trend_all, mwv_var_all = calculate_trend_and_variability2(mwv_allsky)

    # Clear-sky
    total_trend_clear, total_var_clear = calculate_trend_and_variability2(total_clearsky)
    fir_trend_clear, fir_var_clear = calculate_trend_and_variability2(fir_clearsky)
    fir_full_trend_clear, fir_full_var_clear = calculate_trend_and_variability2(fir_full_clearsky)
    fir_med_trend_clear, fir_med_var_clear = calculate_trend_and_variability2(fir_med_clearsky)
    aw_trend_clear, aw_var_clear = calculate_trend_and_variability2(aw_clearsky)
    mwv_trend_clear, mwv_var_clear = calculate_trend_and_variability2(mwv_clearsky)

    # Cloud Radiative Effect (CRE)
    total_trend_cre, total_var_cre = calculate_trend_and_variability2(total_cre)
    fir_trend_cre, fir_var_cre = calculate_trend_and_variability2(fir_cre)
    fir_full_trend_cre, fir_full_var_cre = calculate_trend_and_variability2(fir_full_cre)
    fir_med_trend_cre, fir_med_var_cre = calculate_trend_and_variability2(fir_med_cre)
    aw_trend_cre, aw_var_cre = calculate_trend_and_variability2(aw_cre)
    mwv_trend_cre, mwv_var_cre = calculate_trend_and_variability2(mwv_cre)

    # Squeeze the variability terms
    total_var_all = total_var_all.sel(duration=20).squeeze()
    fir_var_all = fir_var_all.sel(duration=20).squeeze()
    fir_full_var_all = fir_full_var_all.sel(duration=20).squeeze()
    fir_med_var_all = fir_med_var_all.sel(duration=20).squeeze()
    aw_var_all = aw_var_all.sel(duration=20).squeeze()
    mwv_var_all = mwv_var_all.sel(duration=20).squeeze()

    total_var_clear = total_var_clear.sel(duration=20).squeeze()
    fir_var_clear = fir_var_clear.sel(duration=20).squeeze()
    fir_full_var_clear = fir_full_var_clear.sel(duration=20).squeeze()
    fir_med_var_clear = fir_med_var_clear.sel(duration=20).squeeze()
    aw_var_clear = aw_var_clear.sel(duration=20).squeeze()
    mwv_var_clear = mwv_var_clear.sel(duration=20).squeeze()

    total_var_cre = total_var_cre.sel(duration=20).squeeze()
    fir_var_cre = fir_var_cre.sel(duration=20).squeeze()
    fir_full_var_cre = fir_full_var_cre.sel(duration=20).squeeze()
    fir_med_var_cre = fir_med_var_cre.sel(duration=20).squeeze()
    aw_var_cre = aw_var_cre.sel(duration=20).squeeze()
    mwv_var_cre = mwv_var_cre.sel(duration=20).squeeze()

    # SNR
    total_snr_all = total_trend_all / total_var_all
    total_snr_clear = total_trend_clear / total_var_clear
    total_snr_cre = total_trend_cre / total_var_cre
    fir_snr_all = fir_trend_all / fir_var_all
    fir_snr_clear = fir_trend_clear / fir_var_clear
    fir_snr_cre = fir_trend_cre / fir_var_cre
    fir_full_snr_all = fir_full_trend_all / fir_full_var_all
    fir_full_snr_clear = fir_full_trend_clear / fir_full_var_clear
    fir_full_snr_cre = fir_full_trend_cre / fir_full_var_cre
    fir_med_snr_all = fir_med_trend_all / fir_med_var_all
    fir_med_snr_clear = fir_med_trend_clear / fir_med_var_clear
    fir_med_snr_cre = fir_med_trend_cre / fir_med_var_cre
    aw_snr_all = aw_trend_all / aw_var_all
    aw_snr_clear = aw_trend_clear / aw_var_clear
    aw_snr_cre = aw_trend_cre / aw_var_cre
    mwv_snr_all = mwv_trend_all / mwv_var_all
    mwv_snr_clear = mwv_trend_clear / mwv_var_clear
    mwv_snr_cre = mwv_trend_cre / mwv_var_cre

    print("\nAll-sky results:")
    print(f"  Broadband OLR trend: {total_trend_all.mean().values:.3f} W/m²/decade")
    print(f"  FIR trend: {fir_trend_all.mean().values:.3f} W/m²/decade")
    print(f"  FIR (full, 10-630cm-1) trend: {fir_full_trend_all.mean().values:.3f} W/m²/decade")
    print(f"  FIR (med, 10-500cm-1) trend: {fir_med_trend_all.mean().values:.3f} W/m²/decade")
    print(f"  AW trend: {aw_trend_all.mean().values:.3f} W/m²/decade")
    print(f"  MWV trend: {mwv_trend_all.mean().values:.3f} W/m²/decade")

    print("\nCRE results:")
    print(f"  Broadband OLR CRE trend: {total_trend_cre.mean().values:.3f} W/m²/decade")
    print(f"  Broadband OLR CRE variability: {total_var_cre.mean().values:.3f} W/m²")

    # %%
    # Create figure
    print("\n" + "="*70)
    print("CREATING 4-COLUMN COMPARISON FIGURE")
    print("="*70)

    fig_width_cm = 40
    fig_height_cm = 45
    inches_per_cm = 1 / 2.54
    fig_width = fig_width_cm * inches_per_cm
    fig_height = fig_height_cm * inches_per_cm

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Grid: 4 rows × 4 columns
    widths = [1, 1, 1, 1]  # Clear-sky, All-sky, CRE, SNR
    heights = [1, 1, 1, 1, 1]   # MWV, AW, FIR full, FIR med, Broadband OLR
    gs = gridspec.GridSpec(len(heights), len(widths), 
                        left=0.07, right=0.98, top=0.90, bottom=0.08,
                        wspace=0.2, hspace=0.09,
                        width_ratios=widths, height_ratios=heights)

    # Color scheme
    cmap = sns.color_palette("colorblind")
    MWV_c = cmap[0]
    AW_c = cmap[1]
    FIR_c = cmap[2]
    FIR2_c = "#2BCFA1"
    BB_c = 'black'  # Broadband OLR

    # Order: MWV, AW, FIR, Broadband OLR
    panel_data = [
        ("MWV (1390-2080 cm⁻¹)",
        mwv_trend_clear, mwv_var_clear,
        mwv_trend_all, mwv_var_all,
        mwv_trend_cre, mwv_var_cre,
        mwv_snr_all, mwv_snr_clear, mwv_snr_cre,
        MWV_c),
        ("AW (820-980 cm⁻¹)",
        aw_trend_clear, aw_var_clear,
        aw_trend_all, aw_var_all,
        aw_trend_cre, aw_var_cre,
        aw_snr_all, aw_snr_clear, aw_snr_cre,
        AW_c),
        ("FIR (10-630 cm⁻¹)",
        fir_full_trend_clear, fir_full_var_clear,
        fir_full_trend_all, fir_full_var_all,
        fir_full_trend_cre, fir_full_var_cre,
        fir_full_snr_all, fir_full_snr_clear, fir_full_snr_cre,
        FIR2_c),
        ("FIR-H$_2$O (10-500 cm⁻¹)",
        fir_med_trend_clear, fir_med_var_clear,
        fir_med_trend_all, fir_med_var_all,
        fir_med_trend_cre, fir_med_var_cre,
        fir_med_snr_all, fir_med_snr_clear, fir_med_snr_cre,
        FIR_c),
    ]

    # Create figure panels
    create_figure_panels(fig, gs, panel_data, lats, fontsize=14)
    axes = fig.get_axes()
    # Align the ratios of the FIR rows so that the trends and variability are visually comparable
    for (fir1_ax, fir2_ax) in zip(axes[8:11], axes[12:15]):  # FIR full and FIR med rows
        ratio1 = fir1_ax.get_xlim()[0] / fir1_ax.get_xlim()[1]
        min2 = ratio1 * fir2_ax.get_xlim()[1]
        fir2_ax.set_xlim(min2, None)

    for ax in axes[:12]:
        ax.set_xticklabels([None for i in range(7)])
        ax.set_xlabel("")

    # %%
    # Save
    output_file_png = os.path.join(output_dir, 'fig4_CESMtrends.png')
    fig.savefig(output_file_png, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved as: {output_file_png}")

    plt.close()

    # %%
    # Compute GMST change for paper
    ds_ts = load_variable_multicase(data_dir, [case_2015_2089], "TS")
    ds_ts_globalmean = ds_ts['mean'].weighted(np.cos(np.deg2rad(ds_ts.lat))).mean(dim="lat")
    ts_change_decades = ds_ts_globalmean.isel(year=slice(-10, None)).mean("year") - ds_ts_globalmean.isel(year=slice(None,10)).mean("year")
    ts_change_years = ds_ts_globalmean.isel(year=-1) - ds_ts_globalmean.isel(year=0)

    print("\nGlobal mean surface temperature change (2089 minus 2015): ", ts_change_years.values, "K")
    print("\nGlobal mean surface temperature change (2079 - 2089 minus 2015 - 2024): ", ts_change_decades.values, "K")

# %%
