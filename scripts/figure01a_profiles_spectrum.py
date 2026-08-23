"""
Create intuition about the behavior of the FIR by using the standard atmospheric profiles and the TIRS2 spectral response function.
"""

# %%
import numpy as np
import netCDF4
import xarray as xr
import matplotlib.pyplot as plt
from scipy import constants
import matplotlib.gridspec as gridspec
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

import seaborn as sns
cmap = sns.color_palette("colorblind")

# %%
def planck_function(wavelength, temperature):

    wavelength = wavelength*1e-6

    # Calculate the Planck function
    step_1 = wavelength*temperature
    radiance_micron = 1e-6*(2 * h * c**2) / (wavelength**5 * (np.exp((h * c) / (step_1 * k)) - 1))
    
    return radiance_micron


def rad_calc(ds):

    wn = ds.wavenum
    bt = ds.btemp
    rad = ds.rad
    
    wn_list = wn.data
    wl_list = 1e4/wn_list
    
    inv_wl_list = 1/wl_list # µm-1
    
    rad_units = rad/1e3*1e4 #mW / (m^2 sr cm^-1) -> (W/m^2 sr µm-1)
    rad_wl = np.array(rad_units*inv_wl_list**2) # (W/m^2 sr µm)
    
    new_wl = np.arange(4,50.03,0.03)
    
    rad_wl_resamp = np.array(np.interp(new_wl, np.flip(wl_list), np.flip(rad_wl)))
    
    return(new_wl,rad_wl_resamp)


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

    save_dir = OUTPUT_ROOT
    load_dir = DATA_ROOT
    SRFfile1 = load_dir / 'PREFIRE_TIRS2_SRF_v13_2024-09-15.nc'

    SRF_f = netCDF4.Dataset(SRFfile1, 'r')

    print("Spectral response function wavelengths:")
    print(SRF_f['channel_mean_wavelen'][5,0])
    print(SRF_f['channel_mean_wavelen'][12,0])
    print(SRF_f['channel_mean_wavelen'][27,0])

    SRF = SRF_f['srf'][:,:,0]
    SRF_wl = SRF_f['wavelen'][:]

    MWV_SRF = SRF[:,5]
    AW_SRF = SRF[:,12]
    FIR_SRF = SRF[:,27]

    h = constants.h  # Planck's constant
    c = constants.speed_of_light  # Speed of light
    k = constants.k  # Boltzmann's constant

    nc = netCDF4.Dataset(str(load_dir / 'Modtran_standard_profiles_PCRTM_levels.nc'),'r')
    a = 0

    pres = nc['pres'][:,a]
    dp = pres[1:] - pres[:-1]

    tlev = nc['temp'][:,a].astype(np.float32)
    # approximate conversion from ppm to q [g/kg]
    q = nc['h2o'][:,a] * 0.622 * 1e-6
    avg_mixing_ratio = (q[1:] + q[:-1]) / 2.0

    pwv = np.sum(avg_mixing_ratio * dp / 9.8)

    t_tropical = tlev[-1]*0.988

    a = 4
    tskin = nc['temp'][-1,a]
    tlev = nc['temp'][:,a].astype(np.float32)
    q = nc['h2o'][:,a] * 0.622 * 1e-6 #* 1e3

    t_polar = tlev[-1]*0.996

    new_wl = np.arange(4,50.03,0.03)
    tropical_planck = planck_function(new_wl, t_tropical)
    polar_planck = planck_function(new_wl, t_polar)

    FIR_wn = np.array([10,630])
    FIR_wl = 1e4/FIR_wn

    FIRH2O_wn = np.array([10, 500])
    FIRH2O_wl = 1e4/FIRH2O_wn

    NIR_wn = np.array([2080,1390])
    NIR_wl = 1e4/NIR_wn

    WIN_wn = np.array([980,820])
    WIN_wl = 1e4/WIN_wn

    MWV_c = cmap[0]
    AW_c = cmap[1]
    FIR_c = cmap[2]

    # %%
    # Create figure
    fig_width_cm = 30
    fig_height_cm = 30
    inches_per_cm = 1 / 2.54                      # Convert cm to inches
    fig_width = fig_width_cm * inches_per_cm         # width in inches
    fig_height = fig_height_cm * inches_per_cm       # height in inches
    fig_size = [fig_width, fig_height]

    plt.rc('text', usetex=False) # so that LaTeX is not needed when creating a PDF with PdfPages later on
    fig = plt.figure()
    fig.set_size_inches(fig_size)

    widths = [1]
    heights = [1,1,1]
    print(np.shape(heights))
    gs = gridspec.GridSpec(len(heights),1, wspace=0.1, hspace=0.5,width_ratios=widths,height_ratios=heights)
    axes = []

    axes.append(fig.add_subplot(gs[0]))

    # 4 mm TPW, -16.9 C
    wl_list,rad_wl = rad_calc(xr.open_dataset(str(load_dir / 'PCRTM_forward_stdatm5_subarc_winter.nc')))

    # In terms of wl (Wm^-2 sr^-1 µm^-1)
    axes[-1].fill_between(new_wl, rad_wl, polar_planck,facecolor='lightgray',zorder=2)

    plt.plot(wl_list,rad_wl,lw=1.8,c='k',zorder=3)
    plt.plot(new_wl,polar_planck,lw=1,c='k',ls=':',zorder=3)

    axes[-1].fill_between(NIR_wl, 0, 10, color=MWV_c, alpha=0.2,zorder=1)
    axes[-1].fill_between(WIN_wl, 0, 10, color=AW_c, alpha=0.2,zorder=1)
    axes[-1].fill_between(FIR_wl, 0, 10, color=FIR_c, alpha=0.12,zorder=1)
    axes[-1].fill_between(FIRH2O_wl, 0, 10, color=FIR_c, alpha=0.20,zorder=1)

    plt.xlim(4,30)
    plt.ylim(0,10)
    axes[-1].set_xticks(np.arange(4,31,1))
    plt.tick_params(axis='x', labelsize=13)
    plt.tick_params(axis='y', labelsize=13)
    plt.grid()


    axes[-1].text(6, 9.9, 'MWV',fontsize=15,c=MWV_c,ha='center',va='top')
    axes[-1].text(11.5, 0.1, 'AW',fontsize=15,c=AW_c,ha='center',va='bottom')
    axes[-1].text(21.1, 0.1, 'FIR-H$_2$O',fontsize=15,c=FIR_c,ha='left',va='bottom')
    axes[-1].text(17.1, 0.1, 'FIR',fontsize=15,c=FIR_c,ha='left',va='bottom')

    axes[-1].text(0, 1.02, '(a)',fontsize=15, transform=axes[-1].transAxes,ha='left',va='bottom', fontweight='bold')
    axes[-1].text(0.5, 1, 'Polar',fontsize=23, transform=axes[-1].transAxes,ha='center',va='bottom')
    axes[-1].text(0.995, 1, 'T = -17 °C, CWV = 4 mm',fontsize=14, transform=axes[-1].transAxes,ha='right',va='bottom')

    plt.xlabel('Wavelength (μm)',size=12)
    plt.ylabel('Radiance (Wm$^{-2}$ sr$^{-1}$ μm$^{-1}$)',size=12)


    secax = axes[-1].twinx()
    secax.yaxis.set_inverted(True)
    plt.plot(SRF_wl,MWV_SRF,c='k',lw=1,ls='--')
    plt.plot(SRF_wl,AW_SRF,c='k',lw=1,ls='--')
    plt.plot(SRF_wl,FIR_SRF,c='k',lw=1,ls='--')
    plt.ylim(0.5,0)
    secax.yaxis.set_visible(False)

    secax.text(5.8, 0.22, '5.9 μm',fontsize=10,c='k',ha='center',va='top')
    secax.text(12, 0.24, '11.8 μm',fontsize=10,c='k',ha='center',va='top')
    secax.text(24.5, 0.28, '24.5 μm',fontsize=10,c='k',ha='center',va='top')
    #######################################################

    axes.append(fig.add_subplot(gs[1]))

    # 51 mm TPW, 23 C
    wl_list,rad_wl = rad_calc(xr.open_dataset(f'{load_dir}/PCRTM_forward_stdatm1_tropic.nc'))

    # In terms of wl (Wm^-2 sr^-1 µm^-1)
    axes[-1].fill_between(new_wl, rad_wl, tropical_planck,facecolor='lightgray',zorder=2)

    plt.plot(wl_list,rad_wl,lw=1.8,c='k',zorder=3)
    plt.plot(new_wl,tropical_planck,lw=1,c='k',ls=':',zorder=3)

    axes[-1].fill_between(NIR_wl, 0, 10, color=MWV_c, alpha=0.2,zorder=1)
    axes[-1].fill_between(WIN_wl, 0, 10, color=AW_c, alpha=0.2,zorder=1)
    axes[-1].fill_between(FIR_wl, 0, 10, color=FIR_c, alpha=0.12,zorder=1)
    axes[-1].fill_between(FIRH2O_wl, 0, 10, color=FIR_c, alpha=0.20,zorder=1)

    plt.xlim(4,30)
    plt.ylim(0,10)
    axes[-1].set_xticks(np.arange(4,31,1))
    plt.tick_params(axis='x', labelsize=13)
    plt.tick_params(axis='y', labelsize=13)
    plt.grid()

    axes[-1].text(6, 9.9, 'MWV',fontsize=15,c=MWV_c,ha='center',va='top')
    axes[-1].text(11.5, 0.1, 'AW',fontsize=15,c=AW_c,ha='center',va='bottom')
    axes[-1].text(21.1, 0.1, 'FIR-H$_2$O',fontsize=15,c=FIR_c,ha='left',va='bottom')
    axes[-1].text(17.1, 0.1, 'FIR',fontsize=15,c=FIR_c,ha='left',va='bottom')

    axes[-1].text(0, 1.02, '(b)',fontsize=15, transform=axes[-1].transAxes,ha='left',va='bottom', fontweight='bold')
    axes[-1].text(0.5, 1, 'Tropical',fontsize=23, transform=axes[-1].transAxes,ha='center',va='bottom')
    axes[-1].text(0.995, 1, 'T = 23 °C, CWV = 50 mm',fontsize=14, transform=axes[-1].transAxes,ha='right',va='bottom')
    plt.xlabel('Wavelength (μm)',size=12)
    plt.ylabel('Radiance (Wm$^{-2}$ sr$^{-1}$ μm$^{-1}$)',size=12)

    # %%
    to_png(fig, 'fig1_upperpanels', loc=save_dir, dpi=300, bbox_inches="tight", ext='png')

# %%