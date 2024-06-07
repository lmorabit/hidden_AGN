import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join, vstack
from helper_functions import *
from astropy.io import fits
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

##############################################################
## Plotting housekeeping

matplotlib.rcParams['legend.frameon'] = False
matplotlib.rcParams['axes.labelsize'] = 'large'
## set up some colours
n = 255
mycols = plt.cm.viridis(np.linspace(0, 1,n))
mycols_m = plt.cm.inferno(np.linspace(0, 1,n))

##############################################################

## set a default spectral index
si = -0.7
zmin = 0.003  ## matches Mauch & Sadler 2007
zmax = 0.3   ## matches Mauch & Sadler 2007
dz = 0.0001  ## matches Cochrane and Kondapally

## read in SIMBA
Simba_SF = Table.read( paths.static / 'RLFS_50MYR_SF.csv', format='csv', delimiter=',' )
Simba_AGN = Table.read( paths.static / 'RLFS_50MYR_AGN.csv', format='csv', delimiter=',' )
## shift using spectral index
Simba_SF['x']  = Simba_SF['x'] + np.log10( np.power( (144./1400.), si ) )
Simba_AGN['x']  = Simba_AGN['x'] + np.log10( np.power( (144./1400.), si ) )


## read in vmaxes
lockman_vmaxes = Table.read( paths.static / 'lockman_vmaxes_zmin0.003_zmax0.3.fits', format='fits' )
elais_vmaxes = Table.read( paths.static / 'en1_vmaxes_zmin0.003_zmax0.3.fits', format='fits' )

## combine these into a single catalogue for plotting
keep_cols = ['Total_flux_dr','Z_BEST','vmax','agn_vmax','sf_vmax','AGN_flux','SF_flux', 'Overall_class','Mass_cons','SFR_cons']
lockman_tmp = lockman_vmaxes[keep_cols]
elais_tmp = elais_vmaxes[keep_cols]
vmaxes = vstack([lockman_tmp,elais_tmp])

lum_bins, lum_func, agn_lum_func, sf_lum_func, gal_agn_lum_func, gal_sf_lum_func = get_RLFs( vmaxes, zmin, zmax, lmin=20.5, lmax=27, dl=0.3, si=si )
lum_bin_cens = lum_bins[0:-1] + 0.5*(lum_bins[1]-lum_bins[0])

plxlims = (20.1,27)
plylims = (-6.5,-1)

## set some colors
kond = 'black'
coch = 'gray'
sfggalc = mycols[80]
agngalc = mycols_m[120]
sfc = mycols[180]
agnc = mycols_m[180]

fig = plt.figure( figsize=(5,5) )
## plot simba - convert from mag to Lum
# mag = -2.5 * Lum + zeropoint
# so Lum = mag / -2.5
plt.plot( Simba_SF['x'], (Simba_SF['Curve1']+np.log10(2.5)), color='green', linewidth=2.5, label='Simba SF' )
plt.plot( Simba_AGN['x'], (Simba_AGN['Curve2']+np.log10(2.5)), color='blue', linewidth=2.5, label='Simba AGN' )

non_zero = np.where( gal_agn_lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], gal_agn_lum_func[non_zero], color=agngalc, label='AGN galaxies', linewidth=3, linestyle='dotted' )
non_zero = np.where( gal_sf_lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], gal_sf_lum_func[non_zero], color=sfggalc, label='SF galaxies', linewidth=3, linestyle='dotted' )

non_zero = np.where( agn_lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], agn_lum_func[non_zero], color=agnc, label='AGN activity', linewidth=3, alpha=0.75 )
non_zero = np.where( sf_lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], sf_lum_func[non_zero], color=sfc, label='SF activity', linewidth=3, alpha=0.75 )
plt.xlim(plxlims)
plt.ylim(plylims)
plt.xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'simba_comparison.png',dpi=300)
fig.clear()


