import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

## set a default spectral index
si = -0.8

## read in Mauch & Sadler Table 5
mauch_sadler = Table.read( paths.static / 'mauch_sadler_table5.csv', format='csv', delimiter=',' )
## shift using spectral index
ms_144MHz = mauch_sadler['log10_P1p4GHz'] + np.log10( np.power( (144./1400.), si ) ) 

## read in lofar data
lotss = Table.read( paths.static / 'lockman_final_cross_match_catalogue-v1.0.fits', format='fits' )
## flux cut
flux_cut = 150e-6 ## 150 uJy
valid = np.where( np.logical_or(lotss['Total_flux'] > flux_cut, lotss['Peak_flux'] > flux_cut) )[0]
lotss = lotss[valid]
## remove things with no redshifts
## for testing
#good_z = np.where( np.isfinite(lotss['Z_BEST']) )[0]
zmax = 0.5
good_z = np.where(lotss['Z_BEST'] < 0.5 )[0]
lotss = lotss[good_z]

## Total area
rms_image = fits.open( paths.static / 'lockman_rms_starmask_optical.fits' )
rms = rms_image[0].data
rms_idx = np.where(rms <= flux_cut)
npix_total = len(rms_idx[0])
pixscale = np.abs(rms_image[0].header['CDELT1'])
pixarea = np.power(pixscale, 2.)
area_deg2 = npix_total * pixarea

## get z-max
lotss_zmax = RLF_calculate_zmax( flux_cut, lotss['Source_Name'], lotss['Total_flux'], lotss['E_Total_flux'], lotss['Z_BEST'], np.repeat(si,len(lotss)), outfile='zmaxes' )

## calculate RLFs
redshift_bins = np.arange( 0.0, zmax+0.001,0.5 )
lum_bins = np.arange( 19.8, 27, 0.4 ) + np.log10( np.power( (144./1400.), si ) )
Lrad = radio_power( lotss['Total_flux'], lotss['Z_BEST'], spectral_index=si )
RLF, RLF_up, RLF_lo, Lmed, Lmed_err, N_obj = RLF_from_zmax( Lrad, lotss['Z_BEST'], lum_bins, redshift_bins, lotss_zmax, area_deg2, area_units='deg2', error_type='rms' )


fig = plt.figure( figsize=(5,5) )
## plot
plt.plot( ms_144MHz, mauch_sadler['log10RLF_all'], color='black', linewidth=2.5, label='All' )
plt.plot( ms_144MHz, mauch_sadler['log10RLF_SF'], color='magenta', linewidth=2.5, label='SFG' )
plt.plot( ms_144MHz, mauch_sadler['log10RLF_RLAGN'], color='orange', linewidth=2.5, label='AGN' )
## plot the lofar data
plt.fill_between( Lmed, RLF_lo[:,1].flatten(), RLF_up[:,1].flatten(), color='purple', alpha=0.5 )
plt.plot( Lmed, RLF[:,1], 'o', color='purple', label='data' )
plt.xlim((20,28))
plt.ylim(-7.5,-2)
plt.xlabel('log('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log('+r'$\Phi$'+' [mag'+r'$^{-1}$'+' Mpc'+r'$^{-3}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'mauch_sadler_RLFs.png',dpi=300)
fig.clear()



