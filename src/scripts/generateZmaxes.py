import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

## test comment

## set a default spectral index
si = -0.8

## read in Mauch & Sadler Table 5
mauch_sadler = Table.read( paths.static / 'mauch_sadler_table5.csv', format='csv', delimiter=',' )
## shift using spectral index
ms_144MHz = mauch_sadler['log10_P1p4GHz'] + np.log10( np.power( (144./1400.), si ) ) 

## read in lofar data
lotss = Table.read( paths.static / 'lockman_final_cross_match_catalogue-v1.0.fits', format='fits' )
## remove things with no redshifts
## for testing
#good_z = np.where(lotss['Z_BEST'] < 0.5 )[0]
good_z = np.where( np.isfinite(lotss['Z_BEST']) )[0]
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
