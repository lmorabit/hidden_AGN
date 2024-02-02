import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

#########################################
## default parameters to adjust
si = -0.8
zmin = 0.003  ## matches Mauch & Sadler 2007
zmax = 0.3   ## matches Mauch & Sadler 2007
flux_cut = 150e-6

## read in lofar data
lotss = Table.read( paths.static / 'lockman_final_cross_match_catalogue-v1.0.fits', format='fits' )
## remove things with no redshifts, i.e. those which are masked
good_z = np.where( np.logical_not(lotss['Z_BEST'].mask) )[0]
lotss = lotss[good_z]
## remove the mask to avoid warnings later on
lotss['Z_BEST'] = np.ma.getdata(lotss['Z_BEST'])

test = True
if test:
    good_z = np.where(lotss['Z_BEST'] < zmax )[0]
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
redshift_bins = np.arange( 0.0, zmax+0.001,zmax )

redshift_bins = np.array([zmin,zmax])
lum_bins = np.arange( 19.8, 27, 0.4 ) + np.log10( np.power( (144./1400.), si ) )
Lrad = radio_power( lotss['Total_flux'], lotss['Z_BEST'], spectral_index=si )
RLF = RLF_from_zmax( Lrad, lotss['Z_BEST'], lum_bins, redshift_bins, lotss_zmax, area_deg2, area_units='deg2', error_type='rms' )

## first  
#outname = 'RLF_{:s}_{:s}.fits'.format(str(redshift_bins[0]).replace('.','p'), str(redshift_bins[1]).replace('.','p') )
outname = 'RLF.fits'
RLF.write( paths.data / outname, format='fits' )


