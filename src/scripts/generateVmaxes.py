import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
import os
from astropy.cosmology import WMAP9 as cosmo
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

#########################################
## default parameters to adjust
si = -0.8
zmin = 0.003  ## matches Mauch & Sadler 2007
zmax = 0.3   ## matches Mauch & Sadler 2007
dz = 0.0001  ## matches Cochrane and Kondapally


## Lockman
infile = paths.static / 'lockman_final_cross_match_catalogue-v1.0_classifications.fits'
field = os.path.basename(infile).split('_')[0]
outfits = field + '_vmax.fits'
if not os.path.exists(outfits):
    lotss = Table.read( infile, format='fits' )
    vmaxes = get_vmax( lotss, field, zmin=zmin, zmax=zmax, dz=dz, si=si )
    vmaxes.write( paths.data / outfits, format='fits', overwrite=True )
else:
    vmaxes = Table.read( outfits )


## calculate RLFs
redshift_bins = np.array([zmin,zmax])
lum_bins = np.arange( 19.5, 27, 0.5 ) + np.log10( np.power( (144./1400.), si ) )
Lrad = radio_power( lotss['Total_flux'], lotss['Z_BEST'], spectral_index=si )



RLF = RLF_from_zmax( Lrad, lotss['Z_BEST'], lum_bins, redshift_bins, lotss_zmax, area_deg2, area_units='deg2', error_type='rms' )

## first  
#outname = 'RLF_{:s}_{:s}.fits'.format(str(redshift_bins[0]).replace('.','p'), str(redshift_bins[1]).replace('.','p') )
outname = 'RLF.fits'
RLF.write( paths.data / outname, format='fits' )

outname = 'RLF_corr.fits'
RLF_corr.write( paths.data / outname, format='fits' )
