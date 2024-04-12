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
import time

#########################################
## default parameters to adjust
si = -0.8
zmin = 0.003  ## matches Mauch & Sadler 2007
zmax = 0.3   ## matches Mauch & Sadler 2007
dz = 0.0001  ## matches Cochrane and Kondapally

## completeness corrections
cochrane = Table.read( paths.static / 'cochrane_2023_tableA1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table1.csv', format='csv', delimiter=',' )

## Lockman
infile = paths.static / 'lockman_final_cross_match_catalogue-v1.0_classifications_catalogue_filtered_full_SNR5_fluxscaled_withoffset_noduplicates_with_lotss_DR1_detectable.fits'
field = os.path.basename(infile).split('_')[0]
rms_image = paths.static / '{:s}_rms_starmask_optical.fits'.format(field)
lotss = Table.read( infile, format='fits' )

## 6 arcsec vmaxes
outfits = field + '_6arcsec_vmax.fits'
vmaxes = get_vmax( lotss, field, col_suffix='_dr1', zmin=zmin, zmax=zmax, dz=dz, si=si, rms_image=rms_image, cochrane=cochrane, kondapally=kondapally )
vmaxes.write( paths.data / outfits, format='fits', overwrite=True )

print('writing file {:s} and then sleeping for 30 sec'.format(paths.data / outfits))
time.sleep(30)


## uncertainties done by bootstrapping. From Kondapally et al.:
## we performed bootstrap sampling (random sampling by replacement) of the catalogue to generate a distribution of 1000 realizations of the luminosity function. The lower and upper 1σ uncertainties on our luminosity functions are then determined from the 16th and 84th percentiles of the bootstrap samples. For the faint luminosity bins, where the samples are large and the uncertainties computed from bootstrapping correspondingly small, the true uncertainties are likely to be dominated by other factors such as the photometric redshift errors and source classiﬁcation uncertainties. Therefore, we set a minimum uncertainty of 0.03 dex in the luminosity functions reported, based on the ∼7 per cent photometric redshift outlier fraction in these ﬁelds
