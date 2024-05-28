import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
import os
## cosmology to match Kondapally and Cochrane
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import time

#########################################
## default parameters to adjust
si = -0.7
#zmin = 0.003  ## matches Mauch & Sadler 2007
#zmax = 0.3   ## matches Mauch & Sadler 2007
dz = 0.0001  ## matches Cochrane and Kondapally
## for T_b
T_e = 1e4
ref_freqs = np.array([0.003,0.01,0.03,0.1,0.3,1,3])*1e9
freqs_GHz = np.arange( 1e-3, 1e2, 1e-3 )

sigma_cut = 5.

## completeness corrections
cochrane = Table.read( paths.static / 'cochrane_2023_tableA1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table1.csv', format='csv', delimiter=',' )

fields = ['lockman', 'en1']

zbin_starts = [ 0.003, 0.5, 1.0, 1.5, 2.0, 0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.5, 3.3, 4.6 ]
zbin_ends = [ 0.3, 1.0, 1.5, 2.0, 2.5, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 3.3, 4.6, 5.7 ]

for field in fields:
    print('Starting with field: {:s}'.format(field))
    infile = paths.static / '{:s}_03_matched_inMOC_inHR.fits'.format(field)
    field = os.path.basename(infile).split('_')[0]
    rms_image = paths.static / '{:s}_DR1_rms_masked.fits'.format(field)
    lotss = Table.read( infile, format='fits' )
    ## add brightness temperature information
    lotss = get_tb_information( lotss, im_weight=0.5, maj_lim=0.4, min_lim=0.3, T_e=T_e, alpha=si, ref_freqs=ref_freqs, freqs_GHz=freqs_GHz, use_z=False )
    lotss = do_SFR_AGN_separation( lotss )
    ## calculate vmaxes
    for i in np.arange(0,len(zbin_ztarts)):
        zmin = zbin_starts[i]
        zmax = zbin_ends[i]
        outfits = field + '_vmaxes_zmin'+str(zmin)+'_zmax'+str(zmax)+'.fits'
        vmaxes = get_vmax( lotss, field, col_suffix='_dr', zmin=zmin, zmax=zmax, dz=dz, si=si, sigma_cut=sigma_cut, rms_image=rms_image, cochrane=cochrane, kondapally=kondapally, test=False )
        vmaxes.write( paths.data / outfits, format='fits', overwrite=True )

    


## uncertainties done by bootstrapping. From Kondapally et al.:
## we performed bootstrap sampling (random sampling by replacement) of the catalogue to generate a distribution of 1000 realizations of the luminosity function. The lower and upper 1σ uncertainties on our luminosity functions are then determined from the 16th and 84th percentiles of the bootstrap samples. For the faint luminosity bins, where the samples are large and the uncertainties computed from bootstrapping correspondingly small, the true uncertainties are likely to be dominated by other factors such as the photometric redshift errors and source classiﬁcation uncertainties. Therefore, we set a minimum uncertainty of 0.03 dex in the luminosity functions reported, based on the ∼7 per cent photometric redshift outlier fraction in these ﬁelds
