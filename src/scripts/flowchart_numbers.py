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


fields = ['lockman', 'en1']

for field in fields:
    print('Starting with field: {:s}'.format(field))
    infile = paths.static / '{:s}_03_matched_inMOC_inHR.fits'.format(field)
    field = os.path.basename(infile).split('_')[0]
    rms_image = paths.static / '{:s}_DR1_rms_masked.fits'.format(field)
    lotss = Table.read( infile, format='fits' )
    ## add brightness temperature information
    lotss = get_tb_information( lotss, im_weight=0.5, maj_lim=0.4, min_lim=0.3, T_e=T_e, alpha=si, ref_freqs=ref_freqs, freqs_GHz=freqs_GHz, use_z=False )
    lotss = do_SFR_AGN_separation( lotss )
    ## get number of detectable sources
    ndet = len( np.where(lotss['Detectability_SNR'] > 5.)[0] )
    with open( paths.output / '{:s}_detectable.txt'.format(field), 'w' ) as f:
        f.write( str(ndet) )





