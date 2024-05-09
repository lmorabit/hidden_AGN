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
si = -0.7
zmin = 0.003  ## matches Mauch & Sadler 2007
zmax = 0.3   ## matches Mauch & Sadler 2007
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

for field in fields:
    print('Starting with field: {:s}'.format(field))
    infile = paths.static / '{:s}_03_matched_inMOC_inHR.fits'.format(field)
    field = os.path.basename(infile).split('_')[0]
    rms_image = paths.static / '{:s}_DR1_rms_masked.fits'.format(field)
    lotss = Table.read( infile, format='fits' )
    ## add brightness temperature information
    compact_flux_per_SA, tb_from = get_tb_information( lotss, im_weight=0.5, maj_lim=0.4, min_lim=0.3, T_e=T_e, alpha=si, ref_freqs=ref_freqs, freqs_GHz=freqs_GHz, use_z=False )

############################################################
## ONLY VALID FOR T_b IDENTIFIED AND UNRESOLVED
    ## separate star formation and AGN luminosities
    lotss_sub_peak = lotss['Total_flux_dr'] - lotss['Peak_flux']
    e_dr_sub_peak = add_sub_error( lotss['E_Total_flux_dr'], lotss['E_Peak_flux'])
    AGN_lum = radio_power(lotss['Peak_flux'],lotss['z_best'])
    e_AGN_lum = radio_power( lotss['E_Peak_flux'], lotss['z_best'])
    SFR_lum = radio_power(lotss_sub_peak,lotss['z_best'])
    e_SFR_lum = radio_power( e_dr_sub_peak, lotss['z_best'] )
    ## calculate the star formation rates ... 
    ## Smith et al 2021:
    ## log10(L150) = (0.9+-0.01)*log10(SFR)+(0.33+-0.04)*log10(M/10^10)+22.22+-0.02
    ## log10(SFR) = ( log10(L150) - (22.22+-0.02) - (0.33+-0.04)*log10(M/10^10) ) / (0.9+-0.01)
    SFR_lum[np.where(SFR_lum < 0)] = np.nan
    nan_idx = np.unique(np.concatenate([np.where(np.isnan(SFR_lum))[0], np.where(np.isnan(lotss['Mass_cons']))[0]]))
    SFR_lum[nan_idx] = 1.
    lotss['Mass_cons'][nan_idx] = 1.
    sfr = ( np.log10(SFR_lum) - 22.22 - 0.33*np.log10(lotss['Mass_cons']) ) / 0.9
    ## re-set the nans
    SFR_lum[nan_idx] = np.nan
    lotss['Mass_cons'][nan_idx] = np.nan
    sfr[nan_idx] = np.nan
##########################################################


    ## 6 arcsec vmaxes
    outfits = field + '_6arcsec_vmaxes.fits'
    vmaxes = get_vmax( lotss, field, col_suffix='_dr', zmin=zmin, zmax=zmax, dz=dz, si=si, sigma_cut=sigma_cut, rms_image=rms_image, cochrane=cochrane, kondapally=kondapally, test=False )
    vmaxes.write( paths.data / outfits, format='fits', overwrite=True )

    unresolved_idx = np.where(lotss['Resolved'] == 'U')
    




#print('writing file {:s} and then sleeping for 30 sec'.format(str(paths.data / outfits)))
#time.sleep(30)


## uncertainties done by bootstrapping. From Kondapally et al.:
## we performed bootstrap sampling (random sampling by replacement) of the catalogue to generate a distribution of 1000 realizations of the luminosity function. The lower and upper 1σ uncertainties on our luminosity functions are then determined from the 16th and 84th percentiles of the bootstrap samples. For the faint luminosity bins, where the samples are large and the uncertainties computed from bootstrapping correspondingly small, the true uncertainties are likely to be dominated by other factors such as the photometric redshift errors and source classiﬁcation uncertainties. Therefore, we set a minimum uncertainty of 0.03 dex in the luminosity functions reported, based on the ∼7 per cent photometric redshift outlier fraction in these ﬁelds
