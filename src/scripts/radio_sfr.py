import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join, vstack
from helper_functions import *
from astropy.io import fits
import glob
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
## for T_b
T_e = 1e4
ref_freqs = np.array([0.003,0.01,0.03,0.1,0.3,1,3])*1e9
freqs_GHz = np.arange( 1e-3, 1e2, 1e-3 )


## lotss data
cochrane = Table.read( paths.static / 'cochrane_2023_table1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table2.csv', format='csv', delimiter=',' )

t = Table.read( paths.static / 'redshift_bins.csv', format='csv' )
zbin_starts = t['zbin_starts']
zbin_ends = t['zbin_ends']

## read in values from catalogues, combine, and get Tb information and do brightness temp separation
lockman = Table.read( paths.static / 'lockman_03_matched_inMOC_inHR.fits', format='fits' )
elais = Table.read( paths.static / 'en1_03_matched_inMOC_inHR.fits', format='fits' )
tmp = vstack([lockman,elais])
tmp = get_tb_information( tmp, im_weight=0.5, maj_lim=0.4, min_lim=0.3, T_e=T_e, alpha=si, ref_freqs=ref_freqs, freqs_GHz=freqs_GHz, use_z=False )
tmp = do_SFR_AGN_separation( tmp )

## combine these into a single catalogue for plotting
keep_cols = ['Total_flux_dr', 'Z_BEST', 'AGN_flux', 'SF_flux', 'Overall_class', 'Mass_cons', 'SFR_cons', 'Resolved']
lotss = tmp[keep_cols]

## get rid of SFR_cons = -99
good_SFR = np.where( lotss['SFR_cons'] > -99 )
lotss = lotss[good_SFR]
## and where mass = 1.0
good_masses = np.where( lotss['Mass_cons'] > 1 )
lotss = lotss[good_masses]

## find an index for unresolved, non radio-excess sources
non_radio_excess = np.where( np.logical_and( lotss['Overall_class'] != 'HERG', lotss['Overall_class'] != 'LERG' ) )[0]
unresolved = np.where( lotss['Resolved'] == 'U' )[0]
combined_idx = np.intersect1d( non_radio_excess, unresolved )
#lotss = lotss[combined_idx]

## use only star forming galaxies
#sfg_idx = np.where( lotss['Overall_class'] == 'SFG' )[0]
#lotss = lotss[sfg_idx]

## limit redshift
#z_idx = np.where( np.logical_and( lotss['Z_BEST'] > 0.003, lotss['Z_BEST'] < 0.3 ) )[0]
#lotss = lotss[z_idx]

pxlims = (-4,5)
pylims = (18,30)

## set some colors
kond = 'black'
coch = 'gray'
sfggalc = mycols[80]
agngalc = mycols_m[120]
sfc = mycols[180]
agnc = mycols_m[180]


Lrad = radio_power( lotss['Total_flux_dr'], lotss['Z_BEST'], spectral_index=si )
agn_Lrad = radio_power( lotss['AGN_flux'], lotss['Z_BEST'], spectral_index=si )
sf_Lrad = radio_power( lotss['SF_flux'], lotss['Z_BEST'], spectral_index=si )
## calculate the SFR
sfrs = calculate_SFR_from_Lrad_mass( sf_Lrad, lotss['Mass_cons'] )
sfrs_smith = calculate_SFR_from_Lrad_mass( Lrad, lotss['Mass_cons'] )

Lrad_corr = sf_Lrad
#Lrad_corr = np.copy(Lrad)
#Lrad_corr[combined_idx] = sf_Lrad[combined_idx]

nbin = 1000
rl_xvals, rl_yvals = find_ridgeline( lotss['SFR_cons'], np.log10(Lrad), min_num=nbin, sfrmin=-0.4, nbins=20 )
#rl_xvals, rl_yvals = find_ridgeline( lotss['SFR_cons'][non_radio_excess], np.log10(Lrad[non_radio_excess]), min_num=nbin, sfrmin=-1.0, nbins=25 )
rl_xvals_corr, rl_yvals_corr = find_ridgeline( lotss['SFR_cons'][combined_idx], np.log10(Lrad_corr[combined_idx]), min_num=nbin, sfrmin=-1.0, nbins=23 )

## SED SFR vs radio luminosity before and after correcting for AGN component
fig = plt.figure( figsize=(5,5) )
plt.scatter( lotss['SFR_cons'], np.log10(Lrad), marker='o', color='gray', edgecolor='none', label='original', alpha=0.75 )
non_zero = np.where( sf_Lrad > 1 )[0]
plt.scatter( lotss['SFR_cons'][non_zero], np.log10(sf_Lrad[non_zero]), marker='o', color='none', edgecolor=mycols_m[20], label='SF only', alpha=0.75 )
x_sfrs = np.arange(-6,5)
plt.plot( x_sfrs, linear(x_sfrs, 1.08, 22.24 ))
plt.plot( rl_xvals, rl_yvals, color=mycols[100], label='Original')
plt.plot( rl_xvals_corr, rl_yvals_corr, color=mycols_m[200], label='SF only')
plt.xlim(pxlims)
plt.ylim(pylims)
plt.legend()
plt.savefig( paths.figures / 'SFR_Lrad.png', dpi=300 )
fig.clear()
plt.close()

## difference between SED SFR and the SFR calculated from the SF luminosity after removing AGN component
fig = plt.figure( figsize=(9,5) )
gs = fig.add_gridspec(1,2,hspace=0)
axs = gs.subplots( sharex=True, sharey=True )
axs[0].plot( lotss['SFR_cons'],sfrs,'.' )
axs[0].plot( (-6,5),(-6,5),color='gray')
axs[0].set_xlim(pxlims)
axs[0].set_ylim(pxlims)
axs[0].set_xlabel('SED SFR')
axs[0].set_ylabel('SFR from SF separation')
axs[1].plot( sfrs_smith, sfrs, '.' )
axs[1].plot( (-6,5),(-6,5),color='gray')
axs[1].set_xlim(pxlims)
axs[1].set_ylim(pxlims)
axs[1].set_xlabel('SFR from total Lrad')
fig.tight_layout()
fig.savefig( paths.figures / 'SFR_vs_SFR.png', dpi=300 )
fig.clear()
plt.close()


