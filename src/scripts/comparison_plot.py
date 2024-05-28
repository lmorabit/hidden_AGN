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


## lotss data
cochrane = Table.read( paths.static / 'cochrane_2023_table1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table2.csv', format='csv', delimiter=',' )

## read in vmaxes
lockman_vmaxes = Table.read( paths.data / 'lockman_vmaxes.fits', format='fits' )
elais_vmaxes = Table.read( paths.data / 'en1_vmaxes.fits', format='fits' )

## combine these into a single catalogue for plotting
keep_cols = ['Total_flux_dr','Z_BEST','vmax','agn_vmax','sf_vmax','AGN_flux','SF_flux', 'Overall_class','Mass_cons','SFR_cons']
lockman_tmp = lockman_vmaxes[keep_cols]
elais_tmp = elais_vmaxes[keep_cols]
vmaxes = vstack([lockman_tmp,elais_tmp])

## get indices for galaxy identifications
SFG_idx = np.where( vmaxes['Overall_class'] == 'SFG' )[0]
Cochrane_SFG_idx = np.where( np.logical_or( vmaxes['Overall_class'] == 'SFG', vmaxes['Overall_class'] == 'RQAGN' ) )[0]

SFG_idx = Cochrane_SFG_idx

RG_idx = np.where( np.logical_or( vmaxes['Overall_class'] == 'HERG', vmaxes['Overall_class'] == 'LERG' ) )[0]
RQ_idx = np.where( vmaxes['Overall_class'] == 'RQAGN' )[0]
AGN_idx = np.union1d(RG_idx, RQ_idx)
AGN_idx = RG_idx

## RLF bins
redshift_bins = np.array([zmin,zmax])
lum_bins = np.arange( 20.5, 27, 0.3 ) # + np.log10( np.power( (144./1400.), si ) )

## get radio luminosities
Lrad = radio_power( vmaxes['Total_flux_dr'], vmaxes['Z_BEST'], spectral_index=si )
agn_Lrad = radio_power( vmaxes['AGN_flux'], vmaxes['Z_BEST'], spectral_index=si )
sf_Lrad = radio_power( vmaxes['SF_flux'], vmaxes['Z_BEST'], spectral_index=si )
## calculate the SFR
sfrs = calculate_SFR( sf_Lrad, vmaxes['Mass_cons'] )

## take the log
log10_Lrad = np.log10(Lrad)
agn_log10_Lrad = log10_when_zeros(agn_Lrad)
sf_log10_Lrad = log10_when_zeros(sf_Lrad)

## by galaxy
gal_sfg_log10_Lrad = log10_Lrad[SFG_idx]
gal_rg_log10_Lrad = log10_Lrad[AGN_idx]
gal_sfg_vmax = vmaxes['vmax'][SFG_idx]
gal_rg_vmax = vmaxes['vmax'][AGN_idx]


## calculate RLFs
rhos = []
agn_rhos = []
sf_rhos = []
gal_sfg_rhos = []
gal_rg_rhos = []
for i in np.arange(1,len(lum_bins)):
    ## change base from 10 to e to match units in Kondapally and Cochrane
    delta_log_L = (lum_bins[i] - lum_bins[i-1]) * np.log(10.)
    lum_idx = np.where(np.logical_and( log10_Lrad >= lum_bins[i-1], log10_Lrad < lum_bins[i] ) )[0]
    agn_lum_idx = np.where(np.logical_and( agn_log10_Lrad >= lum_bins[i-1], agn_log10_Lrad < lum_bins[i] ) )[0]
    sf_lum_idx = np.where(np.logical_and( sf_log10_Lrad >= lum_bins[i-1], sf_log10_Lrad < lum_bins[i] ) )[0]
    gal_sfg_idx = np.where( np.logical_and( gal_sfg_log10_Lrad >= lum_bins[i-1], gal_sfg_log10_Lrad < lum_bins[i] ) )[0]
    gal_rg_idx = np.where( np.logical_and( gal_rg_log10_Lrad >= lum_bins[i-1], gal_rg_log10_Lrad < lum_bins[i] ) )[0]
    if len(lum_idx) > 0:
        rho = np.log10( np.sum( 1. / vmaxes['vmax'][lum_idx] ) / delta_log_L ) 
    else:
        rho = 0
    if len(agn_lum_idx) > 0:
        agn_rho = np.log10( np.sum( 1. / vmaxes['agn_vmax'][agn_lum_idx] ) / delta_log_L ) 
    else:
        agn_rho = 0
    if len(sf_lum_idx) > 0:
        sf_rho = np.log10( np.sum( 1. / vmaxes['sf_vmax'][sf_lum_idx] ) / delta_log_L ) 
    else:
        sf_rho = 0
    if len(gal_sfg_idx) > 0:
        gal_sf_rho = np.log10( np.sum( 1. / gal_sfg_vmax[gal_sfg_idx] ) / delta_log_L )
    else:
        gal_sf_rho = 0
    if len(gal_rg_idx) > 0:
        gal_rg_rho = np.log10( np.sum( 1. / gal_rg_vmax[gal_rg_idx] ) / delta_log_L )
    else:
        gal_rg_rho = 0

    rhos.append(rho)
    agn_rhos.append(agn_rho)
    sf_rhos.append(sf_rho)
    gal_sfg_rhos.append(gal_sf_rho)
    gal_rg_rhos.append(gal_rg_rho)

lum_func = np.asarray(rhos)
agn_lum_func = np.asarray(agn_rhos)
sf_lum_func = np.asarray(sf_rhos)
gal_agn_lum_func = np.asarray(gal_rg_rhos)
gal_sf_lum_func = np.asarray(gal_sfg_rhos)

lum_bin_cens = lum_bins[0:-1] + 0.5*(lum_bins[1]-lum_bins[0])


fsizex = 14
fsizey = 5
sbsizex = 0.8
sbsizey = 0.8
plxlims = (20.1,27)
plylims = (-6.5,-1)

## set some colors
kond = 'black'
coch = 'gray'
sfggalc = mycols[80]
agngalc = mycols_m[120]
sfc = mycols[180]
agnc = mycols_m[180]

fig = plt.figure( figsize=(fsizex,fsizey) )
## Left panel: Showing that we reproduce the RLFs for galaxies even when using smaller area
p1 = plt.axes([0.05,0.1,sbsizex*fsizey/fsizex,sbsizey])
## plot the previous data
p1.plot( cochrane['logL150'], cochrane['logPhi'], color=coch, label='Cochrane et al. (2023)', linewidth=2 )
p1.plot( kondapally['logL150'], kondapally['logPhi'], color=kond, label='Kondapally et al. (2022)', linewidth=2)
## plot the new data, filtering zeros
non_zero = np.where( gal_agn_lum_func != 0.0 )[0]
p1.plot( lum_bin_cens[non_zero], gal_agn_lum_func[non_zero], color=agngalc, label='AGN galaxies', linewidth=3, linestyle='dotted' )
non_zero = np.where( gal_sf_lum_func != 0.0 )[0]
p1.plot( lum_bin_cens[non_zero], gal_sf_lum_func[non_zero], color=sfggalc, label='SF galaxies', linewidth=3, linestyle='dotted' )
p1.axes.set_xlim(plxlims)
p1.axes.set_ylim(plylims)
p1.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p1.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

handles, labels = p1.get_legend_handles_labels()
order = [3,2,0,1]
p1.legend([handles[idx] for idx in order],[labels[idx] for idx in order])


## Middle panel: same as left panel but new data by activity now
p2 = plt.axes([0.05+sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,sbsizey])
## plot the previous data
p2.plot( cochrane['logL150'], cochrane['logPhi'], color=coch, label='Cochrane et al. (2023)', linewidth=2)
p2.plot( kondapally['logL150'], kondapally['logPhi'], color=kond, label='Kondapally et al. (2022)', linewidth=2)
## plot the lofar data, filtering zeros
#non_zero = np.where( lum_func != 0.0 )[0]
#p2.plot( lum_bin_cens[non_zero], lum_func[non_zero], color='black', label='Total' )
non_zero = np.where( agn_lum_func != 0.0 )[0]
p2.plot( lum_bin_cens[non_zero], agn_lum_func[non_zero], color=agnc, label='AGN activity', linewidth=3, alpha=0.75 )
non_zero = np.where( sf_lum_func != 0.0 )[0]
p2.plot( lum_bin_cens[non_zero], sf_lum_func[non_zero], color=sfc, label='SF activity', linewidth=3, alpha=0.75 )
p2.axes.set_xlim(plxlims)
p2.axes.set_ylim(plylims)
p2.yaxis.set_visible(False)
p2.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
#p2.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
handles, labels = p2.get_legend_handles_labels()
order = [3,2,0,1]
p2.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

## Right panel (top): RLFs using new data, galaxies and activity
p3 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.42,sbsizex*fsizey/fsizex,0.6*sbsizey])
non_zero = np.where( gal_agn_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], gal_agn_lum_func[non_zero], color=agngalc, label='AGN galaxies', linewidth=3, alpha=0.75, linestyle='dotted' )
non_zero = np.where( gal_sf_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], gal_sf_lum_func[non_zero], color=sfggalc, label='SF galaxies', linewidth=3, alpha=0.75, linestyle='dotted' )
non_zero = np.where( agn_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], agn_lum_func[non_zero], color=agnc, label='AGN activity', linewidth=3, alpha=0.75 )
non_zero = np.where( sf_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], sf_lum_func[non_zero], color=sfc, label='SF activity', linewidth=3, alpha=0.75 )
p3.axes.set_xlim(plxlims)
p3.axes.set_ylim(plylims)
p3.xaxis.set_visible(False)
p3.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p3.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
handles, labels = p3.get_legend_handles_labels()
order = [3,1,2,0]
p3.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

## Right panel (bottom): ratio of the RLFs by galaxies and activity
p4 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,0.4*sbsizey])
p4.plot( (19,27), (1,1), color='gray', linestyle='dashed', linewidth=1.5 )
non_zero_idx = np.where( np.logical_and( agn_lum_func != 0, gal_agn_lum_func != 0 ) )[0]
ratio = np.power(10., agn_lum_func[non_zero_idx]) / np.power( 10., gal_agn_lum_func[non_zero_idx]) 
p4.plot( lum_bin_cens[non_zero_idx], ratio, color=agnc, label='AGN', linewidth=3 )
non_zero_idx = np.where( np.logical_and( sf_lum_func != 0, gal_sf_lum_func != 0 ) )[0]
ratio = np.power(10., sf_lum_func[non_zero_idx] ) / np.power( 10., gal_sf_lum_func[non_zero_idx] )
p4.plot( lum_bin_cens[non_zero_idx], ratio, color=sfc, label='SF', linewidth=3 )
p4.axes.set_xlim(plxlims)
#p4.axes.set_ylim((0.8,1.1))
p4.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p4.set_ylabel('Activity / Galaxy')
p4.legend()
fig.savefig(paths.figures / 'deep_fields_RLFs.png',dpi=300)
fig.clear()
plt.close()

'''
Redshift bins

Rohit:
0.5 - 1.0
1.0 - 1.5
1.5 - 2.0
2.0 - 2.5


Rachael:
0.1 - 0.4 
0.4 - 0.6 
0.6 - 0.8
0.8 - 1.0 
1.0 - 1.3
1.3 - 1.6
1.6 - 2.0 
2.0 - 2.5
2.5 - 3.3
3.3 - 4.6
4.6 - 5.7

'''


## read in SIMBA
Simba_SF = Table.read( paths.static / 'RLFS_50MYR_SF.csv', format='csv', delimiter=',' )
Simba_AGN = Table.read( paths.static / 'RLFS_50MYR_AGN.csv', format='csv', delimiter=',' )
## shift using spectral index
Simba_SF['x']  = Simba_SF['x'] + np.log10( np.power( (144./1400.), si ) )
Simba_AGN['x']  = Simba_AGN['x'] + np.log10( np.power( (144./1400.), si ) )

fig = plt.figure( figsize=(5,5) )
## plot simba
#plt.plot( Simba_SF['x'], Simba_SF['Curve1']+np.log10(np.log(10)), color='green', linewidth=2.5, label='Simba SF' )
#plt.plot( Simba_AGN['x'], Simba_AGN['Curve2']+np.log10(np.log(10)), color='blue', linewidth=2.5, label='Simba AGN' )

plt.plot( Simba_SF['x'], Simba_SF['Curve1'], color='green', linewidth=2.5, label='Simba SF' )
plt.plot( Simba_AGN['x'], Simba_AGN['Curve2'], color='blue', linewidth=2.5, label='Simba AGN' )


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
plt.xlabel('log('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log('+r'$\Phi$'+' [mag'+r'$^{-1}$'+' Mpc'+r'$^{-3}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'test.png',dpi=300)
fig.clear()





## SED SFR vs radio luminosity before and after correcting for AGN component
fig = plt.figure( figsize=(5,5) )
plt.plot( vmaxes['SFR_cons'][SFG_idx], np.log10(Lrad[SFG_idx]), '.', color='gray', label='original' )
plt.plot( vmaxes['SFR_cons'][SFG_idx], np.log10(sf_Lrad[SFG_idx]), '.', color=mycols_m[20], label='corrected' )
plt.xlim((-6,4))
plt.legend()
plt.savefig( paths.figures / 'SFR_Lrad.png', dpi=300 )
fig.clear()
plt.close()

## difference between SED SFR and the SFR calculated from the SF luminosity after removing AGN component
fig = plt.figure( figsize=(5,5) )
plt.plot( (-6,4),(-6,4),color='gray')
plt.plot( vmaxes['SFR_cons'][SFG_idx],sfrs[SFG_idx],'.' )
plt.xlim((-6,4))
plt.ylim((-6,4))
plt.savefig( paths.figures / 'SFR_vs_SFR.png', dpi=300 )
fig.clear()
plt.close()

## true radio SFR (i.e., AGN component removed) vs. total radio luminosity relationship.
## do this for SFGs only?
fig = plt.figure( figsize=(5,5) )
plt.plot( sfrs[SFG_idx], np.log10(Lrad[SFG_idx]),'.')
plt.xlim((-6,4))
plt.savefig( paths.figures / 'radioSFR_vs_Lrad.png', dpi=300 )
fig.clear()
plt.close()
