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
mycols_m = plt.cm.magma(np.linspace(0, 1,n))

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
lockman_vmaxes = Table.read( paths.data / 'lockman_6arcsec_vmaxes.fits', format='fits' )
elais_vmaxes = Table.read( paths.data / 'en1_6arcsec_vmaxes.fits', format='fits' )

keep_cols = ['Total_flux_dr','Z_BEST','vmax']

lockman_tmp = lockman_vmaxes[keep_cols]
elais_tmp = elais_vmaxes[keep_cols]

vmaxes = vstack([lockman_tmp,elais_tmp])


## calculate RLFs
redshift_bins = np.array([zmin,zmax])
lum_bins = np.arange( 19.5, 27, 0.5 ) # + np.log10( np.power( (144./1400.), si ) )
Lrad = radio_power( vmaxes['Total_flux_dr'], vmaxes['Z_BEST'], spectral_index=si )
agn_Lrad = radio_power( vmaxes['AGN_flux'], vmaxes['Z_BEST'], spectral_index=si )
sf_Lrad = radio_power( vmaxes['SF_flux'], vmaxes['Z_BEST'], spectral_index=si )
log10_Lrad = np.log10(Lrad)
agn_log10_Lrad = np.log10(agn_Lrad)
sf_log10_Lrad = np.log10(sf_Lrad)

rhos = []
agn_rhos = []
sf_rhos = []
for i in np.arange(1,len(lum_bins)):
    delta_log_L = lum_bins[i] - lum_bins[i-1]
    lum_idx = np.where(np.logical_and( log10_Lrad >= lum_bins[i-1], log10_Lrad < lum_bins[i] ) )[0]
    agn_lum_idx = np.where(np.logical_and( agn_log10_Lrad >= lum_bins[i-1], agn_log10_Lrad < lum_bins[i] ) )[0]
    sf_lum_idx = np.where(np.logical_and( sf_log10_Lrad >= lum_bins[i-1], sf_log10_Lrad < lum_bins[i] ) )[0]
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
    rhos.append(rho)
    agn_rhos.append(agn_rho)
    sf_rhos.append(sf_rho)

lum_func = np.asarray(rhos)
agn_lum_func = np.asarray(agn_rhos)
sf_lum_func = np.asarray(sf_rhos)

lum_bin_cens = lum_bins + 0.5*(lum_bins[1]-lum_bins[0])

fig = plt.figure( figsize=(5,5) )
## plot
#plt.plot( ms_144MHz, mauch_sadler['log10RLF_all'], color='black', linewidth=2.5, label='MS07 All' )
#plt.plot( ms_144MHz, mauch_sadler['log10RLF_SF'], color='magenta', linewidth=2.5, label='MS07 SFG' )
#plt.plot( ms_144MHz, mauch_sadler['log10RLF_RLAGN'], color='orange', linewidth=2.5, label='MS07 AGN' )
plt.plot( cochrane['logL150'], cochrane['logPhi'], color='blue', label='Cochrane et al. 2023, SFGs')
plt.plot( kondapally['logL150'], kondapally['logPhi'], color='red', label='Kondapally et al. 2022, RLAGNs')
## plot the lofar data, filtering zeros
non_zero = np.where( lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], lum_func[non_zero], 'o', color='red', label='data_corr' )
non_zero = np.where( agn_lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], agn_lum_func[non_zero], 'o', color='purple', label='AGN' )
non_zero = np.where( sf_lum_func != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], sf_lum_func[non_zero], 'o', color='pink', label='SF' )

#plt.fill_between( RLF['Lmedian'], RLF['RLF_lo'], RLF['RLF_up'], color='green', alpha=0.5 )
#plt.plot( RLF['Lmedian'], RLF['RLF'], 'o', color='green', label='data' )
#plt.plot( RLF_corr['Lmedian'], RLF_corr['RLF'], 'o', color='red', label='data_corr' )
plt.xlim((20,28))
plt.ylim(-7.5,-1)
plt.xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'mauch_sadler_RLFs.png',dpi=300)
fig.clear()
plt.close()



