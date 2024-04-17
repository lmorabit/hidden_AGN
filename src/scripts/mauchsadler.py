import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
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
si = -0.8
zmin = 0.003  ## matches Mauch & Sadler 2007
zmax = 0.3   ## matches Mauch & Sadler 2007
dz = 0.0001  ## matches Cochrane and Kondapally

## read in Mauch & Sadler Table 5
mauch_sadler = Table.read( paths.static / 'mauch_sadler_table5.csv', format='csv', delimiter=',' )
## shift using spectral index
ms_144MHz = mauch_sadler['log10_P1p4GHz'] + np.log10( np.power( (144./1400.), si ) ) 

## lotss data
cochrane = Table.read( paths.static / 'cochrane_2023_table1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table2.csv', format='csv', delimiter=',' )

## read in vmaxes
vmaxes = Table.read( paths.data / 'lockman_6arcsec_vmaxes.fits', format='fits' )


## calculate RLFs
redshift_bins = np.array([zmin,zmax])
lum_bins = np.arange( 19.5, 27, 0.5 ) + np.log10( np.power( (144./1400.), si ) )
Lrad = radio_power( vmaxes['Total_flux'], vmaxes['Z_BEST'], spectral_index=si )
log10_Lrad = np.log10(Lrad)

phis = []
for i in np.arange(1,len(lum_bins)):
    delta_log_L = lum_bins[i] - lum_bins[i-1]
    lum_idx = np.logical_and( log10_Lrad >= lum_bins[i-1], log10_Lrad < lum_bins[i] )
    phi = ( 1. / np.power( 10., delta_log_L ) ) * np.sum( 1. / vmaxes['vmax'][lum_idx] )
    phis.append(phi)

## convert back to "mag" 
phis = np.asarray(phis)/1e7

lum_bin_cens = lum_bins + 0.5*(lum_bins[1]-lum_bins[0])



fig = plt.figure( figsize=(5,5) )
## plot
#plt.plot( ms_144MHz, mauch_sadler['log10RLF_all'], color='black', linewidth=2.5, label='MS07 All' )
#plt.plot( ms_144MHz, mauch_sadler['log10RLF_SF'], color='magenta', linewidth=2.5, label='MS07 SFG' )
#plt.plot( ms_144MHz, mauch_sadler['log10RLF_RLAGN'], color='orange', linewidth=2.5, label='MS07 AGN' )
plt.plot( cochrane['logL150'], cochrane['logPhi'], color='blue', label='Cochrane et al. 2023, SFGs')
plt.plot( kondapally['logL150'], kondapally['logPhi'], color='red', label='Kondapally et al. 2022, RLAGNs')
## plot the lofar data, filtering zeros
non_zero = np.where( phis != 0.0 )[0]
plt.plot( lum_bin_cens[non_zero], np.log10(phis[non_zero]), 'o', color='red', label='data_corr' )
#plt.fill_between( RLF['Lmedian'], RLF['RLF_lo'], RLF['RLF_up'], color='green', alpha=0.5 )
#plt.plot( RLF['Lmedian'], RLF['RLF'], 'o', color='green', label='data' )
#plt.plot( RLF_corr['Lmedian'], RLF_corr['RLF'], 'o', color='red', label='data_corr' )
plt.xlim((20,28))
plt.ylim(-7.5,-2)
plt.xlabel('log('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log('+r'$\Phi$'+' [mag'+r'$^{-1}$'+' Mpc'+r'$^{-3}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'mauch_sadler_RLFs.png',dpi=300)
fig.clear()



