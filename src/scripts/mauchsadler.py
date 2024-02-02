import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

## set a default spectral index
si = -0.8

## read in Mauch & Sadler Table 5
mauch_sadler = Table.read( paths.static / 'mauch_sadler_table5.csv', format='csv', delimiter=',' )
## shift using spectral index
ms_144MHz = mauch_sadler['log10_P1p4GHz'] + np.log10( np.power( (144./1400.), si ) ) 

## read in RLF
RLF = Table.read( paths.data / 'RLF.fits' )

fig = plt.figure( figsize=(5,5) )
## plot
plt.plot( ms_144MHz, mauch_sadler['log10RLF_all'], color='black', linewidth=2.5, label='All' )
plt.plot( ms_144MHz, mauch_sadler['log10RLF_SF'], color='magenta', linewidth=2.5, label='SFG' )
plt.plot( ms_144MHz, mauch_sadler['log10RLF_RLAGN'], color='orange', linewidth=2.5, label='AGN' )
## plot the lofar data, filtering zeros
non_zero = np.where( RLF['RLF'] != 0.0 )[0]
RLF = RLF[non_zero]
plt.fill_between( RLF['Lmedian'], RLF['RLF_lo'], RLF['RLF_up'], color='green', alpha=0.5 )
plt.plot( RLF['Lmedian'], RLF['RLF'], 'o', color='green', label='data' )
plt.xlim((20,28))
plt.ylim(-7.5,-2)
plt.xlabel('log('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log('+r'$\Phi$'+' [mag'+r'$^{-1}$'+' Mpc'+r'$^{-3}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'mauch_sadler_RLFs.png',dpi=300)
fig.clear()



