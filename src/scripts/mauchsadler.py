import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

si = -0.8

## read in SIMBA
mauch_sadler = Table.read( paths.static / 'mauch_sadler_table5.csv', format='csv', delimiter=',' )
## shift using spectral index
Simba_SF['x']  = Simba_SF['x'] + np.log10( np.power( (144./1400.), si ) )
ms_144MHz = mauch_sadler['log10_P1p4GHz'] + np.log10( np.power( (144./1400.), si ) ) 


fig = plt.figure( figsize=(5,5) )
## plot
plt.plot( ms_144MHz, mauch_sadler['log10RLF_all'], color='black', linewidth=2.5, label='All' )
plt.plot( ms_144MHz, mauch_sadler['log10RLF_SF'], color='magenta', linewidth=2.5, label='SFG' )
plt.plot( ms_144MHz, mauch_sadler['log10RLF_AGN'], color='orange', linewidth=2.5, label='AGN' )
plt.xlim((20,28))
plt.ylim(-7.5,-2)
plt.xlabel('log('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log('+r'$\Phi$'+' [mag'+r'$^{-1}$'+' Mpc'+r'$^{-3}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'mauch_sadler_RLFs.png',dpi=300)
fig.clear()



