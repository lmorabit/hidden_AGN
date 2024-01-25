import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
## WMAP9 is Hinshaw et al. 2013, H_0=69.3, Omega=0.287

si = -0.8

## read in SIMBA
Simba_SF = Table.read( paths.data / 'RLFS_50MYR_SF.csv', format='csv', delimiter=',' )
Simba_AGN = Table.read( paths.data / 'RLFS_50MYR_AGN.csv', format='csv', delimiter=',' )
## shift using spectral index
Simba_SF['x']  = Simba_SF['x'] + np.log10( np.power( (144./1400.), si ) )
Simba_AGN['x']  = Simba_AGN['x'] + np.log10( np.power( (144./1400.), si ) )

fig = plt.figure( figsize=(5,5) )
## plot simba
plt.plot( Simba_SF['x'], Simba_SF['Curve1'], color='red', linewidth=2.5, label='Simba SF' )
plt.plot( Simba_AGN['x'], Simba_AGN['Curve2'], color='blue', linewidth=2.5, label='Simba AGN' )
plt.xlim((20,28))
plt.ylim(-7.5,-2)
plt.xlabel('log('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
plt.ylabel('log('+r'$\Phi$'+' [mag'+r'$^{-1}$'+' Mpc'+r'$^{-3}$'+'])')
plt.legend()
plt.savefig(paths.figures / 'test.png',dpi=300)
fig.clear()



