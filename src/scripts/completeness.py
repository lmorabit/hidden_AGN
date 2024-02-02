import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
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


## read in completeness
cochrane = Table.read( paths.static / 'cochrane_2023_tableA1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table1.csv', format='csv', delimiter=',' )

fig = plt.figure( figsize=(5,5) )
## plot simba
plt.plot( cochrane['FluxDensity_mJy'], cochrane['Lockman'], color=mycols_m[50], linewidth=2.5, label='SFGs, Cochrane' )
plt.plot( kondapally['FluxDensity_mJy'], kondapally['Lockman'], color=mycols[20], linewidth=2.5, label='RLAGN, Kondapally' )
plt.xlabel(r'$S_{\mathrm{144 MHz}}$'+' mJy')
plt.ylabel('Completeness')
plt.ylim((0.8,1.0))
plt.legend()
plt.savefig(paths.figures / 'completeness.png',dpi=300)
fig.clear()



