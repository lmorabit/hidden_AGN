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

fig = plt.figure( figsize=(9,5) )
gs = fig.add_gridspec(1,2,hspace=0,wspace=0)
axs = gs.subplots(sharex=True,sharey=True)
## plot Lockman
axs[0].plot( cochrane['FluxDensity_mJy'], cochrane['Lockman'], color=mycols_m[250], linewidth=2.5, label='SFGs, Cochrane' )
axs[0].plot( kondapally['FluxDensity_mJy'], kondapally['Lockman'], color=mycols[20], linewidth=2.5, label='RLAGN, Kondapally' )
axs[0].set_xlabel(r'$S_{\mathrm{144 MHz}}$'+' mJy')
axs[0].set_ylabel('Completeness')
axs[0].set_ylim((0.8,1.0))
axs[0].legend()
axs[0].text(0.05,0.9,'Lockman',transform=axs[0].transAxes)
## plot Elais
axs[1].plot( cochrane['FluxDensity_mJy'], cochrane['Elais'], color=mycols_m[250], linewidth=2.5, label='SFGs, Cochrane' )
axs[1].plot( kondapally['FluxDensity_mJy'], kondapally['Elais'], color=mycols[20], linewidth=2.5, label='RLAGN, Kondapally' )
axs[1].set_xlabel(r'$S_{\mathrm{144 MHz}}$'+' mJy')
#axs[1].set_ylabel('Completeness')
axs[1].set_ylim((0.8,1.0))
axs[1].legend()
axs[1].text(0.05,0.9,'ELAIS',transform=axs[1].transAxes)
fig.tight_layout()
fig.savefig(paths.figures / 'completeness.png',dpi=300)
fig.clear()



## write out the average completeness ... 
## could weight by fraction of SFG to AGN from Philip's paper
