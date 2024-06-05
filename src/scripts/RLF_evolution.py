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


t = Table.read( paths.static / 'redshift_bins.csv', format='csv' )
zbin_starts = t['zbin_starts']
zbin_ends = t['zbin_ends']

fields = ['lockman'] # ,'en1']

z_lum_bins = []
z_lum_func = []
z_agn_lum_func = []
z_sf_lum_func = []
z_gal_agn_lum_func = []
z_gal_sf_lum_func = []
for i in np.arange(0,5):
    zmin = zbin_starts[i]
    zmax = zbin_ends[i]
    keep_cols = ['Total_flux_dr','Z_BEST','vmax','agn_vmax','sf_vmax','AGN_flux','SF_flux', 'Overall_class','Mass_cons','SFR_cons']
    vmaxes = Table()
    for field in fields:
        infile = '{:s}_vmaxes_zmin{:s}_zmax{:s}.fits'.format(field,str(zmin),str(zmax))
        tmp = Table.read( paths.data / infile, format='fits' )
        vmaxes = vstack([vmaxes,tmp[keep_cols]])
    lum_bins, lum_func, agn_lum_func, sf_lum_func, gal_agn_lum_func, gal_sf_lum_func = get_RLFs( vmaxes, zmin, zmax, lmin=20.5, lmax=27, dl=0.3, si=si )
    z_lum_bins.append(lum_bins)
    z_lum_func.append(lum_func)
    z_agn_lum_func.append(agn_lum_func)
    z_sf_lum_func.append(sf_lum_func)
    z_gal_agn_lum_func.append(gal_agn_lum_func)
    z_gal_sf_lum_func.append(gal_sf_lum_func)


lum_bin_cens = z_lum_bins[0][0:-1] + 0.5*(z_lum_bins[0][1]-z_lum_bins[0][0])


fsizex = 14
fsizey = 5
sbsizex = 0.8
sbsizey = 0.8
plxlims = (20.1,27)
plylims = (-7.5,-1)


##########################################################
## Star formation

zcols_sf = mycols[np.arange(0,len(z_lum_bins))*int(n/len(z_lum_bins))]

fig = plt.figure( figsize=(fsizex,fsizey) )
## Left panel: Galxies
p1 = plt.axes([0.05,0.1,sbsizex*fsizey/fsizex,sbsizey])
for i in np.arange(0,len(z_lum_bins)):
    non_zero = np.where( np.logical_and( z_gal_sf_lum_func[i] != 0.0, np.isfinite(z_gal_sf_lum_func[i]) ) )[0]
    p1.plot( lum_bin_cens[non_zero], z_gal_sf_lum_func[i][non_zero], color=zcols_sf[i], linewidth=2,linestyle='dotted' )
p1.set_title('SF Galaxies')
p1.axes.set_xlim(plxlims)
p1.axes.set_ylim(plylims)
p1.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p1.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Middle panel: Activity
p2 = plt.axes([0.05+sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,sbsizey])
## plot the lofar data, filtering zeros
for i in np.arange(0,len(z_lum_bins)):
    non_zero = np.where( z_sf_lum_func[i] != 0.0 )[0]
    p2.plot( lum_bin_cens[non_zero], z_sf_lum_func[i][non_zero], color=zcols_sf[i], label='{:s} < z < {:s}'.format(str(zbin_starts[i]),str(zbin_ends[i])), linewidth=2, alpha=0.75 )
p2.legend()
p2.axes.set_xlim(plxlims)
p2.axes.set_ylim(plylims)
p2.yaxis.set_visible(False)
p2.set_title('SF Activity')
p2.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
#p2.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Right panel (top): galaxies and activity together
p3 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.42,sbsizex*fsizey/fsizex,0.6*sbsizey])
for i in np.arange(0,len(z_lum_bins)):
    non_zero = np.where( z_gal_sf_lum_func[i] != 0.0 )[0]
    p3.plot( lum_bin_cens[non_zero], z_gal_sf_lum_func[i][non_zero], color=zcols_sf[i], linewidth=3, alpha=0.75, linestyle='dotted' )
    non_zero = np.where( z_sf_lum_func[i] != 0.0 )[0]
    p3.plot( lum_bin_cens[non_zero], z_sf_lum_func[i][non_zero], color=zcols_sf[i], label='SF activity', linewidth=3 )
p3.axes.set_xlim(plxlims)
p3.axes.set_ylim(plylims)
p3.xaxis.set_visible(False)
p3.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p3.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Right panel (bottom): ratio of the RLFs by galaxies and activity
p4 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,0.4*sbsizey])
p4.plot( (19,27), (1,1), color='gray', linestyle='dashed', linewidth=1.5 )
for i in np.arange(0,len(z_lum_bins)):
    non_zero_idx = np.where( np.logical_and( z_sf_lum_func[i] != 0, z_gal_sf_lum_func[i] != 0 ) )[0]
    ratio = np.power(10., z_sf_lum_func[i][non_zero_idx] ) / np.power( 10., z_gal_sf_lum_func[i][non_zero_idx] )
    p4.plot( lum_bin_cens[non_zero_idx], ratio, color=zcols_sf[i], linewidth=3 )
p4.axes.set_xlim(plxlims)
p4.axes.set_ylim((0.5,1.1))
p4.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p4.set_ylabel('Activity / Galaxy')
fig.savefig(paths.figures / 'RLF_evolution_SF.png',dpi=300)
fig.clear()
plt.close()

##########################################################
## AGN

zcols_agn = mycols_m[np.arange(0,len(z_lum_bins))*int(n/len(z_lum_bins))]

fig = plt.figure( figsize=(fsizex,fsizey) )
## Left panel: Galxies
p1 = plt.axes([0.05,0.1,sbsizex*fsizey/fsizex,sbsizey])
for i in np.arange(0,len(z_lum_bins)):
    non_zero = np.where( np.logical_and( z_gal_agn_lum_func[i] != 0.0, np.isfinite(z_gal_agn_lum_func[i]) ) )[0]
    p1.plot( lum_bin_cens[non_zero], z_gal_agn_lum_func[i][non_zero], color=zcols_agn[i], linewidth=2,linestyle='dotted' )
p1.set_title('AGN Galaxies')
p1.axes.set_xlim(plxlims)
p1.axes.set_ylim(plylims)
p1.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p1.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Middle panel: Activity
p2 = plt.axes([0.05+sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,sbsizey])
## plot the lofar data, filtering zeros
for i in np.arange(0,len(z_lum_bins)):
    non_zero = np.where( z_agn_lum_func[i] != 0.0 )[0]
    p2.plot( lum_bin_cens[non_zero], z_agn_lum_func[i][non_zero], color=zcols_agn[i], label='{:s} < z < {:s}'.format(str(zbin_starts[i]),str(zbin_ends[i])), linewidth=2, alpha=0.75 )
p2.legend()
p2.axes.set_xlim(plxlims)
p2.axes.set_ylim(plylims)
p2.yaxis.set_visible(False)
p2.set_title('AGN Activity')
p2.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
#p2.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Right panel (top): galaxies and activity together
p3 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.42,sbsizex*fsizey/fsizex,0.6*sbsizey])
for i in np.arange(0,len(z_lum_bins)):
    non_zero = np.where( z_gal_agn_lum_func[i] != 0.0 )[0]
    p3.plot( lum_bin_cens[non_zero], z_gal_agn_lum_func[i][non_zero], color=zcols_agn[i], linewidth=3, alpha=0.75, linestyle='dotted' )
    non_zero = np.where( z_agn_lum_func[i] != 0.0 )[0]
    p3.plot( lum_bin_cens[non_zero], z_agn_lum_func[i][non_zero], color=zcols_agn[i], label='agn activity', linewidth=3 )
p3.axes.set_xlim(plxlims)
p3.axes.set_ylim(plylims)
p3.xaxis.set_visible(False)
p3.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p3.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Right panel (bottom): ratio of the RLFs by galaxies and activity
p4 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,0.4*sbsizey])
p4.plot( (19,27), (1,1), color='gray', linestyle='dashed', linewidth=1.5 )
for i in np.arange(0,len(z_lum_bins)):
    non_zero_idx = np.where( np.logical_and( z_agn_lum_func[i] != 0, z_gal_agn_lum_func[i] != 0 ) )[0]
    ratio = np.power(10., z_agn_lum_func[i][non_zero_idx] ) / np.power( 10., z_gal_agn_lum_func[i][non_zero_idx] )
    p4.plot( lum_bin_cens[non_zero_idx], ratio, color=zcols_agn[i], linewidth=3 )
p4.axes.set_xlim(plxlims)
p4.axes.set_ylim((0.7,2.5))
p4.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' W Hz'+r'$^{-1}$'+'])')
p4.set_ylabel('Activity / Galaxy')
fig.savefig(paths.figures / 'RLF_evolution_AGN.png',dpi=300)
fig.clear()
plt.close()






'''
Redshift bins

Rohit (AGN)
0.5 - 1.0
1.0 - 1.5
1.5 - 2.0
2.0 - 2.5


Rachael (SFGs)
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

