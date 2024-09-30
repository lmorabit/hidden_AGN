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
lmin = 20.5
lmax = 27
dl = 0.3


## lotss data
cochrane = Table.read( paths.static / 'cochrane_2023_table1.csv', format='csv', delimiter=',' )
kondapally = Table.read( paths.static / 'kondapally_2022_table2.csv', format='csv', delimiter=',' )

## read in vmaxes
lockman_vmaxes = Table.read( paths.static / 'lockman_vmaxes_zmin0.003_zmax0.3.fits', format='fits' )
elais_vmaxes = Table.read( paths.static / 'en1_vmaxes_zmin0.003_zmax0.3.fits', format='fits' )

## combine these into a single catalogue for plotting
keep_cols = ['Total_flux_dr','Z_BEST','vmax','agn_vmax','sf_vmax','AGN_flux','SF_flux', 'Overall_class','Mass_cons','SFR_cons']
lockman_tmp = lockman_vmaxes[keep_cols]
elais_tmp = elais_vmaxes[keep_cols]
vmaxes = vstack([lockman_tmp,elais_tmp])

lum_bins, lum_func, agn_lum_func, sf_lum_func, gal_agn_lum_func, gal_sf_lum_func = get_RLFs( vmaxes, zmin, zmax, lmin=lmin, lmax=lmax, dl=dl, si=si )

e_agn_lum_func, e_sf_lum_func, e_gal_agn_lum_func, e_gal_sf_lum_func = random_resample( agn_lum_func, sf_lum_func, gal_agn_lum_func, gal_sf_lum_func, vmaxes, zmin, zmax, lmin=lmin, lmax=lmax, dl=dl, si=si, nsamp=1000 )

'''
## randomly re-sample 1000 times
print('randomly resampling to get uncertainties')
nsamp = 1000
agn_lfs = np.zeros((len(lum_func),nsamp))
sf_lfs = np.zeros((len(lum_func),nsamp))
g_agn_lfs = np.zeros((len(lum_func),nsamp))
g_sf_lfs = np.zeros((len(lum_func),nsamp))
fracs = np.zeros(nsamp)

for i in np.arange(0,nsamp):
    np.random.seed(i)
    idx = np.random.randint(low=0,high=len(vmaxes),size=len(vmaxes))
    fracs[i] = float(len(np.unique(idx))) / float(len(vmaxes))
    lb, lf, agn_lf, sf_lf, g_agn_lf, g_sf_lf = get_RLFs( vmaxes[idx], zmin, zmax, lmin=lmin, lmax=lmax, dl=dl, si=si )
    agn_lfs[:,i] = agn_lf
    sf_lfs[:,i] = sf_lf
    g_agn_lfs[:,i] = g_agn_lf
    g_sf_lfs[:,i] = g_sf_lf


## as I'm sampling with replacement, need to scale by the size of the resulting random sample
tmp = np.tile( agn_lum_func, (nsamp,1) ).transpose()
e_agn_lum_func = np.std( np.sqrt(fracs)*(agn_lfs-tmp),axis=1)
tmp = np.tile( sf_lum_func, (nsamp,1) ).transpose()
e_sf_lum_func = np.std( np.sqrt(fracs)*(sf_lfs-tmp),axis=1)
tmp = np.tile( gal_agn_lum_func, (nsamp,1) ).transpose()
e_gal_agn_lum_func = np.std( np.sqrt(fracs)*(g_agn_lfs-tmp), axis=1 )
tmp = np.tile( gal_sf_lum_func, (nsamp,1) ).transpose()
e_gal_sf_lum_func = np.std(np.sqrt(fracs)*(g_sf_lfs-tmp),axis=1)

## replace small values with 0.03
e_agn_lum_func[np.where(e_agn_lum_func < 0.03)[0]] = 0.03
e_sf_lum_func[np.where(e_sf_lum_func < 0.03)[0]] = 0.03
e_gal_agn_lum_func[np.where(e_gal_agn_lum_func < 0.03)[0]] = 0.03
e_gal_sf_lum_func[np.where(e_gal_sf_lum_func < 0.03)[0]] = 0.03
'''


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
x, y, dy, idx1, idx2 = get_values( lum_bin_cens, gal_agn_lum_func, e_gal_agn_lum_func )
p1.plot( x, y, color=agngalc, label='AGN galaxies', linewidth=3, linestyle='dotted' )
p1.fill_between(x[idx1], y[idx1]-dy[idx1], y[idx1]+dy[idx1], color=agngalc, alpha=0.4, ec=None)
if len(idx2) > 0:
    for idx22 in idx2:
        p1.fill_between(x[idx22], y[idx22]-dy[idx22], y[idx22]+dy[idx22], color=agngalc, alpha=0.1, hatch='xxx', ec=None)
x, y, dy, idx1, idx2 = get_values( lum_bin_cens, gal_sf_lum_func, e_gal_sf_lum_func )
p1.plot( x, y, color=sfggalc, label='SF galaxies', linewidth=3, linestyle='dotted' )
p1.fill_between( x[idx1], y[idx1]-dy[idx1], y[idx1]+dy[idx1] , color=sfggalc, alpha=0.4, ec=None)
if len(idx2) > 0:
    for idx22 in idx2:
        p1.fill_between( x[idx22], y[idx22]-dy[idx22], y[idx22]+dy[idx22] , color=sfggalc, alpha=0.1, hatch='xxx', ec=None)
p1.text( 20.5, -6.2, '0.003 < z < 0.3', fontsize=14 )
p1.axes.set_xlim(plxlims)
p1.axes.set_ylim(plylims)
p1.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p1.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
handles, labels = p1.get_legend_handles_labels()
order = [3,2,0,1]
p1.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

## Middle panel: same as left panel but new data by process now
p2 = plt.axes([0.05+sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,sbsizey])
## plot the previous data
p2.plot( cochrane['logL150'], cochrane['logPhi'], color=coch, label='Cochrane et al. (2023)', linewidth=2)
p2.plot( kondapally['logL150'], kondapally['logPhi'], color=kond, label='Kondapally et al. (2022)', linewidth=2)
## plot the lofar data, filtering zeros
x, y, dy, idx1, idx2 = get_values( lum_bin_cens, agn_lum_func, e_agn_lum_func )
p2.plot( x, y, color=agnc, label='AGN process', linewidth=3, alpha=0.75 )
p2.fill_between( x[idx1], y[idx1]-dy[idx1], y[idx1]+dy[idx1], color=agnc, alpha=0.4, ec=None)
if len(idx2) > 0:
    for idx22 in idx2:
        p2.fill_between( x[idx22], y[idx22]-dy[idx22], y[idx22]+dy[idx22], color=agnc, alpha=0.1, hatch='xxx', ec=None)

x, y, dy, idx1, idx2 = get_values( lum_bin_cens, sf_lum_func, e_sf_lum_func )
p2.plot( x, y, color=sfc, label='SF process', linewidth=3, alpha=0.75 )
p2.fill_between( x[idx1], y[idx1]-dy[idx1], y[idx1]+dy[idx1], color=sfc, alpha=0.4, ec=None)
if len(idx2) > 0:
    for idx22 in idx2:
        p2.fill_between( x[idx22], y[idx22]-dy[idx22], y[idx22]+dy[idx22], color=sfc, alpha=0.1, hatch='xxx', ec=None)
p2.text( 20.5, -6.2, '0.003 < z < 0.3', fontsize=14 )
p2.axes.set_xlim(plxlims)
p2.axes.set_ylim(plylims)
p2.yaxis.set_visible(False)
p2.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
#p2.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
handles, labels = p2.get_legend_handles_labels()
order = [3,2,0,1]
p2.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

## Right panel (top): RLFs using new data, galaxies and process
p3 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.42,sbsizex*fsizey/fsizex,0.6*sbsizey])
non_zero = np.where( gal_agn_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], gal_agn_lum_func[non_zero], color=agngalc, label='AGN galaxies', linewidth=3, alpha=0.75, linestyle='dotted' )
non_zero = np.where( gal_sf_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], gal_sf_lum_func[non_zero], color=sfggalc, label='SF galaxies', linewidth=3, alpha=0.75, linestyle='dotted' )
non_zero = np.where( agn_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], agn_lum_func[non_zero], color=agnc, label='AGN process', linewidth=3 )
non_zero = np.where( sf_lum_func != 0.0 )[0]
p3.plot( lum_bin_cens[non_zero], sf_lum_func[non_zero], color=sfc, label='SF process', linewidth=3 )
p3.axes.set_xlim(plxlims)
p3.axes.set_ylim(plylims)
p3.xaxis.set_visible(False)
p3.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p3.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')
handles, labels = p3.get_legend_handles_labels()
order = [3,1,2,0]
p3.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

## Right panel (bottom): ratio of the RLFs by galaxies and process
p4 = plt.axes([0.12+2*sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,0.4*sbsizey])
p4.plot( (19,27), (1,1), color='gray', linestyle='dashed', linewidth=1.5 )
### SF
non_zero_idx = np.where( np.logical_and( sf_lum_func != 0, gal_sf_lum_func != 0 ) )[0]
ratio = np.power(10., sf_lum_func[non_zero_idx] ) / np.power( 10., gal_sf_lum_func[non_zero_idx] )
e_ratio = np.sqrt( np.power( e_sf_lum_func[non_zero_idx], 2. ) + np.power( e_gal_sf_lum_func[non_zero_idx], 2. ) )
x, y, dy, idx1, idx2 = get_values( lum_bin_cens[non_zero_idx], ratio, e_ratio, useidx=False )
p4.plot( x, y, color=sfc, label='SF', linewidth=3 )
p4.fill_between( x[idx1], y[idx1]-dy[idx1], y[idx1]+dy[idx1], alpha=0.4, color=sfc, ec=None )
if len(idx2) > 0:
    for idx22 in idx2:
        p4.fill_between( x[idx22], y[idx22]-dy[idx22], y[idx22]+dy[idx22], alpha=0.1, hatch='xxx', color=sfc, ec=None )
### AGN
non_zero_idx = np.where( np.logical_and( agn_lum_func != 0, gal_agn_lum_func != 0 ) )[0]
ratio = np.power(10., agn_lum_func[non_zero_idx]) / np.power( 10., gal_agn_lum_func[non_zero_idx]) 
e_ratio = np.sqrt( np.power( e_agn_lum_func[non_zero_idx], 2. ) + np.power( e_gal_agn_lum_func[non_zero_idx], 2. ) )
x, y, dy, idx1, idx2 = get_values( lum_bin_cens[non_zero_idx], ratio, e_ratio, useidx=False )
p4.plot( x, y, color=agnc, label='AGN', linewidth=3 )
p4.fill_between( x[idx1], y[idx1]-dy[idx1], y[idx1]+dy[idx1], alpha=0.4, color=agnc, ec=None )
if len(idx2) > 0:
    for idx22 in idx2:
        p4.fill_between( x[idx22], y[idx22]-dy[idx22], y[idx22]+dy[idx22], alpha=0.1, hatch='xxx', color=agnc, ec=None )
p4.axes.set_xlim(plxlims)
p4.axes.set_ylim((0.45,1.9))
p4.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p4.set_ylabel(r'$\Delta$RLF')
p4.legend()
fig.savefig(paths.figures / 'deep_fields_RLFs.png',dpi=300)
fig.clear()
plt.close()

