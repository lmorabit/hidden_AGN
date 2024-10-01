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

## read in the data 
t = Table.read( paths.static / 'redshift_bins.csv', format='csv' )
zbin_starts = t['zbin_starts']
zbin_ends = t['zbin_ends']

fields = ['lockman','en1']

z_lum_bins = []
z_lum_func = []
z_agn_lum_func = []
z_sf_lum_func = []
z_gal_agn_lum_func = []
z_gal_sf_lum_func = []

e_z_agn_lum_func = []
e_z_sf_lum_func = []
e_z_gal_agn_lum_func = []
e_z_gal_sf_lum_func = []


for i in np.arange(0,5):
    zmin = zbin_starts[i]
    zmax = zbin_ends[i]
    keep_cols = ['Total_flux_dr','Z_BEST','vmax','agn_vmax','sf_vmax','AGN_flux','SF_flux', 'Overall_class','Mass_cons','SFR_cons']
    vmaxes = Table()
    for field in fields:
        infile = '{:s}_vmaxes_zmin{:s}_zmax{:s}.fits'.format(field,str(zmin),str(zmax))
        tmp = Table.read( paths.static / infile, format='fits' )
        vmaxes = vstack([vmaxes,tmp[keep_cols]],metadata_conflicts='silent')
    lum_bins, lum_func, agn_lum_func, sf_lum_func, gal_agn_lum_func, gal_sf_lum_func = get_RLFs( vmaxes, zmin, zmax, lmin=lmin, lmax=lmax, dl=dl, si=si )
    
    e_agn_lum_func, e_sf_lum_func, e_gal_agn_lum_func, e_gal_sf_lum_func = random_resample( agn_lum_func, sf_lum_func, gal_agn_lum_func, gal_sf_lum_func, vmaxes, zmin, zmax, lmin=lmin, lmax=lmax, dl=dl, si=si, nsamp=1000 )

    ## filter out unconstrained bins - AGN
    filter_idx = np.where( np.abs(agn_lum_func - lum_func) > 4. )[0]
    agn_lum_func[filter_idx] = 0. 
    e_agn_lum_func[filter_idx] = 0. 
    filter_idx = np.where( np.abs(gal_agn_lum_func - lum_func) > 4. )[0]
    gal_agn_lum_func[filter_idx] = 0.
    e_gal_agn_lum_func[filter_idx] = 0.  
    ## filter out unconstrained bins - SF
    filter_idx = np.where( np.logical_or( np.abs(sf_lum_func - lum_func) > 4., sf_lum_func > -1.8 ) )[0]
    sf_lum_func[filter_idx] = 0.
    e_sf_lum_func[filter_idx] = 0.  
    filter_idx = np.where( np.logical_or( np.abs(gal_sf_lum_func - lum_func) > 4., gal_sf_lum_func > -1.8 ) )[0]
    gal_sf_lum_func[filter_idx] = 0. 
    e_gal_sf_lum_func[filter_idx] = 0. 
    ## add to redshift list
    z_lum_bins.append(lum_bins)
    z_lum_func.append(lum_func)
    z_agn_lum_func.append(agn_lum_func)
    z_sf_lum_func.append(sf_lum_func)
    z_gal_agn_lum_func.append(gal_agn_lum_func)
    z_gal_sf_lum_func.append(gal_sf_lum_func)

    e_z_agn_lum_func.append(e_agn_lum_func)
    e_z_sf_lum_func.append(e_sf_lum_func)
    e_z_gal_agn_lum_func.append(e_gal_agn_lum_func)
    e_z_gal_sf_lum_func.append(e_gal_sf_lum_func)


## get the luminosity bin centres
lum_bin_cens = z_lum_bins[0][0:-1] + 0.5*(z_lum_bins[0][1]-z_lum_bins[0][0])

##########################################################
## Redshift evolution plot


## figure configuration
fsizex = 10 
fsizey = 5
sbsizex = 0.8
sbsizey = 0.8
plxlims = (20.1,27)
plylims = (-7.5,-1)

## colours
zcols_sf = mycols[np.arange(0,len(z_lum_bins))*int(n/len(z_lum_bins))]
zcols_agn = mycols_m[np.arange(0,len(z_lum_bins))*int(n/len(z_lum_bins))]

## start the figure
fig = plt.figure( figsize=(fsizex,fsizey) )

print('SF')
sf_delta_int = []
e_sf_delta_int = []

## STAR FORMATION
p1 = plt.axes([0.07,0.42,sbsizex*fsizey/fsizex,0.6*sbsizey])
for i in np.arange(0,len(z_lum_bins)):
    ## Galaxies
    x, y, dy, idx1, idx2 = get_values( lum_bin_cens, z_gal_sf_lum_func[i], e_z_gal_sf_lum_func[i] )
    p1.plot( x, y, color=zcols_sf[i], linewidth=3, alpha=0.75, linestyle='dotted' )
    # use rectangular integration 
    gal_trapz = np.sum( np.power( 10., y ) * dl )
    e_gal_trapz = np.sqrt( np.sum( np.power( dy * np.log( 10. ) * np.power( 10., y ), 2. ) ) ) * dl 
    ## Process
    x, y, dy, idx1, idx2 = get_values( lum_bin_cens, z_sf_lum_func[i], e_z_sf_lum_func[i] )
    p1.plot( x, y, color=zcols_sf[i], label='{:s} < z < {:s}'.format(str(zbin_starts[i]),str(zbin_ends[i])), linewidth=3 )
    # use rectangular integration
    trapz = np.sum( np.power( 10., y ) * dl )
    e_trapz = np.sqrt( np.sum( np.power( dy * np.log( 10. ) * np.power( 10., y ), 2. ) ) ) * dl
    ## ratio and errors
    sf_delta_int.append( trapz / gal_trapz )
    e_sf_delta_int.append( mult_div_error( trapz/gal_trapz, np.asarray([trapz,gal_trapz]), np.asarray([e_trapz,e_gal_trapz]) ) )
    
l1 = p1.legend()
dline = matplotlib.lines.Line2D([0],[0], color='black', linewidth=3 )
lline = matplotlib.lines.Line2D([0],[0], color='black', linewidth=3, linestyle='dotted' )
l2 = p1.legend((lline,dline),('Galaxies', 'Process'),loc='upper right')
p1.add_artist(l1)
p1.axes.set_xlim(plxlims)
p1.axes.set_ylim(plylims)
p1.xaxis.set_visible(False)
p1.set_title('Star Formation',fontsize=20)
p1.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p1.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Bottom panel: ratio of the RLFs by galaxies and process
p2 = plt.axes([0.07,0.1,sbsizex*fsizey/fsizex,0.4*sbsizey])
p2.plot( (19,27), (1,1), color='gray', linestyle='dashed', linewidth=1.5 )
p2.plot( (19,27), (0.5,0.5), color='gray', linewidth=1, alpha=0.25 )
for i in np.arange(0,len(z_lum_bins)):
    non_zero_idx = np.where( np.logical_and( z_sf_lum_func[i] < 0, z_gal_sf_lum_func[i] < 0 ) )[0]
    ratio = np.power(10., z_sf_lum_func[i][non_zero_idx] ) / np.power( 10., z_gal_sf_lum_func[i][non_zero_idx] )
#    e_ratio = 
#    x, y, dy, idx1, idx2 = get_values( lum_bin_cens


    p2.plot( lum_bin_cens[non_zero_idx], ratio, color=zcols_sf[i], linewidth=3 )    
p2.axes.set_xlim(plxlims)
p2.axes.set_ylim((0.4,1.1))
p2.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p2.set_ylabel(r'$\Delta$'+'RLF')

## integrate the ratio

print('AGN')
agn_delta_int = []
e_agn_delta_int = []

## ACTIVE GALACTIC NUCLEI
p3 = plt.axes([0.14+sbsizex*fsizey/fsizex,0.42,sbsizex*fsizey/fsizex,0.6*sbsizey])
for i in np.arange(0,len(z_lum_bins)):
    ## Galaxies
    x, y, dy, idx1, idx2 = get_values( lum_bin_cens, z_gal_agn_lum_func[i], e_z_gal_agn_lum_func[i] )
    non_zero = np.where( z_gal_agn_lum_func[i] != 0.0 )[0]
    print( x )
    print( lum_bin_cens[non_zero] )
    print( y )
    print( z_gal_agn_lum_func[i][non_zero] )
    p3.plot( lum_bin_cens[non_zero], z_gal_agn_lum_func[i][non_zero], color=zcols_agn[i], linewidth=3, alpha=0.75, linestyle='dotted' )
    # rectangular integration
    gal_trapz = np.sum( np.power( 10., y ) * dl )
    e_gal_trapz = np.sqrt( np.sum( np.power( dy * np.log( 10. ) * np.power( 10., y ), 2. ) ) ) * dl 
    non_zero = np.where( z_agn_lum_func[i] != 0.0 )[0]
    p3.plot( lum_bin_cens[non_zero], z_agn_lum_func[i][non_zero], color=zcols_agn[i], label='{:s} < z < {:s}'.format(str(zbin_starts[i]),str(zbin_ends[i])), linewidth=3 )
    trapz = np.trapz( np.power( 10., z_agn_lum_func[i][non_zero] ), lum_bin_cens[non_zero] )
    agn_delta_int.append(trapz / gal_trapz)
l1 = p3.legend()
dline = matplotlib.lines.Line2D([0],[0], color='black', linewidth=3 )
lline = matplotlib.lines.Line2D([0],[0], color='black', linewidth=3, linestyle='dotted' )
l2 = p3.legend((lline,dline),('Galaxies', 'Process'),loc='upper right')
p3.add_artist(l1)
p3.axes.set_xlim(plxlims)
p3.axes.set_ylim(plylims)
p3.xaxis.set_visible(False)
p3.set_title('Active Galactic Nuclei',fontsize=20)
p3.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p3.set_ylabel('log'+r'$_{10}$'+'('+r'$\rho$'+' [Mpc'+r'$^{-3}$'+' log'+r'$L^{-1}$'+'])')

## Right panel (bottom): ratio of the RLFs by galaxies and process
p4 = plt.axes([0.14+sbsizex*fsizey/fsizex,0.1,sbsizex*fsizey/fsizex,0.4*sbsizey])
p4.plot( (19,27), (1,1), color='gray', linestyle='dashed', linewidth=1.5 )
p4.plot( (19,27), (2,2), color='gray', linewidth=1, alpha=0.25 )
for i in np.arange(0,len(z_lum_bins)):
    non_zero_idx = np.where( np.logical_and( z_agn_lum_func[i] != 0, z_gal_agn_lum_func[i] != 0 ) )[0]
    ratio = np.power(10., z_agn_lum_func[i][non_zero_idx] ) / np.power( 10., z_gal_agn_lum_func[i][non_zero_idx] )
    p4.plot( lum_bin_cens[non_zero_idx], ratio, color=zcols_agn[i], linewidth=3 )
p4.axes.set_xlim(plxlims)
p4.axes.set_ylim((0.7,2.5))
p4.set_xlabel('log'+r'$_{10}$'+'('+r'$L_{\mathrm{144 MHz}}$'+' [W Hz'+r'$^{-1}$'+'])')
p4.set_ylabel(r'$\Delta$'+'RLF')

fig.savefig(paths.figures / 'RLF_evolution.png',dpi=300)
fig.clear()
plt.close()

## write an output table
with open( paths.output / 'integrated_differences.txt', 'w' ) as f:
    f.write( '\\begin{table}\n' )
    f.write( '    \\centering\n' )
    f.write( '    \\begin{tabular}{lccccc}\n' )
    zmin = ' $z_{min}$ '
    zmax = ' $z_{max}$ '
    for i in np.arange(0,len(z_lum_bins)):
        zmin = zmin + ' & {:s}'.format(str(zbin_starts[i]))
        zmax = zmax + ' & {:s}'.format(str(zbin_ends[i]))
    zmin = zmin + ' \\\\ \n'
    zmax = zmax + ' \\\\ \hline \n'
    f.write( zmin ) 
    f.write( zmax )
    sfstr = 'SF ' 
    for i in np.arange(0,len(sf_delta_int)):
        sfstr = sfstr + ' & {:1.2f}'.format(sf_delta_int[i])
        sfstr = sfstr + '$\\pm$' + '{:1.1f}'.format(e_sf_delta_int[i]) 
    sfstr = sfstr + ' \\\\ \n' 
    f.write(sfstr)
    agnstr = 'AGN ' 
    for i in np.arange(0,len(agn_delta_int)):
        agnstr = agnstr + ' & {:1.2f}'.format(agn_delta_int[i])
        agnstr = agnstr + '$\\pm$' + '{:1.1f}'.format(e_agn_delta_int[i]) 
    agnstr = agnstr + ' \\\\ \n' 
    f.write(agnstr)
    f.write( '    \\end{tabular}\n' )
    f.write( '    \\caption{Integrated $\\Delta$RLF, calculated as the ratio of areas under the RLF curve by process to the RLF curve by galaxy classification.}\n' )
    f.write( '    \\label{tab:intvalues}\n' )
    f.write( '\\end{table}' )



