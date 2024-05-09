import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table, join
import astropy.units as u
from progress.bar import Bar
import paths
from astropy.stats import median_absolute_deviation as apy_mad
import os
from astropy.io import fits

## error propagation functions
def add_sub_error( eX, eY ):
    result = np.sqrt( np.power( eX, 2. ) + np.power( eY, 2. ) )
    return( result )


def radio_power( obs_flux_Jy, redshift, spectral_index=-0.8):
    flux_cgs = obs_flux_Jy * 1e-23
    DL = cosmo.luminosity_distance(redshift).value * 3.086e24 ## convert to cm
    power_cgs = 4. * np.pi * np.power(DL,2.) * flux_cgs / np.power( 1+redshift, 1.+spectral_index )
    ## convert to W/Hz
    power = power_cgs * 1e-7
    return( power )

def trapezoid_rule( yvals, dx ):
    first_term = dx / 2. * yvals[0]
    last_term = dx / 2. * yvals[-1]
    middle_terms = dx * np.nansum(yvals[1:-2])
    add_terms = first_term + middle_terms + last_term
    return( add_terms )

def midpoint_rule( yvals, dx ):
    sum_vals = np.nansum( yvals * dx )
    return( sum_vals )

######################################################################
## Functions for calculating brightness temp

def get_tb_information( mycat, im_weight=0.5, maj_lim=0.4, min_lim=0.3, T_e=1e4, alpha=-0.7, ref_freqs=np.array([0.003,0.01,0.03,0.1,0.3,1,3])*1e9, freqs_GHz=np.arange( 1e-3, 1e2, 1e-3 ), use_z=True ):
    freqs = freqs_GHz*1e9  ## Hz

    ## step 1: sizes
    beam_solid_angle = get_beam_sizes( mycat, im_weight=im_weight, maj_lim=maj_lim, min_lim=min_lim )
    ## step 2: flux density per solid angle in units of mJy per arcsec squared
    compact_flux_per_SA_peak = mycat['Peak_flux']*1e3 / beam_solid_angle
    compact_flux_per_SA_total = mycat['Total_flux']*1e3 / beam_solid_angle

    if use_z:
        print('Using redshift information to find T_b-identified AGN')
        ## use redshift information redshift information to find high T_b sources
        peak_agn = find_agn_withz( compact_flux_per_SA_peak, mycat['Z_BEST'], T_e=T_e, rf_array=ref_freqs, freq=144e6, alpha=alpha)
        total_agn = find_agn_withz( compact_flux_per_SA_total, mycat['Z_BEST'], T_e=T_e, rf_array=ref_freqs, freq=144e6, alpha=alpha)
    else:
        ## High T_b sources, assuming nu_0 = 3 GHz
        peak_agn = find_agn( compact_flux_per_SA_peak, T_e=T_e, rf_array=ref_freqs, freq=144e6, alpha=alpha)
        total_agn = find_agn( compact_flux_per_SA_total, T_e=T_e, rf_array=ref_freqs, freq=144e6, alpha=alpha)

    agn_final = np.unique( np.concatenate((peak_agn,total_agn)))
    print('Number of identified AGN: ', len(agn_final), ' which is ', len(agn_final)/len(mycat)*100, ' percent of sources.')
    print(len(peak_agn), 'are identified via peak brightness')
    print(len(total_agn), 'are identified via total brightness')

    ## add some information to the table
    ## column for peak (=1) or total (=2) T_b identification
    tb_from = np.zeros(len(mycat))
    tb_from[total_agn] = 2.
    tb_from[peak_agn] = 1. ## overwrites total
    ## add single column for flux density per solid angle 
    compact_flux_per_SA = np.copy(compact_flux_per_SA_total)
    compact_flux_per_SA[peak_agn] = compact_flux_per_SA_peak[peak_agn]

    ## return the flux density per solid angle and whether it was identified using peak or total
    return( compact_flux_per_SA, tb_from )

def get_beam_solid_angle( theta1, theta2 ):
    ## for a 2D Gaussian: https://www.cv.nrao.edu/~sransom/web/Ch3.html eqn 3.118
    result =  np.pi / (4. * np.log(2.)) * theta1 * theta2 
    return( result )

def get_beam_sizes( mycat, im_weight=0.5, maj_lim=0.4, min_lim=0.3 ):
    ## where the deconvolved resolution is smaller than the limiting resolution, use the limiting resolution
    ## See Radcliffe et al. (2018) Eqn. 6 (also Lobanov 2005)
    lim_factor = np.power( 2., 2-0.5*im_weight ) * np.sqrt( np.log(2.)/np.pi * np.log( mycat['Peak_flux']/mycat['Isl_rms'] / ( mycat['Peak_flux']/mycat['Isl_rms'] -1. ) ) )
    maj_lim_theory = maj_lim * lim_factor / 60. / 60. ## convert to degrees
    min_lim_theory = min_lim * lim_factor / 60. / 60. ## convert to degrees
    maj_idx = np.where( mycat['DC_Maj'] < maj_lim_theory )[0]
    min_idx = np.where( mycat['DC_Min'] < min_lim_theory )[0]
    beam_solid_angle = get_beam_solid_angle( mycat['DC_Maj']*60.*60., mycat['DC_Min']*60.*60. )
    return(beam_solid_angle)

def find_agn( flux_per_SA, T_e=1e4, rf_array=np.array([0.003,0.01,0.03,0.1,0.3,1,3])*1e9, freq=144e6, alpha=-0.8 ):
    ## predict brightness temps at frequency
    flux_per_SA_rfs = []
    for rf in rf_array:        
        tau_norm = 1. 
        opt_depth_term = 1. - np.exp( - tau_norm * ( freq / rf )**-2.1 )
        tb = T_e * opt_depth_term * ( 1. + 10. * (freq/1e9)**(0.1+alpha) ) 
        mJyperasec2 = tb * (freq)**2. / 1.38e24 * 1e3
        flux_per_SA_rfs.append(mJyperasec2)
    maxval = np.max(np.array(flux_per_SA_rfs)) 
    idx = np.where(flux_per_SA >= maxval)[0]
    return( idx )

def find_agn_withz( flux_per_SA_vec, redshift_vec, T_e=1e4, rf_array=np.array([0.003,0.01,0.03,0.1,0.3,1,3])*1e9, freq=144e6, alpha=-0.8 ):
    idx = []
    for i in np.arange(0,len(flux_per_SA_vec)):
        flux_per_SA = flux_per_SA_vec[i]
        if np.isfinite(flux_per_SA):
            redshift = redshift_vec[i]    
            ## predict brightness temps at frequency
            flux_per_SA_rfs = []
            for rf in rf_array:        
                tau_norm = 1. 
                opt_depth_term = 1. - np.exp( - tau_norm * ( freq / rf )**-2.1 )
                tb = T_e * opt_depth_term * ( 1. + 10. * (freq/1e9)**(0.1+alpha) ) 
                mJyperasec2 = tb * (freq/(1.+redshift))**2. / 1.38e24 * 1e3
                flux_per_SA_rfs.append(mJyperasec2)
            maxval = np.max(np.array(flux_per_SA_rfs))
            if flux_per_SA >= maxval:
                idx.append(i)
    return( idx )

######################################################################
## Functions for calculating vmaxes

def effective_area( rms_file ):
    tmp = os.path.basename( rms_file )
    outfile = paths.data / tmp.replace( '.fits', '_effective_areas.fits' ) 
    print(outfile)
    if os.path.exists( outfile ):
        t = Table.read( outfile, format='fits' )
        noise_bins = t['noise_bin']
        eff_area_deg2 = t['eff_area_deg2']
    else:
        ## Effective areas as function of flux density 
        rms_image = fits.open( rms_file )
        rms = rms_image[0].data.flatten()
        ## pixel scales -- this is in degrees
        pixscale = np.abs(rms_image[0].header['CDELT1'])
        pixarea = np.power(pixscale, 2.)
        ## first get where this is finite
        rms_idx = np.where(np.isfinite(rms))
        ## get cumulative distribution of effective areas
        rms = rms[rms_idx]  ## cuts down on size of rms to make it faster
        noise_min = np.min(rms)
        noise_max = np.max(rms)
        binsize = 1. * np.power( 10., -(len(str(1./noise_min).split('.')[0])+1) )
        noise_bins = np.arange(noise_min, noise_max, binsize)
        npix = []
        for nb in noise_bins:
            tmp_idx = np.where(rms <= nb)
            npix.append(len(tmp_idx[0]))
        eff_area_deg2 = np.asarray(npix) * pixarea
        t = Table()
        t.add_column( noise_bins, name='noise_bin' )
        t.add_column( eff_area_deg2, name='eff_area_deg2' )
        t.write( outfile, format='fits' )
    return( noise_bins, eff_area_deg2 )

def get_source_vmax( power, power_to_flux_factor, V_DLs, completeness, fluxdens, solid_angle, noise_bins, noise_min, sigma_cut = 5. ):
    ## find the flux at all possible redshifts
    flux_at_z = power_to_flux_factor * power  
    ## find the maximum redshift at which it can be observed       
    noise_diff = np.where( flux_at_z >= sigma_cut*noise_min )[0]
    ## if this is in the middle of the redshift bin, use an index to limit to max redshift
    if len(noise_diff) < len(flux_at_z):
        flux_at_z = power_to_flux_factor[noise_diff] * power
        V_DLs = V_DLs[noise_diff]
    ## initialise some empty lists
    solid_angle_at_z = []
    completeness_at_z = []
    ## loop over the fluxes at different redshift
    for j in np.arange(0,len(flux_at_z)):
        ## find the solid angle for which this S(z) is observable
        idx = np.where(flux_at_z[j] >= sigma_cut*noise_bins)[0]
        solid_angle_at_z.append(solid_angle[np.max(idx)])
        ## find the completeness correction for this S(z)
        idx = np.where( np.abs(flux_at_z[j]-fluxdens) == np.min( np.abs(flux_at_z[j]-fluxdens) ) )[0]
        completeness_at_z.append(completeness[idx[0]])
    ## calculate theta_z
    theta_sz = np.asarray(solid_angle_at_z) * np.asarray(completeness_at_z)
    integrand = V_DLs * theta_sz
    ## trapezoid rule it?? 
    #vmax = trapezoid_rule( integrand, 0.0001 )
    vmax = np.nansum( integrand )
    return(vmax)

def get_vmax( lotss, field, col_suffix='', zmin=0.003, zmax=0.3, dz=0.0001, si=-0.8, sigma_cut=5., rms_image='', cochrane='', kondapally='', test=False ):
    
    nsrc = len(lotss)
    print('Starting with {:s} sources in the area of coverage'.format(str(nsrc)))

    ## to select the rested flux density columns
    fluxcol = 'Total_flux{:s}'.format(col_suffix)
    fluxerrcol = 'E_Total_flux{:s}'.format(col_suffix)

    ## remove things with no redshifts, i.e. those which are masked
    good_z = np.where( np.logical_not(lotss['Z_BEST'].mask) )[0]
    lotss = lotss[good_z]
    ## remove the mask to avoid warnings later on
    lotss['Z_BEST'] = np.ma.getdata(lotss['Z_BEST'])

    nsrc = len(lotss)
    print('Continuing with {:s} sources with unmasked redshifts.'.format(str(nsrc)))

    ## limit to the redshift range requested
    z_idx = np.where( np.logical_and( lotss['Z_BEST'] >= zmin, lotss['Z_BEST'] <= zmax  ) )[0]
    lotss = lotss[z_idx]

    nsrc = len(lotss)
    print('Continuing with {:s} sources in the redshift range of {:s} to {:s}'.format(str(nsrc),str(zmin),str(zmax)))

    if test:
        ## randomly select 200 sources
        test_idx = np.random.choice(np.arange(0,len(lotss)),200)
        lotss = lotss[test_idx]
        nsrc = len(lotss)
        print('Test is true, continuing with {:s} sources.'.format(str(nsrc)))        

    ## get solid angle
    noise_bins, eff_area_deg2 = effective_area( rms_image )
    eff_area_sr = eff_area_deg2 * np.power( ( np.pi / 180. ), 2. )
    solid_angle = eff_area_sr / ( 4. * np.pi )
    ## convert the noise bins from Jy to W/Hz/m^2
    noise_bins = noise_bins * 1e-26
    noise_min = np.min(noise_bins)

    ## calculate intrinsic powers: first convert from Jy to W/Hz/m^2
    fluxes = lotss[fluxcol] * 1e-26
    flux_errors = lotss[fluxerrcol] * 1e-26
    redshifts = lotss['Z_BEST']
    ## calculate the power
    tmp = cosmo.luminosity_distance(redshifts)
    ## convert to metres so units are consistent
    DL = tmp.to(u.m).value
    power = 4. * np.pi * np.power( DL, 2. ) * fluxes / np.power( (1.+redshifts), (1+si) )
    power_error = 4. * np.pi * np.power( DL, 2. ) * flux_errors / np.power( (1.+redshifts), (1+si) )

    ## calculate distances and volumes 
    zbins = np.arange(zmin,zmax+dz,dz)
    tmp = cosmo.comoving_volume(zbins)  ## these are in Mpc3
    V_DLs = tmp[1:] - tmp[0:-1] ## length of one per bin
    zbin_cens = np.arange(zmin+0.5*dz,zmax,dz)
    tmp = cosmo.luminosity_distance(zbin_cens)
    DLs = tmp.to(u.m).value

    ## keep some calculations out of the loop
    DL_fac = 4. * np.pi * np.power( DLs, 2. )
    redshift_kcorr = np.power( (1.+zbin_cens), (1+si) ) 
    power_to_flux_factor = redshift_kcorr / DL_fac 

    ## now loop over every source to get the vmax
    vmaxes = []
    #vmaxes_lo = []
    #vmaxes_up = []
    tmp = field.capitalize()
    pb = Bar('Processing', max=len(lotss))
    for i in np.arange(0,len(lotss)):
        if lotss['Overall_class'][i] == 'SFG':
            field_idx = [ i for i,cn in enumerate(cochrane.colnames) if tmp[0] in cn ]
            colname = np.asarray(cochrane.colnames)[field_idx][0]
            completeness = cochrane[colname]
            completeness_fluxDens = cochrane['FluxDensity_mJy']*1e-29
        else:
            field_idx = [ i for i,cn in enumerate(kondapally.colnames) if tmp[0] in cn ]
            colname = np.asarray(kondapally.colnames)[field_idx][0]
            completeness = kondapally[colname]
            completeness_fluxDens = kondapally['FluxDensity_mJy']*1e-29
        try:
            completeness = completeness.filled(0)
        except:
            completeness = completeness

        vmax = get_source_vmax( power[i], power_to_flux_factor, V_DLs, completeness, completeness_fluxDens, solid_angle, noise_bins, noise_min, sigma_cut = sigma_cut )
        #vmax_lo = get_source_vmax( power[i]-power_error[i], power_to_flux_factor, V_DLs, completeness, completeness_fluxDens, solid_angle, noise_bins, noise_min, sigma_cut = sigma_cut )
        #vmax_up = get_source_vmax( power[i]+power_error[i], power_to_flux_factor, V_DLs, completeness, completeness_fluxDens, solid_angle, noise_bins, noise_min, sigma_cut = sigma_cut )
        vmaxes.append(vmax)
        #vmaxes_lo.append(vmax_lo)
        #vmaxes_up.append(vmax_up)
        pb.next()
    pb.finish()

    lotss.add_column( vmaxes, name='vmax' )
    #lotss.add_column( vmaxes_lo, name='vmax_lo' )
    #lotss.add_column( vmaxes_up, name='vmax_up' )
    return( lotss )
