import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table, join
import astropy.units as u
from progress.bar import Bar
import paths
from astropy.stats import median_absolute_deviation as apy_mad
import os
from astropy.io import fits

def radio_power( obs_flux_Jy, redshift, spectral_index=-0.8):
    flux_cgs = obs_flux_Jy * 1e-23
    DL = cosmo.luminosity_distance(redshift).value * 3.086e24 ## convert to cm
    power_cgs = 4. * np.pi * np.power(DL,2.) * flux_cgs / np.power( 1+redshift, 1.+spectral_index )
    ## convert to W/Hz
    power = power_cgs * 1e-7
    return( power )

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
    flux_at_z = power_to_flux_factor * power         
    noise_diff = np.where( flux_at_z >= sigma_cut*noise_min )[0]
    if len(noise_diff) < len(flux_at_z):
        flux_at_z = power_to_flux_factor[noise_diff] * power
        V_DLs = V_DLs[noise_diff]
    completeness = completeness.filled(0)
    solid_angle_at_z = []
    completeness_at_z = []
    for j in np.arange(0,len(flux_at_z)):
        ## first check when/if the source is above the noise at a given redshift
        idx = np.where(flux_at_z[j] >= sigma_cut*noise_bins)[0]
        ## sources will become too faint to be observed above a redshift and then this fails
        sa_idx = np.min(idx)
        solid_angle_at_z.append(solid_angle[sa_idx])
        ## completeness
        idx = np.where( np.abs(flux_at_z[j]-fluxdens) == np.min( np.abs(flux_at_z[j]-fluxdens) ) )[0]
        completeness_at_z.append(completeness[idx[0]])
    theta_sz = np.asarray(solid_angle_at_z) * np.asarray(completeness_at_z)
    vmax = np.nansum( V_DLs * theta_sz )
    return(vmax)

def get_vmax( lotss, field, zmin=0.003, zmax=0.3, dz=0.0001, si=-0.8 ):
    ## remove things with no redshifts, i.e. those which are masked
    good_z = np.where( np.logical_not(lotss['Z_BEST'].mask) )[0]
    lotss = lotss[good_z]
    ## remove the mask to avoid warnings later on
    lotss['Z_BEST'] = np.ma.getdata(lotss['Z_BEST'])

    test = True
    if test:
        good_z = np.where(lotss['Z_BEST'] < zmax )[0]
        lotss = lotss[good_z]

    ## get solid angle
    rms_image = paths.static / '{:s}_rms_starmask_optical.fits'.format(field)
    print(rms_image)
    noise_bins, eff_area_deg2 = effective_area( rms_image )
    eff_area_sr = eff_area_deg2 * np.power( ( np.pi / 180. ), 2. )
    solid_angle = eff_area_sr / ( 4. * np.pi )
    ## convert the noise bins from Jy to W/Hz/m^2
    noise_bins = noise_bins * 1e-26
    noise_min = np.min(noise_bins)

    ## calculate intrinsic powers: first convert from Jy to W/Hz/m^2
    fluxes = lotss['Total_flux'] * 1e-26
    flux_errors = lotss['E_Total_flux'] * 1e-26
    redshifts = lotss['Z_BEST']
    ## calculate the power
    tmp = cosmo.luminosity_distance(redshifts)
    DL = tmp.to(u.m).value
    power = 4. * np.pi * np.power( DL, 2. ) * fluxes / np.power( (1.+redshifts), (1+si) )
    power_error = 4. * np.pi * np.power( DL, 2. ) * flux_errors / np.power( (1.+redshifts), (1+si) )

    ## calculate distances and volumes 
    zbins = np.arange(zmin,zmax+dz,dz)
    tmp = cosmo.comoving_volume(zbins)
    ## these are in Mpc3
    V_DLs = tmp[1:] - tmp[0:-1]
    zbin_cens = np.arange(zmin+0.5*dz,zmax,dz)
    tmp = cosmo.luminosity_distance(zbin_cens)
    DLs = tmp.to(u.m).value

    ## sigma cutoff
    sigma_cut = 5.

    ## completeness corrections
    cochrane = Table.read( paths.static / 'cochrane_2023_tableA1.csv', format='csv', delimiter=',' )
    kondapally = Table.read( paths.static / 'kondapally_2022_table1.csv', format='csv', delimiter=',' )

    ## keep some calculations out of the loop
    DL_fac = 4. * np.pi * np.power( DLs, 2. )
    redshift_kcorr = np.power( (1.+zbin_cens), (1+si) ) 
    power_to_flux_factor = redshift_kcorr / DL_fac 

    ## now loop over every source to get the vmax
    vmaxes = []
    vmaxes_lo = []
    vmaxes_up = []
    pb = Bar('Processing', max=len(lotss))
    for i in np.arange(0,len(lotss)):
    #for i in np.arange(0,10):
        if lotss['Overall_class'][i] == 'SFG':
            completeness = cochrane[field.capitalize()]
            completeness_fluxDens = cochrane['FluxDensity_mJy']*1e-29
        else:
            completeness = kondapally[field.capitalize()]
            completeness_fluxDens = kondapally['FluxDensity_mJy']*1e-29
        vmax = get_source_vmax( power[i], power_to_flux_factor, V_DLs, completeness, completeness_fluxDens, solid_angle, noise_bins, noise_min, sigma_cut = sigma_cut )
        vmax_lo = get_source_vmax( power[i]-power_error[i], power_to_flux_factor, V_DLs, completeness, completeness_fluxDens, solid_angle, noise_bins, noise_min, sigma_cut = sigma_cut )
        vmax_up = get_source_vmax( power[i]+power_error[i], power_to_flux_factor, V_DLs, completeness, completeness_fluxDens, solid_angle, noise_bins, noise_min, sigma_cut = sigma_cut )
        vmaxes.append(vmax)
        vmaxes_lo.append(vmax_lo)
        vmaxes_up.append(vmax_up)
        pb.next()
    pb.finish()

    lotss.add_column( vmaxes, name='vmax' )
    lotss.add_column( vmaxes_lo, name='vmax_lo' )
    lotss.add_column( vmaxes_up, name='vmax_up' )
    return( lotss )

def RLF_from_zmax( Lum, zz, lum_bins, redshift_bins, myzmax, area_cov, area_units='deg2', error_type='rms', completeness=[] ):
    ## error_type can be "data" or "rms"
    ## "data" will use measured uncertainties to calculate upper and lower values and propagate the errors
    ## "rms" will use the estimator from Marshall 1985 and is equal to
    ## 1/(delta log L) * sqrt( sum( 1/vmax^2 )

    ## get log10 of the values
    log10_L_total = np.log10( Lum )

    ## create an array for the RLF
    lum_func = np.zeros( (len(lum_bins),len(redshift_bins)-1) )
    lum_func_up = np.copy(lum_func)
    lum_func_lo = np.copy(lum_func)
    ## median luminosities
    med_lum = np.copy( lum_func )
    med_lum_err = np.copy( lum_func )
    ## number of objects
    N_obj = np.copy( lum_func )
    
    if area_units == 'deg2':
        ## area is in sqare degrees, convert to sr
        area_sr = area_cov * np.power( ( np.pi / 180. ), 2. )
        area_val = area_sr / ( 4. * np. pi )
        ## to avoid evaluating if statements every loop
        area_scaling = np.repeat( area_val, len(redshift_bins) )
        area_scaling_zmax = np.repeat( 1., len(Lum))
    elif area_units == 'Mpc2':
        ## area in Mpc^2
        ## this is redshift invariant because H_0 is invariant
        angular_distance = cosmo.angular_diameter_distance( redshift_bins )
        area_scaling = area_cov / ( 4. * np.pi * np.power( angular_distance, 2. ) )
        angular_distance_zmax = cosmo.agular_diameter_distance( myzmax['z_max'] )
        area_scaling_zmax = area_cov / ( 4. * np.pi * np.power( angular_distance_zmax, 2. ) )

    pb = Bar('Processing', max=len(lum_bins)*len(redshift_bins))
    for rr in np.arange( 1,len(lum_bins) ):
        for ss in np.arange( 1, len(redshift_bins) ):
            ## volumes for the redshift bin
            vmin_z = cosmo.comoving_volume( redshift_bins[ss-1] ) * area_scaling[ss-1]
            vmax_z = cosmo.comoving_volume( redshift_bins[ss] ) * area_scaling[ss]
            lum_idx = np.where( np.logical_and( log10_L_total >= lum_bins[rr-1], log10_L_total < lum_bins[rr] ) )
            redshift_idx = np.where( np.logical_and( zz >= redshift_bins[ss-1], zz < redshift_bins[ss] ) )
            bin_idx = np.intersect1d( lum_idx, redshift_idx )
            N_obj[rr-1,ss-1] = len(bin_idx)
            if len(bin_idx) > 0:
                vmax = cosmo.comoving_volume( myzmax['z_max'][bin_idx] )
                greater_vmax = np.where( vmax > vmax_z )
                vmax[greater_vmax] = vmax_z
                vmax = vmax - vmin_z
                if len(completeness) > 0:
                    bin_completeness = np.asarray(completeness)[bin_idx]
                    tmp_idx = np.where(vmax + vmin_z == vmax_z)[0]
                    vmax = vmax / bin_completeness
                ## sum 1/vmax to get lum func
                lum_func[rr,ss-1] = np.log10( np.nansum( 1./vmax.value ) / ( lum_bins[rr] - lum_bins[rr-1] ) )
                ## now ... error propagation
                if error_type == 'data':
                    ## intuitively this seems better, but then I end up with such small errors that vmax_up = vmax or vmax_lo = vmax
                    vmax_up = cosmo.comoving_volume( myzmax['z_max_up'][ss-1] ) * area_scaling[ss-1]
                    vmax_up = vmax_up - vmin_z
                    greater_vmax = np.where( vmax_up > vmax_z )
                    vmax_up[greater_vmax] = vmax_up[greater_vmax] - vmax_z
                    vmax_lo = cosmo.comoving_volume( myzmax['z_max_lo'][ss-1] ) * area_scaling[ss-1]
                    vmax_lo = vmax_lo - vmin_z
                    greater_vmax = np.where( vmax_lo > vmax_z )
                    vmax_lo[greater_vmax] = vmax_lo[greater_vmax] - vmax_z
                    lum_func_up[rr,ss-1] = np.log10( np.nansum( 1./vmax_up.value ) )
                    lum_func_lo[rr,ss-1] = np.log10( np.nansum( 1./vmax_lo.value ) )
                elif error_type == 'rms':
                    err_term = np.sqrt( np.nansum( 1. / np.power(vmax.value, 2.) ) ) / ( lum_bins[rr] - lum_bins[rr-1] )
                    lum_func_up[rr,ss-1] = np.log10( np.power( 10., lum_func[rr,ss-1] ) + err_term )
                    test = np.power( 10., lum_func[rr,ss-1] ) - err_term
                    if test == 0:
                        ## set to zero
                        lum_func_lo[rr,ss-1] = 0.
                    else:
                        lum_func_lo[rr,ss-1] = np.log10( np.power( 10., lum_func[rr,ss-1] ) - err_term )
                med_lum[rr,ss-1] = np.nanmedian( np.log10(Lum[bin_idx]) )
                med_lum_err[rr,ss-1] = apy_mad( Lum[bin_idx] )
            else:
                print( 'Skipping empty bin.' )
                ## set to zero
                lum_func[rr,ss-1] = 0
                lum_func_up[rr,ss-1] = 0
                lum_func_lo[rr,ss-1] = 0
                med_lum[rr,ss-1] = 0
                med_lum_err[rr,ss-1] = 0
            pb.next()
    pb.finish()

    ## create a final array to return
    RLF = Table()
    RLF.add_column( lum_func.flatten(), name = 'RLF' )
    RLF.add_column( lum_func_up.flatten(), name = 'RLF_up' )
    RLF.add_column( lum_func_lo.flatten(), name = 'RLF_lo' )
    RLF.add_column( med_lum.flatten(), name = 'Lmedian' )
    RLF.add_column( med_lum_err.flatten(), name = 'Lmedian_err' )
    RLF.add_column( N_obj, name='N_obj' )

    return( RLF )

