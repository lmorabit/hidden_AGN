import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
import astropy.units as u
from progress.bar import Bar
import paths
from astropy.stats import median_absolute_deviation as apy_mad

def radio_power( obs_flux_Jy, redshift, spectral_index=-0.8):
    flux_cgs = obs_flux_Jy * 1e-23
    DL = cosmo.luminosity_distance(redshift).value * 3.086e24 ## convert to cm
    power_cgs = 4. * np.pi * np.power(DL,2.) * flux_cgs / np.power( 1+redshift, 1.+spectral_index )
    ## convert to W/Hz
    power = power_cgs * 1e-7
    return( power )

def RLF_calculate_zmax( flux_limit, source_id, fluxes, flux_errors, redshifts, spectral_indices, outfile='' ):

	##########################################
	## this function has as input:
	## flux_limit -- survey flux limit in Jy
	## source_id -- a unique source identifier
	## fluxes -- source flux densities in Jy
	## flux_errors -- source flux density errors in Jy
	## redshifts -- source redshifts
	## spectral_indices -- source spectral indices
	## outfile -- the outfile name
	##
	## the output is the maximum redshift to which the source could be seen (with errors)

    ## convert from Jy to W/Hz/m^2
    flux_limit = flux_limit * 1e-26
    fluxes = fluxes * 1e-26
    flux_errors = flux_errors * 1e-26

    ## convert to intrinsic power 
    tmp = cosmo.luminosity_distance(redshifts)
    DL = tmp.to(u.m).value
    power = 4. * np.pi * np.power( DL, 2. ) * fluxes / np.power( (1.+redshifts), (1+spectral_indices) )
    power_error = 4. * np.pi * np.power( DL, 2. ) * flux_errors / np.power( (1.+redshifts), (1+spectral_indices) )

    ## limiting distance will be when 
    ## P_rest = 4*pi*flux_limit * D_L(z_limit)^2 / ( 1 + z_limit )
    ## or when D_L(z_limit)^2 / ( 1 + z_limit ) = P_rest / ( 4 * pi * flux_limit )
    galaxy_values = power / flux_limit / ( 4. * np.pi )
    galaxy_values_err = power_error / flux_limit / ( 4. * np.pi )
    redshift_range = np.arange( 0.001, 10.001, 0.001 )
    redshift_range2 = np.arange( 9, 30.001, 0.001 )
    tmp = cosmo.luminosity_distance(redshift_range)
    V_DL = tmp.to(u.m).value
    tmp = cosmo.luminosity_distance(redshift_range2)
    V_DL2 = tmp.to(u.m).value
    
    z_limits = []
    z_limits_up = []
    z_limits_lo = []

    pb = Bar('Processing', max=len(galaxy_values))
    for i in np.arange(0,len(galaxy_values)):
        distance_values = np.power( V_DL, 2. ) / np.power( (1. + redshift_range), (1. + spectral_indices[i]) )
        diff_vals = galaxy_values[i] - distance_values
        limit_idx = np.where( np.abs(diff_vals) == np.min(np.abs(diff_vals)) )[0]
        if limit_idx == len(distance_values):
            ## limiting redshift is larger than z = 10, use distance_values2
            distance_values = np.power( V_DL2, 2. ) / np.power( (1. + redshift_range2), (1. + spectral_indices[i]) )
            diff_vals = galaxy_values[i] - distance_values
            limit_idx = np.where( np.abs(diff_vals) == np.min(np.abs(diff_vals)) )[0]
            if limit_idx == len(distance_values2):
                ## extreme (probably unphysical) correction - get rid of it
                z_limits.append(0)
                z_limits_up.append(0)
                z_limits_lo.append(0)
            else:
                z_limits.append( redshift_range2[limit_idx] )
                diff_vals_lo = (galaxy_values[i] - galaxy_values_err[i]) - distance_values
                limits_err_idx = np.where( np.abs(diff_vals_lo) == np.min(np.abs(diff_vals_lo)) )
                z_limits_lo.append(redshift_range2[limits_err_idx])
                diff_vals_up = (galaxy_values[i] + galaxy_values_err[i]) - distance_values
                limits_err_idx = np.where( np.abs(diff_vals_up) == np.min(np.abs(diff_vals_up)) )
                z_limits_up.append( redshift_range2[limits_err_idx] )
        else:
            z_limits.append( redshift_range[limit_idx] )
            diff_vals_lo = (galaxy_values[i] - galaxy_values_err[i]) - distance_values
            limits_err_idx = np.where( np.abs(diff_vals_lo) == np.min(np.abs(diff_vals_lo)) )
            z_limits_lo.append(redshift_range[limits_err_idx])
            diff_vals_up = (galaxy_values[i] + galaxy_values_err[i]) - distance_values
            limits_err_idx = np.where( np.abs(diff_vals_up) == np.min(np.abs(diff_vals_up)) )
            z_limits_up.append( redshift_range[limits_err_idx] )
        pb.next()
    pb.finish()

    z_limits = np.array(z_limits)
    z_limits_lo = np.array(z_limits_lo)
    z_limits_up = np.array(z_limits_up)
    t = Table()
    t.add_column( source_id, name='Source_id' )
    t.add_column( z_limits, name='z_max' )
    t.add_column( z_limits_lo, name='z_max_lo' )
    t.add_column( z_limits_up, name='z_max_up' )

    outfits = outfile + '.fits'

    t.write( paths.data / outfits , format='fits', overwrite=True )
    return( t )

def RLF_from_zmax( Lum, zz, lum_bins, redshift_bins, myzmax, area_cov, area_units='deg2', error_type='rms' ):
    ## error_type can be "data" or "rms"
    ## "data" will use measured uncertainties to calculate upper and lower values and propagate the errors
    ## "rms" will use the estimator from Marshall 1985 and is equal to
    ## 1/(delta log L) * sqrt( sum( 1/vmax^2 )

    ## get log10 of the values
    log10_L_total = np.log10( Lum )

    ## create an array for the RLF
    lum_func = np.zeros( (len(lum_bins),len(redshift_bins)) )
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
                ## sum 1/vmax to get lum func
                lum_func[rr,ss] = np.log10( np.nansum( 1./vmax.value ) / ( lum_bins[rr] - lum_bins[rr-1] ) )
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
                    lum_func_up[rr,ss] = np.log10( np.nansum( 1./vmax_up.value ) )
                    lum_func_lo[rr,ss] = np.log10( np.nansum( 1./vmax_lo.value ) )
                elif error_type == 'rms':
                    err_term = np.sqrt( np.nansum( 1. / np.power(vmax.value, 2.) ) ) / ( lum_bins[rr] - lum_bins[rr-1] )
                    lum_func_up[rr,ss] = np.log10( np.power( 10., lum_func[rr,ss] ) + err_term )
                    test = np.power( 10., lum_func[rr,ss] ) - err_term
                    if test == 0:
                        ## set to zero
                        lum_func_lo[rr,ss] = 0.
                    else:
                        lum_func_lo[rr,ss] = np.log10( np.power( 10., lum_func[rr,ss] ) - err_term )
                med_lum[rr,ss] = np.nanmedian( np.log10(Lum[bin_idx]) )
                med_lum_err[rr,ss] = apy_mad( Lum[bin_idx] )
            else:
                print( 'Skipping empty bin.' )
                ## set to zero
                lum_func[rr,ss] = 0
                lum_func_up[rr,ss] = 0
                lum_func_lo[rr,ss] = 0
                med_lum[rr,ss] = 0
                med_lum_err[rr,ss] = 0
            pb.next()
    pb.finish()

    return( lum_func, lum_func_up, lum_func_lo, med_lum, med_lum_err, N_obj )

