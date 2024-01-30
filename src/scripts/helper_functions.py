import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
import astropy.units as u
from progress.bar import Bar
import paths

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

    t.write( paths.static / outfits , format='fits', overwrite=True )
    return( t )

