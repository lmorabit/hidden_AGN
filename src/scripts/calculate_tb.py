#!/usr/bin/python3

import paths
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.table import Table, Column, join
from helper_functions import *
from astropy.io import fits
import os
import time
## cosmology to match Kondapally and Cochrane
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

## read in catalogue

field = 'lotss'
infile = paths.static / '{:s}_03_matched_inMOC_inHR.fits'.format(field)
lotss = Table.read( infile, format='fits' )

## get an index for where there is high resolution information

hr_idx = np.where( lotss['Total_flux'] > 0 )[0]

tmpcat = lotss[hr_idx]

compact_flux_per_SA, tb_from = add_tb_information( tmpcat )

tmp_compcat_flux_per_SA = np.zeros(len(lotss))
tmp_tb_from = np.zeros(len(lotss))

tmp_compact_flux_per_SA[hr_idx] = compact_flux_per_SA
tmp_tb_from[hr_idx] = tb_from






