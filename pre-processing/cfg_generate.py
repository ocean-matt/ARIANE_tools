#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
cfg_generate.py
Created on Mon Jan 13 14:44:04 2020

@author: mpc2g13
"""

import sys
sys.path.append('/home/users/mpc2g13/tools/ARIANE_tools/pre-processing/')

from ariane_convert import cfg_convert
hgr = '/gws/nopw/j04/nemo_vol1/ORCA0083-N006/domain/mesh_hgr.nc'
zgr = '/gws/nopw/j04/nemo_vol1/ORCA0083-N006/domain/mesh_zgr.nc'
bat = '/gws/nopw/j04/nemo_vol1/ORCA0083-N006/domain/bathymetry_ORCA12_V3.3.nc'

cfg_convert(hgr, zgr, bat, ln_zco=0, 
                ln_zps=0, ln_sco=0, ln_isfcav=0)
