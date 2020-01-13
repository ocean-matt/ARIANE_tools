#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 14:15:04 2020

@author: mpc2g13 (using script from jdha)
"""

from netCDF4 import Dataset
import numpy as np

def cfg_convert(hgr, zgr, bat, ln_zco=0, 
                ln_zps=0, ln_sco=0, ln_isfcav=0):
    '''
    Uses mesh mask files to produce a domain_cfg.nc file.
    
    Args:
            
    Returns:
        ''' 
    # Open pointers to file and get dimensions 
    nc_hgr = Dataset(hgr)
    nc_zgr = Dataset(zgr)
    nc_bat = Dataset(bat)
    
    bat = nc_bat.variables['Bathymetry'][:]
    
    nc_bat.close()
    
    nz, ny, nx = nc_zgr.variables['e3t_0'][0,:,:,:].shape
    
    # Create a blank domain_cfg.nc 
    
    dataset = Dataset('domain_cfg.nc', 'w', format='NETCDF4_CLASSIC')

    dataset.createDimension('x', nx)
    dataset.createDimension('y', ny)
    dataset.createDimension('z', nz)
    
    # create and populate
    nav_lon = dataset.createVariable('nav_lon', np.float32, ('y', 'x'))
    nav_lat = dataset.createVariable('nav_lat', np.float32, ('y', 'x'))
    nav_lev = dataset.createVariable('nav_lev', np.float32, 'z')
    
    nav_lon[:, :] = nc_hgr.variables['nav_lon'][:, :]
    nav_lat[:, :] = nc_hgr.variables['nav_lat'][:, :]
    nav_lev[:]    = nc_hgr.variables['nav_lev'][:]
    
    giglo = dataset.createVariable('jpiglo', "i4")
    gjglo = dataset.createVariable('jpjglo', "i4")
    gkglo = dataset.createVariable('jpkglo', "i4")
    
    giglo[:] = nx
    gjglo[:] = ny
    gkglo[:] = nz
    
    gperio = dataset.createVariable('jperio', "i4")
    
    gperio[:] = 0
    
    gzco = dataset.createVariable('ln_zco', "i4")
    gzps = dataset.createVariable('ln_zps', "i4")
    gsco = dataset.createVariable('ln_sco', "i4")
    gcav = dataset.createVariable('ln_isfcav', "i4")
    
    gzco[:] = ln_zco
    gzps[:] = ln_zps
    gsco[:] = ln_sco
    gcav[:] = ln_isfcav

    ge3t1d = dataset.createVariable('e3t_1d', np.float64, 'z')
    ge3w1d = dataset.createVariable('e3w_1d', np.float64, 'z')
    
    ge3t1d[:] = nc_zgr.variables['e3t_1d'][0, :]
    ge3w1d[:] = nc_zgr.variables['e3w_1d'][0, :]
    
    gitop = dataset.createVariable('top_level', "i4", ('y', 'x'))
    gibot = dataset.createVariable('bottom_level', "i4", ('y', 'x'))
    
    gitop[:, :] = nc_zgr.variables['mbathy'][0, :, :] * 0 + 1
    gibot[:, :] = nc_zgr.variables['mbathy'][0, :, :]
    
    gbat = dataset.createVariable('Bathymetry', np.float64, ('y', 'x'))
    gbat[:, :] = bat
    
    glamt = dataset.createVariable('glamt', np.float64, ('y', 'x'))
    glamu = dataset.createVariable('glamu', np.float64, ('y', 'x'))
    glamv = dataset.createVariable('glamv', np.float64, ('y', 'x'))
    glamf = dataset.createVariable('glamf', np.float64, ('y', 'x'))
    gphit = dataset.createVariable('gphit', np.float64, ('y', 'x'))
    gphiu = dataset.createVariable('gphiu', np.float64, ('y', 'x'))
    gphiv = dataset.createVariable('gphiv', np.float64, ('y', 'x'))
    gphif = dataset.createVariable('gphif', np.float64, ('y', 'x'))
    
    glamt[:, :] = nc_hgr.variables['glamt'][0, :, :]
    glamu[:, :] = nc_hgr.variables['glamu'][0, :, :]
    glamv[:, :] = nc_hgr.variables['glamv'][0, :, :]
    glamf[:, :] = nc_hgr.variables['glamf'][0, :, :]
    gphit[:, :] = nc_hgr.variables['gphit'][0, :, :]
    gphiu[:, :] = nc_hgr.variables['gphiu'][0, :, :]
    gphiv[:, :] = nc_hgr.variables['gphiv'][0, :, :]
    gphif[:, :] = nc_hgr.variables['gphif'][0, :, :]
    
    ge1t = dataset.createVariable('e1t', np.float64, ('y', 'x'))
    ge1u = dataset.createVariable('e1u', np.float64, ('y', 'x'))
    ge1v = dataset.createVariable('e1v', np.float64, ('y', 'x'))
    ge1f = dataset.createVariable('e1f', np.float64, ('y', 'x'))
    ge2t = dataset.createVariable('e2t', np.float64, ('y', 'x'))
    ge2u = dataset.createVariable('e2u', np.float64, ('y', 'x'))
    ge2v = dataset.createVariable('e2v', np.float64, ('y', 'x'))
    ge2f = dataset.createVariable('e2f', np.float64, ('y', 'x'))
    
    ge1t[:, :] = nc_hgr.variables['e1t'][0, :, :]
    ge1u[:, :] = nc_hgr.variables['e1u'][0, :, :]
    ge1v[:, :] = nc_hgr.variables['e1v'][0, :, :]
    ge1f[:, :] = nc_hgr.variables['e1f'][0, :, :]
    ge2t[:, :] = nc_hgr.variables['e2t'][0, :, :]
    ge2u[:, :] = nc_hgr.variables['e2u'][0, :, :]
    ge2v[:, :] = nc_hgr.variables['e2v'][0, :, :]
    ge2f[:, :] = nc_hgr.variables['e2f'][0, :, :]
    
    gfff = dataset.createVariable('ff_f', np.float64, ('y', 'x'))
    gfft = dataset.createVariable('ff_t', np.float64, ('y', 'x'))
    ge3t = dataset.createVariable('e3t_0', np.float64, ('z', 'y', 'x'))
    ge3w = dataset.createVariable('e3w_0', np.float64, ('z', 'y', 'x'))
    ge3u = dataset.createVariable('e3u_0', np.float64, ('z', 'y', 'x'))
    ge3v = dataset.createVariable('e3v_0', np.float64, ('z', 'y', 'x'))
    ge3f = dataset.createVariable('e3f_0', np.float64, ('z', 'y', 'x'))
    ge3uw = dataset.createVariable('e3uw_0', np.float64, ('z', 'y', 'x'))
    ge3vw = dataset.createVariable('e3vw_0', np.float64, ('z', 'y', 'x'))
    
    gfff[:, :] = nc_hgr.variables['ff'][0, :, :]
    gfft[:, :] = nc_hgr.variables['ff'][0, :, :] 
    ge3t[:, :, :]  = nc_zgr.variables['e3t_0'][0, :, :, :]
    ge3w[:, :, :]  = nc_zgr.variables['e3w_0'][0, :, :, :]
    ge3u[:, :, :]  = nc_zgr.variables['e3u_0'][0, :, :, :]
    ge3v[:, :, :]  = nc_zgr.variables['e3v_0'][0, :, :, :]
    ge3f[:, :, :]  = nc_zgr.variables['e3t_0'][0, :, :, :]  #TODO: need to calc
    ge3uw[:, :, :]  = nc_zgr.variables['e3t_0'][0, :, :, :] #TODO: need to calc
    ge3vw[:, :, :]  = nc_zgr.variables['e3t_0'][0, :, :, :] #TODO: need to calc
    
    nav_lon.units, nav_lon.long_name = 'km', 'X'
    nav_lat.units, nav_lat.long_name = 'km', 'Y'
    
    nc_hgr.close()
    nc_zgr.close()
    dataset.close()
