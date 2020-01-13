# ===================================================================
# The contents of this file are dedicated to the public domain.  To
# the extent that dedication to the public domain is not available,
# everyone is granted a worldwide, perpetual, royalty-free,
# non-exclusive license to exercise all rights associated with the
# contents of this file for any purpose whatsoever.
# No rights are reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ===================================================================
'''
Created on Thur Aug 22 07:39:56 2019

A tool to generate indices for ingesting into ARIANE

Example usage:
    # first create a logical 2D mask where particles are to be released
    import ariane_setup as ar
    ar_ind=ar.ariane_indices(mask,'domain_cfg.nc')
    ar_ind.get_dir_list('path_to_nemo_data','*U.nc','*V.nc')
    ar_ind.gen_ind()
    ar_ind.plot_map(0,0)
    ar_ind.write_file('initial_positions.txt',1)

@author James Harle

$Last commit on:$
'''

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from matplotlib.colors import ListedColormap
from os import listdir, path
import fnmatch
#from netCDF4 import netcdftime
import cftime

class ariane_indices(object):
    """
    An object that stores the information required to generate an ARIANE 
    initial positions file
    """
    def __init__(self, mask, cfgfile, pvol=None):
        """ 
        Initialise IceCav object.
        
        Args:
            cfgfile    (str): configuration filename

        Returns:

        """
        self.mask     = mask    # numpy 2D logical array of active cells
        self.cfgfile  = cfgfile # domain_cfg filename
        self.pvol     = None    # 
        self.flist    = None    # 
        self.a_ind    = np.zeros((0,5))    # 
        self.datetime = None    #
        
        # Check dimensions are correct
        nc = Dataset(cfgfile)
        idim = nc.dimensions['x'].size
        jdim = nc.dimensions['y'].size
        kdim = nc.dimensions['z'].size
        tdim = nc.dimensions['t'].size
        
        j_msk, i_msk = mask.shape
        
        if (i_msk != idim) | (j_msk != jdim):
            print('Dimension mismatch') # TODO: Throw an exception here
            
        # Find the subdomain for ARIANE release        
        k_ind = (0, kdim) # TODO: add option for selective depths
        t_ind = (0, tdim) # TODO: add option for multiple time slices
        
        j_ind, i_ind = self._find_sub_inds(mask)
        
        # Open pointer to cfgfile
        nc    = Dataset(cfgfile)
        
        # TODO: Clean up adding 1 to indices as we may already be reading in
        #       the whole array
        tl = nc.variables['top_level'][   0, 
                                          j_ind[0]:j_ind[1],
                                          i_ind[0]:i_ind[1]]
        bl = nc.variables['bottom_level'][0, 
                                          j_ind[0]:j_ind[1],
                                          i_ind[0]:i_ind[1]]
        
        self.tmask = self._gen_mask(tl, bl, kdim, 'T')[:, 1:-1, 1:-1]
        self.umask = self._gen_mask(tl, bl, kdim, 'U')[:, 1:-1, :   ]
        self.vmask = self._gen_mask(tl, bl, kdim, 'V')[:, :   , 1:-1]
        
        self.e2u   = nc.variables['e2u'][0, 
                                         j_ind[0]+1:j_ind[1]-1,
                                         i_ind[0]  :i_ind[1]-1]
        self.e1v   = nc.variables['e1v'][0, 
                                         j_ind[0]  :j_ind[1]-1,
                                         i_ind[0]+1:i_ind[1]-1]
        self.e3u   = nc.variables['e3u_0'][0, :,
                                           j_ind[0]+1:j_ind[1]-1,
                                           i_ind[0]  :i_ind[1]-1]
        self.e3v   = nc.variables['e3v_0'][0, :,
                                           j_ind[0]  :j_ind[1]-1,
                                           i_ind[0]+1:i_ind[1]-1]
        
        # Close pointer to cfgfile
        nc.close()
        
        # Assign indices to self
        self.t_ind = t_ind
        self.k_ind = k_ind
        self.j_ind = j_ind
        self.i_ind = i_ind
        
    def gen_ind(self):
        """ 
        Generate the location indices for an ARIANE experiment.
    
        Args:
            
    
        Returns:
        """
        # Shorten indices
        i = self.i_ind
        j = self.j_ind
        
        # Set up array to hold partical info
        mask_local = self.mask[j[0]+1:j[1]-1,
                               i[0]+1:i[1]-1] * self.tmask > 0
        k_loc, j_loc, i_loc = np.where(mask_local==1)
        k_loc += 0.5
        active_cells = np.sum(mask_local, dtype=int)
        P = np.zeros((active_cells,len(self.flist)))
        U = np.zeros_like(P)
        V = np.zeros_like(P)
        
        # Check if flist is empty - if so ask for list of directory to populate
        # list
        
        if self.flist == None:
            print('No file list to work from') # TODO: Throw an exception
        
        # TODO: select the files required based on the t_indices input
           
        count = 0
        
        # TODO: if pvol is negative asign constant volume flux (therefore 
        # constant number of particles) to each cell.
        
        for fn in self.flist:
            
            nc  = Dataset(fn[0])
            u   = nc.variables['vozocrtx'][0, :,
                                       j[0]+1:j[1]-1,
                                       i[0]  :i[1]-1]
            nc.close()
            nc  = Dataset(fn[1])
            v   = nc.variables['vomecrty'][0, :,
                                       j[0]  :j[1]-1,
                                       i[0]+1:i[1]-1]
            nc.close()
            
            U[:,count] = ( ( 
                    u[:,:,:-1] * self.e2u[:,:-1] * 
                    self.e3u[:,:,:-1] * self.umask[:,:,:-1] +
                    u[:,:,1: ] * self.e2u[:,1: ] * 
                    self.e3u[:,:,1: ] * self.umask[:,:,1: ] 
                    ) * 0.5 )[mask_local]
            
            V[:,count] = ( (
                    v[:,:-1,:] * self.e1v[:-1,:] * 
                    self.e3v[:,:-1,:] * self.vmask[:,:-1,:] +
                    v[:,1: ,:] * self.e1v[1: ,:] * 
                    self.e3v[:,1: ,:] * self.vmask[:,1: ,:] 
                    ) * 0.5 )[mask_local]
            
            count += 1
        
        if self.pvol == None:
            
            # Derive a pvol number based on the max flow through a cell, let us
            # set the maximum number of particles to 64 as a first estimate
            
            self.pvol = np.amax(U+V)/64.
        
        pvol_r = 1/self.pvol
            
            #TODO: shall we set all values >0 & <1 =1?
            
        # Work out the decomposition of particles per cell    
        nx, ny = self._decomp(U*pvol_r, V*pvol_r)
        
        self.P = nx*ny
        
        for t in range(count):
            for n in range(active_cells):
            
                jj, ii = np.meshgrid(np.linspace(0, 1, ny[n,t]*2+1), 
                                     np.linspace(0, 1, nx[n,t]*2+1))
            
                ii = ii[1:-1,1:-1].flatten()[:,np.newaxis]
                jj = jj[1:-1,1:-1].flatten()[:,np.newaxis]
                
                self.a_ind = np.vstack( (self.a_ind,
                                         np.hstack((ii+i_loc[n]+i[0], 
                                                    jj+j_loc[n]+j[0], 
                                                    np.ones_like(ii)+k_loc[n], 
                                                    np.ones_like(ii)+t, 
                                                    np.ones_like(ii)))
                                         ) )
                                           
    def write_file(self, filename, t_ind):
        """ 
        Writes file

        Args:

        Returns:
            
        """
        
        # Make sure that self.a_ind is 
        np.savetxt(filename, self.a_ind, fmt='%6.2f', delimiter=' ')
        
    def plot_map(self, tslice, kslice):
        """ 
        Creates a file list

        Args:


        Returns:

        """
        
        # Shorten indices
        i = self.i_ind
        j = self.j_ind
        mask_lcl = self.mask[j[0]+1:j[1]-1,
                               i[0]+1:i[1]-1] * self.tmask > 0
        data_gbl = np.zeros_like(self.mask)
        data_lcl = np.zeros_like(self.tmask)
        data_lcl[mask_lcl] = self.P[:,tslice]
        data_gbl[j[0]+1:j[1]-1, i[0]+1:i[1]-1] = np.squeeze(
                                                 data_lcl[kslice,:,:])
        
        ax = sns.heatmap(data_gbl)
        ax.invert_yaxis()
        return ax
    
    def _decomp(self, nu_int, nv_int): 
        """ 
        Work out the closest to square decomposition for particles within a 
        grid cell. At the moment this is only done in the 2D horizontal plane. 
        Perhaps this could be expanded to the vertical plane too.
    
        Args:
            nu_int    (int): total number of E-W particles
            nv_int    (int): total number of S-N particles
    
        Returns:
            nx        (int): E-W dimension
            ny        (int): S-N dimension
        """
        # Make sure we're dealing with floats
        nu = nu_int.astype(float)
        nv = nv_int.astype(float)
        nu = np.where(nu<0.5, 0.5, nu)
        nv = np.where(nv<0.5, 0.5, nv)
        
        # Total number of particles
        N = nu + nv
        
        # Set min(N)=1 if there is a flow? 
        
        # nx * ny = N
        nx = ( nu * N / nv ) ** 0.5 
        ny = N/nx

        # NB nx and ny are now floats. We convert back to integer values
        nx = nx.astype(int)
        ny = ny.astype(int)
        
        # Check which square value is closest to N
        solution_0 = np.abs( ( nx  * ny ) - N.astype(int) )
        solution_1 = np.abs( (( nx + 1 ) * ny ) - N.astype(int) )
        solution_2 = np.abs( ( nx * ( ny + 1 )) - N.astype(int) )
    
        # Which solution is the closest    
        ny[ (solution_0 != 0) & (solution_1 >  solution_2) ] += 1
        nx[ (solution_0 != 0) & (solution_1 <= solution_2) ] += 1
        
        return nx, ny
    
    def get_dir_list(self, data_dir, regexU, regexV):
        """
        Create a filename list for the dataset
    
        Args:
            data_dir  (str): path to the data
            regex     (str): search string for files
    
        Returns:
            self
        """
        
        # Create an empty list to populate
        flistU = []
        flistV = []
        
        # Get directory listing and check against regex
        for file in listdir(data_dir):
            if fnmatch.fnmatch(file, regexU):
                flistU.append(path.join(data_dir, file))
            if fnmatch.fnmatch(file, regexV):
                flistV.append(path.join(data_dir, file))
                
        self.flist = zip(flistU, flistV)
    
        # TODO: Get the datetime information from the files
        # _get_source_timedata()
    
    def _get_source_timedata(self):
        """
        Get time information for the dataset
    
        Args:
    
        Returns:
            self
        """
        
        dir_list = self._get_dir_list(grid)
        group = GridGroup()
        group.data_list = []
        group.time_counter = []        
        group.date_counter = []
        for filename in dir_list:   
            nc = Dataset(filename, 'r')
            varid = nc.variables['time_counter'] 
            for index in range(0,len(varid)):
                x = [filename, index]
                group.data_list.append(x)
                group.time_counter.append(varid[index]+t_adjust)
                group.date_counter.append(cftime.utime(varid.units,
                                                           varid.calendar).num2date(varid[index]+t_adjust))
            group.units = varid.units
            group.calendar = varid.calendar
            nc.close()
        tmp_data_list = copy.deepcopy(group.data_list)
        tmp_time_counter = copy.deepcopy(group.time_counter)
        for index in range(len(group.time_counter)):
            tmp_data_list[index] = group.data_list[index]
            tmp_time_counter[index] = group.time_counter[index]
        group.data_list = tmp_data_list
        group.time_counter = tmp_time_counter
        return group
    
    def _gen_mask(self, top, bot, kdim, grd):
        """
        Generate 3D mask from teh top_level and bottom_level provided in the
        domain_cfg.nc file.
    
        Args:
            top       (numpy.ndarray): top grid cell index array
            bot       (numpy.ndarray): bottom grid cell index array
            kdim      (int)          : length of the depth dimension
            grd       (str)          : T/U/V indicating grid type
        
        Returns:
            msk       (numpy.ndarray): 3D mask array (1==wet)
        """
        
        # Collect dimensions
        jdim, idim = top.shape
        
        # Create 3D arrays of vertical index and bottom level to avoid looping
        wrk = np.tile(bot[np.newaxis,:,:],(kdim,1,1))
        kin = np.tile(np.arange(1, kdim+1)[:, np.newaxis, np.newaxis], 
                      (1, jdim, idim))
        
        # Define mask array
        msk = np.zeros((kdim, jdim, idim))
        
        # Determine where active wet cells are from bottom level
        msk = np.where(kin<=wrk, 1, msk)
        
        # Now adjust for top level (this would only apply for the isf case)
        wrk = np.tile(top[np.newaxis,:,:],(kdim,1,1))
        msk = np.where(kin<wrk, 0, msk)
        
        # Adjust for grid type
        
        if grd=='U':
            msk = msk[:, :, :-1] * msk[:, :, 1:]
        elif grd=='V':
            msk = msk[:, :-1, :] * msk[:, 1:, :]
            
        return msk
    
    def _find_sub_inds(self, mask):
        """
        Find the i/j extent of a 2D mask and return as a tuple. We add a 1 grid
        cell halo here to allow umask/vmask generation.
    
        Args:
            mask      (numpy.ndarray): 2D array 1==active domain, 0 elsewhere
        
        Returns:
            j_ind     (tuple)        : min/max j extent
            i_ind     (tuple)        : min/max i extent
        """      
        
        j, i = np.where(mask>0)
        
        j_ind = ((np.amin(j)-1, np.amax(j)+2))
        i_ind = ((np.amin(i)-1, np.amax(i)+2))
        
        return j_ind, i_ind


