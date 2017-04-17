#! /usr/bin/env python

import numpy as np
import xarray as xr
import progressbar
import pickle
import timeseries_tools as tstools
from scipy.io import savemat
from salishsea_tools import nc_tools

# Timerange
timerange = ['2016 Jan 1 00:00', '2016 Dec 31 23:00']
timeslice = slice(timerange[0], timerange[1])

# Get mask and grid files
mask_NEMO = xr.open_dataset(
    '/ocean/bmoorema/research/MEOPAR/NEMO-forcing/grid/mesh_mask_downbyone2.nc')

GRIDMASK = {}

# Mask
tmask = mask_NEMO.tmask.isel(t=0).values.astype(bool)
tmask[:, 750:, :] = 0
tmask[:, :350, :] = 0
tmask[:, :, :100] = 0

# Grid and depth
deptht, gridY, gridX = np.meshgrid(
    mask_NEMO.gdept_1d.isel(t=0), mask_NEMO.y, mask_NEMO.x, indexing='ij')
deptht = mask_NEMO.gdept_0.isel(t=0).values

# Nowcast domain slice parameters
indices = [0, 0, 20, 450, 520, 680]
spacing = [1, 5, 5, 1, 1, 1]
dims = [0, 0, 0, 1, 1, 1]

# Create mask and grid arrays for Nowcast domain slices
for index, space, dim in zip(indices, spacing, dims):
    
    # Depth 20 is index 18
    if index == 20:
        ii = 18
    else:
        ii = index
        
    # Store mask and grid arrays in dict
    GRIDMASK[f'spc{space}_{index}'] = tstools.reshape_grid(
        tmask, deptht, gridY, gridX, index=ii, dim=dim, spacing=space)

# dz
deptht = deptht.reshape(
    -1, GRIDMASK['spc1_0']['ngrid'])[:, GRIDMASK['spc1_0']['tmask']]
dz = np.diff(deptht, axis=0)

# Build filename lists
filenames_T = nc_tools.make_filename_list(
    timerange, 'T', model='nowcast-green', resolution='h')
filenames_U = nc_tools.make_filename_list(
    timerange, 'U', model='nowcast-green', resolution='h')
filenames_V = nc_tools.make_filename_list(
    timerange, 'V', model='nowcast-green', resolution='h')

# Predefine storage dictionaries
NEMO, FULL, TRIM = {}, {}, {}

# -------- Main loop ---------------
with progressbar.ProgressBar(max_value=len(filenames_T)) as bar:
    for i, filename in bar(enumerate(zip(filenames_U, filenames_V, filenames_T))):
        
        # Open NEMO results into numpy arrays
        U = xr.open_dataset(filename[0]).vozocrtx
        V = xr.open_dataset(filename[1]).vomecrty
        T = xr.open_dataset(filename[2]).votemper
        S = xr.open_dataset(filename[2]).vosaline
        
        # Calculate density
        rho = tstools.calc_rho(S, T, mask_NEMO.gdept_0.isel(t=0)).values
        rho = rho.reshape((-1, rho.shape[1], GRIDMASK['spc1_0']['ngrid']))
        TRIM['pycnocline'] = np.zeros((S.shape[0], GRIDMASK['spc1_0']['ngrid_water']))
        
        # Reshape Parameters
        indices = [0, 20, 450, 520, 680]
        spacing = [5, 5, 1, 1, 1]
        dims = ['depth', 'depth', 'y', 'y', 'y']
        tracers = [False, False, True, True, True]
        
        # Slice and reshape NEMO results and allocate trimmed output arrays
        for index, space, dim, tracer in zip(indices, spacing, dims, tracers):
            
            # Hash slice dims
            if dim is 'depth':
                dimu, dimv, dimt = 'depthu', 'depthv', 'deptht'
            else:
                dimu, dimv, dimt = dim, dim, dim
            
            # Reshape tracers
            if tracer:
                FULL[f'T{index}'], TRIM[f'T{index}'] = tstools.reshape_data(
                    T, GRIDMASK[f'spc{space}_{index}'], dimt, index=index, spacing=space)
                FULL[f'S{index}'], TRIM[f'S{index}'] = tstools.reshape_data(
                    S, GRIDMASK[f'spc{space}_{index}'], dimt, index=index, spacing=space)
            
            # Reshape velocity
            if space > 1:
                dimx, dimy = 'x', 'y'
            else:
                dimx, dimy = None, None
            FULL[f'U{index}'], TRIM[f'U{index}'] = tstools.reshape_data(
                U, GRIDMASK[f'spc{space}_{index}'], dimu,
                index=index, spacing=space, unstagger_dim=dimx)
            FULL[f'V{index}'], TRIM[f'V{index}'] = tstools.reshape_data(
                V, GRIDMASK[f'spc{space}_{index}'], dimv,
                index=index, spacing=space, unstagger_dim=dimy)
        
        # Trim Land Points and calculate halocline depth
        for tindex, timestamp in enumerate(U.time_counter):
            
            # Loop through slice parameters
            for index, space, tracer in zip(indices, spacing, tracers):
                
                # Trim tracers
                if tracer:
                    TRIM[f'T{index}'][tindex, :] = FULL[f'T{index}'][tindex, :][
                        GRIDMASK[f'spc{space}_{index}']['tmask']]
                    TRIM[f'S{index}'][tindex, :] = FULL[f'S{index}'][tindex, :][
                        GRIDMASK[f'spc{space}_{index}']['tmask']]
                
                # Trim velocity
                TRIM[f'U{index}'][tindex, :] = FULL[f'U{index}'][tindex, :][
                    GRIDMASK[f'spc{space}_{index}']['tmask']]
                TRIM[f'V{index}'][tindex, :] = FULL[f'V{index}'][tindex, :][
                    GRIDMASK[f'spc{space}_{index}']['tmask']]
            
            # Pycnocline depth
            idrhodz = (np.diff(rho[tindex][:, GRIDMASK['spc1_0']['tmask']],
                               axis=0) / dz).argmax(axis=0)
            TRIM['pycnocline'][tindex, :] = np.array([depth[dindex]
                                for dindex, depth in zip(idrhodz, deptht.T)])
        
        # Concatenate into storage arrays
        for index, tracer in zip(indices, tracers):
            
            # Assign if first file
            if i == 0:
                if tracer:
                    NEMO[f'T{index}'] = TRIM[f'T{index}']
                    NEMO[f'S{index}'] = TRIM[f'S{index}']
                
                NEMO[f'U{index}'] = TRIM[f'U{index}']
                NEMO[f'V{index}'] = TRIM[f'V{index}']
            
            # Otherwise concatenate
            else:
                if tracer:
                    NEMO[f'T{index}'] = np.concatenate(
                        [NEMO[f'T{index}'], TRIM[f'T{index}']], axis=0)
                    NEMO[f'S{index}'] = np.concatenate(
                        [NEMO[f'S{index}'], TRIM[f'S{index}']], axis=0)
                
                NEMO[f'U{index}'] = np.concatenate(
                    [NEMO[f'U{index}'], TRIM[f'U{index}']], axis=0)
                NEMO[f'V{index}'] = np.concatenate(
                    [NEMO[f'V{index}'], TRIM[f'V{index}']], axis=0)
        
        # Concatenate pycnocline depth
        NEMO['pycnocline'] = TRIM['pycnocline']
        NEMO['pycnocline'] = np.concatenate(
            [NEMO['pycnocline'], TRIM['pycnocline']], axis=0)
        
        # Update progress bar
        bar.update(i)

# Save NEMO and GRIDMASK
savemat('/ocean/bmoorema/research/MEOPAR/analysis-ben/data/NEMO_2016.mat', NEMO)
fid = open('/ocean/bmoorema/research/MEOPAR/analysis-ben/data/GRIDMASK_2016', 'wb') 
pickle.dump(GRIDMASK, fid)   
fid.close()
