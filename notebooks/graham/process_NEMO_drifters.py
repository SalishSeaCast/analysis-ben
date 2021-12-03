#!/usr/bin/env python

import numpy as np
import xarray as xr
import pandas as pd
import sqlite3 as sql
from datetime import datetime, timedelta
from dateutil.parser import parse
from scipy import io, interpolate
import sys


def mtime_to_datetime(mtime):
    """Convert from Matlab time to datetime.datetime
    """
    
    return datetime.fromordinal(int(mtime) - 366) + timedelta(days=mtime%1)


def get_interpolation_args(data, time, lon, lat, gridvars, gridref, size=1):
    """Return a bounding box of coordinate, value pairs around a lonlat point
    """
    
    n, m = gridvars['lonlat'][0].shape
    j, i = [gridref[var].sel(lats=lat, lons=lon, method='nearest').item() for var in ('jj', 'ii')]
    if any([j < size, i < size, j >= n - size, i >= m - size]): return np.nan, np.nan
    jslc, islc = [slice(coord - size, coord + size + 1) for coord in (j, i)]
    points = np.vstack([coord[jslc, islc].ravel() for coord in gridvars['lonlat']]).T
    dataarray = data.isel(y=jslc, x=islc).interp(time_counter=time).where(gridvars['mask'][..., jslc, islc])
    
    return points, dataarray


def get_grid_variables(variables=['u', 'v'], path='/home/bmoorema/MEOPAR/grid/'):
    """Return grid variables dict GRIDVARS with LONLAT and MASK fields for U and V.
    Also return the GRIDREF xarray object for looking up ji from lonlat
    """

    # Master files
    coords = xr.open_dataset(path + 'coordinates_seagrid_SalishSea201702.nc', decode_times=False)
    mask = xr.open_dataset(path + 'mesh_mask202108.nc')
    gridref = xr.open_dataset(path + 'grid_from_lat_lon_mask999.nc')

    # Grid variables
    gridvars = {}
    for var in variables:
        gridvars[var] = {}
        gridvars[var]['lonlat'] = [coords[key][0, ...].values for key in (f'glam{var}', f'gphi{var}')]
        gridvars[var]['mask'] = mask[f'{var}mask'][0, ...].values
    
    return gridvars, gridref


def get_drifter_IDs(daterange, drifters):
    """Return drifter IDs and unique observation times on the half-hour
    within DATERANGE
    """
    
    
    # Get IDs list
    start, end = [parse(date) for date in daterange]
    IDs, IDout = list(np.unique(drifters['id'].astype(int))), []
    for ID in IDs:
        dindex = np.where(drifters['id'].astype(int) == ID)[0][0]
        if drifters['lon'][dindex][0] < -124.5: IDout.append(ID)
        if not all(start <= mtime_to_datetime(t) <= end for t in drifters['mtime'][dindex][:, 0]): IDout.append(ID)
    IDout = list(np.unique(IDout))
    for ID in IDout: IDs.remove(ID)

    # Get inclusive times on the half-hour
    times, kwargs = [], {'minute': 30, 'second': 0, 'microsecond': 0}
    for ID in IDs:
        dindex = np.where(drifters['id'].astype(int) == ID)[0][0]
        for t in drifters['mtime'][dindex][:, 0]:
            for hour in range(-1, 2):
                times.append((mtime_to_datetime(t) + timedelta(hours=hour)).replace(**kwargs))
    
    return IDs, np.unique(times)


def get_NEMO_velocities(ID, drifters, NEMO, gridvars, gridref):
    """Interpolate NEMO velocites at times and positions based on drifter ID
    and return dict of arrays
    """
    
    # Load drifters
    dindex = np.where(drifters['id'].astype(int) == ID)[0][0]
    times = np.array([mtime_to_datetime(t) for t in drifters['mtime'][dindex][:, 0]])
    lons, lats = [drifters[key][dindex][:, 0] for key in ('lon', 'lat')]
    
    # Load NEMO
    vel = {'u': [], 'v': []}
    for time, lon, lat in zip(times, lons, lats):
        for var in ['u', 'v']:
            points, mask, values = get_interpolation_args(NEMO[var], time, lon, lat, gridvars[var], gridref)
            values = np.ma.masked_where(mask == 0, values).ravel()
            vel[var].append(float(interpolate.griddata(points, values, (lon, lat))))
    
    return {'time': times, 'longitude': lons, 'latitude': lats, 'u': np.array(vel['u']), 'v': np.array(vel['v'])}


def main_routine(
    runID, startdate, enddate,
    path_NEMO='/scratch/bmoorema/Results/Currents',
    path_obs='/home/bmoorema/project/SalishSea/obs',
    path_save='/scratch/bmoorema/Drifters',

):
    """Process NEMO velocities from runID onto drifter trajectories for STARTDATE, ENDDATE
    (STARTDATE, ENDDATE must match the results file path)
    """

    # Initialize
    daterange = [startdate, enddate]
    datestr = '_'.join(daterange)
    sys.stdout = open(f'{path_save}/stdout/SSCdrifters_{runID}_{datestr}.stdout', 'w')
    rundir = f'SalishSeaCast_currenttuning_{runID}_{datestr}'
    print('Processing results from ' + rundir)
    sys.stdout.flush()
    
    # Load drifters
    drifters = io.loadmat(f'{path_obs}/Salish_L3_20190728T103529.mat')['drift'][0]
    IDs, times = get_drifter_IDs(daterange, drifters)

    # Load NEMO
    NEMO = {}
    for var, name in zip(['u', 'v'], ['vozocrtx', 'vomecrty']):
        print(f'Loading {var} ...')
        sys.stdout.flush()
        fn = f'{path_NEMO}/{rundir}/SalishSea_1h_{datestr}_grid_{var.upper()}.nc'
        NEMO[var] = xr.open_dataset(fn)[name].isel({f'depth{var}': 0}).sel(time_counter=times).load()

    # Get velocities and save to database
    con = sql.connect(f'{path_save}/SSCdrifters_{runID}_{datestr}.sqlite')
    gridvars, gridref = get_grid_variables()
    for n, ID in enumerate(IDs):
        print(f'Loading drifter {ID:03d} ({n+1} of {len(IDs)}) ...')
        sys.stdout.flush()
        velocities = get_NEMO_velocities(ID, drifters, NEMO, gridvars, gridref)
        pd.DataFrame(velocities).to_sql(f'drifter{ID:03d}', con, index=False)
    con.close()
    print('Done!')
    sys.stdout.close()


if __name__ == "__main__":
    main_routine(*[str(arg) for arg in sys.argv[1:]])
