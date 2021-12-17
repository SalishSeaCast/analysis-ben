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
    if any([j < size, i < size, j >= n - size, i >= m - size]): return None, None
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


def get_NEMO_tracers(daterange, dataframe, NEMO, gridvars, gridref, tol_time=600, tol_lonlat=0.001):
    """
    """
    
    # Extract fields
    start, end = [parse(date) for date in daterange]
    data = {'time': dataframe.index.values.astype('datetime64[s]').astype(datetime)}
    variables = ['longitude', 'latitude', 'depth', 'temperature', 'salinity']
    names = ['long_deg', 'lat_deg', 'pressure_dbar', 'temp_degc', 'salinity_teos10']
    for var, name in zip(variables, names):
        data[var] = dataframe[name].values.astype(float)
        
    # Remove nan and out-of-range values
    index = np.logical_and.reduce([~np.isnan(data[var]) for var in variables])
    index = np.logical_and(index, [start <= t <= end for t in data['time']])
    for var in data: data[var] = data[var][index]
    
    # Interpolate NEMO
    for var in ['votemper', 'vosaline']: data[var] = []
    time_last, lon_last, lat_last = datetime(2016, 1, 1), 0, 0
    dataarray = {}
    counter = 0
    for time, lon, lat, depth in zip(*[data[var] for var in ('time', 'longitude', 'latitude', 'depth')]):
        for var in ['votemper', 'vosaline']:

            # Check if next station
            next_station = np.logical_or.reduce((
                abs(time - time_last) > timedelta(seconds=tol_time),
                abs(lon - lon_last) > tol_lonlat,
                abs(lat - lat_last) > tol_lonlat,
            ))
            if next_station:
                points, dataarray[var] = get_interpolation_args(NEMO[var], time, lon, lat, gridvars, gridref)
                dataarray[var] = dataarray[var].interpolate_na(dim='deptht', limit=2, fill_value='extrapolate')

            # Interpolate NEMO
            values = dataarray[var].interp(deptht=depth, kwargs={'fill_value': 'extrapolate'}).values.ravel()
            data[var].append(float(interpolate.griddata(points, values, (lon, lat))))

        # Update previous values
        time_last, lon_last, lat_last = time, lon, lat
        counter += 1
        if counter == 10000:
            print('Next 10000 complete ...')
            sys.stdout.flush()
            counter = 0
    
    # Lists to arrays
    for var in ['votemper', 'vosaline']: data[var] = np.array(data[var])
    
    return data


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
            points, dataarray = get_interpolation_args(NEMO[var], time, lon, lat, gridvars[var], gridref)
            value = float(interpolate.griddata(points, dataarray.values.ravel(), (lon, lat))) if points is not None else np.nan
            vel[var].append(value)
    
    return {'time': times, 'longitude': lons, 'latitude': lats, 'u': np.array(vel['u']), 'v': np.array(vel['v'])}


def main_routine_drifters(
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


def main_routine(
    runID, startdate, enddate,
    path_NEMO='/scratch/bmoorema/Results/Currents',
    path_obs='/home/bmoorema/project/SalishSea/obs',
    path_save='/scratch/bmoorema/Tracers',

):
    """Process NEMO tracers for STARTDATE, ENDDATE
    (STARTDATE, ENDDATE must match the results file path)
    """

    # Initialize
    daterange = [startdate, enddate]
    datestr = '_'.join(daterange)
    sys.stdout = open(f'{path_save}/stdout/SSCtracers_{runID}_{datestr}.stdout', 'w')
    rundir = f'SalishSeaCast_currenttuning_{runID}_{datestr}'
    print('Processing results from ' + rundir)
    sys.stdout.flush()
    
    # Load PSF CitSci data
    dataframe = pd.read_csv(f'{path_obs}/CitSci2015_2019.csv', parse_dates=True, index_col='datetime_utc')

    # Load NEMO
    NEMO = {}
    fn = f'{path_NEMO}/{rundir}/SalishSea_1h_{datestr}_grid_T.nc'
    for var in ['votemper', 'vosaline']: NEMO[var] = xr.open_dataset(fn)[var]

    # Get tracers and save to database
    con = sql.connect(f'{path_save}/SSCtracers_{runID}_{datestr}.sqlite')
    gridvars, gridref = get_grid_variables(variables=['t'])
    gridvars = gridvars['t']
    data = get_NEMO_tracers(daterange, dataframe, NEMO, gridvars, gridref)
    pd.DataFrame(data).to_sql('PSF', con, index=False)
    con.close()
    print('Done!')
    sys.stdout.close()


if __name__ == "__main__":
    main_routine(*[str(arg) for arg in sys.argv[1:]])
