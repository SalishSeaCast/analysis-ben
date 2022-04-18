#!/usr/bin/env python

import numpy as np
import xarray as xr
import pandas as pd
import sqlite3 as sql
import os
from datetime import datetime, timedelta
from dateutil.parser import parse
from scipy import io, interpolate
import sys


def initialize(
    runID, category, startdate, enddate,
    path_save='/scratch/bmoorema',
):
    """
    """
    
    # Define strings and paths
    daterange = [parse(d) for d in (startdate, enddate)]
    datestr = f'{startdate}_{enddate}'
    fn_sql = f'SSC{category}_{runID}_{datestr}'
    rundir = f'SalishSeaCast_currenttuning_{runID}_{datestr}'
    path_run = os.path.join(rundir, f'SalishSea_1h_{datestr}')
    
    # Open sql and stdout files and print status
    con = sql.connect(os.path.join(path_save, category, f'{fn_sql}.sqlite'))
    sys.stdout = open(os.path.join(path_save, category, 'stdout', f'{fn_sql}.stdout'), 'w')
    print('Processing results from ' + rundir)
    sys.stdout.flush()
    
    return con, path_run, daterange


def loopstatus(k, n):
    """
    """
    
    nscale = 100 / n
    percent = k * nscale
    if percent % 10 < nscale:
        print(f'{int(percent)}% complete ...')
        sys.stdout.flush()


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


def get_drifter_coords(ID, drifters):
    """
    """
    
    index = np.where(drifters['id'].astype(int) == ID)[0][0]
    times, lons, lats = [drifters[var][index].squeeze() for var in ('mtime', 'lon', 'lat')]
    times = np.array([datetime.fromordinal(int(t) - 366) + timedelta(days=t%1) for t in times])
    
    return times, lons, lats


def get_bbox_index(coords, bbox):
    """
    """
    
    index = True
    for coord, bounds in zip(coords, bbox):
        for bound, func in zip(bounds, ['greater_equal', 'less_equal']):
            if bound is not None:
                index = np.logical_and(index, getattr(np, func)(coord, bound))
    
    return index


def rotate_vel(u_in, v_in, theta=29):
    """Rotate velocities from map coordinates to grid coordinates
    """
    
    theta_rad = math.radians(theta)
    u_out = u_in * np.cos(theta_rad) + v_in * np.sin(theta_rad)
    v_out = -u_in * np.sin(theta_rad) + v_in * np.cos(theta_rad)
    
    return u_out, v_out


def interpolate_NEMO(data, time, lon, lat, gridvars, gridref, depths=None, size=1):
    """Interpolate NEMO results at a time,lon,lat point and in the vertical if
    depths is provided. Uses Xarray dimensional interpolation for time and depth,
    and Scipy griddata for lon and lat. The size variable defines a bounding box
    to minimize the cost of griddata.
    """
    
    # Find j,i from lon,lat
    n, m = gridvars['lonlat'][0].shape
    j, i = [gridref[var].sel(lats=lat, lons=lon, method='nearest').item() for var in ('jj', 'ii')]
    
    # Return NaN if j,i is outside of domain
    if any([j < size, i < size, j >= n - size, i >= m - size]):
        return np.nan

    # Extract coordinate-results pairs inside bbox
    jslc, islc = [slice(coord - size, coord + size + 1) for coord in (j, i)]
    points = np.vstack([coord[jslc, islc].ravel() for coord in gridvars['lonlat']]).T
    dataarray = data.isel(y=jslc, x=islc).interp(time_counter=time).where(gridvars['mask'][..., jslc, islc])
    
    # 3D interpolation (e.g. CTD data)
    if depths is not None:
        dataarray = dataarray.interpolate_na(dim='deptht', limit=2, fill_value='extrapolate')
        interpolated = []
        for depth in depths:
            values = dataarray.interp(deptht=depth, kwargs={'fill_value': 'extrapolate'}).values.ravel()
            interpolated.append(float(interpolate.griddata(points, values, (lon, lat))))
        interpolated = np.array(interpolated)
    
    # 2D interpolation (e.g. drifters, ferry data)
    else:
        interpolated = float(interpolate.griddata(points, dataarray.values.ravel(), (lon, lat)))
    
    return interpolated


def process_drifters(
    daterange, NEMO, gridvars, gridref,
    bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='Salish_L3_20190728T103529.mat',

):
    """Interpolate NEMO velocites at times and positions based on drifter ID
    and return dict of arrays
    """
    
    # Load drifters file
    drifters = io.loadmat(os.path.join(path, file))['drift'][0]
    
    # Initialize data dict
    variables = ['time', 'longitude', 'latitude', 'u_obs', 'v_obs', 'u_nemo', 'v_nemo']
    data = {var: [] for var in variables}
    
    # Get drifter IDs within daterange-bbox
    # and times on the half-hour for loading NEMO results into memory
    IDs, times_NEMO = list(np.unique(drifters['id'].astype(int))), []
    for ID in IDs:
        times, lons, lats = get_drifter_coords(ID, drifters)
        index = get_bbox_index([times, lons, lats], [daterange, bbox[:2], bbox[2:]])
        if all(index):  # -- ID inside range: keep and get times for NEMO
            for time in times:
                t = time.replace(minute=30, second=0, microsecond=0)
                times_NEMO.append([t + timedelta(hours=hour) for hour in (-1, 0, 1)])
        else:  # ----------- ID outside range: remove
            IDs.remove(ID)
    
    # Load NEMO time into memory
    times_NEMO = np.unique(times_NEMO)
    for var in ['u', 'v']:
        NEMO[var].sel(time_counter=times_NEMO).load()

    # Loop through drifter IDs
    n = len(IDs)
    for k, ID in zip(range(n), IDs):
    
        # Get drifter coordinates
        times, lons, lats = get_drifter_coords(ID, drifters)
        
        # Calculate drifter velocities based on geolocation
        deltat, uv = np.array([dt.total_seconds() for dt in np.diff(times)]), []
        for var, coord, theta in zip(['u', 'v'], [lons, lats], [np.deg2rad(lats[:-1]), 0]):
            uv.append(np.diff(coord) / deltat * 111000 * np.cos(theta))
        u, v = rotate_vel(*uv)
        
        # Interpolate coords to midpoint and append drifter to data dict
        times, lons, lats = [(var[:-1] + var[1:]) / 2 for var in (times, lons, lats)]
        for var, values in zip(variables, [times, lons, lats, u, v]):
            data[var].append(values)

        # Get interpolated NEMO results
        for var in ['u', 'v']:
            for time, lon, lat in zip([times, lons, lats]):
                value = interpolate_NEMO(NEMO[var], time, lon, lat, gridvars[var], gridref)
                data[var + '_nemo'].append(value)
        
        # Print status
        loopstatus(k, n)
    
    # Concatenate
    for var in variables:
        data[var] = np.hstack(data[var])
        
    return data


def process_PSF(
    daterange, NEMO, gridvars, gridref,
    bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='CitSci2015_2019.csv',
):
    """
    """
    
    # Load data
    dataframe = pd.read_csv(os.path.join(path, file), parse_dates=[3])
    
    # Define naming conventions and initialize data dict
    variables = ['time', 'longitude', 'latitude', 'depth', 'temperature', 'salinity']
    names = ['datetime_utc', 'long_deg', 'lat_deg', 'pressure_dbar', 'temp_degc', 'salinity_teos10']

    # Extract variables from dataframe and build NaN index
    data, index = {}, True
    for var, name in zip(variables, names):
        data[var] = dataframe[name].values
        if var == 'time':
            data[var] = data[var].astype('datetime64[s]').astype(datetime)
        else:
            data[var] = data[var].astype(float)
            index = np.logical_and(index, ~np.isnan(data[var]))

    # Add daterange-bbox criteria to index and extract variables
    index = get_bbox_index([data[var] for var in variables[:3]], [daterange, bbox[:2], bbox[2:]])
    for var in variables:
        data[var] = data[var][index]

    # Get ID and station arrays
    IDs, stations = [dataframe[name][index].values.astype(str) for name in ('ID', 'station')]

    # Initialize NEMO to NaN
    n = sum(index)
    for var in ['votemper', 'vosaline']:
        data[var] = np.empty(n)
        data[var][:] = np.nan

    # Nested loop over cruise IDs and stations
    IDs_unique = np.unique(IDs)
    n = len(IDs_unique)
    for k, ID in zip(range(n), IDs_unique):
        index_ID = (IDs == ID)
        for station in np.unique(stations[index_ID]):
            
            # Get cast index and interpolation coordinates
            index = np.logical_and(index_ID, stations == station)
            times, lons, lats, depths = [data[var][index] for var in variables[:4]]
            
            # Get interpolated NEMO results
            for var in ['votemper', 'vosaline']:
                
                # coords are constant over cast, so use zero index
                values = interpolate_NEMO(NEMO[var], times[0], lons[0], lats[0], gridvars, gridref, depths=depths)
                data[var][index] = np.array(values)
        
        # Print status
        loopstatus(k, n)
        
    return data


def process_DFO(
    daterange, NEMO, gridvars, gridref,
    bbox=[-125.5, -122.7, 48, 50.6],
    path='/home/bmoorema/project/SalishSea/obs',
    file='DFO_CTD.sqlite',
):
    """
    """
    
    # Load dataframes from sql
    con = sql.connect(os.path.join(path, file))
    df_coords = pd.read_sql_query('SELECT * from StationTBL', con)
    dataframe = pd.read_sql_query('SELECT * from CalcsTBL', con)
    con.close()
    
    # Define variable names, datanames and initialize storage dict
    variables = ['depth', 'temperature', 'salinity']
    datanames = [
        ['Z'],                                                                    # depth
        ['Temperature_Primary_CT', 'Temperature_Secondary_CT', 'Temperature_CT'], # temperature
        ['Salinity_T0_C0_SA', 'Salinity_T1_C1_SA', 'Salinity_SA'],                # salinity
    ]
    data = {var: [] for var in variables}
    
    # Get IDs, dates, lons, lats and remove IDs outside daterange and bbox
    coordnames = ['ID', 'Lon', 'Lat', 'StartYear', 'StartMonth', 'StartDay', 'StartHour']
    IDs, lons, lats, y, m, d, h = [df_coords[name].values for name in coordnames]
    y, m, d = [var.astype(int) for var in (y, m, d)]
    times = np.array([datetime(*args[:3]) + timedelta(hours=args[3]) for args in zip(y, m, d, h)])
    index = get_bbox_index([times, lons, lats], [daterange, bbox[:2], bbox[2:]])
    IDs, times, lons, lats = IDs[index], times[index], lons[index], lats[index]
    
    # Loop through casts and interpolate NEMO
    n = len(IDs)
    for k, ID, time, lon, lat in zip(range(n), IDs, times, lons, lats):
        
        # Extract fields
        index_ID, index = dataframe['StationTBLID'] == ID, True
        ncast = sum(index_ID)
        for var, names in zip(variables, datanames):
            for name in names:  # Try all field names since they vary between separate casts
                if any(dataframe[name][index_ID].notna()):
                    profile = dataframe[name][index_ID]
                    break
            else:
                profile = np.ones(ncast) * np.nan
            index = np.logical_and(index, ~np.isnan(profile))
            data[var].append(profile)

        # Trim NaNs
        for var in variables:
            data[var][-1] = data[var][-1][index]
        
        # Add coordinate arrays to cast
        ncast = sum(index)
        for var, value in zip(['time', 'longitude', 'latitude'], [time, lon, lat]):
            data[var].append(np.ones(ncast) * value)

        # Interpolate NEMO
        for var in ['votemper', 'vosaline']:
            
            values = interpolate_NEMO(NEMO[var], time, lon, lat, gridvars, gridref, depths=data['depth'][-1])
            data[var].append(values)
        
        # Print status
        loopstatus(k, n)
    
    # Concatenate casts
    for var in data:
        data[var] = np.hstack(data[var])
    
    return data


def main_routine(
    runID, category, startdate, enddate,
    path_NEMO='/project/def-allen/bmoorema/results/Currents',
):
    """Process NEMO tracers for STARTDATE, ENDDATE
    (STARTDATE, ENDDATE must match the results file path)
    """

    # Initialize
    con, path_run, daterange = initialize(runID, category, startdate, enddate)
    
    # Process NEMO tracers
    if category in ['DFO', 'PSF']:
        gridvars, gridref = get_grid_variables(variables=['t'])
        gridvars = gridvars['t']
        NEMO = xr.open_dataset(os.path.join(path_NEMO, f'{path_run}_grid_T.nc'))
        
        # Process DFO
        if category == 'DFO':
            data = process_DFO(daterange, NEMO, gridvars, gridref)
        
        # Process PSF
        elif category == 'PSF':
            data = process_PSF(daterange, NEMO, gridvars, gridref)
    
    # Process NEMO velocities
    elif category == 'drifters':
        gridvars, gridref = get_grid_variables()
        NEMO = {}
        for var, name in zip(['u', 'v'], ['vozocrtx', 'vomecrty']):
            print(f'Loading {var} ...')
            sys.stdout.flush()
            fn = os.path.join(path_NEMO, f'{path_run}_grid_{var.upper()}.nc')
            NEMO[var] = xr.open_dataset(fn)[name].isel({f'depth{var}': 0})
        
        # Process drifters
        data = process_drifters(daterange, NEMO, gridvars, gridref)
        
    else:
        raise ValueError(f'Unrecognized category: {category}')

    # Close out
    pd.DataFrame(data).to_sql(category, con, index=False)
    con.close()
    print('Done!')
    sys.stdout.close()


if __name__ == "__main__":
    main_routine(*[str(arg) for arg in sys.argv[1:]])
