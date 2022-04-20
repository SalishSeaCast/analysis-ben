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


def loopstatus(k, n):
    """
    """
    
    nscale = 100 / n
    percent = k * nscale
    if percent % 10 < nscale:
        print(f'{int(percent)}% complete ...')
        sys.stdout.flush()


def get_bbox_index(coords, bbox):
    """
    """
    
    index = True
    for coord, bounds in zip(coords, bbox):
        for bound, func in zip(bounds, ['greater_equal', 'less_equal']):
            if bound is not None:
                index = np.logical_and(index, getattr(np, func)(coord, bound))
    
    return index


def get_drifter_coords(ID, data):
    """
    """
    
    index = np.where(data['id'].astype(int) == ID)[0][0]
    times, lons, lats = [data[var][index].squeeze() for var in ('mtime', 'lon', 'lat')]
    times = np.array([datetime.fromordinal(int(t) - 366) + timedelta(days=t%1) for t in times])
    
    return times, lons, lats


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


def load_drifter_data(
    daterange, bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='Salish_L3_20190728T103529.mat',
):
    """
    """
    
    # Load data from matfile
    data = io.loadmat(os.path.join(path, file))['drift'][0]
    
    # Get drifter IDs within daterange-bbox
    # and times on the half-hour for loading NEMO results into memory
    IDs, times = list(np.unique(data['id'].astype(int))), []
    IDs = IDs[:10]
    for ID in IDs:
        coords = get_drifter_coords(ID, data)
        index = get_bbox_index(coords, [daterange, bbox[:2], bbox[2:]])
        if all(index):  # -- ID inside range: keep and get times for NEMO
            for time in coords[0]:
                t = time.replace(minute=30, second=0, microsecond=0)
                times.append([t + timedelta(hours=hour) for hour in (-1, 0, 1)])
        else:  # ----------- ID outside range: remove
            IDs.remove(ID)
    
    return data, IDs, np.unique(times)


def load_PSF_data(
    daterange, bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='CitSci2015_2019.csv',
):
    """
    """
    
    # Load data from csv
    data = pd.read_csv(os.path.join(path, file), parse_dates=[3])
    
    # Get coords
    names, dtypes = ['datetime_utc', 'long_deg', 'lat_deg'], ['datetime64[s]', float, float]
    times, lons, lats = [data[name].values.astype(dtype) for name, dtype in zip(names, dtypes)]
    times = times.astype(datetime)

    # Apply daterange-bbox criteria and filter NaNs from lon and lat
    index = get_bbox_index([times, lons, lats], [daterange, bbox[:2], bbox[2:]])
    index = np.logical_and.reduce((index, ~np.isnan(lons), ~np.isnan(lats)))
    data = data[index]

    # Append station string to ID column to make unique identifiers
    data['ID'] = data['ID'].astype(str) + data['station']
    IDs = np.unique(data['ID'])
    IDs = IDs[:10]
    
    return data, IDs


def load_DFO_data(
    daterange, bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='DFO_CTD.sqlite',
):
    """
    """
    
    # Load data from sql
    con = sql.connect(os.path.join(path, file))
    coords = pd.read_sql_query('SELECT * from StationTBL', con)
    data = pd.read_sql_query('SELECT * from CalcsTBL', con)
    con.close()
    
    # Get coords
    coordnames = ['ID', 'Lon', 'Lat', 'StartYear', 'StartMonth', 'StartDay', 'StartHour']
    IDs, lons, lats, y, m, d, h = [coords[name].values for name in coordnames]
    y, m, d = [var.astype(int) for var in (y, m, d)]
    times = np.array([datetime(*args[:3]) + timedelta(hours=args[3]) for args in zip(y, m, d, h)])
    
    # Apply daterange-bbox criteria and filter NaNs using predetermined IDs
    IDs_nan = [790, 1147, 1148, 1999, 6226] + list(range(2847, 2853))
    index_nan = np.logical_or.reduce([IDs == ID for ID in IDs_nan])
    index = get_bbox_index([times, lons, lats], [daterange, bbox[:2], bbox[2:]])
    index = np.logical_and(index, ~index_nan)
    IDs, times, lons, lats = IDs[index], times[index], lons[index], lats[index]
    IDs, times, lons, lats = IDs[:10], times[:10], lons[:10], lats[:10]
    
    return data, IDs, times, lons, lats


def get_drifter_track(ID, data, theta=29):
    """
    """
    
    # Get drifter coordinates, deltas and midpoint values
    coords = get_drifter_coords(ID, data)
    deltas = [np.diff(coord) for coord in coords]
    times, lons, lats = [coord[:-1] + delta / 2 for coord, delta in zip(coords, deltas)]
    
    # Estimate velocities
    dt = np.array([delta.total_seconds() for delta in deltas[0]])
    scales = [np.cos(np.deg2rad(lats)), 1]
    u, v = [dx * 111000 * scale / dt for dx, scale in zip(deltas[1:], scales)]
    
    # Rotate velocities to NEMO grid
    theta_rad = math.radians(theta)
    cos, sin = np.cos(theta_rad), np.sin(theta_rad)
    R = np.array([[cos, sin], [-sin, cos]])
    u, v = R.dot([u, v])
    
    return times, lons, lats, u, v


def get_PSF_cast(ID, data):
    """
    """
    
    fieldnames = ['datetime_utc', 'long_deg', 'lat_deg', 'pressure_dbar', 'temp_degc', 'salinity_teos10']
    index, variables = data['ID'] == ID, []
    for name in fieldnames:
        dtype = 'datetime64[s]' if name == 'datetime_utc' else float
        variables.append(data[name][index].values.astype(dtype))
    variables[0] = variables[0].astype(datetime)
    
    return variables


def get_DFO_cast(ID, data):
    """
    """
    
    fieldnames = [
        ['Z'],                                                                    # depth
        ['Temperature_Primary_CT', 'Temperature_Secondary_CT', 'Temperature_CT'], # temperature
        ['Salinity_T0_C0_SA', 'Salinity_T1_C1_SA', 'Salinity_SA'],                # salinity
    ]
    index, index_nan, variables = data['StationTBLID'] == ID, True, []
    for names in fieldnames:
        for name in names:  # Try all field names since they vary between separate casts
            notnan = data[name][index].notna()
            if any(notnan):
                profile = data[name][index].values.astype(float)
                break
        index_nan = np.logical_and(index_nan, notnan)
        variables.append(profile)

    # Trim NaNs
    variables = [var[index_nan] for var in variables]
    
    return variables


def interpolate_NEMO(data, coords, gridvars, gridref, depths=None, size=1):
    """Interpolate NEMO results at a time,lon,lat point and in the vertical if
    depths is provided. Uses Xarray dimensional interpolation for time and depth,
    and Scipy griddata for lon and lat. The size variable defines a bounding box
    to minimize the cost of griddata.
    """
    
    # Dimensions
    n, m = gridvars['lonlat'][0].shape
    
    # Loop through coords
    values = []
    for time, lon, lat in zip(*coords):
        
        # Find j,i from lon,lat
        j, i = [gridref[var].sel(lats=lat, lons=lon, method='nearest').item() for var in ('jj', 'ii')]

        # Return NaN if j,i is outside of domain
        if any([j < size, i < size, j >= n - size, i >= m - size]):
            interpolated = np.ones(len(depths)) * np.nan if depths is not None else np.nan
        
        else:
            # Extract coordinate-results pairs inside bbox
            jslc, islc = [slice(coord - size, coord + size + 1) for coord in (j, i)]
            points = np.vstack([coord[jslc, islc].ravel() for coord in gridvars['lonlat']]).T
            dataarray = data.isel(y=jslc, x=islc).interp(time_counter=time).where(gridvars['mask'][..., jslc, islc])

            # 3D interpolation (e.g. CTD data)
            if depths is not None:
                dataarray = dataarray.interpolate_na(dim='deptht', limit=2, fill_value='extrapolate')
                interpolated = []
                for depth in depths:
                    try:
                        raw = dataarray.interp(deptht=depth, kwargs={'fill_value': 'extrapolate'}).values.ravel()
                        interpolated.append(float(interpolate.griddata(points, raw, (lon, lat))))
                    except:
                        interpolated.append(np.nan)

            # 2D interpolation (e.g. drifters, ferry data)
            else:
                raw = dataarray.values.ravel()
                interpolated = float(interpolate.griddata(points, raw, (lon, lat)))
        
        # Append
        values.append(interpolated)
    
    return np.hstack(values)


def process_NEMO_obs_matching(
    runID, category, startdate, enddate,
    bbox=[-125.4, -122.5, 48.1, 50.5],
    root_nemo='/project/def-allen/bmoorema/results/Currents',
    root_save='/scratch/bmoorema/evaluation',
):
    """Process NEMO-obs matching for CATEGORY, STARTDATE, ENDDATE
    (STARTDATE, ENDDATE must match the results file path)
    """
    
     # Define paths
    daterange = [parse(d) for d in (startdate, enddate)]
    datestr = f'{startdate}_{enddate}'
    rundir = f'SalishSeaCast_currenttuning_{runID}_{datestr}'
    path_nemo = os.path.join(root_nemo, rundir, f'SalishSea_1h_{datestr}')
    file_save = f'SSC{category}_{runID}_{datestr}'
    path_save = os.path.join(root_save, category, f'{file_save}.csv')
    
    # Open stdout files and print status
    sys.stdout = open(os.path.join(root_save, category, 'stdout', f'{file_save}.stdout'), 'w')
    print('Processing results from ' + rundir)
    sys.stdout.flush()
    
    # Load drifter data
    if category == 'drifters':
        names = ['time', 'longitude', 'latitude', 'u_obs', 'v_obs', 'u_nemo', 'v_nemo']
        data_obj, IDs, times = load_drifter_data(daterange, bbox=bbox)
        gridvars, gridref = get_grid_variables()
        NEMO = {}
        for var, name in zip(['u', 'v'], ['vozocrtx', 'vomecrty']):
            print(f'Loading {var} ...')
            sys.stdout.flush()
            fn = f'{path_nemo}_grid_{var.upper()}.nc'
            NEMO[var] = xr.open_dataset(fn)[name].isel({f'depth{var}': 0}).sel(time_counter=times).load()
    
    # Load cruise data
    elif category in ['PSF', 'DFO']:
        names = ['time', 'longitude', 'latitude', 'depth', 'T_obs', 'S_obs', 'T_nemo', 'S_nemo']
        if category == 'PSF':
            data_obj, IDs = load_PSF_data(daterange, bbox=bbox)
        elif category == 'DFO':
            data_obj, IDs, times, lons, lats = load_DFO_data(daterange, bbox=bbox)
        gridvars, gridref = get_grid_variables(variables=['t'])
        gridvars = gridvars['t']
        NEMO = xr.open_dataset(f'{path_nemo}_grid_T.nc')
    
    # Unknown data category
    else:
        raise ValueError(f'Unrecognized category: {category}')
    
    # Loop through sampling IDs
    data, n = {name: [] for name in names}, len(IDs)
    coordinate_list = [times, lons, lats] if category == 'DFO' else [range(n)]
    for k, ID, coordinates in zip(range(n), IDs, zip(*coordinate_list)):
    
        # Drifters
        if category == 'drifters':
            
            # Get drifter track observations
            times, lons, lats, u_obs, v_obs = get_drifter_track(ID, data_obj)
            u_nemo, v_nemo = [interpolate_NEMO(NEMO[name], (times, lons, lats), gridvars[name], gridref) for var in ('u', 'v')]
            
            # Append values
            values = [times, lons, lats, u_obs, v_obs, u_nemo, v_nemo]
            for name, value in zip(names, values):
                data[name].append(value)
        
        # Cruises
        elif category in ['PSF', 'DFO']:
            
            # PSF
            if category == 'PSF':
                times, lons, lats, depths, T_obs, S_obs = get_PSF_cast(ID, data_obj)
            
            # DFO
            elif category == 'DFO':
                depths, T_obs, S_obs = get_DFO_cast(ID, data_obj)
                ncast = len(depths)
                times, lons, lats = [np.repeat(coord, ncast) for coord in coordinates]
                
            # Get interpolated NEMO results
            coords = [[var[0]] for var in (times, lons, lats)]
            T_nemo, S_nemo = [interpolate_NEMO(NEMO[var], coords, gridvars, gridref, depths=depths) for var in ('votemper', 'vosaline')]
            
            # Append values
            values = [times, lons, lats, depths, T_obs, S_obs, T_nemo, S_nemo]
            for name, value in zip(names, values):
                data[name].append(value)
        
        # Unknown data category
        else:
            raise ValueError(f'Unrecognized category: {category}')
        
        # Print status
        loopstatus(k, n)
        
    # Concatenate
    for name in names:
        data[name] = np.hstack(data[name])
    
    # Close out
    pd.DataFrame(data).to_csv(path_save, index=False)
    print('Done!')
    sys.stdout.close()


if __name__ == "__main__":
    process_NEMO_obs_matching(*[str(arg) for arg in sys.argv[1:]])
