#!/usr/bin/env python

import numpy as np
import xarray as xr
import pandas as pd
import sqlite3 as sql
import math
import os
import sys
from datetime import datetime, timedelta
from dateutil.parser import parse
from scipy import io, interpolate


def loopstatus(k, n, interval=10):
    """Print loop progress at percentage intervals given an iteration k
    and a total number of iterations n
    """
    
    nscale = 100 / n
    percent = k * nscale
    if percent % interval < nscale:
        print(f'{int(percent)}% complete ...')
        sys.stdout.flush()


def get_bbox_index(coords, bbox):
    """Get the indices for a list of given coordinate arrays that fall
    within the corresponding boundary pairs specified by the bbox list
    """
    
    index = True
    for coord, bounds in zip(coords, bbox):
        for bound, func in zip(bounds, ['greater_equal', 'less_equal']):
            if bound is not None:
                index = np.logical_and(index, getattr(np, func)(coord, bound))
    
    return index


def get_drifter_coords(ID, data):
    """Get the times, lons, lats coordinate arrays from the drifter data object
    for a given drifter track specified by ID
    """
    
    index = np.where(data['id'].astype(int) == ID)[0][0]
    times, lats, lons = [data[var][index].squeeze() for var in ('mtime', 'lat', 'lon')]
    times = np.array([datetime.fromordinal(int(t) - 366) + timedelta(days=t%1) for t in times])
    
    return times, lats, lons


def get_grid_variables(gridtypes=['u', 'v'], depth=False, path='/home/bmoorema/MEOPAR/grid/'):
    """Return grid variables dict GRIDVARS with LONLAT and MASK fields for U and V.
    Also return the GRIDREF xarray object for looking up ji from lonlat.
    If depth=True, also include the 3D depth field at rest from the meshmask file.
    Extra variables necessary to account for VVL are included in comments, to be
    incorporated later...
    """

    # Master files
    coords = xr.open_dataset(path + 'coordinates_seagrid_SalishSea201702.nc', decode_times=False)
    mask = xr.open_dataset(path + 'mesh_mask202108.nc')
    gridref = xr.open_dataset(path + 'grid_from_lat_lon_mask999.nc')

    # Grid variables
    gridvars = {}
    for gridtype in gridtypes:
        
        # Extract latlon coordinates
        latlon = [coords[name + gridtype].values[0, ...] for name in ('gphi', 'glam')]
        
        # Extract depth-dependent grid variables
        if depth:
            names = [f'{gridtype}mask', f'gdep{gridtype}'] #, f'e3{gridtype}_0']
            if gridtype == 't': names[1] = names[1] + '_0'
            meshmask, depth = [mask[name][0, ...].values for name in names] #, e30
            # H = np.sum(e30 * meshmask, axis=0)    needed if using vvl calculation
            gridvars[gridtype] = {'latlon': latlon, 'mask': meshmask, 'depth': depth}
        
        # Only use latlon and surface meshmask
        else:
            meshmask = mask[f'{gridtype}mask'][0, 0, ...].values
            gridvars[gridtype] = {'latlon': latlon, 'mask': meshmask}
    
    # Return inner dict if only a single gridtype was requested
    if len(gridvars) == 1:
        gridvars = next(iter(gridvars.values()))
    
    return gridvars, gridref


def load_drifter_data(
    daterange, bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='Salish_L3_20190728T103529.mat',
):
    """Load the drifter dataset from a matfile and construct a list
    of drifter IDs that fall inside the daterange and bbox criteria
    specified.
    """
    
    # Load data from matfile
    data = io.loadmat(os.path.join(path, file))['drift'][0]
    
    # Get drifter IDs within daterange-bbox
    # and times on the half-hour for loading NEMO results into memory
    IDs, times = list(np.unique(data['id'].astype(int))), []
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
    """Load the PSF citizen science dataset from .csv, removing values
    that fall outside the daterange and bbox criteria, and also removing
    rows that have entirely nan coordinates.
    """
    
    # Load data from csv
    data = pd.read_csv(os.path.join(path, file), parse_dates=[3])
    
    # Get coords
    names, dtypes = ['datetime_utc', 'lat_deg', 'long_deg'], ['datetime64[s]', float, float]
    times, lats, lons = [data[name].values.astype(dtype) for name, dtype in zip(names, dtypes)]
    times = times.astype(datetime)

    # Apply daterange-bbox criteria and filter NaNs from lon and lat
    index = get_bbox_index([times, lats, lons], [daterange, bbox[:2], bbox[2:]])
    index = np.logical_and.reduce((index, ~np.isnan(lats), ~np.isnan(lons)))
    data = data[index]

    # Append station string to ID column to make unique identifiers
    data['ID'] = data['ID'].astype(str) + data['station']
    IDs = np.unique(data['ID'])
    
    return data, IDs


def load_DFO_data(
    daterange, bbox=[None, None, None, None],
    path='/home/bmoorema/project/SalishSea/obs',
    file='DFO_CTD.sqlite',
):
    """Load the DFO SoG cruise data record from SQL database and construct
    a list of cast IDs that fall inside the daterange and bbox criteria
    specified. Also remove IDs predetermined to contain entirely nan data.
    """
    
    # Load data from sql
    con = sql.connect(os.path.join(path, file))
    coords = pd.read_sql_query('SELECT * from StationTBL', con)
    data = pd.read_sql_query('SELECT * from CalcsTBL', con)
    con.close()
    
    # Get coords
    names = ['ID', 'Lat', 'Lon', 'StartYear', 'StartMonth', 'StartDay', 'StartHour']
    IDs, lats, lons, y, m, d, h = [coords[name].values for name in names]
    y, m, d = [var.astype(int) for var in (y, m, d)]
    times = np.array([datetime(*args[:3]) + timedelta(hours=args[3]) for args in zip(y, m, d, h)])
    
    # Apply daterange-bbox criteria and filter all NaN casts using predetermined IDs
    IDs_nan = [790, 1147, 1148, 1999, 6226] + list(range(2847, 2853))
    index_nan = np.logical_or.reduce([IDs == ID for ID in IDs_nan])
    index = get_bbox_index([times, lats, lons], [daterange, bbox[:2], bbox[2:]])
    index = np.logical_and(index, ~index_nan)
    IDs, times, lats, lons = IDs[index], times[index], lats[index], lons[index]
    
    return data, IDs, times, lats, lons


def get_drifter_track(ID, data, theta=29):
    """Get the GPS time and position arrays for drifter ID and
    estimate the Eulerian velocity using a centered difference
    between two points.
    """
    
    # Get drifter coordinates, deltas and midpoint values
    coords = get_drifter_coords(ID, data)
    deltas = [np.diff(coord) for coord in coords]
    times, lats, lons = [coord[:-1] + delta / 2 for coord, delta in zip(coords, deltas)]
    
    # Estimate velocities
    dt = np.array([delta.total_seconds() for delta in deltas[0]])
    scales = [1, np.cos(np.deg2rad(lats))]
    v, u = [dx * 111000 * scale / dt for dx, scale in zip(deltas[1:], scales)]
    
    # Rotate velocities to NEMO grid
    theta_rad = math.radians(theta)
    cos, sin = np.cos(theta_rad), np.sin(theta_rad)
    R = np.array([[cos, sin], [-sin, cos]])
    u, v = R.dot([u, v])
    
    # Return variables as list
    variables = [times, lats, lons, u, v]
    
    return variables


def get_PSF_cast(ID, data):
    """Get the coordinate, temperature and salinity arrays for cast ID
    from the PSF dataset.
    """
    
    # Define fieldnames and extract cast from data using ID
    fieldnames = ['datetime_utc', 'lat_deg', 'long_deg', 'pressure_dbar', 'temp_degc', 'salinity_teos10']
    index, variables = data['ID'] == ID, []
    for name in fieldnames:
        dtype = 'datetime64[s]' if name == 'datetime_utc' else float
        variables.append(data[name][index].values.astype(dtype))
    variables[0] = variables[0].astype(datetime)
    
    return variables


def get_DFO_cast(ID, data):
    """Get the coordinate, temperature and salinity arrays for cast ID
    from the DFO dataset.
    """
    
    # Define fieldnames and extract cast from data using ID
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
    depths is provided. Uses Xarray dimensional interpolation for time,
    and Scipy griddata for lat, lon and depth. The size variable defines a
    bounding box to minimize the cost of griddata. Depth is interpolated
    using the 3D gridded depth field from the meshmask file, VVL can be
    accounted for by calculating the time-dependent depth scale from SSH
    and e30 (to be implimented ...)
    """
    
    # Find j,i from gridref lookup and create bbox slices
    ji = gridref.sel(lats=coords[1], lons=coords[2], method='nearest')
    ji = [ji[name].item() for name in ('jj', 'ii')]
    
    # Failure output to nan
    nan = np.ones(len(depths)) * np.nan if depths is not None else np.nan
    
    # Return nan if latlon is outside model points
    if -999 in ji:
        return nan
    
    # Prepare interpolation points on common grid
    jslc, islc = [slice(coord - size, coord + size + 1) for coord in ji]
    mask = gridvars['mask'][..., jslc, islc].astype(bool).ravel()
    points = [coord[jslc, islc] for coord in gridvars['latlon']]
    xi = coords[1:]
    
    # Add depths
    if depths is not None:
        depth = gridvars['depth'][:, jslc, islc] # multiply by (1 + ssh/H for vvl)
        points = [depth] + [np.broadcast_to(point, depth.shape) for point in points]
        xi = [depths] + [np.broadcast_to(x, len(depths)) for x in xi]
    
    # Define points, values and xi
    points, xi = tuple(point.ravel()[mask] for point in points), tuple(xi)
    values = data.isel(y=jslc, x=islc).interp(time_counter=coords[0]).values.ravel()[mask]
    
    # Interpolate using `scipy.griddata`
    try:
        interpolated = interpolate.griddata(points, values, xi)
    except:
        interpolated = nan
    
    return interpolated


def process_NEMO_obs_matching(
    runID, category, startdate, enddate,
    bbox=[48.1, 50.5, -125.4, -122.5],
    root_nemo='/project/def-allen/bmoorema/results/Currents',
    root_save='/scratch/bmoorema/evaluation',
):
    """Process NEMO-obs matching for RUNID, CATEGORY, STARTDATE, ENDDATE
    (STARTDATE, ENDDATE must match the results file path)
    """
    
     # Define strings and paths
    daterange = [parse(d) for d in (startdate, enddate)]
    datestr = f'{startdate}_{enddate}'
    rundir = f'SalishSeaCast_currenttuning_{runID}_{datestr}'
    evalID = f'SSC{category}_{runID}_{datestr}'
    path_nemo = os.path.join(root_nemo, rundir, f'SalishSea_1h_{datestr}')
    path_save = os.path.join(root_save, category, f'{evalID}.csv')
    path_stdout = os.path.join(root_save, category, 'stdout', f'{evalID}.stdout')
    
    # Open stdout files and print status
    sys.stdout = open(path_stdout, 'w')
    print('Processing results from ' + rundir)
    sys.stdout.flush()
    
    # Load drifter data
    # Loads u and v into memory at times returned by `load_drifter_data`
    if category == 'drifters':
        names = ['time', 'longitude', 'latitude', 'u_obs', 'v_obs', 'u_nemo', 'v_nemo']
        data_obj, IDs, times = load_drifter_data(daterange, bbox=bbox)
        gridvars, gridref = get_grid_variables(gridtypes=['u', 'v'], depth=False)
        NEMO, tslc = {}, {'time_counter': times}
        for gridtype, varname in zip(['u', 'v'], ['vozocrtx', 'vomecrty']):
            print(f'Loading {gridtype} ...')
            sys.stdout.flush()
            fn = f'{path_nemo}_grid_{gridtype.upper()}.nc'
            zlsc = {f'depth{gridtype}': 0}
            NEMO[gridtype] = xr.open_dataset(fn)[varname].isel(zlsc).sel(tslc).load()
    
    # Load cruise data
    # Does not load tracers into memory
    elif category in ['PSF', 'DFO']:
        names = ['time', 'longitude', 'latitude', 'depth', 'T_obs', 'S_obs', 'T_nemo', 'S_nemo']
        if category == 'PSF':
            data_obj, IDs = load_PSF_data(daterange, bbox=bbox)
        elif category == 'DFO':
            data_obj, IDs, times, lats, lons = load_DFO_data(daterange, bbox=bbox)
        gridvars, gridref = get_grid_variables(gridtypes=['t'], depth=True)
        NEMO = xr.open_dataset(f'{path_nemo}_grid_T.nc')
    
    # Unknown data category
    else:
        raise ValueError(f'Unrecognized category: {category}')
    
    # Loop through sampling IDs
    data, n = {name: [] for name in names}, len(IDs)
    coordinate_list = [times, lats, lons] if category == 'DFO' else [range(n)]
    for k, ID, coordinates in zip(range(n), IDs, zip(*coordinate_list)):
    
        # Drifters
        if category == 'drifters':
            
            # Get drifter track observations
            variables = get_drifter_track(ID, data_obj)
            
            # Get interpolated NEMO results for each drifter point
            for coords in zip(*variables[:3]):
                for name, gridtype in zip(['u_nemo', 'v_nemo'], ['u', 'v']):
                    interpolated = interpolate_NEMO(NEMO[gridtype], coords, gridvars[gridtype], gridref, size=2)
                    data[name].append(interpolated)
        
        # Cruises
        elif category in ['PSF', 'DFO']:
            
            # PSF
            if category == 'PSF':
                variables = get_PSF_cast(ID, data_obj)
                coords = [var[0] for var in variables[:3]]
            
            # DFO
            elif category == 'DFO':
                variables = get_DFO_cast(ID, data_obj)
                coords, ncast = coordinates, len(variables[0])
                variables = [np.broadcast_to(coord, ncast) for coord in coords] + variables
                
            # Get interpolated NEMO results for whole cast
            for name, varname in zip(['T_nemo', 'S_nemo'], ['votemper', 'vosaline']):
                interpolated = interpolate_NEMO(NEMO[varname], coords, gridvars, gridref, depths=variables[3])
                data[name].append(interpolated)
        
        # Unknown data category
        else:
            raise ValueError(f'Unrecognized category: {category}')
        
        # Append obs variables
        for name, var in zip(names, variables):
            data[name].append(var)
        
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