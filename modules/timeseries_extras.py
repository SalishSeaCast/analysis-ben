from salishsea_tools import timeseries_tools, utilities
from scipy.io import savemat
from datetime import timedelta
from dateutil.parser import parse
import numpy as np
import xarray as xr
import os


def reshape_GEM(data, grid, mask):
    """
    """

    # Build GEM coordinates
    data_flat = {'raw': {'grid': {}}}
    mask_flt, coords, ngrid, ngrid_water = timeseries_tools.reshape_coords_GEM(
        grid, mask,
    )

    # Reshape, flatten, and trim
    data_flat['raw']['grid']['u'] = timeseries_tools.reshape_to_ts(
        data.u_wind.values, mask_flt.astype(bool), ngrid, ngrid_water,
    )
    data_flat['raw']['grid']['v'] = timeseries_tools.reshape_to_ts(
        data.v_wind.values, mask_flt.astype(bool), ngrid, ngrid_water,
    )

    # Assign parameters to output dict
    data_flat['mask'] = mask_flt
    data_flat['coords'] = coords
    data_flat['ngrid'] = ngrid
    data_flat['ngrid_water'] = ngrid_water

    return data_flat


def iterate_NEMO_timeseries(
        timerange, variables, mask, dims, indices, datadir='nowcast2016',
        datapath='/ocean/bmoorema/research/MEOPAR/analysis-ben/data',
):
    """
    """

    # Loop through results slices
    for var in variables:  # --- Loop through variables

        for dim, index in zip(dims, indices):  # --- Loop through slices

            # Initialize dictionaries
            coords = {}
            data = {}

            # Determine filenames and unstaggering from variable names
            dim_in, unstagger = dim, None
            if var is 'vozocrtx':
                file = 'U'
                if dim is 'depth':
                    dim_in, unstagger = 'depthu', 'x'
            elif var is 'vomecrty':
                file = 'V'
                if dim is 'depth':
                    dim_in, unstagger = 'depthv', 'y'
            else:
                file = 'T'
                if dim is 'depth':
                    dim_in = 'deptht'

            # Determine spacing from variable names and slicing
            if (var is 'vozocrtx' or var is 'vomecrty') and dim is 'depth':
                spacing = 1
            else:
                spacing = 1

            # Load Results
            filenames = timeseries_tools.make_filename_list(
                timerange, file, model='nowcast-green', resolution='h',
            )
            data, coords = timeseries_tools.load_NEMO_timeseries(
                filenames, mask, var, dim_in, index=index,
                spacing=spacing, shape='flat', unstagger_dim=unstagger,
            )

            # Export timerange
            coords['timerange'] = timerange

            # Save model results
            savemat(os.path.join(
                datapath, datadir, var, f'{var}_{dim}{index}.mat'), data)
            savemat(os.path.join(
                datapath, datadir, var, f'coords_{dim}{index}.mat'), coords)


def load_hindcast_timeseries_location(
    varlist, locs, daterange, loadpath, writepath,
    ftype='grid_T', res='h', date_cutoff=None, loadpath_cutoff=None,
):
    """Loads a list of NEMO hindcast variables from the daily file
    record at specified point locations over a specified date range
    and writes the resulting timeseries to a netCDF file.
    """

    # Parse daterange
    dates = [parse(date) for date in daterange]

    # Create dict to hold timeseries at locs
    data = {'time': np.empty(0, dtype='datetime64')}
    for key in locs:
        data[key] = {}
        for var in varlist:
            data[key][var] = np.empty(0)

    # Loop through all hindcast hourly files
    bar = utilities.statusbar('Loading NEMO record ...')
    for day in bar(range(np.diff(dates)[0].days)):

        # Parse date info
        date = dates[0] + timedelta(days=day)
        datestr = date.strftime('%Y%m%d')

        # Consider date cutoff for path switch
        path = loadpath
        if date_cutoff is not None:
            if date >= parse(date_cutoff):
                path = loadpath_cutoff

        # Load timeseries at points
        fn = os.path.join(
            path, date.strftime('%d%b%y').lower(),
            f'SalishSea_1{res}_{datestr}_{datestr}_{ftype}.nc',
        )
        with xr.open_dataset(fn) as ds:
            data['time'] = np.concatenate(
                [data['time'], ds.time_counter.values],
            )
            for key in locs:
                for var in varlist:
                    data[key][var] = np.concatenate([
                        data[key][var],
                        ds[var][(slice(None),) + locs[key][::-1]].values,
                    ])

    # Save to xarray dataset
    prefix = os.path.split(loadpath)[1].replace('.', '')
    dates = [date.strftime('%Y%m%d') for date in dates]
    for key in locs:
        for var in varlist:
            data[key][var] = ('time', data[key][var])
        fn = os.path.join(
            writepath, '_'.join([prefix] + dates + [key] + varlist) + '.nc',
        )
        xr.Dataset(data[key], coords={'time': data['time']}).to_netcdf(fn)

    return
