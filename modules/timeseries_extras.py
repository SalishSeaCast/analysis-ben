from salishsea_tools import timeseries_tools
from scipy.io import savemat
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
