import numpy as np
import xarray as xr
import gsw
import os
from salishsea_tools import grid_tools, viz_tools


def load_results(t, rundir, prefix, maskfile, x_range=[None], y_range=[None]):
    """
    """

    # Load netCDF files
    T = xr.open_dataset(os.path.join(rundir, f'{prefix}_grid_T.nc'))
    U = xr.open_dataset(os.path.join(rundir, f'{prefix}_grid_U.nc'))
    V = xr.open_dataset(os.path.join(rundir, f'{prefix}_grid_V.nc'))
    mask = xr.open_dataset(maskfile)

    # Define slices
    xslice = slice(*x_range)
    yslice = slice(*y_range)

    # Define parameters
    params = {
        'x': mask.x[xslice].values,
        'y': mask.y[yslice].values,
        'u': U.vozocrtx[t, :, yslice, xslice].values,
        'v': V.vomecrty[t, :, yslice, xslice].values,
        'eta': T.sossheig[t, yslice, xslice].values,
        'tmask': mask.tmask[0, :, yslice, xslice].values,
        'gdept_0': mask.gdept_0[0, :, yslice, xslice].values,
        'e1t': mask.e1t[0, yslice, xslice].values,
        'e2t': mask.e2t[0, yslice, xslice].values,
        'e3t_0': mask.e3t_0[0, :, yslice, xslice].values,
    }

    # Obtain the time dependent grid parameters
    VVL = grid_tools.calculate_time_dependent_grid(
        params['e3t_0'], params['tmask'], params['eta'][np.newaxis, ...],
        {
            'e3t_t': params['e3t_0'][np.newaxis, ...],
            'gdept_t': params['gdept_0'][np.newaxis, ...],
        },
    )

    # Calculate rho
    rho = gsw.rho(
        T.vosaline[t, :, yslice, xslice],
        T.votemper[t, :, yslice, xslice],
        VVL['gdept_t'][0, ...],
    )

    # Finalize output dict
    params.update({
        'rho': rho,
        'e3t_t': VVL['e3t_t'][0, ...],
        'gdept_t': VVL['gdept_t'][0, ...],
    })

    return params


def calc_geostrophic_velocities(z, params):
    """Calculate pressure gradients and geostrophic velocities at idepth
    """

    # Geostrophic parameters dict
    GEO = {}

    # Define constants
    g = 9.81
    f = 1e-4

    # Calculate z surface displacement xi
    stretching = (
        params['e3t_t'][:z, ...].sum(axis=0) -
        params['e3t_0'][:z, ...].sum(axis=0)
    )
    xi = params['eta'] - stretching

    # Define rho at z surface
    rho_bot = params['rho'][z, ...]
    rho_bot[xi < 0] = params['rho'][z - 1, ...][xi < 0]

    # Calculate pressure
    GEO['pressure'] = g * (
        (params['rho'][:z, ...] * params['e3t_t'][:z, ...]).sum(axis=0) +
        rho_bot * xi
    )

    # Extract and unstagger the model velocity fields
    GEO['u'] = params['u'][z, ...]
    GEO['v'] = params['v'][z, ...]
    GEO['u'][1:, 1:], GEO['v'][1:, 1:] = viz_tools.unstagger(
        GEO['u'], GEO['v'],
    )

    # Calculate the pressure gradient
    dpdy, dpdx = np.gradient(GEO['pressure'], axis=(0, 1))

    # Calculate the geostrophic velocities
    GEO['u_g'] = -1 / (f * params['rho'][z, ...]) * dpdy / params['e2t']
    GEO['v_g'] = 1 / (f * params['rho'][z, ...]) * dpdx / params['e1t']

    # Calculate ageostrophic velocities
    GEO['u_a'] = GEO['u'] - GEO['u_g']
    GEO['v_a'] = GEO['v'] - GEO['v_g']

    return GEO
