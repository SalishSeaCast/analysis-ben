# --- Postprocessing code for SalishSeaCast NEMO upwelling results

import numpy as np
import xarray as xr
import os
import gsw
from datetime import datetime, timedelta
from dateutil.parser import parse
from scipy import signal
from scipy.interpolate import interp1d
from tqdm import tqdm_notebook as tqdm
from salishsea_tools import geo_tools
from dynmodes import dynmodes
from warnings import simplefilter

# Filter invalid warnings
simplefilter('ignore')


def calc_rho(data, depth, tmask):
    """Calculate the density, rho
    """

    tracers = ['vosaline', 'votemper']
    rho = gsw.rho(
        *[np.ma.masked_where(tmask == 0, data[k]) for k in tracers], depth,
    )

    return rho


def make_prefix(date, path, res='h'):
    """Construct path prefix for local SalishSeaCast results given date object
    and paths dict. e.g.,
    /results/SalishSea/hindcast.201812/ddmmmyy/SalishSea_1h_yyyymmdd_yyyymmdd
    """

    datestr = '_'.join(np.repeat(date.strftime('%Y%m%d'), 2))
    prefix = os.path.join(
        path, date.strftime('%d%b%y').lower(), f'SalishSea_1{res}_{datestr}',
    )

    return prefix


def calc_coastline_indices(mask, bathy, e1t=440):
    """Find the NEMO-SalishSeaCast j, i pairs along the western coastline.
    Also finds the coastline angle CW from north-south, and the average cross-
    shore bottom slope between zero and 200 m depth. I would prefer this
    routine to be flexible to various choices of jmin, jmax and jsub but in
    reality it takes a bit of tuning to avoid errors.
    """

    # Indexing parameters
    imin, jmin, jmax, jsub = 115, 365, 745, 5
    winlen = 15

    # Initialize index array and seed deep landpoint index
    index = np.empty(0, dtype=int)
    i_d_prev = 200

    # Iterate through jindex
    for j in range(jmin-winlen, jmax+winlen):

        # Extract surface and deep mask rows
        i_s = np.where(abs(np.diff(mask.tmask[0, 0, j, imin:])) > 0)[0]
        i_d = np.where(abs(np.diff(mask.tmask[0, 29, j, imin:])) > 0)[0]

        # Remove discontinuous deep land points (i.e., Gulf Island basins)
        while abs(i_d[0] - i_d_prev) > 30: i_d = i_d[1:]
        i_d_prev = i_d[0]

        # Append coastline point to array
        index = np.append(index, i_s[i_s < i_d[0]][-1] + imin)

    # Find coastline angle from smoothed coastline index
    window = signal.get_window('blackman', 2*winlen+1)
    angle = np.arctan(
        np.diff(signal.convolve(index, window / sum(window), mode='same')),
    )

    # Populate sections dict
    sections = {
        'n': int((jmax - jmin) / jsub),
        'ji': [np.arange(jmin, jmax, jsub), index[winlen:-winlen:jsub]],
        'angle': angle[winlen:-winlen:jsub],
    }

    # Calculate bottom slope
    bottom = np.empty((0, 50))
    for j, i in zip(*sections['ji']):
        bottom = np.concatenate(
            (bottom, bathy[j, i:i+50].values[np.newaxis, :]),
        )
    bottom[bottom > 200] = 200
    sections['slope'] = np.nanmax(bottom, axis=1) / \
      (np.cos(sections['angle']) * np.nanargmax(bottom, axis=1) * e1t)

    return sections


def calc_coastline_indices_HRDPS(sections, mask, grid_HRDPS):
    """Find the HRDPS j, i pairs along the western coastline given the
    NEMO-SalishSeaCast j, i pairs.
    """

    # Initialize index arrays and define search tolerances
    sections['ji_HRDPS'] = [np.empty(0, dtype=int), np.empty(0, dtype=int)]
    tols = {
        'NEMO'  : {'tol_lon': 0.010, 'tol_lat': 0.003},
        'GEM2.5': {'tol_lon': 0.017, 'tol_lat': 0.017},
    }

    # Loop through NEMO coastline indices
    for j, i in zip(*sections['ji']):

        # Search for nearest lon/lat neighbors
        jj, ii = geo_tools.find_closest_model_point(
            mask.glamt[0, j, i].values, mask.gphit[0, j, i].values,
            grid_HRDPS.longitude-360, grid_HRDPS.latitude,
            grid='GEM2.5', tols=tols,
        )

        # Append new values to indices
        for dim, val in zip([0, 1], [jj, ii]):
            sections['ji_HRDPS'][dim] = np.append(
                sections['ji_HRDPS'][dim], val,
            )

    return sections


def calc_stratification_parameters(
    rho, mask, h_ref=None, mode='p', L=50, const={'g': 9.81, 'rho_0': 1e3},
):
    """
    """

    # Extract deptht array
    deptht = mask.gdept_1d[0, :].values
    depthw = mask.gdepw_1d[0, :].values
    e3t = mask.e3t_1d[0, :].values

    # Find average N2 profile on deptht grid within L of coast
    interp_rho = interp1d(deptht, rho, axis=0, fill_value='extrapolate')
    N2 = const['g'] / const['rho_0'] * np.diff(interp_rho(depthw), axis=0) / \
      np.expand_dims(e3t[:-1], axis=1)
    N2 = np.nanmedian(N2[:, :L], axis=1)
    nanindex = np.isnan(N2)
    N2 = N2[~nanindex]

    # Calculate pycnocline depth
    if h_ref is None:

        # Calculate pycnocline depth using vertical mode calculator dynmodes
        zlim = min(len(N2), 24)
        _, pmode, rmode, _ = dynmodes(N2[:zlim], deptht[:zlim], 1)
        if mode is 'r':
            zindex = int(abs(rmode[0, :]).argmax())
        elif mode is 'p':
            zindex = int(np.where(np.diff(np.signbit(pmode[0, :])))[0])
        else:
            raise ValueError(f"Can't interpret mode type: {mode}")
        h_s = deptht[zindex]

    else:

        # Use z_ref as pycnocline depth
        zindex = abs(deptht - h_ref).argmin(axis=0)
        h_s = h_ref

    rho_s = (rho[:zindex, :L].mean(axis=1) * e3t[:zindex]).sum() / \
      depthw[zindex]
    rho_ref = rho[zindex, :L].mean()
    N_int = (np.sqrt(N2) * e3t[:-1][~nanindex])[zindex:31].sum()
    N_s = (np.sqrt(N2)[:zindex] * e3t[:zindex]).sum() / depthw[zindex]

    return h_s, rho_s, rho_ref, N_int, N_s


def calc_upwelling_parameters(rho, rho_t0, rho_ref, deptht, e1t=440, angle=0, L=50):
    """
    """

    # Calculate upwelling parameters
    rho_max = rho[:L].max()
    h_u = deptht[int(np.median(abs(rho_t0[:, :L] - rho_max).argmin(axis=0)))]
    x_u = np.cos(angle) * (rho[:L].compressed() >= rho_ref).sum() * e1t

    return h_u, x_u


def process_idealized_results(
    param, sections, mask, subgrid, hour=24, L=50,
    keys=[], runs={}, vals={}, strat='2layer',
    path='/data/bmoorema/results/Lake/S4d_2layer',
    fn='SalishSeaIdeal_1h_20170701_20170706_grid_T.nc',
):
    """
    """

    # Define subgrid boundaries and deptht array
    isub, jsub, deptht = subgrid[0], subgrid[2], mask.gdept_1d[0, :].values

    # Raise ValueError if no keys list
    if not keys:
        raise ValueError('Must provide list of upwelling parameter keys!')

    # Recursively loop through run parameters
    for val in tqdm(param[keys[0]], leave=False, desc=keys[0]):
        vals[keys[0]] = val

        # Continue recursion and define nested dict fields
        if len(keys) > 1:
            runs[val] = {}
            runs[val] = process_idealized_results(
                param, sections, mask, subgrid, hour=hour, L=L,
                keys=keys[1:], runs=runs[val], vals=vals,
                strat=strat, path=path, fn=fn,
            )

        # Bottom level: stop recursion and process NEMO results
        else:
            runs[val] = {'N_int': [], 'N_s': [], 'h_u': [], 'x_u': []}

            # Define and load results file
            if strat is '2layer':
                runID = f"SalishSeaPond_S4d{vals['u_wind']:02d}ms_" + \
                  f"halocline{vals['h_s']:2d}m_rhosurf{vals['rho_s']:4d}"
                h_ref = vals['h_s']
            elif strat is 'const':
                runID = f"SalishSeaPond_S4d{vals['u_wind']:02d}ms_" + \
                  f"N{vals['N']*1e4:04.0f}s"
                h_ref = 10

            file = os.path.join(path, runID, fn)
            if os.path.exists(file):
                with xr.open_dataset(file) as data:

                    # Loop through sections
                    for j, i, a in zip(*sections['ji'], sections['angle']):

                        # Calculate rho at t=0
                        rho_t0 = calc_rho(
                            data.isel(time_counter=0, y=j-jsub, x=slice(i-isub, None)),
                            np.expand_dims(deptht, axis=1), mask.tmask[0, :, j, i:],
                        )

                        # Calculate stratification parameters
                        h_s, rho_s, rho_ref, N_int, N_s = calc_stratification_parameters(
                            rho_t0, mask, h_ref=h_ref, L=L,
                        )

                        # Calculate rho at t=hour
                        rho = calc_rho(
                            data.isel(time_counter=hour, deptht=0, y=j-jsub, x=slice(i-isub, None)),
                            0, mask.tmask[0, 0, j, i:],
                        )

                        # Calculate upwelling metrics
                        h_u, x_u = calc_upwelling_parameters(
                            rho, rho_t0, rho_ref, deptht, angle=a, L=L,
                        )

                        # Append values to arrays
                        runs[val]['N_int'].append(N_int)
                        runs[val]['N_s'].append(N_s)
                        runs[val]['h_u'].append(h_u)
                        runs[val]['x_u'].append(x_u)

    return runs


def process_hindcast_results(
    events, sections, mask, HRDPS, hour=24, L=50,
    path='/results2/SalishSea/nowcast-green.201905',
):
    """
    """

    # Initialize variables storage, parse dates, and define deptht array
    varnames = ['rho_s', 'z_h', 'N_int', 'N_s', 'tau', 'h_u', 'x_u']
    outnames = ['h_s', 'rho_s', 'rho_ref', 'N_int', 'N_s']
    variables = {}
    for var in varnames: variables[var] = np.empty((0, sections['n']))
    dates = np.array([parse(date) for date in events])
    deptht = mask.gdept_1d[0, :].values

    # Loop through upwelling event dates
    for date in tqdm(dates):

        # Define loop variables
        datestr = date.strftime('%Y%b%d')
        varlist = {'rho_t0': [], 'rho_ref': []}
        for var in varnames: varlist[var] = []

        # Loop through hours
        for hour in [0, hour]:
            time = date + timedelta(hours=hour)

            # Open hindcast record
            with xr.open_dataset(make_prefix(time, path) + '_grid_T.nc') as data:

                # Loop through sections
                for n, j, i, jj, ii, a in zip(
                    range(sections['n']), *sections['ji'],
                    *sections['ji_HRDPS'], sections['angle'],
                ):

                    # Calculations at t=0
                    if hour == 0:

                        # Calculate density
                        rho_t0 = calc_rho(
                            data.sel(time_counter=time, method='nearest').isel(y=j, x=slice(i, None)),
                            np.expand_dims(deptht, axis=1), mask.tmask[0, :, j, i:],
                        )

                        # Calculate stratification parameters
                        outputs = calc_stratification_parameters(rho_t0, mask, L=L)

                        # Append to dict lists
                        varlist['rho_t0'].append(rho_t0)
                        for var, out in zip(outnames, outputs):
                            varlist[var].append(out)

                    # Calculations at hours
                    else:

                        # Calculate tau
                        u, v = [
                            HRDPS.sel(time=slice(*timerange))[k].values[:, jj, slice(ii, ii+10)]
                            for k in ['u_wind', 'v_wind']
                        ]
                        jtau = np.sin(np.arctan2(v, u) - np.pi * (22 / 180) + angle)
                        tau = sum(1.225e-3 * (jtau * np.sqrt(u**2 + v**2)).mean(axis=1)**2) * 3600

                        # Calculate rho at t=hour
                        rho = calc_rho(
                            data.sel(
                                time_counter=timerange[1], method='nearest',
                            ).isel(deptht=0, y=j, x=slice(i, None)),
                            0, mask.tmask[0, 0, j, i:],
                        )

                        # Calculate upwelling parameters
                        h_u, x_u = calc_upwelling_parameters(
                            rho, rho_t0, rho_ref, deptht, angle=a, L=L,
                        )

                        # Append to dict lists
                        varlist['tau'].append(tau)
                        varlist['h_u'].append(h_u)
                        varlist['x_u'].append(x_u)

        # Concatenate lists onto storage arrays
        for var in varnames:
            variables[var] = np.concatenate((
                variables[var], np.expand_dims(varlist[var], axis=0),
            ))

    # Compile storage arrays into xarray dataset
    for var in varnames:
        variables[var] = (['date', 'y'], variables[var])
    ds = xr.Dataset(variables, coords={'date': dates, 'y': sections['ji'][0]})

    return ds
