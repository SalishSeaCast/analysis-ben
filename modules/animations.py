# Copyright 2016 The Salish Sea NEMO Project and
# The University of British Columbia

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Functions for Salish Sea NEMO postprocessing and visualization
"""

import matplotlib           as mpl
mpl.use('Agg')
import numpy                as np
import xarray               as xr
import pandas               as pd
import datetime             as dtm
import pytz
import dateutil.parser      as dparser
import matplotlib.pyplot    as plt
import matplotlib.animation as animation
import scipy.io             as sio
import os
import collections
from salishsea_tools import viz_tools, tidetools
from eofs.xarray     import Eof

mpl.rcParams.update({'font.size': 12})
mpl.rcParams["axes.formatter.useoffset"] = False


def load_results(
        timerange, spacing={'NEMO': 5, 'GEM': 5},
        plot_opts={'wind'    : True, 'currents': True,
                   'salinity': True, 'filtered': True},
        window={'x_NEMO': slice(0, 400), 'x_GEM': slice(0, 637500),
                'y_NEMO': slice(0, 900), 'y_GEM': slice(0, 662500)}
):
    """
    """
    
    # Load Nowcast and GEM results from ERDDAP using xarray
    grid, qty = load_ERDDAP(timerange=timerange, window=window,
                            plot_opts=plot_opts)
    grid['spc_NEMO'] = spacing['NEMO']
    grid['spc_GEM' ] = spacing['GEM' ]
    
    # Initialize filtered qty as empty dict
    qty_f = {}
    
    # Process wind
    if plot_opts['wind']:
        # Subsample wind fields
        qty['u_GEM'] = qty['u_GEM'][:, ::spacing['GEM' ], ::spacing['GEM' ]]
        qty['v_GEM'] = qty['v_GEM'][:, ::spacing['GEM' ], ::spacing['GEM' ]]
        
        if plot_opts['filtered']:
            # Apply Doodson filter to GEM winds
            qty_f['u_GEM'] = tidetools.filter_timeseries(
                             qty['u_GEM'], method='doodson')
            qty_f['v_GEM'] = tidetools.filter_timeseries(
                             qty['v_GEM'], method='doodson')
    
    # Process currents
    if plot_opts['currents']:
        # Unstagger NEMO currents
        qty['u_NEMO'], qty['v_NEMO'] = viz_tools.unstagger(
                                       qty['u_NEMO'], qty['v_NEMO'])
        qty['v_NEMO'] = qty['v_NEMO'].reindex_like(qty['u_NEMO'])
        
        # Subsample velocity fields
        qty['u_NEMO'] = qty['u_NEMO'][:, ::spacing['NEMO'], ::spacing['NEMO']]
        qty['v_NEMO'] = qty['v_NEMO'][:, ::spacing['NEMO'], ::spacing['NEMO']]
        
        # Rotate NEMO currents
        qty['u_NEMO'], qty['v_NEMO'] = viz_tools.rotate_vel(
                              qty['u_NEMO'], qty['v_NEMO'], origin='grid')
        
        if plot_opts['filtered']:
            # Apply Doodson filter to NEMO currents
            qty_f['u_NEMO'] = tidetools.filter_timeseries(
                              qty['u_NEMO'], method='doodson')
            qty_f['v_NEMO'] = tidetools.filter_timeseries(
                              qty['v_NEMO'], method='doodson')
    
    # Process tracers
    if plot_opts['salinity']:
        if plot_opts['filtered']:
            # Apply Doodson filter to NEMO tracers
            qty_f['S_NEMO'] = tidetools.filter_timeseries(
                              qty['S_NEMO'], method='doodson')
    
    return grid, qty, qty_f


def plot_tracers(time_ind, ax, cmap, clim, qty, NEMO, zorder=0):
    """
    """
    
    # NEMO horizontal tracers
    C = ax.contourf(NEMO['lon'], NEMO['lat'],
        np.ma.masked_values(NEMO[qty].sel(time=time_ind, method='nearest'), 0),
        range(clim[0], clim[1], clim[2]), cmap=cmap, zorder=zorder)
    
    return C


def plot_currents(time_ind, ax, spacing, NEMO, zorder=5):
    """
    """
    
    # NEMO horizontal currents
    Q = ax.quiver(
        NEMO['lon'][1::spacing, 1::spacing],
        NEMO['lat'][1::spacing, 1::spacing],
        np.ma.masked_values(
            NEMO['u'].sel(time=time_ind, method='nearest'), 0),
        np.ma.masked_values(
            NEMO['v'].sel(time=time_ind, method='nearest'), 0),
        scale=10, zorder=zorder)
    
    return Q


def plot_wind(time_ind, ax, spacing, GEM, zorder=10, processed=False,
              color='gray'):
    """
    """
    
    # Determine whether to space vectors
    spc = spacing
    if processed: spc = 1
    
    # GEM winds
    Q = ax.quiver(
        GEM['lon'][::spacing, ::spacing],
        GEM['lat'][::spacing, ::spacing],
        GEM['u_wind'].sel(time=time_ind, method='nearest')[::spc, ::spc],
        GEM['v_wind'].sel(time=time_ind, method='nearest')[::spc, ::spc],
        color=color, edgecolor='k', scale=40,
        linewidth=0.5, headwidth=4, zorder=zorder)
    
    return Q


def plot_drifters(time_ind, ax, drifters, zorder=15):
    """
    """
    
    # Define color palette
    palette = ['blue', 'teal', 'cyan', 'green', 'lime', 'darkred', 'red',
               'orange', 'magenta', 'purple', 'black', 'dimgray', 'saddlebrown',
               'blue', 'teal', 'cyan', 'green', 'lime', 'darkred', 'red',
               'orange', 'magenta', 'purple', 'black', 'dimgray', 'saddlebrown']
    
    # Plot drifters
    L = collections.OrderedDict()
    P = collections.OrderedDict()
    for i, drifter in enumerate(drifters.keys()):
        # Compare timestep with available drifter data
        dtime = [pd.Timestamp(t.to_pandas()).to_datetime() - time_ind
                 for t in drifters[drifter].time[[0, -1]]]
        # Show drifter track if data is within time threshold
        if (dtime[0].total_seconds() < 3600 and   # 1 hour before deployment
            dtime[1].total_seconds() > -86400):   # 24 hours after failure
            L[drifter] = ax.plot(
                drifters[drifter].lon.sel(time=time_ind, method='nearest'),
                drifters[drifter].lat.sel(time=time_ind, method='nearest'),
                '-', linewidth=2, color=palette[i], zorder=zorder)
            P[drifter] = ax.plot(
                drifters[drifter].lon.sel(time=time_ind, method='nearest'),
                drifters[drifter].lat.sel(time=time_ind, method='nearest'),
                'o', color=palette[i], zorder=zorder+1)
        else: # Hide if outside time threshold
            L[drifter] = ax.plot(
                [], [], '-', linewidth=2, color=palette[i], zorder=zorder)
            P[drifter] = ax.plot(
                [], [], 'o', color=palette[i], zorder=zorder+1)
    
    return L, P



def plot_currents_old(
        ax, fig, init_time, grid, qty, drifters,
        map_bounds=[-124, -122.7, 48.3, 49.7],
        cmap_type='jet', landmask='burlywood',
        plot_opts={'wind'    : True, 'currents': True,
                   'salinity': True, 'drifters': True},
        bathy_path='/ocean/bmoorema/research/MEOPAR/NEMO-forcing/grid',
        bathy_file='bathy_meter_SalishSea2.nc'
):
    """
    """
    
    # Define color palette
    palette = ['blue', 'teal', 'cyan', 'green', 'lime', 'darkred', 'red',
               'orange', 'magenta', 'purple', 'black', 'dimgray', 'saddlebrown',
               'blue', 'teal', 'cyan', 'green', 'lime', 'darkred', 'red',
               'orange', 'magenta', 'purple', 'black', 'dimgray', 'saddlebrown']
    
    # Initialize plot objects dict
    plot_objs = {}
    
    # Bathymetry path
    bathy = os.path.join(bathy_path, bathy_file)
    
    # Plot land mask and coastline
    viz_tools.plot_land_mask(ax, bathy, coords='map', color=landmask, zorder=1)
    viz_tools.plot_coastline(ax, bathy, coords='map', zorder=2)
    
    # Conditional plots
    if plot_opts['salinity']:
        # NEMO surface salinity
        cmap = plt.get_cmap(cmap_type)
        plot_objs['C_NEMO'] = ax.contourf(grid['lon_NEMO'], grid['lat_NEMO'],
                 qty['S_NEMO'].sel(time=init_time, method='nearest'),
                 range(32), cmap=cmap, zorder=0)
        cbar = fig.colorbar(plot_objs['C_NEMO'], ax=ax)
        cbar.set_label('Practical Salinity')
    
    if plot_opts['drifters']:
        # Plot drifters
        plot_objs['L_DRIFTERS'] = collections.OrderedDict()
        plot_objs['P_DRIFTERS'] = collections.OrderedDict()
        for i, drifter in enumerate(drifters.keys()):
            # Compare timestep with available drifter data
            dtime = [pd.Timestamp(t.to_pandas()).to_datetime() - init_time
                     for t in drifters[drifter]['lon'].time[[0, -1]]]
            # Show drifter track if data is within time threshold
            if (dtime[0].total_seconds() < 3600 and   # 1 hour before deployment
                dtime[1].total_seconds() > -86400):   # 24 hours after failure
                plot_objs['L_DRIFTERS'][drifter] = ax.plot(
                    drifters[drifter]['lon'].sel(
                        time=init_time, method='nearest'),
                    drifters[drifter]['lat'].sel(
                        time=init_time, method='nearest'),
                    '-', linewidth=2, color=palette[i], zorder=3)
                plot_objs['P_DRIFTERS'][drifter] = ax.plot(
                    drifters[drifter]['lon'].sel(
                        time=init_time, method='nearest'),
                    drifters[drifter]['lat'].sel(
                        time=init_time, method='nearest'),
                    'o', color=palette[i], zorder=4)
            else: # Hide if outside time threshold
                plot_objs['L_DRIFTERS'][drifter] = ax.plot(
                    [], [], '-', linewidth=2, color=palette[i], zorder=3)
                plot_objs['P_DRIFTERS'][drifter] = ax.plot(
                    [], [], 'o', color=palette[i], zorder=4)

    if plot_opts['currents']:
        # NEMO surface currents
        plot_objs['Q_NEMO'] = ax.quiver(
                grid['lon_NEMO'][1::grid['spc_NEMO'], 1::grid['spc_NEMO']],
                grid['lat_NEMO'][1::grid['spc_NEMO'], 1::grid['spc_NEMO']],
                np.ma.masked_values(
                     qty['u_NEMO'].sel(time=init_time, method='nearest'), 0),
                np.ma.masked_values(
                     qty['v_NEMO'].sel(time=init_time, method='nearest'), 0),
                scale=10, zorder=5)
        Qkey_NEMO = plt.quiverkey(plot_objs['Q_NEMO'], 0.87, 0.93, 1, '1 m/s',
                coordinates='axes')
        Qkey_NEMO.set_zorder(8)
    
    if plot_opts['wind']:
        # GEM winds
        plot_objs['Q_GEM'] = ax.quiver(
                grid['lon_GEM'][::grid['spc_GEM'], ::grid['spc_GEM']],
                grid['lat_GEM'][::grid['spc_GEM'], ::grid['spc_GEM']],
                qty['u_GEM'].sel(time=init_time, method='nearest'),
                qty['v_GEM'].sel(time=init_time, method='nearest'),
                color='gray', edgecolor='k', scale=40,
                linewidth=0.5, headwidth=4, zorder=6)
        Qkey_GEM = plt.quiverkey(plot_objs['Q_GEM'],  0.85, 0.87, 5, '5 m/s', 
                coordinates='axes')
        Qkey_GEM.set_zorder(8)

    # Timestamp
    plot_objs['time_text'] = ax.text(
                0.02, 1.02, init_time.strftime('%a %Y-%m-%d %H:%M:%S %Z'),
                transform=ax.transAxes, zorder=8)
    
    # Prettying
    viz_tools.set_aspect(ax)
    ax.set_xlim(map_bounds[0:2])
    ax.set_ylim(map_bounds[2:4])
    ax.grid()
    
    if plot_opts['currents'] or plot_opts['winds']:
        # Quiver Key
        lbox = ax.add_patch(mpl.patches.Rectangle((0.78, 0.85), 0.22, 0.15,
                facecolor='white', transform=ax.transAxes, zorder=7))
    
    return plot_objs


def update_currents(
        ax, time, init_time, plot_objs, grid, qty, drifters,
        cmap_type='jet', landmask='burlywood',
        plot_opts = {'wind'    : True, 'currents': True,
                     'salinity': True, 'drifters': True}
):
    """
    """
    
    if plot_opts['salinity']:
        # Update NEMO salinity contours
        cmap = plt.get_cmap(cmap_type)
        for C in plot_objs['C_NEMO'].collections: C.remove()
        plot_objs['C_NEMO'] = ax.contourf(
                              grid['lon_NEMO'], grid['lat_NEMO'],
                              qty['S_NEMO'].sel(time=time, method='nearest'),
                              range(32), cmap=cmap, zorder=0)
    
    if plot_opts['drifters']:
        # Update drifter positions
        for drifter in plot_objs['L_DRIFTERS'].keys():
            # Compare timestep with available drifter data
            dtime = [pd.Timestamp(t.to_pandas()).to_datetime() - time
                     for t in drifters[drifter]['lon'].time[[0, -1]]]
            # Show drifter track if data is within time threshold
            if (dtime[0].total_seconds() < 3600 and   # 1 hour before deployment
                dtime[1].total_seconds() > -86400):   # 24 hours after failure
                plot_objs['L_DRIFTERS'][drifter][0].set_data(
                    drifters[drifter]['lon'].sel(time=slice(init_time, time)),
                    drifters[drifter]['lat'].sel(time=slice(init_time, time)))
                plot_objs['P_DRIFTERS'][drifter][0].set_data(
                    drifters[drifter]['lon'].sel(time=time, method='nearest'),
                    drifters[drifter]['lat'].sel(time=time, method='nearest'))
            else: # Hide if outside time threshold
                plot_objs['L_DRIFTERS'][drifter][0].set_data([], [])
                plot_objs['P_DRIFTERS'][drifter][0].set_data([], [])
                
    if plot_opts['currents']:
        # Update NEMO current vectors
        plot_objs['Q_NEMO'].set_UVC(
            np.ma.masked_values(
                qty['u_NEMO'].sel(time=time, method='nearest'), 0),
            np.ma.masked_values(
                qty['v_NEMO'].sel(time=time, method='nearest'), 0))
    
    if plot_opts['wind']:
        # Update GEM wind vectors
        plot_objs['Q_GEM' ].set_UVC(
                qty['u_GEM'].sel(time=time, method='nearest'),
                qty['v_GEM'].sel(time=time, method='nearest'))
    
    # Update timestamp
    plot_objs['time_text'].set_text(time.strftime('%a %Y-%m-%d %H:%M:%S %Z'))
    
    return plot_objs


def animate_currents_winds(
        timerange, depth=0, figsize=[10, 10], spacing={'NEMO': 5, 'GEM': 5},
        framerate=12, filetag='', map_bounds=[-124, -122.7, 48.3, 49.7],
        plot_opts = {'wind'    : True, 'currents': True, 'salinity'  : True,
                     'drifters': True, 'filtered': True, 'unfiltered': True},
        writepath='/ocean/bmoorema/research/MEOPAR/analysis-ben/visualization',
        window={'x_NEMO': slice(100, 398), 'x_GEM': slice(300000, 425000),
                'y_NEMO': slice(230, 570), 'y_GEM': slice(237500, 412500)}
):
    """
    """
    
    print('Initializing ...')
    
    # Parse timerange
    starttime, endtime = map(dparser.parse, timerange)
    frames = (endtime - starttime).days * 24
    winlen = dtm.timedelta(hours=20)
    timerange_full = [starttime - winlen, endtime + winlen]
    
    # Subplots
    if plot_opts['unfiltered'] and plot_opts['filtered']:  nplots = 2
    elif plot_opts['unfiltered'] or plot_opts['filtered']: nplots = 1
    
    # Write file path
    writefile = 'NEMO{}_{}m_{}to{}.mp4'.format(filetag, depth,
                starttime.strftime('%Y%b%dT%H'), endtime.strftime('%Y%b%dT%H'))
    filepath = os.path.join(writepath, writefile)
    
    # Load difters
    drifters = load_drifters()
    
    # Load results into xarray dataset dictionaries
    print('Loading model datasets ...')
    grid, qty, qty_f = load_results(
        timerange_full, window=window, spacing=spacing, plot_opts=plot_opts)
    
    # Generate figure panels
    print('Initializing plots ...')
    fig = plt.figure(figsize=figsize)
    
    # Initialize plot objects dictionaries
    plot_objs, plot_objs_f = [], []
    
    # Plot initial vector fields
    if plot_opts['unfiltered']:
        ax1 = fig.add_subplot(1, nplots, 1)
        plot_objs = plot_currents(
                    ax1, fig, starttime, grid, qty, drifters,
                    plot_opts=plot_opts, map_bounds=map_bounds)
    if plot_opts['filtered']:
        ax2 = fig.add_subplot(1, nplots, nplots)
        plot_objs_f = plot_currents(
                    ax2, fig, starttime, grid, qty_f, drifters,
                    plot_opts=plot_opts, map_bounds=map_bounds)
    
    # Next frame definition
    def update_plot(t, plot_objs, plot_objs_f):
        time = (starttime + dtm.timedelta(hours=t))
        if plot_opts['unfiltered']:
            plot_objs = update_currents(
                    ax1, time, starttime, plot_objs, grid, qty, drifters,
                    plot_opts=plot_opts)
        if plot_opts['filtered']:
            plot_objs_f = update_currents(
                    ax2, time, starttime, plot_objs_f, grid, qty_f, drifters,
                    plot_opts=plot_opts)
        return plot_objs, plot_objs_f
    
    # Animate
    print('Compiling animation ...')
    mywriter = animation.FFMpegWriter(fps=framerate, bitrate=10000)
    ani = animation.FuncAnimation(fig, update_plot, frames=frames,
                    fargs=(plot_objs, plot_objs_f), blit=False)
    ani.save(filepath, writer=mywriter)
    print('Done!')


def plot_EOF(
        timerange, depth=0, plot_bounds=[0, 397, 0, 897],
        plot_opts={'wind': True, 'currents': True, 'salinity': True},
        window={'x_NEMO': slice(0, 400), 'x_GEM': slice(0, 637500),
                'y_NEMO': slice(0, 900), 'y_GEM': slice(0, 662500)}
):
    """ Principle component analysis
    Adapted from http://ajdawson.github.io/eofs/examples
    pc1  = solver.pcs(npcs=1, pcscaling=1)
    """
    
    print('Initializing ...')

    # Bathymetry file path
    bathy_path = '/ocean/bmoorema/research/MEOPAR/NEMO-forcing/grid'
    bathy_file = 'bathy_meter_SalishSea2.nc'
    bathy = os.path.join(bathy_path, bathy_file)
    
    # Load ERDDAP results
    grid, qty = load_ERDDAP(timerange, depth=depth, plot_opts=plot_opts,
                            window=window)
    
    # Calculate EOFs
    print('Calculating cross-strait component ...')
    u_eof   = Eof(qty['u_NEMO'])
    u_corr  = u_eof.eofsAsCorrelation(neofs=3)
    u_field = u_eof.eofs(neofs=3)
    print('Calculating along-strait component ...')
    v_eof   = Eof(qty['v_NEMO'])
    v_corr  = v_eof.eofsAsCorrelation(neofs=3)
    v_field = v_eof.eofs(neofs=3)

    # Unstagger principle component vectors
    print('Compiling figure ...')
    U, V = viz_tools.unstagger(u_field.sel(mode=0), v_field.sel(mode=0))
    V = V.reindex_like(U)
    
    # Make figure
    fig, ax = plt.subplots(1, 1, figsize=(8, 15))
    viz_tools.plot_land_mask(ax, bathy, coords='grid', color='burlywood')
    viz_tools.plot_coastline(ax, bathy, coords='grid')
    Q_NEMO = ax.quiver(
        U.gridX[::5], U.gridY[::5], U[::5, ::5], V[::5, ::5], scale=50)
    ax.set_xlim(plot_bounds[0:2])
    ax.set_ylim(plot_bounds[2:4])
    viz_tools.set_aspect(ax)
    
    # Save figure
    fig_path = '/ocean/bmoorema/research/MEOPAR/analysis-ben/visualization/'
    fig_name = 'EOF.eps'
    fig.savefig(os.path.join(fig_path, fig_name), format='eps')
