#! /usr/bin/env python

from salishsea_tools import (
    nc_tools,
    data_tools,
    visualisations,
    viz_tools)
from collections     import OrderedDict
from matplotlib      import pyplot, animation, rcParams, patches
from dateutil        import parser
import datetime

rcParams.update({'font.size': 12})
rcParams["axes.formatter.useoffset"] = False

# Timerange
timerange = ['2016 Jun 28 00:00', '2016 Jul 12 00:00']
starttime, endtime = map(parser.parse, timerange)
frames = int((endtime - starttime).total_seconds() / 3600)

# Load and process results
print('Loading results ...')
GEM  = nc_tools.load_GEM_from_erddap( timerange)
NEMO = nc_tools.load_NEMO_from_erddap(timerange, depth=[0, 1]).isel(depth=0)
DRIFTERS = data_tools.load_drifters()
NEMO['u_vel'] = viz_tools.unstagger_xarray(NEMO.u_vel, 'gridX')
NEMO['v_vel'] = viz_tools.unstagger_xarray(NEMO.v_vel, 'gridY')

# ------------- FIGURE CODE --------------------
print('Making figure ...')
fig, ax   = pyplot.subplots(1, 3, figsize=(20, 10))
L_drift   = OrderedDict()
palette   = ['blue', 'teal', 'cyan', 'green', 'lime', 'darkred', 'red',
             'orange', 'magenta', 'purple', 'black', 'dimgray', 'saddlebrown']

# Plot Drifters
visualisations.create_figure(ax[0], NEMO.isel(time=0))
for index, drifter in enumerate(DRIFTERS['deployment7'].items()):
    L_drift[drifter[0]] = visualisations.plot_drifters(
        ax[0], drifter[1].sel(time=slice(None, starttime)), color=palette[index])

# Plot surface currents
visualisations.create_figure(ax[1], NEMO.isel(time=0))
Q_vel = visualisations.plot_velocity(ax[1], 'NEMO', NEMO.sel(time=starttime, method='nearest'))
Qkey = pyplot.quiverkey(Q_vel, 0.88, 0.94, 1, '1 m/s', coordinates='axes')

# Plot Wind
visualisations.create_figure(ax[2], NEMO.isel(time=0))
Q_wind = visualisations.plot_velocity(ax[2], 'GEM', GEM.sel(time=starttime, method='nearest'),
                    color='gray', scale=60, linewidth=0.5, headwidth=5, mask=False)
lbox = ax[2].add_patch(patches.Rectangle((0.82, 0.88), 0.18, 0.12,
                    facecolor='white', transform=ax[2].transAxes, zorder=10))
Qkey = pyplot.quiverkey(Q_wind, 0.88, 0.94, 5, '5 m/s', coordinates='axes').set_zorder(11)

# Add timestamp
TXT_time = ax[0].text(0.02, 1.02, starttime.strftime('%a %Y-%m-%d %H:%M:%S'), transform=ax[0].transAxes)

# ------------ ANIMATION CODE ------------------

# Create dict of objects to be modified with each timestep
PLOT_OBJS = {'Q_vel': Q_vel, 'Q_wind': Q_wind, 'L_drift': L_drift, 'TXT_time': TXT_time}

# Create local function that updates these objects (iterates over time integer t)
def next_frame(t, PLOT_OBJS):
    
    # Step time index forward
    time_ind = starttime + datetime.timedelta(hours=t)
    
    # Update drifter tracks (one at a time)
    for index, drifter in enumerate(DRIFTERS['deployment7'].items()):
        PLOT_OBJS['L_drift'][drifter[0]] = visualisations.plot_drifters(
            ax[0], drifter[1].sel(time=slice(None, time_ind)),
            DRIFT_OBJS=PLOT_OBJS['L_drift'][drifter[0]], color=palette[index])
    
    # Update surface currents
    PLOT_OBJS['Q_vel']  = visualisations.plot_velocity(
        ax[1], 'NEMO', NEMO.sel(time=time_ind, method='nearest'), Q=PLOT_OBJS['Q_vel'])
    
    PLOT_OBJS['Q_wind'] = visualisations.plot_velocity(
        ax[2], 'GEM',  GEM.sel( time=time_ind, method='nearest'), Q=PLOT_OBJS['Q_wind'], mask=False)
    
    # Update timestamp
    PLOT_OBJS['TXT_time'].set_text(time_ind.strftime('%a %Y-%m-%d %H:%M:%S'))
    
    return PLOT_OBJS

print('Animating ...')
# Call the animation function (create animation object)
ANI = animation.FuncAnimation(fig, next_frame, fargs=[PLOT_OBJS], frames=frames)

# Save the animation
ANI.save('/ocean/bmoorema/research/MEOPAR/analysis-ben/visualization/test_movie.mp4',
         writer=animation.FFMpegWriter(fps=12, bitrate=10000))
