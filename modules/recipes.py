import numpy                as np
from graphviz             import Source


def calc_rho(Sal, TempC, P):
    """ Calculate rho: Based on SOG code
    """
    
    # Calculate the square root of the salinities
    sqrSal = np.sqrt(Sal)

    # Calculate the density profile at the grid point depths
    # Pure water density at atmospheric pressure
    # (Bigg P.H., (1967) Br. J. Applied Physics 8 pp 521-537)
    R1 = ((((6.536332e-9 * TempC - 1.120083e-6) * TempC + 1.001685e-4)
           * TempC - 9.095290e-3) * TempC + 6.793952e-2) * TempC - 28.263737
    R2 = (((5.3875e-9 * TempC - 8.2467e-7) * TempC + 7.6438e-5)
          * TempC - 4.0899e-3) * TempC + 8.24493e-1
    R3 = (-1.6546e-6 * TempC + 1.0227e-4) * TempC - 5.72466e-3

    # International one-atmosphere equation of state of seawater
    SIG = (4.8314e-4 * Sal + R3 * sqrSal + R2) * Sal + R1

    # Specific volume at atmospheric pressure
    V350P = 1.0 / 1028.1063
    SVA   = -SIG * V350P / (1028.1063 + SIG)

    # Density anomoly at atmospheric pressure
    rho = 28.106331 - SVA / (V350P * (V350P + SVA)) + 1000
    
    return rho


def flow_diagram():
    """
    """
    
    diagram = Source('digraph G { ' +
            'rankdir=TB ' +
            'node[shape=circle]; a1 b1 c1 d1 ' +
            'node[shape=box] ' +
            'a1 -> b1 -> c1 -> d1 ' +
            'b1 -> b2 -> b1 ' +
            'c1 -> c2 -> c3 -> c2 -> c1 ' +
            'a1 [label="Load Data"] ' +
            'b1 [label="Make Figure"] ' +
            'c1 [label="Animate"] ' +
            'd1 [label="Save"] ' +
            'b2 [label="def make_figure( ):"] ' +
            'c2 [label="def next_frame( t ):"] ' +
            'c3 [label="def update_figure( ):"] ' +
            '{rank=same; b1 b2} ' +
            '{rank=same; c1 c2 c3}}')
    
    return diagram


def load_GEM_from_erddap(
        timerange, window=[None, None, None, None],
        fields=['u_wind', 'v_wind'],
        gridpath=(
            'https://salishsea.eos.ubc.ca/erddap/'
            'griddap/ubcSSaAtmosphereGridV1'),
        datapath=(
            'https://salishsea.eos.ubc.ca/erddap/'
            'griddap/ubcSSaSurfaceAtmosphereFieldsV1'),
):
    """Returns surface atmospheric variables from the Environment
    Canada GEM 2.5 km HRDPS atmospheric model, accessed through the
    ERDDAP server.

    :arg timerange: Start and end datetimes for the requested data range.
                    (ex. ['yyyy mmm dd HH:MM', 'yyyy mmm dd HH:MM'])
    :type timerange: list or tuple of str

    :arg window: Model domain slice bounds.
                 (ex. [x_min, x_max, y_min, y_max])
    :type window: list or tuple of integers

    :arg fields: Requested variables (*must match netCDF fields in file*)
    :type fields: list or tuple of str

    :arg str gridpath: Location of grid information on ERDDAP server

    :arg str datapath: Location of data on ERDDAP server

    :returns: :py:class:`xarray.Dataset` of lat/lon coordinates and
              model variables
    :rtype: :py:class:`xarray.Dataset`
    """

    # Create slices
    starttime, endtime = map(dparser.parse, timerange)
    timeslice = slice(starttime, endtime)
    xslice = slice(window[0], window[1])
    yslice = slice(window[2], window[3])

    # Load GEM grid
    grid = xr.open_dataset(gridpath)
    data = xr.open_dataset(datapath)
    GEM = xr.Dataset({
              'longitude': grid.longitude.sel(gridX=xslice, gridY=yslice)-360,
              'latitude': grid.latitude.sel(gridX=xslice, gridY=yslice)})

    # Load GEM variables
    for field in fields:
        GEM = GEM.merge({field: data[field].sel(
                        time=timeslice, gridX=xslice, gridY=yslice)})

    return GEM


def load_GEM_from_path(
        timerange, window=[None, None, None, None],
        fields=['u_wind', 'v_wind'], model='operational',
        path='/results/forcing/atmospheric/GEM2.5',
):
    """Returns surface atmospheric variables from the Environment
    Canada GEM 2.5 km HRDPS atmospheric model, accessed through the
    local filesystem.
    
    :arg timerange: Start and end datetimes for the requested data range.
        (ex. ['yyyy mmm dd HH:MM', 'yyyy mmm dd HH:MM'])
    :type timerange: list or tuple of str
    
    :arg window: Model domain slice bounds.
        (ex. [x_min, x_max, y_min, y_max])
    :type window: list or tuple of integers
    
    :arg fields: Requested variables (*must match netCDF fields in file*)
    :type fields: list or tuple of str
    
    :arg model: Model type
        (one of 'operational', 'research')
    :type model: str
    
    :arg path: Location of local results server
    :type path: str
    
    :returns: :py:class:`xarray.Dataset` of lat/lon coordinates and
        model variables
    :rtype: :py:class:`xarray.Dataset`
    """
    
    # Define filename prefix
    if model is 'operational': prefix='ops'
    elif model is 'research' : prefix='res'
    else: raise ValueError('Unknown model type: {}'.format(model))
    
    # Create slices
    starttime, endtime = map(dparser.parse, timerange)
    timeslice = slice(starttime, endtime)
    xslice    = slice(window[0], window[1])
    yslice    = slice(window[2], window[3])
        
    # Create list of sequential filenames to load
    date = starttime
    filenames = []
    while date < endtime:
        filename = '{}_{}.nc'.format(prefix, date.strftime('y%Ym%md%d'))
        filenames.append(os.path.join(path, model, filename))
        date = date + timedelta(days=1)
    
    # Concatenate files using NCO and load using xarray
    data = xr.open_mfdataset(filenames).rename({
              'time_counter': 'time' ,
              'x'           : 'gridX',
              'y'           : 'gridY'})
    GEM  = xr.Dataset({
              'longitude': data.nav_lon.sel(gridX=xslice, gridY=yslice)-360,
              'latitude' : data.nav_lat.sel(gridX=xslice, gridY=yslice)})
    
    # Load remaining variables into dataset
    for field in fields:
        GEM = GEM.merge({field: data[field].sel(
                        time=timeslice, gridX=xslice, gridY=yslice)})
    
    return GEM


def load_NEMO_from_erddap(
        timerange, depth=[None, None], window=[None, None, None, None],
        fields=['salinity', 'temperature', 'u_vel', 'v_vel'],
        path='https://salishsea.eos.ubc.ca/erddap/griddap',
        bathy_dataset='ubcSSnBathymetry2V1',
):
    """Returns vector and tracer variables from the Salish Sea
    NEMO model, accessed through the ERDDAP server.

    :arg timerange: Start and end datetimes for the requested data range.
        (ex. ['yyyy mmm dd HH:MM', 'yyyy mmm dd HH:MM'])
    :type timerange: list or tuple of str

    :arg depth: Horizontal depth slice.
    :type depth: integer

    :arg window: Model domain slice bounds.
        (ex. [x_min, x_max, y_min, y_max])
    :type window: list or tuple of integers

    :arg fields: Requested variables
        (one of 'u_vel', 'v_vel', 'w_vel', 'salinity', 'temperature')
    :type fields: list or tuple of str

    :arg path: Location of ERDDAP server
    :type path: str

    :arg path: Bathymetry Dataset on ERDDAP server
    :type path: str

    :returns: :py:class:`xarray.Dataset` of lat/lon coordinates and
        model variables
    :rtype: :py:class:`xarray.Dataset`
    """

    # Create slices
    starttime, endtime = map(dparser.parse, timerange)
    timeslice = slice(starttime, endtime)
    depthslice = slice(depth[0], depth[1])
    xslice = slice(window[0], window[1])
    yslice = slice(window[2], window[3])

    # Load NEMO grid variables
    grid = xr.open_dataset(os.path.join(path, bathy_dataset),
                           mask_and_scale=False)
    NEMO = xr.Dataset({
              'longitude': grid.longitude.sel(gridX=xslice, gridY=yslice),
              'latitude': grid.latitude.sel(gridX=xslice, gridY=yslice),
              'bathymetry': grid.bathymetry.sel(gridX=xslice, gridY=yslice)})

    # Load u velocity
    if 'u_vel' in fields:
        u = xr.open_dataset(os.path.join(path, 'ubcSSn3DuVelocity1hV1'))
        NEMO = NEMO.merge({'u_vel': u.uVelocity.sel(
                 time=timeslice, depth=depthslice,
                 gridX=xslice, gridY=yslice)})

    # Load v velocity
    if 'v_vel' in fields:
        v = xr.open_dataset(os.path.join(path, 'ubcSSn3DvVelocity1hV1'))
        NEMO = NEMO.merge({'v_vel': v.vVelocity.sel(
                 time=timeslice, depth=depthslice,
                 gridX=xslice, gridY=yslice)})

    # Load w velocity (reindex to centered depths, unstagger before use)
    if 'w_vel' in fields:
        w = xr.open_dataset(os.path.join(path, 'ubcSSn3DwVelocity1hV1'))
        NEMO = NEMO.merge({'w_vel': xr.DataArray(w.wVelocity.sel(
                 time=timeslice, depth=depthslice, gridX=xslice, gridY=yslice),
                 {'time': NEMO.time,  'depth': NEMO.depth,
                  'gridY': NEMO.gridY, 'gridX': NEMO.gridX})})

    # Load tracers
    if 'salinity' in fields or 'temperature' in fields:
        trc = xr.open_dataset(os.path.join(path, 'ubcSSn3DTracerFields1hV1'))

        # Salinity
        if 'salinity' in fields:
            NEMO = NEMO.merge({'salinity': trc.salinity.sel(
                 time=timeslice, depth=depthslice,
                 gridX=xslice, gridY=yslice)})

        # Temperature
        if 'temperature' in fields:
            NEMO = NEMO.merge({'temperature': trc.temperature.sel(
                 time=timeslice, depth=depthslice,
                 gridX=xslice, gridY=yslice)})

    # Load and merge mesh mask (replace gridZ with NEMO.depth)
    mask = xr.open_dataset(os.path.join(path, 'ubcSSn3DMeshMask2V1'))
    NEMO = NEMO.merge({
        'mask': mask.merge(
            {'depth': ('gridZ', mask.gdept.isel(time=0, gridX=0, gridY=0))}
                ).swap_dims({'gridZ': 'depth'}).isel(time=0).drop(
                    ('time', 'gridZ')).sel(
                        depth=depthslice, gridX=xslice, gridY=yslice).tmask})

    return NEMO


def load_NEMO_from_path(
        timerange, depth=[None, None], window=[None, None, None, None],
        model='nowcast', resolution='h', path='/results/SalishSea',
        bathy_meter='/results/SalishSea/nowcast/01aug16/bathy_meter.nc',
        fields=['salinity', 'temperature', 'u_vel', 'v_vel']
):
    """Returns vector and tracer variables from the Salish Sea
    NEMO model, accessed through the local filesystem.

    This function uses the Python NCO package, and is super slow. If we can
    speed this up somehow that would be awesome!

    :arg timerange: Start and end datetimes for the requested data range.
        (ex. ['yyyy mmm dd HH:MM', 'yyyy mmm dd HH:MM'])
    :type timerange: list or tuple of str

    :arg depth: Horizontal depth slice.
    :type depth: integer

    :arg window: Model domain slice bounds.
        (ex. [x_min, x_max, y_min, y_max])
    :type window: list or tuple of integers

    :arg model: Model run type
        (one of 'forecast', 'forecast2', 'nowcast', 'nowcast-green')
    :type model: str

    :arg resolution: Time resolution ('h' for hourly, 'd' for daily avg)
    :type resolution: str

    :arg path: Location of local results server
    :type path: str

    :arg path: Location of bathymetry file
    :type path: str

    :arg fields: Requested variables
        (one of 'u_vel', 'v_vel', w_vel, 'salinity', 'temperature')
    :type fields: list or tuple of str

    :returns: :py:class:`xarray.Dataset` of lat/lon coordinates and
        model variables
    :rtype: :py:class:`xarray.Dataset`
    """

    # Create slices
    starttime, endtime = map(dparser.parse, timerange)
    timeslice = slice(starttime, endtime)
    depthslice = slice(depth[0], depth[1])
    xslice = slice(window[0], window[1])
    yslice = slice(window[2], window[3])

    if model == 'nowcast-green':
        extra_inds = ('nav_lon', 'nav_lat', 'time_centered',
                      'axis_nbounds', 'nvertex')
    else:
        extra_inds = ('nav_lon', 'nav_lat', 'tbnds')

    # Create lists of sequential filenames to load
    date = starttime
    filenames_T = []
    filenames_U = []
    filenames_V = []
    filenames_W = []
    while date < endtime:
        datestr = date.strftime('%d%b%y').lower()
        prefix = 'SalishSea_1{}_{}_{}_grid_'.format(
                resolution, date.strftime('%Y%m%d'), date.strftime('%Y%m%d'))
        filenames_T.append(os.path.join(
                path, model, datestr, '{}T.nc'.format(prefix)))
        filenames_U.append(os.path.join(
                path, model, datestr, '{}U.nc'.format(prefix)))
        filenames_V.append(os.path.join(
                path, model, datestr, '{}V.nc'.format(prefix)))
        filenames_W.append(os.path.join(
                path, model, datestr, '{}W.nc'.format(prefix)))
        date = date + timedelta(days=1)

    # Load NEMO grid variables
    grid = xr.open_dataset(
                 bathy_meter,
                 mask_and_scale=False)
    NEMO = xr.Dataset({
              'longitude': grid.nav_lon.sel(x=xslice, y=yslice),
              'latitude': grid.nav_lat.sel(x=xslice, y=yslice),
              'bathymetry': grid.Bathymetry.sel(x=xslice, y=yslice)})

    # Load u velocity
    if 'u_vel' in fields:
        u = xr.open_mfdataset(filenames_U)
        u = u.drop(extra_inds).rename({'depthu': 'depth'})
        NEMO = NEMO.merge({'u_vel': u.vozocrtx.sel(
                 time_counter=timeslice, depth=depthslice, x=xslice, y=yslice)})

    # Load v velocity
    if 'v_vel' in fields:
        v = xr.open_mfdataset(filenames_V).drop(extra_inds).rename({'depthv': 'depth'})
        NEMO = NEMO.merge({'v_vel': v.vomecrty.sel(
                 time_counter=timeslice, depth=depthslice, x=xslice, y=yslice)})

    # Load w velocity
    if 'w_vel' in fields:
        w = xr.open_mfdataset(filenames_W).drop(extra_inds).rename({'depthw': 'depth'})
        NEMO = NEMO.merge({'w_vel': w.vovecrtz.sel(
                 time_counter=timeslice, depth=depthslice, x=xslice, y=yslice)})

    # Load tracers
    if 'salinity' in fields or 'temperature' in fields:
        trc = xr.open_mfdataset(filenames_T).drop(extra_inds).rename({'deptht': 'depth'})

        # Salinity
        if 'salinity' in fields:
            NEMO = NEMO.merge({'salinity': trc.vosaline.sel(
                 time_counter=timeslice, depth=depthslice, x=xslice, y=yslice)})

        # Temperature
        if 'temperature' in fields:
            NEMO = NEMO.merge({'temperature': trc.votemper.sel(
                 time_counter=timeslice, depth=depthslice, x=xslice, y=yslice)})

    # Rename indices to ERDDAP conventions
    NEMO = NEMO.rename({'time_counter': 'time',
                        'x': 'gridX',
                        'y': 'gridY'})

    # Load and merge mesh mask (replace gridZ with NEMO.depth)
    mask = xr.open_dataset('https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSn3DMeshMask2V1')
    NEMO = NEMO.merge({
        'mask': mask.merge(
            {'depth': ('gridZ', mask.gdept.isel(time=0, gridX=0, gridY=0))}
                ).swap_dims({'gridZ': 'depth'}).isel(time=0).drop(
                    ('time', 'gridZ')).sel(
                        depth=depthslice, gridX=xslice, gridY=yslice).tmask})

    return NEMO
