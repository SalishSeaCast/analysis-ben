# Extract wind from ERDDAP, average along SoG sections and save to netCDF

import numpy as np
import xarray as xr
import yaml
from tqdm import tqdm
from datetime import datetime, timedelta
from salishsea_tools import places

# Paths
f = 'HRDPS_sections.nc'
fi='https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSaSurfaceAtmosphereFieldsV1'
fo='/ocean/bmoorema/research/MEOPAR/analysis-ben/data/SalishSeaCast/' + f

# Load sections from YAML
with open('parameters.yaml') as f: _, sections, _, _ = yaml.safe_load_all(f)
for sec in sections:
    for key in ['xw', 'yw']:
        sec[key] = np.round(np.linspace(*sec[key], 14)).astype(int)

# Daterange
daterange = [datetime(2014, 12, 15), datetime(2019, 1, 15)]

# Load HRDPS record from ERDDAP
ds = xr.open_dataset(fi)

# Preallocate wind storage array
ws = np.empty((np.diff(daterange)[0].days * 24, 0))

# Iterate through sections
for sec in sections:

    # Preallocate chunk storage array
    data = np.empty(0)

    # Construct box containing wind section
    n = np.prod([sec[dim][-1] - sec[dim][0] + 1 for dim in ['xw', 'yw']])
    ii = xr.DataArray(
        np.full((ds.gridY.size, ds.gridX.size), False),
        dims=('gridY', 'gridX'),
    )
    for i, j in zip(sec['xw'], sec['yw']): ii[j, i] = True

    # Process wind data in ~6 month chunks
    hours = int(np.diff(daterange)[0].total_seconds() / 3600)
    chunk = int(hours / 8)
    for hour in tqdm(range(0, hours, chunk)):

        # Date slice for current chunk
        dslc = slice(
            *[daterange[0] + timedelta(hours=h) for h in (hour, hour+chunk-1)]
        )

        # Extract u and v along sections
        u = ds.sel(time=dslc).u_wind.where(ii, drop=True).values.reshape(-1, n)
        v = ds.sel(time=dslc).v_wind.where(ii, drop=True).values.reshape(-1, n)
        u, v = u[:, ~np.isnan(u[0, :])], v[:, ~np.isnan(v[0, :])]

        # Find along-axis windspeed averaged along section and concatenate
        fac = np.cos(np.pi * (1 + 22 / 180) - np.arctan2(v, u) - np.arctan(2))
        data = np.concatenate((data, (fac * np.sqrt(u**2+v**2)).mean(axis=1)))

    # Concatenate to wind storage array
    ws = np.concatenate((ws, data[:, np.newaxis]), axis=1)

# Save to netCDF
dims = {
    'time': ds.time.sel(time=slice(*daterange))[:-1],
    'section': [sec['name'] for sec in sections],
}
xr.Dataset({'windspeed': (['time', 'section'], ws)}, coords=dims).to_netcdf(fo)
