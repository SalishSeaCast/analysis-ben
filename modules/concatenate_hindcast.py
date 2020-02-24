#!/usr/bin/env python3

# Import libraries
# --------------------------------------------------------------------------
import numpy as np
from time import sleep
from glob import glob
from os import remove, mkdir
from os.path import join
from datetime import timedelta
from dateutil.parser import parse
from subprocess import call, Popen
from salishsea_tools import utilities


# Specify parameters
# --------------------------------------------------------------------------
daterange = ['2010 Jan 1', '2020 Jan 1']
varlist = ['nitrate']
indices = {'x': [115, 360], 'y': [310, 788], 'deptht': [0]}
paths = {
    'loadpath': '/results/SalishSea/hindcast.201812',
    'loadpath_cutoff': '/results2/SalishSea/hindcast.201812_annex',
    'savepath': '/ocean/bmoorema/research/MEOPAR/analysis-ben/data/SalishSeaCast',
    'date_cutoff': '2016 Nov 21',
}


# Functions
# --------------------------------------------------------------------------
def clear_procs(procs, njobs):
    """
    """

    while len(procs) >= njobs:
        sleep(1)
        for n, proc in enumerate(procs):
            if proc.poll() is not None:
                if proc.poll() == 0:
                    proc.communicate()
                    del procs[n]
                else:
                    msg = 'Improper termination of process resulting in'
                    msg = f'{msg} return code: {proc.returncode}'
                    raise ValueError(msg)


def concat_files(path, exp, dates, fstr):
    """
    """

    files = sorted(glob(join(path, exp)))
    fn = '_'.join(fstr[:1] + [d.strftime('%Y%m%d') for d in dates] + fstr[1:])
    call(['ncrcat'] + files + [join(path, fn)])
    for file in files:
        remove(file)


def main_routine(
    daterange, varlist, indices, paths,
    res='h', ftype='ptrc_T', maxjobs=4, chunksize=20,
):
    """
    """

    # Parse dates, make temp directory, define file strings, init process list
    dates = [parse(d) for d in daterange]
    temppath = join(paths['savepath'], 'temp')
    mkdir(temppath)
    filestr = [f'SalishSea_1{res}', f'{ftype}.nc']
    procs = []

    # Loop through filenames
    bar = utilities.statusbar('Running ncks ...')
    for chunk in bar(range(0, np.diff(dates)[0].days, chunksize)):

        # Extract in chunks
        for day in range(chunksize):

            # Parse date info, define filenames and paths
            date = dates[0] + timedelta(days=chunk+day)
            datestr = date.strftime('%Y%m%d')
            fn = '_'.join([filestr[0], datestr, datestr, filestr[1]])
            fpath = paths['loadpath']
            if paths['date_cutoff'] is not None:
                if date >= parse(paths['date_cutoff']):
                    fpath = paths['loadpath_cutoff']
            fpath = join(fpath, date.strftime('%d%b%y').lower(), fn)

            # Extract and slice variables and save to temp dir (ncks)
            cmd = ['ncks', '-v', ','.join(varlist)]
            for dim in indices:
                dimslice = ','.join(str(index) for index in indices[dim])
                cmd = cmd + ['-d', ','.join([dim, dimslice])]
            procs.append(
                Popen(cmd + [fpath, join(temppath, f'chunkfile_{day:02d}.nc')])
            )

            # Wait for jobs to clear if maxjobs is reached
            clear_procs(procs, maxjobs)

        # Clear remaining jobs, concatenate chunk and delete daily files
        clear_procs(procs, 1)
        concat_files(
            temppath, 'chunkfile*',
            [date-timedelta(days=chunksize), date], filestr,
        )

    # Concat chunked files (ncrcat), move up 1 directory, then delete temp dir
    # concat_files(temppath, f'{filestr[0]}*', dates, filestr, cleardir=True)


# Parent function
# --------------------------------------------------------------------------
if __name__ == '__main__':
    main_routine(daterange, varlist, indices, paths)
