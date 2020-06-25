"""
Correction of the long-term evolution of 100m wind speeds in CERA20C.

The approach consists of

    a) calculation of 3 hourly s_100 from u_100 and v_100 in CERA20C and 20CRv3
    b) regridding of 20CRv3 to the CERA20C grid
    c) calculation of low-pass filtered timeseries (Lanczos)
    d) Computing s_100_CERA_corrected = s_100_CERA + <s_100_20CRv3> - <s_100_CERA20C>

This script does c).
"""

import numpy as np
import xarray as xr
from oceans.filters import lanc
import glob
import sys


def apply_lanczos(df, years=4.5):
    """
    :param df: Input data as pandas DataFrame. Has to include df['s100'] wind speeds at 100m
    :param years: filter window width equals years * timesteps/year
    :return: pandas DataFrame containing low-pass filtered time series
    """
    freq = 1. / (8 * 365 * years)
    N = int(1 / freq)
    wt = lanc(N, freq)
    res = np.convolve(wt,
                      df['s100'],
                      mode='same')  # mode='same': length of filtered signal identical to unfiltered
    res[0:N] = np.nan  # decided to keep all values first (simpler for mapping to timesteps) and
    res[-N:] = np.nan  # set to nan here
    df['s100_low'] = res
    return df['s100_low']


def get_filelist(cera_index, noaa_index, data_path):
    """
    calculates list of files needed for the CERA ensemble or the specific noaa 20cr ensemble member
    :param cera_index: index of cera ensemble member
    :param noaa_index: index of noaa 20Cr. Members chosen from a randomly created list indexed with noaa_index
    :param data_path: path to the data directory
    :return:
    """
    if cera_index:
        filelist = glob.glob(data_path + 'CERA20C/' + '*instant.nc')
    elif noaa_index:
        noaa_member = ['001', '011', '014', '020', '027', '037', '045', '046', '069', '075'][noaa_index]
        filelist = glob.glob(data_path + '20CRv3/wspd/daily/' + noaa_member +'/*.nc')  # todo this is just a placeholder. Has to point to regridded 100m 20CRv3 wind speeds later
    return filelist


def main(cera_index, noaa_index):
    data_path = '/cluster/work/apatt/wojan/data/'  # todo update path
    out_path = '/cluster/work/apatt/wojan/data/'  # todo update path
    name = ('CERA_' + str(cera_index) if cera_index else '20CR_' + str(noaa_index)) + '_lanczos_54months.nc'# todo should this be the index [0,1,2,3,4] or the ensemble member used?
    try:
        xr.open_dataset(out_path + name)
    except FileNotFoundError:
        filelist = get_filelist(cera_index, noaa_index, data_path)
        ds = xr.open_mfdataset(filelist,
                               chunks={'latitude': 5, 'longitude': 5},
                               combine='by_coords')
        if cera_index:
            ds = ds.sel({'number': cera_index}, drop=True)  # todo check if number correct name
        # todo from here on very handwavy
        dfs = []
        for loc in ds.loc:
            df = ds.to_dataframe()
            apply_lanczos(df)
            dfs.append(df)
        dfs.merge()
        ds_out = dfs.todataset()
        # todo add 1.5h offset

        ds_out.to_netcdf(out_path + name)


cera_index, noaa_index = int(sys.argv[1])-1, int(sys.argv[2])-1
main(cera_index, noaa_index)