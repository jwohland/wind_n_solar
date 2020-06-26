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
import pandas as pd


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
    :param noaa_index: index of noaa 20Cr
    :param data_path: path to the data directory
    :return:
    """
    if cera_index >= 0:
        filelist = glob.glob(data_path + 'CERA20C/*.nc')
    elif noaa_index >= 0:
        filelist = glob.glob(data_path + '20CRv3/*.nc')
    return sorted(filelist)


def main(cera_index, noaa_index, lat_index):
    """
    :param cera_index: number of the CERA ensemble member to be used. If <0, noaa 20CR is used.
    :param noaa_index: number of the 20CRv3 ensemble member to be used. If <0, CERA20C is used.
    :param lat_index: controls which latitude is to be computed
    :return:
    """
    data_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/'
    out_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/'
    name = ('CERA_' + str(cera_index) if cera_index >= 0 else '20CRv3_' + str(noaa_index)) + \
           '_latindex_' + str(lat_index) + '_lanczos_54months.nc'
    try:
        xr.open_dataset(out_path + name)
    except FileNotFoundError:
        filelist = get_filelist(cera_index, noaa_index, data_path)
        ds = xr.open_mfdataset(filelist,
                               combine='by_coords')
        ds = prepare_data(ds)
        dfs = []
        lat = ds.latitude.values[lat_index]
        for lon in ds.longitude.values:
            print(lon)
            data_tmp = ds.sel({'latitude': lat, 'longitude': lon}).load()
            df = data_tmp.to_dataframe()
            apply_lanczos(df)
            df = df.drop('s100', axis=1)
            dfs.append(pd.pivot_table(df, index=['time', 'number', 'latitude', 'longitude']).to_xarray())
        ds_out = xr.combine_by_coords(dfs)
        ds_out.to_netcdf(out_path + name)


def prepare_data(ds):
    ds = ds['s100']
    if cera_index >= 0:
        ds = ds.sel({'number': cera_index})
    elif noaa_index >= 0:
        ds = ds.sel({'number': noaa_index})
    ds = ds.chunk({'latitude': 1, 'longitude': 1, 'time': 318496})
    return ds


# todo this currently can not work for 20CR

cera_index, noaa_index, lat_index = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
print('cera index: ' + str(cera_index))
print('noaa index: ' + str(noaa_index))
print('lat index: ' + str(lat_index))
if cera_index >= 0 and noaa_index >= 0:
    print('This case is not allowed')
else:
    main(cera_index, noaa_index, lat_index)