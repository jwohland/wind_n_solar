"""
Correction of the long-term evolution of 100m wind speeds in CERA20C.

The approach consists of

    a) calculation of 3 hourly s_100 from u_100 and v_100 in CERA20C and 20CRv3
    b) regridding of 20CRv3 to the CERA20C grid
    c) calculation of low-pass filtered timeseries (Lanczos)
    d) Computing s_100_CERA_corrected = s_100_CERA + <s_100_20CRv3> - <s_100_CERA20C>

This script does d).
"""
import glob
import xarray as xr
import pandas as pd
import numpy as np


def build_filedic(data_path, lanczos_path):
    """
    builds a directory with paths to raw data ['CERA'] and the lanczos-filtered versions ['lanczos(X)']
    :param data_path: path to raw data
    :param lanczos_path: path to lanczos filtered data
    :return:
    """
    filedic = {'CERA': sorted(glob.glob(data_path + 'CERA20C/*.nc')),
               'lanczos(CERA)': sorted(glob.glob(lanczos_path + 'CERA_7*.nc')),
               'lanczos(20CR)': sorted(glob.glob(lanczos_path + '20CRv3_5*.nc'))}
    return filedic


def drop_non_wind(ds):
    """
    Preprocessing function to drop all variables that are not wind.
    :param ds: CERA20C dataset
    :return:
    """
    return ds.drop(['fdir', 't2m', 'ssrd'])


def calculate_correction(filedic):
    """
    calculates low-frequency correction term
    :param filedic: dictionary including lanczos filtered versions of CERA20C and 20CRv3
    :return:
    """
    lanczos_cera = xr.open_mfdataset(filedic['lanczos(CERA)'], combine='by_coords')
    lanczos_noaa = xr.open_mfdataset(filedic['lanczos(20CR)'], combine='by_coords')
    return lanczos_noaa.drop('number').squeeze() - lanczos_cera.drop('number').squeeze()


def main():
    data_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/'
    lanczos_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/lanczos/'
    out_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/trend_corrected/'

    filedic = build_filedic(data_path, lanczos_path)
    diff = calculate_correction(filedic)
    diff.load()
    # calculate for individual years for memory reasons
    for year_index in range(len(filedic['CERA'])):
        cera_raw = xr.open_mfdataset(filedic['CERA'][year_index], combine='by_coords', preprocess=drop_non_wind)
        cera_corrected = (cera_raw['s100'] + diff['s100_low']).to_dataset(name='s100')
        cera_corrected.compute()
        year = pd.to_datetime(cera_corrected.time[0].values).year
        print(year)
        name = 'CERA20C_correcteds100_' + str(year) + '_representative.nc'
        # some years consist of nans only (because of filter boundary effects)
        if not np.isnan(cera_corrected['s100'].values).all():
            cera_corrected.to_netcdf(out_path + name)


main()
