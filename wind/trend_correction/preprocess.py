"""
Correction of the long-term evolution of 100m wind speeds in CERA20C.

The approach consists of

    a) calculation of 3 hourly s_100 from u_100 and v_100 in CERA20C and 20CRv3
    b) regridding of 20CRv3 to the CERA20C grid
    c) calculation of low-pass filtered timeseries (Lanczos)
    d) Computing s_100_CERA_corrected = s_100_CERA + <s_100_20CRv3> - <s_100_CERA20C>

This script does a).
"""
import glob
import xarray as xr


def calc_wind_speeds(ds):
    """
    calculates 100m wind speeds from 100mwind components
    :param ds: dataset containing at least 'u100' and 'v100'
    :return:
    """
    with xr.set_options(keep_attrs=True):
        ds['s100'] = (ds['u100'] ** 2 + ds['v100'] ** 2) ** (1. / 2)
    ds['s100'].attrs['long_name'] = '100 metre wind speed'
    ds = ds.drop(['u100', 'v100'])
    return ds


def preprocess_CERA_wind(file_list, out_path):
    for filepath in file_list:
        ds = xr.open_dataset(filepath, drop_variables=['sro', 'fdir', 'smlt', 'sf', 'u10', 'v10', 't2m', 'ssrd', 'tp'])
        ds = calc_wind_speeds(ds)
        name = filepath.split('/')[-1].split('_')
        name = name[0] + '_' + name[1] + '_' + 's100' + '_' + name[3] + '.nc'
        ds.to_netcdf(out_path + name)
        print(name)  # just to keep track of progress


def preprocess_CERA_wind(file_list, out_path):
    for filepath in file_list:
        ds = xr.open_dataset(filepath)
        ds = calc_wind_speeds(ds)  # todo I think variable names are different in 20CR u100 is UGRD100 or so
        name = filepath.split('/')[-1].split('_')  # todo this also has to be adjusted
        name = name[0] + '_' + name[1] + '_' + 's100' + '_' + name[3] + '.nc'
        ds.to_netcdf(out_path + name)
        print(name)  # just to keep track of progress



CERA_path = "/cluster/work/apatt/wojan/data/CERA20C/"  # folder containing CERA data as downloaded from ECMWF
CERA_list = sorted(glob.glob(CERA_path+'*instant.nc'))  # instant: accumulation period 1 timestep (wind: not important)
NOAA_path = "/cluster/work/apatt/wojan/data/20CRv3/"
NOAA_list = sorted(glob.glob(CERA_path+'.nc'))

out_path = '/cluster/work/apatt/wojan/renewable_generation/wind/data/'

preprocess_CERA_wind(CERA_list, out_path)
preprocess_NOAA_wind(NOAA_list, out_path)

