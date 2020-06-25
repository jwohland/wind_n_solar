# Calculates 100m wind speeds for 20CRv3 data and stores all ensemble member in one file per year
import xarray as xr
import glob


def comp_windspeed(ds):
    """
    compute wind speeds from wind components
    :param ds:
    :return:
    """
    ds['s100'] = (ds['u100']**2 + ds['v100']**2)**(1./2)
    ds = ds.drop(['u100', 'v100'])
    ds['s100'].attrs['standard_name'] = '100m_wind_speed'
    ds['s100'].attrs['long_name'] = '100m wind speed'
    ds['s100'].attrs['units'] = 'm s**-1'
    return ds


data_dir = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/CERA20C/"
filelist = sorted(glob.glob("/cluster/work/apatt/wojan/data/CERA20C/" + '*instant.nc'))

drop_variables = ['sro',
                  'smlt',
                  'sf',
                  'u10', 'v10',
                  'tp']

for filename in filelist:
    year = filename.split('_')[-2]
    ds = xr.open_dataset(filename, drop_variables=drop_variables)
    ds = comp_windspeed(ds)
    ds.to_netcdf(data_dir + 'CERA20C_multiple_' + str(year) + '.nc')
