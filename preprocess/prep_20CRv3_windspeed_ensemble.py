# Calculates 100m wind speeds for 20CRv3 data and stores all ensemble member in one file per year
import xarray as xr
import glob


def comp_windspeed(ds):
    """
    compute wind speeds from wind components
    :param ds:
    :return:
    """
    ds = ds.squeeze('height', drop=True)  # height is 100m always. No need to keep this dim
    ds['s100'] = (ds['UGRD']**2 + ds['VGRD']**2)**(1./2)
    ds = ds.drop(['UGRD', 'VGRD'])
    ds['s100'].attrs['standard_name'] = '100m_wind_speed'
    ds['s100'].attrs['long_name'] = '100m wind speed'
    ds['s100'].attrs['units'] = 'm s**-1'
    return ds


data_dir = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/20CRv3/"

for year in range(1901, 2010):
    print(year)
    try:
        xr.open_dataset(data_dir + '20CRv3_s100_' + str(year) + '.nc')
    except FileNotFoundError:
        ds_list = []
        for mem in range(1, 9):
            filelist = glob.glob(data_dir + '*' + str(year) + '_mem00' + str(mem) + '*.nc')
            ds = xr.open_mfdataset(filelist, combine='by_coords')
            ds = comp_windspeed(ds)
            ds['number'] = mem
            ds = ds.set_coords('number')
            ds_list.append(ds)
        joined_ds = xr.concat(ds_list, 'number')
        joined_ds.to_netcdf(data_dir + '20CRv3_s100_' + str(year) + '.nc')



