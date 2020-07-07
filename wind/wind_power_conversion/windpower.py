import pandas as pd
import xarray as xr
import glob
import numpy as np
import time
import sys


class Power:
    def __init__(self, turbine_index):
        self.power_curve = pd.read_pickle(sorted(glob.glob(power_curve_path + '*.p'))[turbine_index])
        self.turbine_name = self.power_curve.keys()[0].replace('/', '_')

    def power_conversion(self, s):
        s = np.round(s, 2)  # only two decimal accuracy in power curve
        if s < self.power_curve.index[0] or s > self.power_curve.index[-1]:
            # below cut_in or above cut_out
            out = 0
        else:
            out = self.power_curve.loc[s].values[0]
        return out


wind_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/trend_corrected/'
power_curve_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/power_curves/'
out_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/wind_power/'
filelist = sorted(glob.glob(wind_path + '*.nc'))

i = int(sys.argv[1])-1  # year index. Allowed values: 0,1,...,100
wind = xr.open_dataset(filelist[i], chunks={'number': 1, 'time': 100}).load()
year = filelist[i].split('_')[-2]
print(year)

# First and last year contain Nans in first/last half of the year. Remove here
if i in [0, 100]:
    wind_dropped = wind.dropna('time')

for turbine_index in range(3):
    try:
        xr.open_dataset(out_path + P.turbine_name + '/Wind_power_' + year + '.nc')
    except FileNotFoundError:
        P = Power(turbine_index)
        print(P.turbine_name)
        t_0 = time.time()
        wind_power = xr.apply_ufunc(P.power_conversion,
                                    wind['s100'],
                                    vectorize=True,
                                    dask='allowed')
        print('wind power conversion took ' + str(time.time() - t_0))
        t_0 = time.time()
        wind_power.to_netcdf(out_path + P.turbine_name + '/Wind_power_' + year + '.nc')
        print('saving took ' + str(int(time.time() - t_0)) + ' s')
