import numpy as np
import xarray as xr
import pandas as pd
from oceans.filters import lanc
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob

mpl.use('Agg')


def essentials_CERA(ds):
    ds = ds.sel({'number': 0}, drop=True)
    ds['s10'] = np.sqrt(ds['u10'] ** 2 + ds['v10'] ** 2)
    ds['s100'] = np.sqrt(ds['u100'] ** 2 + ds['v100'] ** 2)
    ds = ds.drop(['u10', 'v10', 'u100', 'v100'])
    return ds


data_path = '/cluster/work/apatt/wojan/data/'
plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind/plots/'

filelist = glob.glob(data_path + 'CERA20C/' + '*instant.nc')
filelist.sort()

for i in [0, 50, 100]:
    ds = xr.open_dataset(filelist[i])
    year = filelist[i].split('_')[-2]
    ds = ds.drop(['sro', 'fdir', 'smlt', 'sf', 't2m', 'ssrd', 'tp'])
    ds = essentials_CERA(ds)

    plt.clf()
    plt.scatter(ds['s10'].values, ds['s100'].values, alpha=.4)
    plt.xlim(xmin=-1, xmax=30)
    plt.ylim(ymin=-1, ymax=40)
    plt.xlabel('s10 3h values in ' + year + ' [m/s]')
    plt.ylabel('s100 3h values in ' + year + ' [m/s]')
    plt.plot([0,25], [0,25*10**(1./7)], label='power law, **1/7', color='black')
    plt.legend()
    plt.savefig(plot_path + 'CERA20C_s100_s10_scatter_' + year + '.jpeg', dpi=300)

ds = xr.open_mfdataset(filelist,
                       preprocess=essentials_CERA,
                       chunks={'latitude': 1, 'longitude': 1},
                       combine='by_coords',
                       drop_variables=['sro', 'fdir', 'smlt', 'sf', 't2m', 'ssrd', 'tp'])
ds.load()
