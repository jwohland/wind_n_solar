import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import glob
import numpy as np
import sys

sys.path.append('../../')
from utils import plot_field
import cartopy.crs as ccrs

wind_power_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/wind_power/'
plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/plots/wind_power/'

filelist = glob.glob(wind_power_path + '*/*1914.nc')  # example data from 1914

f, ax = plt.subplots(ncols=3, figsize=(12, 5), subplot_kw={'projection': ccrs.PlateCarree()})
cbar_ax = f.add_axes([0.2, 0.1, 0.6, 0.05])

for i, filename in enumerate(filelist):
    wind_power = xr.open_dataset(filename)
    levels = np.arange(0.1, 0.61, .025)
    cmap = mpl.cm.get_cmap('Greens')
    mean_windpower = wind_power['s100'].mean(dim=['time', 'number'])  # todo s100 -> wind_power after test run
    if i != 1:
        plot_field(mean_windpower,
                   ax[i],
                   filename.split('/')[-2],
                   cmap=cmap,
                   levels=levels,
                   add_colorbar=False)
    else:
        plot_field(mean_windpower,
                   ax[i],
                   filename.split('/')[-2],
                   cmap=cmap,
                   levels=levels,
                   cbar_ax=cbar_ax,
                   cbar_kwargs={'orientation': 'horizontal'})
    mean_windpower.plot.contourf(ax=ax[i],
                                 levels=[0, 0.1], hatches=["...", ""],
                                 alpha=0,
                                 add_colorbar=False,
                                 add_labels=False)

plt.subplots_adjust(0.03, 0.15, 0.97, 0.95)
plt.savefig(plot_path + 'Mean_wind_power_1914_test.pdf')
