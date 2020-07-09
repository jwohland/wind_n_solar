import xarray as xr
from eofs.xarray import Eof
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys
sys.path.append('../../')
from utils import plot_field
import cartopy.crs as ccrs
import numpy as np
import pandas as pd

base_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/'
data_paths = ['output/trend_corrected/annual/', '/data/CERA20C/annual/', 'data/20CRv3/annual/']
plot_paths = ['plots/analysis/eof/annual/trend_corrected/',
              'plots/analysis/eof/annual/CERA20C/',
              'plots/analysis/eof/annual/20CRv3/']

for i, data_path in enumerate(data_paths):
    plot_path = plot_paths[i]
    data_name = plot_path.split('/')[-2]

    filelist = sorted(glob.glob(base_path + data_path + '*.nc'))
    all_wind = xr.open_mfdataset(filelist, combine='by_coords')

    for number in all_wind.number.values:
        wind = all_wind.sel({'number': number})
        solver = Eof(wind['s100'])
        N = 5
        eofs = solver.eofs(neofs=N)
        pcs = solver.pcs(npcs=N)
        variance_fraction = solver.varianceFraction()

        cmap = mpl.cm.get_cmap('RdBu_r')
        plt.close('all')
        f, ax = plt.subplots(ncols=N, nrows=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(12, 5))
        for i in range(5):
            vmax = np.round(np.max([eofs.values.max(), -eofs.values.min()])*1.1, 1)
            plot_field(eofs.sel({'mode': i}, drop=True), ax=ax[0, i],
                       title= str(np.round(variance_fraction.sel({'mode': i}).values*100, 1)) + '% variance ',
                       vmin=-vmax, vmax=vmax, add_colorbar=True if i == N-1 else '',
                       cmap=cmap)
            ax[1, i] = plt.subplot(2, N, N+1+i)  # override the GeoAxes object
            pc = pcs.sel({'mode': i}, drop=True)
            pc = pd.Series(data=pc.values, index=pc.time.values)
            pc.plot(ax=ax[1, i])
            pc.rolling(window=5, center=True).mean().plot(ax=ax[1, i], ls='--', color='black', lw=2)
        plt.subplots_adjust(left=0.05, right=0.92)
        plt.suptitle(data_name)
        plt.savefig(base_path + plot_path + 's100_eofs_' + str(int(number)) + '.png')