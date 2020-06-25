import numpy as np
import xarray as xr
import pandas as pd
from oceans.filters import lanc
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
mpl.use('Agg')


def essentials(ds):
    ds = ds.sel({'number': 0, 'latitude': 54.75, 'longitude': 5.25}, drop=True)
    ds['s10'] = np.sqrt(ds['u10']**2 + ds['v10']**2)
    #ds = ds.drop(['sro', 'fdir', 'smlt', 'sf', 't2m', 'ssrd', 'tp', 'u10', 'v10', 'u100', 'v100'])
    ds = ds.drop(['u10', 'v10'])
    return ds


data_path = '/cluster/work/apatt/wojan/data/CERA20C/'
filelist = glob.glob(data_path + '*instant.nc')
filelist.sort()

ds = xr.open_mfdataset(filelist,
                       preprocess=essentials,
                       chunks={'latitude': 1, 'longitude': 1},
                       combine='by_coords',
                       drop_variables=['sro', 'fdir', 'smlt', 'sf', 't2m', 'ssrd', 'tp', 'u100', 'v100'])
ds.load()

var = 's10'
df = ds[var].to_dataframe()

plt.clf()
f, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 6), sharex=True)
#df[var].plot(ax=ax, alpha=.6)

for i, years in enumerate([1, 5, 10, 20]):
    freq = 1. / (8 * 365 * years)
    N = int(1/freq)
    wt = lanc(N, freq)
    res = np.convolve(wt, df[var], mode='same')  # was macht mode='same'?
    res[0:N] = np.nan
    res[-N:] = np.nan
    df['low'] = res
    df['high'] = df[var] - df['low']
    df['low'].plot(ax=ax.flatten()[i], label='Lanczos ' + str(years) + ' y')
    df[var].rolling(window=N, min_periods=N, center=True).mean().plot(ax=ax.flatten()[i], label='rolling mean'+ str(years) + ' y')
    ax.flatten()[i].set_ylim(ymin=6, ymax=10)
    ax.flatten()[i].set_ylabel('Wind speeds [m/s]')
    ax.flatten()[i].legend()
    ax.flatten()[i].set_title('North Sea, lat=54.75N, lon=5.25E, CERA20C')
plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind/plots/'
plt.savefig(plot_path + 'lanc_test_110y_2x2.jpeg')







plt.clf()
f, ax = plt.subplots(ncols=4, sharey=True, figsize=(12,5))
for plot_i, n in enumerate([100, 200, 300, 400]):
    ax[plot_i].set_title('n = '+str(n))
    for f_c in [1./(8*7), 1./(8*7*2), 1./(8*30), 1./(8*90)]:
        wt = lanc(n, f_c)
        ax[plot_i].plot(wt, label='f_c = '+str(np.round(f_c, 4)))

ax[plot_i].legend()
plt.savefig(plot_path + 'lanc_test_parameters.jpeg')

def lanczos_weights(f_c, n):
    """
    Calculation of the Lanczos filter following Eq. 9 in Duchon 1979.

    At k = 0, the computation of the ideal filter weights yields nan. Using
    the series definition of the sine, it can easily be shown that 2f_c is the
    correct value.


    Reference:
    Duchon, C. E. Lanczos Filtering in One and Two Dimensions.
    Journal of Applied Meteorology 18, 1016–1022 (1979).

    :param f_c: cutoff frequency, given per timestep
    :param n: n determines the number of weights that are used (2n +1)
    :return:
    """
    from scipy.special import sinc
    import warnings
    k = np.arange(-n, n + 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ideal_filter_weights = np.sin(2*np.pi*f_c*k)/(np.pi*k)
    ideal_filter_weights[n] = 2*f_c  # limes considerations yield this value for k=0
    sigma = sinc(np.pi*k/n)
    return ideal_filter_weights * sigma #/ (ideal_filter_weights * sigma).sum()

def lanczos_response(f, f_c, n):
    """
    Calculates the response funtion of a Lanczos filter with
    finite considered frequencies

    Implements Eq. 7 in Duchon 1979

    :param f: frequencies at which response function is to be evaluated
    :param f_c: cutoff frequency, given per timestep
    :param n:  n determines the number of weights that are used (2n +1)
    :return:
    """
    w_bar = lanczos_weights(f_c, n)
    w_bar_0 = w_bar[n+1]


plt.clf()
f, ax = plt.subplots(ncols=4, sharey=True, figsize=(12,5))
for plot_i, n in enumerate([100, 200, 300, 400]):
    ax[plot_i].set_title('n = '+str(n))
    for f_c in [1./(8*7), 1./(8*7*2), 1./(8*30), 1./(8*90)]:
        wt = lanczos_weights(f_c, n)
        ax[plot_i].plot(wt, label='f_c = '+str(np.round(f_c, 4)))

ax[plot_i].legend()
plt.savefig(plot_path + 'lanc_test_parameters_self.jpeg')

plt.clf()
f, ax = plt.subplots(ncols=4, sharey=True, figsize=(12,5))

f_c = 3/(5*365*24)  # corresponds to 1/5y
for plot_i, n in enumerate([15000, 20000, 30000, 40000]):
    ax[plot_i].set_title('n = '+str(n))
    wt = lanczos_weights(f_c, n)
    ax[plot_i].plot(wt, label='f_c = '+str(np.round(f_c, 4)))

ax[plot_i].legend()
plt.savefig(plot_path + 'lanc_test_parameters_self_1over5y.jpeg')