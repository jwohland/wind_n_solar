import numpy as np
import xarray as xr
import pandas as pd
from oceans.filters import lanc
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

mpl.use('Agg')


def essentials_CERA(ds):
    ds = ds.sel({'number': 0, 'latitude': 54.75, 'longitude': 5.25}, drop=True)
    ds['s10'] = np.sqrt(ds['u10'] ** 2 + ds['v10'] ** 2)
    ds = ds.drop(['u10', 'v10'])
    return ds


def essentials_20CRv2c(ds):
    ds = ds.sel({'lat': 54.2846, 'lon': 5.625}, drop=True)
    ds['s10'] = np.sqrt(ds['u10'] ** 2 + ds['v10'] ** 2)
    ds = ds.drop(['u10', 'v10'])
    return ds


def essentials_20CRv3(ds):
    ds = ds.sel({'lat': 55.087563, 'lon': 4.921, 'height': 10}, drop=True)
    ds = ds.rename({'WSPD10m': 's10'})
    ds = ds.drop(['time_bnds'])
    return ds


def filter_and_plot(df, rea_name):
    plt.clf()
    f, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 6), sharex=True)
    var = 's10'
    for i, years in enumerate([1, 5, 10, 20]):
        freq = 1. / (8 * 365 * years)
        if rea_name == '20CRv3':  # monthly data
            freq = 1. / (12 * years)
        N = int(1 / freq)
        wt = lanc(N, freq)
        res = np.convolve(wt, df[var], mode='same')  # mode='same': length of filtered signal identical to unfiltered
        res[0:N] = np.nan  # decided to keep all values first (simpler for mapping to timesteps) and set to nan here
        res[-N:] = np.nan  #
        df['low'] = res
        df['low'].plot(ax=ax.flatten()[i], label='Lanczos ' + str(years) + ' y')
        df[var].rolling(window=N,
                        min_periods=N,
                        center=True).mean().plot(ax=ax.flatten()[i], label='rolling mean' + str(years) + ' y')
        ax.flatten()[i].set_ylim(ymin=6, ymax=10)
        ax.flatten()[i].set_ylabel('10m Wind speeds [m/s]')
        ax.flatten()[i].legend()
        ax.flatten()[i].set_title('North Sea, lat=54.75N, lon=5.25E, ' + rea_name)
    plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind/plots/'
    plt.savefig(plot_path + rea_name + '_lanc_test_110y_2x2.jpeg')


def do(rea_name):
    if rea_name == 'CERA20C':
        try:
            ds = xr.open_dataset(out_path + 'CERA20C_testlocation_northsea.nc')
        except FileNotFoundError:
            filelist = glob.glob(data_path + 'CERA20C/' + '*instant.nc')
            filelist.sort()

            ds = xr.open_mfdataset(filelist,
                                   preprocess=essentials_CERA,
                                   chunks={'latitude': 1, 'longitude': 1},
                                   combine='by_coords',
                                   drop_variables=['sro', 'fdir', 'smlt', 'sf', 't2m', 'ssrd', 'tp', 'u100', 'v100'])
            ds.load()
            ds.to_netcdf(out_path + 'CERA20C_testlocation_northsea.nc')

        var = 's10'
        df = ds[var].to_dataframe()
        filter_and_plot(df, 'CERA20C')
    elif rea_name == '20CRv2c':
        """
        20CRv2c
        !!! Not sure about the wind components, they are ensemble means, I think, and could show large trends due to 
        decreased ensemble spread !!!
        """
        try:
            ds = xr.open_dataset(out_path + '20CRv2c_testlocation_northsea.nc')
        except FileNotFoundError:
            filelist = glob.glob(data_path + '20CRv2c/' + '*.nc')
            filelist.sort()

            ds = xr.open_mfdataset(filelist,
                                   preprocess=essentials_20CRv2c,
                                   chunks={'lat': 10, 'lon': 10},
                                   combine='by_coords',
                                   drop_variables=['t2m', 't2m_spread',
                                                   'ssrd', 'ssrd_spread',
                                                   'tp', 'tp_spread',
                                                   'sro', 'sro_spread',
                                                   'u10_spread', 'v10_spread'])
            ds.load()
            ds.to_netcdf(out_path + '20CRv2c_testlocation_northsea.nc')

        var = 's10'
        df = ds[var].to_dataframe()
        df = df['1901':'2009']
        filter_and_plot(df, '20CRv2c')

    elif rea_name == '20CRv3':
        """
        20CRv3
        !!! Not sure about the wind components, they are ensemble means, I think and could show large trends due to 
        decreased ensemble spread !!!
        """
        try:
            ds = xr.open_dataset(out_path + '20CRv3_testlocation_northsea.nc')
        except FileNotFoundError:
            ds = xr.open_dataset(
                data_path + '20CRv3/wspd/monthly/001/WSPD.mnmean_mem001.nc')  # todo this is one ensemble member
            ds = essentials_20CRv3(ds)
            ds.load()
            ds.to_netcdf(out_path + '20CRv3_testlocation_northsea.nc')
        var = 's10'
        df = ds[var].to_dataframe()
        filter_and_plot(df, '20CRv3')


    elif rea_name == 'both':
        try:
            ds_CERA = xr.open_dataset(out_path + 'CERA20C_testlocation_northsea.nc')
            ds_20CRv2c = xr.open_dataset(out_path + '20CRv2c_testlocation_northsea.nc')
            ds_20CRv3 = xr.open_dataset(out_path + '20CRv3_testlocation_northsea.nc')
            plt.clf()
            f, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 6), sharex=True, sharey=True)
            var = 's10'
            for j, df in enumerate([ds_CERA[var].to_dataframe(),
                                    ds_20CRv2c[var].to_dataframe(),
                                    ds_20CRv3[var].to_dataframe()]):
                #df = df['1901':'2009']
                rea_name = ['CERA20C', '20CRv2c', '20CRv3'][j]
                for i, years in enumerate([1, 5, 10, 20]):
                    freq = 1. / (8 * 365 * years)
                    if rea_name == '20CRv3':
                        freq = 1. / (12 * years)
                    N = int(1 / freq)
                    wt = lanc(N, freq)
                    res = np.convolve(wt, df[var],
                                      mode='same')  # mode='same': length of filtered signal identical to unfiltered
                    res[0:N] = np.nan  # decided to keep all values first (simpler for mapping to timesteps) and
                    res[-N:] = np.nan  # set to nan here
                    df['low'] = res
                    ax.flatten()[i].plot(df.index, df['low'], label=rea_name)
                    #df['low'].plot(x=df.index, ax=ax.flatten()[i], label=rea_name)
                    ax.flatten()[i].set_xlim(xmin='1900', xmax='2009')
                    ax.flatten()[i].set_ylim(ymin=6, ymax=10)
                    ax.flatten()[i].set_ylabel('10m wind speeds [m/s]')
                    ax.flatten()[i].set_title('Lanczos ' + str(years) + ' years')
                ax.flatten()[1].legend()
            plt.suptitle('Single location middle of North Sea')
            plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind/plots/'
            plt.savefig(plot_path + 'comparison_lanc_test_110y_2x2.jpeg')
        except FileNotFoundError:
            print('data not available')

    elif rea_name == 'difference':
        try:
            ds_CERA = xr.open_dataset(out_path + 'CERA20C_testlocation_northsea.nc')
            ds_20CRv2c = xr.open_dataset(out_path + '20CRv2c_testlocation_northsea.nc')
            ds_20CRv3 = xr.open_dataset(out_path + '20CRv3_testlocation_northsea.nc')
            var = 's10'
            filtered_dic = {}
            for j, df in enumerate([ds_CERA[var].to_dataframe(),
                                    ds_20CRv2c[var].to_dataframe(),
                                    ds_20CRv3[var].to_dataframe()]):

                rea_name = ['CERA20C', '20CRv2c', '20CRv3'][j]
                years = 5
                freq = 1. / (8 * 365 * years)
                if rea_name == '20CRv3':
                    freq = 1. / (12 * years)
                N = int(1 / freq)
                wt = lanc(N, freq)
                res = np.convolve(wt, df[var],
                                  mode='same')  # mode='same': length of filtered signal identical to unfiltered
                res[0:N] = np.nan  # decided to keep all values first (simpler for mapping to timesteps) and
                res[-N:] = np.nan  # set to nan here
                df['low'] = res
                filtered_dic[rea_name] = df

            plt.clf()
            f, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
            #for crea in rea_name[1:]:
            if 1 == 1:
                crea = '20CRv2c'
                diff = filtered_dic[crea]['low'] - filtered_dic['CERA20C']['low']
                diff = diff['1906':'2004']
                corrected = filtered_dic['CERA20C']['s10']['1906':'2004'] + diff
                ax[0].plot(diff.index, diff, label=crea)
                ax[0].set_ylabel('wind speed label - CERA20C')
                ax[0].set_title('Lanczos ' + str(years) + ' years')
                ax[0].legend()
                ax[1].scatter(filtered_dic['CERA20C']['s10']['1906':'2004'], corrected, s=1, alpha=.3)
                ax[1].set_xlabel('Original CERA20C 3h wind speeds [m/s]')
                ax[1].set_ylabel('20CRv2c corrected CERA20C 3h wind speeds [m/s]')
                ax[1].plot([0, 25], [0, 25], '--', color='black')
                plt.suptitle('Single location middle of North Sea')
                plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind/plots/'
                plt.savefig(plot_path + 'diff_lanc_test.jpeg')
        except FileNotFoundError:
            print('data not available')


data_path = '/cluster/work/apatt/wojan/data/'
out_path = '/cluster/work/apatt/wojan/renewable_generation/wind/output/'

do('both')
