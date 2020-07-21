# utils.py
import xarray as xr
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=ImportWarning)
    import cartopy
    import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['axes.spines.left'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = False


def load_data(datapath, dataset, chunks):
    if '20CR' in dataset:
        if 'N47' in dataset:
            data = xr.open_mfdataset(datapath + dataset + "/*.nc", combine='by_coords',
                                     chunks=chunks)
        else:
            data = xr.open_mfdataset(datapath + dataset + "/*.nc", combine='by_coords',
                                     chunks=chunks, preprocess=european_20CR)
            try:
                memberdata = xr.open_mfdataset(datapath + dataset + "/members/s*.nc", combine='by_coords',
                                     chunks=chunks, preprocess=european_20CR)
            except OSError:
                print('member data not available')
            else:
                data = data.merge(memberdata, join='override')
    elif 'ERA20CM' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=ensemble_mean)
    elif 'CERA20C' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*instant*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=ensemble_mean)
    elif 'ERA20C' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*instant*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=add_windspeeds)
    else:
        data = xr.open_mfdataset(datapath + dataset + "/*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=add_windspeeds)
    return data


def load_data_10y(datapath, dataset, chunks):
    """
    short version of load_data, for trend testing only
    :param datapath:
    :param dataset:
    :param chunks:
    :return:
    """
    if '20CR' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*197*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=european_20CR)
    elif 'ERA20CM' in dataset or 'CERA20C' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*197*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=ensemble_mean)
    else:
        data = xr.open_mfdataset(datapath + dataset + "/*197*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=add_windspeeds)
    return data


def load_data_1y(datapath, dataset, chunks):
    """
    correlation testing only, always uses N47 data
    :param datapath:
    :param dataset:
    :param chunks:
    :return:
    """
    if '20CR' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*1979*.nc", combine='by_coords',
                                     chunks=chunks)
    elif 'ERA20CM' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*1979*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=ensemble_mean)
    elif 'CERA20C' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*1979*instant*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=ensemble_mean)
    elif 'ERA20C' in dataset:
        data = xr.open_mfdataset(datapath + dataset + "/*1979*instant*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=add_windspeeds)
    else:
        data = xr.open_mfdataset(datapath + dataset + "/*1979*.nc", combine='by_coords',
                                 chunks=chunks, preprocess=add_windspeeds)
    return data


def windspeeds(u, v):
    """
    calculates wind speeds from the wind components
    :param u: east-west component
    :param v: north-south component
    :return: wind speeds in same units as imput data
    """
    with xr.set_options(keep_attrs=True):
        s = (u ** 2 + v ** 2) ** (1. / 2)
    return s


def add_windspeeds(inxarray):
    """
    adds 10m and 100m wind speeds to the xarray
    :param inxarray:
    :return:
    """
    for height in ['10', '100']:
        try:
            inxarray['s' + height] = windspeeds(inxarray['u' + height], inxarray['v' + height])
            update_attrs(inxarray['s' + height], 'm s**-1', 'ensemble mean ' + height + ' m wind speed')
        except KeyError:
            pass  # 100 m not available in all datasets
    return inxarray


def update_attrs(inxarray, unitname, varname):
    """
    Updates metadata of an xarray, for example after trend calculation
    :param inxarray:
    :param unitname: new unit
    :param varname: new name of the variable
    :return:
    """
    inxarray.attrs['units'] = unitname
    inxarray.attrs['long_name'] = varname


def european_20CR(inxarray):
    """
    removes all grid boxes outside of the European domain
    :param inxarray: 20CR global dataset
    :return: 20CR European dataset
    """
    inxarray = inxarray.assign_coords(lon=(((inxarray.lon + 180) % 360) - 180))
    inxarray = inxarray.sortby(inxarray.lon)
    inxarray = inxarray.sortby(inxarray.lat)
    lat_bnds, lon_bnds = [30, 75], [-15, 35]
    inxarray = inxarray.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
    #add_windspeeds(inxarray)
    print('Preprocessing done')
    return inxarray


def ensemble_mean(inxarray):
    """
    calculates the ensemble mean and modifies the variable names to account for this
    :param inxarray:
    :return:
    """
    try:
        add_windspeeds(inxarray)
    except KeyError:
        print('Wind not available')
    out = inxarray.mean('number', keep_attrs=True)
    for var in out.var():
        update_attrs(out[var], unitname=out[var].attrs['units'], varname='ensemble mean ' + out[var].attrs['long_name'])
    return out


def plot_field(data, ax=None, title=None, **kwargs):
    """
    plots maps
    :param data: data xarray that shall be plotted
    :param ax: axes
    :param title: figure title
    :return:
    """
    ax = ax or plt.axes(projection=ccrs.PlateCarree())
    data.plot(ax=ax, **kwargs)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))
    ax.set_title(title)
    return ax


def add_letters(ax, x=-.1, y=1.03, fs=10):
    """
    adds bold letters a,b,c,... to the upper left corner of subplots
    :param ax: axis
    :param x: x location of text
    :param y: ylocation of text
    :param fs: fontsize
    :return:
    """
    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
    for il, tmp_ax in enumerate(ax.flat):
        tmp_ax.text(x, y, letters[il], weight='bold', transform=tmp_ax.transAxes, fontsize=fs)


def load_annual_solar_ensemble(base_path, data_path):
    """
    Loads the annual solar pv generation ensemble which is stored as individual files per year and ensemble member
    :param base_path:
    :param data_path:
    :return:
    """
    all_power_list = []
    for number in range(10):
        filelist = sorted(
            glob.glob(base_path + data_path + "annual/*number_" + str(number) + ".nc")
        )
        all_power_list.append(xr.open_mfdataset(filelist, combine="by_coords"))
    return xr.concat(all_power_list, pd.Index(range(10), name="number"))
