# utils.py
import xarray as xr
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy
    import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
import pandas as pd
import matplotlib as mpl
from scipy.stats import spearmanr, pearsonr

mpl.use("Agg")
mpl.rcParams["axes.spines.left"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.bottom"] = False


def load_data(datapath, dataset, chunks):
    if "20CR" in dataset:
        if "N47" in dataset:
            data = xr.open_mfdataset(
                datapath + dataset + "/*.nc", combine="by_coords", chunks=chunks
            )
        else:
            data = xr.open_mfdataset(
                datapath + dataset + "/*.nc",
                combine="by_coords",
                chunks=chunks,
                preprocess=european_20CR,
            )
            try:
                memberdata = xr.open_mfdataset(
                    datapath + dataset + "/members/s*.nc",
                    combine="by_coords",
                    chunks=chunks,
                    preprocess=european_20CR,
                )
            except OSError:
                print("member data not available")
            else:
                data = data.merge(memberdata, join="override")
    elif "ERA20CM" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=ensemble_mean,
        )
    elif "CERA20C" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*instant*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=ensemble_mean,
        )
    elif "ERA20C" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*instant*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=add_windspeeds,
        )
    else:
        data = xr.open_mfdataset(
            datapath + dataset + "/*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=add_windspeeds,
        )
    return data


def load_data_10y(datapath, dataset, chunks):
    """
    short version of load_data, for trend testing only
    :param datapath:
    :param dataset:
    :param chunks:
    :return:
    """
    if "20CR" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*197*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=european_20CR,
        )
    elif "ERA20CM" in dataset or "CERA20C" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*197*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=ensemble_mean,
        )
    else:
        data = xr.open_mfdataset(
            datapath + dataset + "/*197*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=add_windspeeds,
        )
    return data


def load_data_1y(datapath, dataset, chunks):
    """
    correlation testing only, always uses N47 data
    :param datapath:
    :param dataset:
    :param chunks:
    :return:
    """
    if "20CR" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*1979*.nc", combine="by_coords", chunks=chunks
        )
    elif "ERA20CM" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*1979*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=ensemble_mean,
        )
    elif "CERA20C" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*1979*instant*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=ensemble_mean,
        )
    elif "ERA20C" in dataset:
        data = xr.open_mfdataset(
            datapath + dataset + "/*1979*instant*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=add_windspeeds,
        )
    else:
        data = xr.open_mfdataset(
            datapath + dataset + "/*1979*.nc",
            combine="by_coords",
            chunks=chunks,
            preprocess=add_windspeeds,
        )
    return data


def windspeeds(u, v):
    """
    calculates wind speeds from the wind components
    :param u: east-west component
    :param v: north-south component
    :return: wind speeds in same units as imput data
    """
    with xr.set_options(keep_attrs=True):
        s = (u ** 2 + v ** 2) ** (1.0 / 2)
    return s


def add_windspeeds(inxarray):
    """
    adds 10m and 100m wind speeds to the xarray
    :param inxarray:
    :return:
    """
    for height in ["10", "100"]:
        try:
            inxarray["s" + height] = windspeeds(
                inxarray["u" + height], inxarray["v" + height]
            )
            update_attrs(
                inxarray["s" + height],
                "m s**-1",
                "ensemble mean " + height + " m wind speed",
            )
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
    inxarray.attrs["units"] = unitname
    inxarray.attrs["long_name"] = varname


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
    # add_windspeeds(inxarray)
    print("Preprocessing done")
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
        print("Wind not available")
    out = inxarray.mean("number", keep_attrs=True)
    for var in out.var():
        update_attrs(
            out[var],
            unitname=out[var].attrs["units"],
            varname="ensemble mean " + out[var].attrs["long_name"],
        )
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
    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"))
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"))
    ax.set_title(title)
    return ax


def add_letters(ax, x=-0.1, y=1.03, fs=10):
    """
    adds bold letters a,b,c,... to the upper left corner of subplots
    :param ax: axis
    :param x: x location of text
    :param y: ylocation of text
    :param fs: fontsize
    :return:
    """
    letters = [
        "a",
        "b",
        "c",
        "d",
        "e",
        "f",
        "g",
        "h",
        "i",
        "j",
        "k",
        "l",
        "m",
        "n",
        "o",
        "p",
    ]
    for il, tmp_ax in enumerate(ax.flat):
        tmp_ax.text(
            x, y, letters[il], weight="bold", transform=tmp_ax.transAxes, fontsize=fs
        )


def load_solar_ensemble(base_path, data_path, annual=True):
    """
    Loads the annual solar pv generation ensemble which is stored as individual files per year and ensemble member
    :param base_path:
    :param data_path:
    :return:
    """
    all_power_list = []
    for number in range(10):
        if annual:
            filelist = sorted(
                glob.glob(
                    base_path + data_path + "annual/*number_" + str(number) + "*.nc"
                )
            )
        else:
            filelist = sorted(
                glob.glob(base_path + data_path + "*198*number_" + str(number) + "*.nc") +
                glob.glob(base_path + data_path + "*199*number_" + str(number) + "*.nc")
            )
        all_power_list.append(xr.open_mfdataset(filelist, combine="by_coords"))
    return xr.concat(all_power_list, pd.Index(range(10), name="number"))


class Generation_type:
    def __init__(self, name, scenarios, data_path, plot_path, base_path):
        """
        Class to handle wind and solar generation considering different turbines or tilt/azimuth angles
        :param name: wind or solar
        :param scenarios: wind turbine name or name of the tilt/azimuth angle combination
        :param data_path: path to data, relative to base_path
        :param plot_path: path to plots, relative to base_path
        :param base_path: path to directory containing code and plots
        """
        self.name = name
        self.scenarios = scenarios
        self.data_path = base_path + data_path
        self.plot_path = base_path + plot_path
        self.base_path = base_path
        self.all_power = {}
        self.monthly_power = {}
        if name == "solar":
            self.var = "PV"
            self.CF_threshold = 0.1
        else:
            self.var = "wind_power"
            self.CF_threshold = 0.15

    def get_filelist(self, scenario, annual=True):
        if annual:
            return sorted(glob.glob(self.data_path + scenario + "/annual/*.nc"))
        else:
            return sorted(
                glob.glob(self.data_path + scenario + "/*198*.nc")
                + glob.glob(self.data_path + scenario + "/*199*.nc")
            )

    def open_data(self, scenario):
        # check if already calculated
        if not scenario in self.all_power:
            if self.name == "wind":
                self.all_power[scenario] = xr.open_mfdataset(
                    self.get_filelist(scenario), combine="by_coords"
                ).rename({"latitude": "lat", "longitude": "lon"})
            else:
                self.all_power[scenario] = load_solar_ensemble(
                    self.data_path, scenario + "/"
                )
        return self.all_power[scenario]

    def open_data_monthly_ensmean(self, scenario):
        """
        Open 3h data from 1980 to 2000, calculate monthly mean and ensemble mean
        :param scenario:
        :return:
        """
        if not scenario in self.monthly_power:
            if self.name == "wind":
                tmp = xr.open_mfdataset(
                    self.get_filelist(scenario, annual=False), combine="by_coords"
                ).rename({"latitude": "lat", "longitude": "lon"})
            else:
                tmp = load_solar_ensemble(self.data_path, scenario + "/", annual=False)
            tmp = tmp.mean(["number"])  # ensemble mean
            self.monthly_power[scenario] = tmp.resample(time="1MS").mean(
                dim="time",
            )  # monthly mean
        return self.monthly_power[scenario]

    def reset_plot_path(self, plot_path):
        self.plot_path = self.base_path + plot_path

    def get_mask(self, scenario):
        mask_path = self.base_path + "output/optimization/masks/"
        mask_name = "mask_" + scenario
        try:
            return xr.open_dataset(mask_path + mask_name + ".nc")
        except FileNotFoundError:
            print("Mask not computed yet.")


def corr_pearson(ts1, ts2):
    """
    calculate Pearson correlation
    :param ts1: first timeseries
    :param ts2: second timeseries
    :return:
    """
    r = pearsonr(ts1, ts2)[0]
    return r


def corr_spearman(ts1, ts2):
    """
    calculate spearman correlation
    :param ts1: first timeseries
    :param ts2: second timeseries
    :return:
    """
    rho = spearmanr(ts1, ts2)[0]
    return rho


def calc_correlations_xarray(ds1, ds2, var1, var2, measure):
    """
    computes a map of correlation between the datasets variables
    :param ds1: first dataset
    :param ds2: second dataset
    :param var1: variable of interest in first dataset
    :param var2: variable of interest in second dataset
    :param measure: type of correlation pearson or spearman
    :return:
    """
    ds1, ds2 = xr.align(ds1, ds2)
    if measure == "pearson":
        res = xr.apply_ufunc(
            corr_pearson,
            ds1[var1][:, :, :],
            ds2[var2][:, :, :],
            input_core_dims=[["time"], ["time"]],
            vectorize=True,
            dask="allowed",
        )
    elif measure == "spearman":
        res = xr.apply_ufunc(
            corr_spearman,
            ds1[var1][:, :, :],
            ds2[var2][:, :, :],
            input_core_dims=[["time"], ["time"]],
            vectorize=True,
            dask="allowed",
        )
    else:
        warnings.warn("correlation type ill defined")
    return res


def add_row_label(ax, text):
    ax.text(
        -0.1,
        0.5,
        text,
        horizontalalignment="center",
        rotation=90,
        verticalalignment="center",
        transform=ax.transAxes,
    )
