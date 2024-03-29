# Calculation of PV generation from CERA20C data
# !! As of July 16th 2020, this needs a branched version of GSEE (fix-geometry-around-sunset) rather than master !!

import xarray as xr
import glob
import time
import sys
import gsee.pv as pv  # branched case
import numpy as np
import os


def open_CERA(year_index, number):
    """
    Opens CERA20C 3 hourly data for one year and one realization and modifies names and units in line with GSEE.
    :param year_index:
    :param number:
    :return:
    """
    data_path = (
        "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/CERA20C/"
    )
    filelist = sorted(glob.glob(data_path + "*.nc"))
    data = xr.open_dataset(filelist[year_index], drop_variables=["s100"],)
    data = data.rename(
        {
            "latitude": "lat",
            "longitude": "lon",
            "ssrd": "global_horizontal",
            "fdir": "direct_radiation",
            "t2m": "temperature",
        }
    )
    data = data.sel({"number": number}, drop=True)
    data = convert_units(data)
    data = provide_GSEE_fields(data)
    return data


def convert_units(data):
    """
    CERA20C data comes in wrong units for GSEE:
           J/m**2 rather than W/m**2 and gives 3 hours cumulative values
           Temperature in Kelvin rather than Celsius
    "This parameter is accumulated over a particular time period which depends on the data extracted.
    The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the
    accumulated values should be divided by the accumulation period expressed in seconds." ECMWF, Parameter database
    https://apps.ecmwf.int/codes/grib/param-db/?id=169
    :param data:
    :return:
    """
    secs_per_3h = 60 * 60 * 3.0
    data["global_horizontal"] /= secs_per_3h
    data["direct_radiation"] /= secs_per_3h
    data["temperature"] -= 273.15
    return data


def provide_GSEE_fields(data):
    """
    GSEE needs diffuse fraction as input which is not directly provided by ECMWF. Calculation based on
        global horizontal = direct_radiation + diffuse_radiation
        direct_fraction = direct_radiation / global_horizontal
        direct_fraction + diffuse_fraction = 1
    :param data:
    :return:
    """
    diff_frac = 1 - (data["direct_radiation"] / data["global_horizontal"]).fillna(0)
    # division by zero leads to lots of nans
    diff_frac = diff_frac.where(diff_frac > 0, 0)
    diff_frac = diff_frac.where(diff_frac < 1, 0)
    # remove inf and -inf that sometimes apeear when direct_radiation approximately zero
    data["diffuse_fraction"] = diff_frac
    data = data.drop("direct_radiation")
    return data


def prep_dataset_gsee(ds):
    """
    Preparation for CERA20C data to be passed to gsee.pv.run_model

    Involves (a) conversion to a pandas dataframe and (b) meaningful shifting/interpolation of the data.

    Regarding b:
        ECMWF radiation data is cumulative. Value at time t corresponds to average from (t - delta t) to (t).
        GSEE using irradiance_type == 'cumulative' expects centered means which is achieved by shifting the time axis
            t ' = t - 01:30

        ECMWF temperature data is instantaneous. To calculate meaningful temperatures at the new time steps t'
        (04:30, 07:30, 10:30 etc) from the old time steps t (03:00, 06:00, 09:00), the values are linearly interpolated
            Example: T ( t = 04:30 )   =  0.5 * [ T ( t = 03:00 ) + T ( t = 06:00 ) ]
                                       =  0.5 * [ T ( t' = 01:30 ) + T ( t' = 04:30 ) ]

    :param ds: dataset containing temperature and radiation (direct and diffuse)
    :return:
    """
    tmp_df = ds.to_dataframe()
    tmp_df.index = tmp_df.index.shift(-1, (tmp_df.index[1] - tmp_df.index[0]) / 2)
    tmp_df["temperature"] = (
        tmp_df["temperature"].rolling(window=2, min_periods=1).mean()
    )
    return tmp_df


def draw_sample(mean, std, shape, seed):
    """
    Draw a sample for the tilt and azimuth angles of the PV panels based on gaussians.
    
    Follows the approach by Pfenninger & Staffell, 2016, Energy (last sentences in Sec 2.2.)

    :param mean: mean value
    :param std: standard deviation
    :param shape: shape of the output
    :return: 
    """
    np.random.seed(seed)
    return np.random.normal(loc=mean, scale=std, size=shape)


def get_parameters(out_path, scenario, shape):
    if scenario == 0:
        out_path += "both_constant/"
        params = {
            "tilt": 25,
            "azim": 180,
            "tracking": 0,
            "capacity": 1,
        }  # tilt & azimuth converted to radians in pv.run_model
    elif scenario == 1:
        out_path += "tilt_constant/"
        params = {
            "tilt": 25,
            "azim_mean": 180,
            "azim_std": 40,
            "tracking": 0,
            "capacity": 1,
        }  # tilt & azimuth converted to radians in pv.run_model
    elif scenario == 2:
        out_path += "neither_constant/"
        params = {
            "tilt_mean": 25,
            "tilt_std": 15,
            "azim_mean": 180,
            "azim_std": 40,
            "tracking": 0,
            "capacity": 1,
        }  # tilt & azimuth converted to radians in pv.run_model
    # draw samples for tilt and azimuth
    if "tilt" in params.keys():
        params["tilt_array"] = np.zeros(shape) + params["tilt"]
    else:
        params["tilt_array"] = draw_sample(
            params["tilt_mean"], params["tilt_std"], shape, seed=0
        )
    if "azim" in params.keys():
        params["azim_array"] = np.zeros(shape) + params["azim"]
    else:
        params["azim_array"] = draw_sample(
            params["azim_mean"], params["azim_std"], shape, seed=1
        )

    return out_path, params


def run_GSEE_CERA(year_index=57, number=0, scenario=0):
    """
    Calls GSEEs pv.run_model() for all latlon pairs in the input data.
    Results are converted into a xarray Dataset and saved.
    :param year_index:
    :param number: CERA20C ensemble member (allowed values: 0 to 9)
    :param scenario: scenario index (0 is constant tilt and azimuth, 1 is constant tilt, 2 is tilt and azimuth non constant)
    :return:
    """
    out_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/solar_power/"

    test_data = open_CERA(year_index, number=number)
    lats, lons = test_data.lat.values, test_data.lon.values
    shape = (lats.size, lons.size)
    out_path, params = get_parameters(out_path, scenario, shape)
    output_name = "Solar_power_" + str(year_index + 1901) + "_number_" + str(number) + "_scenario_" + str(scenario) + ".nc"

    if not os.path.isfile(out_path + output_name):
        # loop over all locations and call GSEE pv.runmodel()
        lat_list = []
        for i_lat, lat in enumerate(lats):
            print(lat)
            start = time.time()
            lon_list = []
            for i_lon, lon in enumerate(lons):
                tmp_data = test_data.sel({"lat": lat, "lon": lon})
                tmp_df = prep_dataset_gsee(tmp_data)

                tmp_df["PV"] = pv.run_model(
                    data=tmp_df,
                    coords=(lat, lon),
                    tilt=params["tilt_array"][i_lat, i_lon],  # 30 degrees tilt angle
                    azim=params["azim_array"][i_lat, i_lon],  # facing towards equator,
                    tracking=params["tracking"],  # fixed - no tracking
                    capacity=params["capacity"],  # 1 W
                    irradiance_type="cumulative",  # CERA data is cumulative
                )

                tmp_df.drop(
                    ["temperature", "global_horizontal", "diffuse_fraction"],
                    axis=1,
                    inplace=True,
                )
                lon_list.append(
                    tmp_df.reset_index()
                    .pivot_table(values="PV", index=["time", "lat", "lon"])
                    .to_xarray()
                )
            lat_list.append(xr.concat(lon_list, dim="lon"))
            end = time.time()
            print(str(int(end - start)) + " s")  # around 25s each
        pv_data = xr.concat(lat_list, dim="lat")
        for var in ["tracking", "capacity"]:
            pv_data[var] = params[var]
        for var in ["tilt_array", "azim_array"]:
            pv_data[var] = xr.DataArray(
                dims=["lat", "lon"], coords=(lats, lons), data=params[var]
            )
        pv_data.to_netcdf(out_path + output_name)


try:
    year_index = (
        int(sys.argv[1]) - 1
    )  # JOB IDs start at 1 and not at 0, allowed values for year_index 0..108
    number = int(sys.argv[2]) - 1  # allowed values 0..9
    scenario = int(sys.argv[3]) - 1  # allowed values 0..2
except IndexError:
    year_index = 0
    number = 0
    scenario = 0

print(
    "year_index = "
    + str(year_index)
    + ", number = "
    + str(number)
    + ", scenario = "
    + str(scenario)
)
run_GSEE_CERA(year_index, number, scenario)
