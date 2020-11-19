# Make essential data available on Zenodo. Storage is limited to 50GB. Consequently, the full dataset can not be
# made available. Instead, we provide
    # a) Ensemble mean PV and wind power generation on a gridbox level
    # b) Country group timeseries
# The required postprocessing is performed below. Data is stored in
    #  /cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/zenodo
# and copied to the zenodo directory.

import xarray as xr
import glob
import pandas as pd

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/"
zen_path = base_path + "zenodo/"

"""
Calculate solar power ensemble mean 
"""

for solar_scenario in ["both_constant", "tilt_constant", "neither_constant"]:
    solar_path = zen_path + "Solar/" + solar_scenario + "/"
    for year in range(1901, 2010):
        file_list = glob.glob(base_path + "solar_power/" + solar_scenario + "/Solar*" + str(year) + "*.nc")
        assert len(file_list) == 10  # 10 ensemble members
        ds_list = []
        for file_name in file_list:
            ds_list.append(xr.open_dataset(file_name))
        ds = xr.concat(ds_list, pd.Index(range(10), name="number")).mean("number")
        # store as Solar_power_YYYY_ensmean.nc in zenodo/Solar/solar_scenario
        ds.to_netcdf(solar_path + file_name.split("/")[-1].split("_number")[0]+"_ensmean.nc")

"""
Calculate wind power ensemble mean
"""

for wind_scenario in ["E-126_7580", "SWT120_3600", "SWT142_3150"]:
    wind_path = zen_path + "Wind/" + wind_scenario + "/"
    for year in range(1905, 2006):
        file_name = glob.glob(base_path + "wind_power/" + wind_scenario + "/Wind*" + str(year) + "*.nc")
        assert len(file_name) == 1  # ensemble members stored in one file
        ds = xr.open_dataset(file_name[0]).mean("number")
        # store as Wind_power_YYYY_ensmean.nc in zenodo/Wind/wind_scenario
        ds.to_netcdf(wind_path + file_name[0].split("/")[-1].split(".nc")[0] + "_ensmean.nc")

"""
Copy country group timeseries
"""

