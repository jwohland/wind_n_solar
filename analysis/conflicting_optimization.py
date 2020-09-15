import salem
from salem.utils import get_demo_file
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("../")
from utils import plot_field, Generation_type

"""
1) select country groups from netcdf data
"""
plot_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/plots/optimization/country_domains"
shape_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/shapefile/EEZ_land_union_v3_202003/"
path = (
    "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/CERA20C/annual/"
)
ds = xr.open_dataset(path + "annual_CERA20C_multiple_1949.nc")

# uncomment if onshore only
# shdf = salem.read_shapefile(get_demo_file("world_borders.shp"))
shdf = salem.read_shapefile(shape_path + "EEZ_Land_v3_202030.shp")

for country_group in [
    ["United Kingdom", "Ireland"],
    ["Portugal", "Spain"],
    ["France", "Belgium", "Netherlands"],
    ["Germany", "Denmark"],
    ["Italy", "Austria", "Switzerland", "Slovenia"],
    ["Sweden", "Norway"],
    ["Poland", "Czech Republic", "Slovakia", "Hungary"],
    ["Lithuania", "Latvia", "Estonia", "Finland"],
    [
        "Greece",
        "Bulgaria",
        "Serbia",
        "Croatia",
        "Bosnia and Herzogovina",
        "Romania",
        "Albania",
    ],
]:
    plt.clf()
    # uncomment if onshore only
    # shdf_tmp = shdf.loc[shdf["CNTRY_NAME"].isin(country_group)]
    shdf_tmp = shdf.loc[shdf["UNION"].isin(country_group)]
    ds_tmp = ds.salem.roi(shape=shdf_tmp, all_touched=True)["s100"].sel({"number": 0})
    N_gridboxes = np.isfinite(ds_tmp.values).sum()
    ds_tmp.salem.quick_map()
    plt.savefig(
        plot_path + "/" + country_group[0] + "_N_" + str(N_gridboxes) + ".jpeg", dpi=300
    )


"""
2) Load technology masks and prepare wind and solar generation timeseries
"""
def combined_generation(Generation, annual=False):
    """
    Calculate combined generation of different wind turbines/solar panels, accounting for location suitability
    as defined in the masks
    :param Generation: wind or solar
    :param annual: if annual=True, annual means for 1901 - 2009 are used. Otherwise, monthly means for 1980- 2000.
    :return:
    """
    for i, scenario in enumerate(Generation.scenarios):
        if annual:
            data = Generation.open_data(scenario).mean(["number"])
        else:
            data = Generation.open_data_monthly_ensmean(scenario)
        mask = Generation.get_mask(scenario)
        mask_name = "mask_" + scenario
        data_masked = data.where(mask[mask_name] == 1, other=0)  # other = 0  needed for the sum
        if i == 0:
            combined = data_masked
        else:
            combined += data_masked
    combined = combined.where(combined[Generation.var] > 0)  # set unsuitable locations to nan again
    return combined


base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
Solar = Generation_type(
    "solar",
    ["both_constant"],
    "output/solar_power/",
    "plots/optimization/",
    base_path,
)
Wind = Generation_type(
    "wind",
    ["E-126_7580", "SWT120_3600", "SWT142_3150"],
    "output/wind_power/",
    "plots/optimization/",
    base_path,
)

# wind

f, ax = plt.subplots(
    ncols=4, figsize=(12, 6), subplot_kw={"projection": ccrs.PlateCarree()}
)
cbar_ax = f.add_axes([0.3, 0.15, 0.4, 0.03])
for i, scenario in enumerate(Wind.scenarios):
    wind_data = Wind.open_data_monthly_ensmean(scenario)
    mask = Wind.get_mask(scenario)
    mask_name = "mask_" + scenario
    wind_data_masked = wind_data.where(mask[mask_name] == 1)
    plot_field(
        wind_data_masked["wind_power"].mean(dim="time"),
        ax[i],
        scenario,
        vmin=0.2,
        vmax=0.55,
        add_colorbar=False,
    )
combined_wind_data = combined_generation(Wind)
plot_field(
    combined_wind_data["wind_power"].mean(dim="time"),
    ax[3],
    "combination",
    vmin=0.2,
    vmax=0.55,
    cbar_ax=cbar_ax,
    cbar_kwargs={
        "orientation": "horizontal",
        "label": "Mean wind capacity factor (1980 - 2000)",
    },
)
plt.savefig(Wind.plot_path + "test_masked_monthly_CF.jpeg", dpi=200)


# solar

combined_solar = combined_generation(Solar)
plt.clf()
f, ax = plt.subplots(
    ncols=1, figsize=(5, 5), subplot_kw={"projection": ccrs.PlateCarree()}
)
plot_field(combined_solar["PV"].mean(dim="time"), ax, "Solar PV")
plt.savefig(Wind.plot_path + "test_masked_solar_monthly_CF.jpeg", dpi=200)