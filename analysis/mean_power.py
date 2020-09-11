# plots mean capacity factors and calculates masks for conflicting_optimization

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys
import xarray as xr
import salem


sys.path.append("../")
from utils import plot_field, Generation_type


def make_mask(mean_CF, scenario, mask_path):
    mask_name = "mask_" + scenario
    try:
        xr.open_dataset(mask_path + mask_name + ".nc")
    except FileNotFoundError:
        thr_dic = {
            "E-126_7580": 0.3,
            "SWT120_3600": 0.25,
            "SWT142_3150": 0.2,
            "both_constant": 0.1,
        }  # determines which turbine can be built where
        threshold = thr_dic[scenario]
        mask = mean_CF.where(
            mean_CF > threshold
        )  # keep only values larger than threshold
        mask /= mask  # set values to one or nan
        mask = mask.rename(mask_name)
        # Ensure that wind turbine masks do not overlap. Rationale: Build biggest suitable turbine
        if "SWT" in scenario:
            print("SWT")
            mask_E126 = xr.open_dataset(mask_path + "mask_E-126_7580.nc")[
                "mask_E-126_7580"
            ]
            mask = mask.where(mask_E126 != 1.0)
        if scenario == "SWT142_3150":
            print("SWT120")
            mask_SWT120 = xr.open_dataset(mask_path + "mask_SWT120_3600.nc")[
                "mask_SWT120_3600"
            ]
            mask = mask.where(mask_SWT120.variable != 1.0)
        if scenario == "both_constant":  # mask offshore for solar
            shape_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/data/shapefile/ne_110m_land/"
            land = salem.read_shapefile(shape_path + "ne_110m_land.shp", cached=True)
            mask = mask.salem.roi(shape=land, all_touched=True)
        mask.to_netcdf(mask_path + mask_name + ".nc")


base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
Solar = Generation_type(
    "solar",
    ["both_constant", "tilt_constant", "neither_constant"],
    "output/solar_power/",
    "plots/analysis/mean_CF/",
    base_path,
)
Wind = Generation_type(
    "wind",
    ["E-126_7580", "SWT120_3600", "SWT142_3150"],
    "output/wind_power/",
    "plots/analysis/mean_CF/",
    base_path,
)

for Generation in [Wind, Solar]:
    f, ax = plt.subplots(
        ncols=3, figsize=(12, 5), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    cbar_ax = f.add_axes([0.3, 0.15, 0.4, 0.03])
    if Generation.name == "solar":
        vmin, vmax = 0.05, 0.25
    else:
        vmin, vmax = 0.1, 0.5
    for i, scenario in enumerate(Generation.scenarios):
        # load data
        mean_CF = (
            Generation.open_data(scenario)
            .load()[Generation.var]
            .mean(dim=["time", "number"])
        )
        # compute and store mask
        if not scenario in ["tilt_constant", "neither_constant"]:
            make_mask(mean_CF, scenario, base_path + "output/optimization/masks/")

        # plotting
        if i == 1:
            plot_field(
                mean_CF,
                ax[i],
                scenario,
                vmin=vmin,
                vmax=vmax,
                add_colorbar=True,
                cbar_ax=cbar_ax,
                cbar_kwargs={
                    "orientation": "horizontal",
                    "label": "Mean " + Generation.name + " capacity factor",
                },
                levels=9,
            )
        else:
            plot_field(
                mean_CF,
                ax[i],
                scenario,
                vmin=vmin,
                vmax=vmax,
                add_colorbar=False,
                levels=9,
            )
    plt.subplots_adjust(0.05, 0.15, 0.95, 0.98)
    plt.savefig(Generation.plot_path + Generation.name + "_mean_CF.pdf")
