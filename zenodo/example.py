# Example script that loads wind and solar generation for all scenarios and plots the long-term mean

import xarray as xr
import glob
import matplotlib.pyplot as plt
from cartopy import feature
import cartopy.crs as ccrs

zen_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/zenodo/"

f, axs = plt.subplots(
    nrows=2, ncols=3, subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(14, 10)
)


# load solar both constant scenario and compute long-term mean
solar_scenarios = [
    "Solar/both_constant/",
    "Solar/tilt_constant/",
    "Solar/neither_constant/",
]
for i, solar_scenario in enumerate(solar_scenarios):
    filelist = glob.glob(zen_path + solar_scenario + "*.nc")
    ds = xr.open_mfdataset(filelist, combine="by_coords")
    ds["PV"].mean("time").plot(
        ax=axs.flatten()[i],
        vmin=0.05,
        vmax=0.25,
        cbar_kwargs={
            "orientation": "horizontal",
            "label": "Mean capacity factor",
            "fraction": 0.056,
        },
        levels=9,
    )
    axs.flatten()[i].set_title(solar_scenario.split("/")[-2])

# same for wind
wind_scenarios = [
    "Wind/E-126_7580/",
    "Wind/SWT120_3600/",
    "Wind/SWT142_3150/",
]
for i, wind_scenario in enumerate(wind_scenarios):
    filelist = glob.glob(zen_path + wind_scenario + "*.nc")
    ds = xr.open_mfdataset(filelist, combine="by_coords")
    ds["wind_power"].mean("time").plot(
        ax=axs.flatten()[i + 3],
        vmin=0.1,
        vmax=0.5,
        cbar_kwargs={
            "orientation": "horizontal",
            "label": "Mean capacity factor",
            "fraction": 0.056,
        },
        levels=9,
    )
    axs.flatten()[i + 3].set_title(wind_scenario.split("/")[-2])

# add country boarders and coastlines
for ax in axs.flatten():
    ax.add_feature(feature.COASTLINE.with_scale("50m"))
    ax.add_feature(feature.BORDERS.with_scale("50m"))

plt.subplots_adjust(left=0.05, right=.95, bottom=0.07)
plt.savefig(zen_path + "example_plot.jpeg", dpi=300)
