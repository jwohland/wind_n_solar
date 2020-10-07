"""
Plotting of solar-wind correlation maps
- using monthly means 1980-2000
- using 20y running means 1905 - 2005
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import sys
import cartopy.crs as ccrs

sys.path.append("../")
from utils import plot_field, Generation_type, calc_correlations_xarray, add_letters

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
Solar = Generation_type(
    "solar",
    ["both_constant", "tilt_constant", "neither_constant"],
    "output/solar_power/",
    "plots/analysis/correlation/",
    base_path,
)
Wind = Generation_type(
    "wind",
    ["E-126_7580", "SWT120_3600", "SWT142_3150"],
    "output/wind_power/",
    "plots/analysis/correlation/",
    base_path,
)
cmap = cm.get_cmap("coolwarm")
for solar_scenario in Solar.scenarios:
    # load 20y running mean data
    solar_power = Solar.open_data(solar_scenario).load().mean(["number"])
    solar_power["time"] += pd.Timedelta(
        90, unit="m"
    )  # Solar power is 1:30h shifted back
    solar_power = solar_power.rolling(time=20, center=True).mean().dropna("time")
    solar_power_monthly = Solar.open_data_monthly_ensmean(solar_scenario)
    for i, wind_scenario in enumerate(Wind.scenarios):
        print(str(i))
        f, ax = plt.subplots(
            ncols=2,
            subplot_kw={"projection": ccrs.PlateCarree()},
            figsize=(8, 4),
        )
        cbar_ax = f.add_axes([0.15, 0.12, 0.7, 0.03])
        wind_power = Wind.open_data(wind_scenario).load().mean(["number"])
        wind_power = wind_power.rolling(time=20, center=True).mean().dropna("time")
        wind_power_monthly = Wind.open_data_monthly_ensmean(wind_scenario)
        # calculate correlations of running mean smoothed data
        correlations_longterm = calc_correlations_xarray(
            solar_power, wind_power, Solar.var, Wind.var, "pearson"
        )
        # calculate correlations of monthly mean data
        correlations_monthly = calc_correlations_xarray(
            solar_power_monthly, wind_power_monthly, Solar.var, Wind.var, "pearson"
        )
        # plotting
        plot_field(
            correlations_longterm,
            ax[0],
            cmap=cmap,
            title="Multidecadal",
            add_colorbar=False,
            vmin=-1,
            vmax=1,
        )
        #ax[i].set_ylabel(wind_scenario)
        plot_field(
            correlations_monthly,
            ax[1],
            title="Seasonal",
            cbar_ax=cbar_ax,
            cbar_kwargs={
                "orientation": "horizontal",
                "label": "Correlation between wind and solar generation",
            },
            cmap=cmap,
            vmin=-1,
            vmax=1,
        )
        #add_row_label(ax[i, 0], wind_scenario)
        add_letters(ax, fs=12)
        plt.subplots_adjust(left=0.05, right=0.97, top=0.94, bottom=0.17)
        plt.savefig(Wind.plot_path + "corr_" + solar_scenario + "_" + wind_scenario +".png")




