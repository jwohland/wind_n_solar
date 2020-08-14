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
from utils import plot_field, Generation_type, calc_correlations_xarray

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
    f, ax = plt.subplots(
        ncols=3,
        nrows=3,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(10, 10),
    )
    cbar_ax = f.add_axes([0.15, 0.05, 0.7, 0.02])
    # load 20y running mean data
    solar_power = Solar.open_data(solar_scenario).load().mean(["number"])
    solar_power["time"] += pd.Timedelta(
        90, unit="m"
    )  # Solar power is 1:30h shifted back
    solar_power = solar_power.rolling(time=20, center=True).mean().dropna("time")
    solar_power_monthly = Solar.open_data_monthly_ensmean(solar_scenario)
    for i, wind_scenario in enumerate(Wind.scenarios):
        print(str(i))
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
        # calculate difference
        correlations_diff = correlations_monthly - correlations_longterm
        # plotting
        plot_field(
            correlations_longterm,
            ax[i, 0],
            cmap=cmap,
            title="longterm" if i == 0 else "",
            add_colorbar=False,
            extend='both',
            vmin=-1,
            vmax=1,
        )
        ax[i, 0].set_ylabel(wind_scenario)
        plot_field(
            correlations_monthly,
            ax[i, 1],
            title="seasonal to interannual" if i == 0 else "",
            add_colorbar=False,
            extend='both',
            cmap=cmap,
            vmin=-1,
            vmax=1,
        )
        if i == 0:
            plot_field(
                correlations_diff,
                ax[i, 2],
                title="difference" if i == 0 else "",
                cbar_ax=cbar_ax,
                cmap=cmap,
                cbar_kwargs={
                    "orientation": "horizontal",
                    "label": "Correlation between wind and solar generation",
                },
                extend='both',
                vmin=-1,
                vmax=1,
            )
        else:
            plot_field(
                correlations_diff,
                ax[i, 2],
                title="difference" if i == 0 else "",
                add_colorbar=False,
                extend='both',
                vmin=-1,
                vmax=1,
                cmap=cmap,
            )
        add_row_label(ax[i, 0], wind_scenario)

    plt.subplots_adjust(left=.08, right=.93, top=.95)
    plt.savefig(Wind.plot_path + "corr_" + solar_scenario + ".png")
