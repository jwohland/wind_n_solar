"""
Identification of hotspots of multidecacal variabilty in solar and wind power
"""
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys

sys.path.append("../")
from utils import plot_field, load_annual_solar_ensemble, Generation_type
import cartopy.crs as ccrs
import numpy as np


def plot_CDF(ds, ax, title, scenario):
    """
    plots CDF of dataset ds with variable wind_power on axis ax
    :param ds:
    :param ax:
    :param title
    :return:
    """
    values = ds.to_array().values.flatten()
    ax.hist(
        values[np.isfinite(values)],  # drop nan values that occur in masked calculation
        cumulative=True,
        bins=3000,
        density=True,
        histtype="step",
        label=scenario,
    )
    ax.set_title(title)


def relative_change(ds):
    """
    Calculate relative change in per cent between maximum and minimum value
    :return:
    """
    return (ds.max(dim="time") - ds.min(dim="time")) * 100.0 / ds.min(dim="time")


def plot_hotspot(
    Generation, delta_power, all_power, rel, scenario, ensemblemean=False,
):
    """

    :param Generation:
    :param delta_power:
    :param all_power:
    :param CF_treshold:
    :return:
    """
    # prepare plotting
    cmap = mpl.cm.get_cmap("coolwarm")
    plt.close("all")
    if ensemblemean:
        f, ax = plt.subplots(
            subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(5, 5),
        )
    else:
        f, ax = plt.subplots(
            ncols=5,
            nrows=2,
            subplot_kw={"projection": ccrs.PlateCarree()},
            figsize=(12, 5),
        )
    cbar_ax = f.add_axes([0.15, 0.1, 0.7, 0.03])
    vmax = np.round(delta_power[Generation.var].values.max() * 1.1, 2)
    label_name = "Max - Min 20y " + Generation.name + " power"
    if rel:
        vmax = 20
        label_name = "(Max - Min)/Min 20y " + Generation.name + "power [%]"
    if ensemblemean:
        plot_field(
            delta_power.mean("number")[Generation.var],
            ax=ax,
            title="ensemble mean",
            cmap=cmap,
            add_colorbar=True,
            cbar_ax=cbar_ax,
            cbar_kwargs={"orientation": "horizontal", "label": label_name},
            vmin=0,
            vmax=vmax,
            levels=7,
        )
        # shade low generation areas
        mean_CF = all_power.mean(["time", "number"])[Generation.var]
        mean_CF.plot.contourf(
            ax=ax,
            levels=[0, Generation.CF_threshold],
            hatches=["...", ""],
            alpha=0,
            add_colorbar=False,
            add_labels=False,
        )
    else:
        for i, number in enumerate(all_power.number.values):
            plot_field(
                delta_power.sel({"number": number}, drop=True)[Generation.var],
                ax=ax.flatten()[i],
                title=str(number),
                cmap=cmap,
                add_colorbar=True,
                cbar_ax=cbar_ax,
                cbar_kwargs={"orientation": "horizontal", "label": label_name},
                vmin=0,
                vmax=vmax,
                levels=9,
            )
            # shade low generation areas
            mean_CF = all_power.sel({"number": number}).mean("time")[Generation.var]
            mean_CF.plot.contourf(
                ax=ax.flatten()[i],
                levels=[0, Generation.CF_threshold],
                hatches=["...", ""],
                alpha=0,
                add_colorbar=False,
                add_labels=False,
            )

    plt.subplots_adjust(left=0.05, right=0.92, bottom=0.18)
    plt.suptitle(scenario)
    plotname = "_hotspots"
    if rel:
        plotname += "_rel"
    if ensemblemean:
        plotname += "_ensmean"
    plt.savefig(
        Generation.plot_path
        + "/"
        + Generation.name
        + "/"
        + scenario
        + plotname
        + ".png"
    )


def plot_CDFs(Generation, cf_mins=[0.3, 0.2]):
    cf_mins.append(Generation.CF_threshold)
    mpl.rcParams["axes.spines.left"] = True
    mpl.rcParams["axes.spines.bottom"] = True
    f, ax = plt.subplots(ncols=3, figsize=(15, 5))
    for i, scenario in enumerate(Generation.scenarios):
        all_power = Generation.open_data(scenario).load()
        # evaluate data
        all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
        for j, cf_min in enumerate(cf_mins):
            if j == 2:  # all locations with CF larger than threshold
                masked_power = all_power.where(
                    (all_power[Generation.var].mean(dim=["time", "number"]) > cf_min)
                )
                title = "CF > " + str(cf_min)
            else:  # locations within a CF band
                masked_power = all_power.where(
                    (all_power[Generation.var].mean(dim=["time", "number"]) > cf_min)
                    & (
                        all_power[Generation.var].mean(dim=["time", "number"])
                        < cf_min + 0.1
                    )
                )
                title = str(np.round(cf_min + 0.1, 2)) + " > CF > " + str(cf_min)
            delta_power = relative_change(masked_power)
            plot_CDF(
                delta_power,
                ax[j],
                title,
                scenario,
            )

    ax[0].set_ylabel("Cumulative density function")
    for i in range(3):
        ax[i].set_xlabel(r"$\frac{P_{max} - P_{min}}{P_{min}}$ [%]")  #
        ax[i].set_xlim(xmin=2, xmax=15)
        ax[i].grid(True)
    ax[0].legend(loc=2)
    plt.tight_layout()
    plt.savefig(Generation.plot_path + "/" + Generation.name + "/CDF_rel.png")


base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
Solar = Generation_type(
    "solar",
    ["default_panel"],
    "output/solar_power/",
    "plots/analysis/hotspots/",
    base_path,
)  # todo update with scenario names once three solar options ready
Wind = Generation_type(
    "wind",
    ["E-126_7580", "SWT120_3600", "SWT142_3150"],
    "output/wind_power/",
    "plots/analysis/hotspots/",
    base_path,
)

for Generation in [Solar, Wind]:
    for rel in [True, False]:
        for scenario in Generation.scenarios:
            # load data
            all_power = Generation.open_data(scenario).load()
            # evaluate data
            all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
            delta_power = all_power.max(dim="time") - all_power.min(dim="time")
            if rel:
                delta_power *= 100.0 / all_power.min(dim="time")
            plot_hotspot(Generation, delta_power, all_power, rel, scenario)
            plot_hotspot(
                Generation, delta_power, all_power, rel, scenario, ensemblemean=True
            )  # todo plot ensemble mean in 3x1 subplot for all scenarios

###
# CDF analysis
###
plot_CDFs(Solar, [0.1, 0.0])
plot_CDFs(Wind, [0.3, 0.2])
