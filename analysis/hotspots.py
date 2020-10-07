"""
Identification of hotspots of multidecacal variabilty in solar and wind power
"""
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys

sys.path.append("../")
from utils import plot_field, load_solar_ensemble, Generation_type, add_letters
import cartopy.crs as ccrs
import numpy as np

TITLE_DICT = {
    "E-126_7580": "Wind turbine E-126_7580",
    "SWT120_3600": "Wind turbine SWT120_3600",
    "SWT142_3150": "Wind turbine SWT142_3150",
    "both_constant": "Solar constant panel geometry",
    "tilt_constant": "Solar variable azimuth",
    "neither_constant": "Solar variable azimuth and tilt",
}

COLORS = {"wind": ["Olive", "forestgreen", "yellowgreen"],
          "solar": ["darkgoldenrod", "darkorange", "darkkhaki"]}

def plot_CDF(ds, ax, title, scenario, color):
    """
    plots CDF of dataset ds with variable wind_power on axis ax
    :param ds:
    :param ax:
    :param title
    :return:
    """
    values = ds.values.flatten()
    ax.hist(
        values[np.isfinite(values)],  # drop nan values that occur in masked calculation
        cumulative=True,
        bins=3000,
        density=True,
        histtype="step",
        label=TITLE_DICT[scenario],
        linewidth=1.7,
        color=color
    )
    ax.set_title(title)


def relative_change(ds):
    """
    Calculate relative change in per cent between maximum and minimum value
    :return:
    """
    return (ds.max(dim="time") - ds.min(dim="time")) * 100.0 / ds.min(dim="time")


def plot_hotspot(
    Generation,
    delta_power,
    all_power,
    rel,
    scenario,
    ensemblemean=False,
    letter_offset=0,
):
    """

    :param Generation:
    :param delta_power:
    :param all_power:
    :param CF_treshold:
    :return:
    """
    # prepare plotting
    cmap = mpl.cm.get_cmap("YlOrBr")
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
    cbar_ax = f.add_axes([0.15, 0.12, 0.7, 0.03])
    vmax = np.round(delta_power[Generation.var].values.max() * 1.1, 2)
    label_name = "Max - Min 20y " + Generation.name + " power"
    if rel:
        vmax = 20
        label_name = r"Amplitude of multidecadal variability   $\frac{G_\mathrm{max} -G_\mathrm{min}}{G_\mathrm{min}}$ [%]"
    if ensemblemean:
        map_plot = plot_field(
            delta_power.mean("number")[Generation.var],
            ax=ax,
            cmap=cmap,
            add_colorbar=True,
            title=TITLE_DICT[scenario],
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
        add_letters(ax, x=-0.04, y=1.01, fs=12, letter_offset=letter_offset)
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

    plt.subplots_adjust(left=0.05, right=0.96, top=0.95, bottom=0.18)

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


def plot_CDFs(Generations, cf_mins):
    """
    Plots cumulative density functionfor capacity factor bands defined in
    :param Generations: list of Generation_type
    :param cf_mins: dictionary with keys Generation type and entries defining cf ranges
    :return:
    """
    mpl.rcParams["axes.spines.left"] = True
    mpl.rcParams["axes.spines.bottom"] = True
    f, ax = plt.subplots(ncols=3, nrows=2, figsize=(13, 8))
    for j, Generation in enumerate(Generations):
        for i, scenario in enumerate(Generation.scenarios):
            all_power = Generation.open_data(scenario).load()
            # evaluate data
            all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
            cf_inc = cf_mins[Generation.name][1] - cf_mins[Generation.name][0]
            for k, cf_min in enumerate(cf_mins[Generation.name]):
                # locations within a CF band
                masked_power = all_power.where(
                    (all_power[Generation.var].mean(dim=["time", "number"]) > cf_min)
                    & (
                        all_power[Generation.var].mean(dim=["time", "number"])
                        < cf_min + cf_inc
                    )
                )
                title = (
                    str(cf_min)
                    + " < CF < "
                    + str(np.round(cf_min + cf_inc, 2))
                )

                delta_power = relative_change(masked_power[Generation.var])
                plot_CDF(
                    delta_power, ax[j, k], title, scenario, color=COLORS[Generation.name][i]
                )

    for j in range(2):
        ax[j, 0].set_ylabel(["Wind", "Solar"][j] + " CDF", fontsize=12)
        for i in range(3):
            ax[j, 1].set_xlabel(
                r"Amplitude of multidecadal variability   $\frac{G_\mathrm{max} -G_\mathrm{min}}{G_\mathrm{min}}$ [%]", fontsize=12
            )  #
            ax[j, i].set_xlim(xmin=0, xmax=14.5)
            ax[j, i].set_xticks(np.arange(2, 15, 2))
            ax[j, i].grid(True)
        ax[j, 1].legend(bbox_to_anchor=(-1, -.32, 3, .09), loc="center", borderaxespad=0., ncol=3, frameon=False)


    plt.subplots_adjust(0.06, 0.11, 0.97, 0.96, None, 0.435)
    add_letters(ax, fs=12)
    plt.savefig(Generation.plot_path + "/CDF_rel.png")


base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
Solar = Generation_type(
    "solar",
    ["both_constant", "tilt_constant", "neither_constant"],
    "output/solar_power/",
    "plots/analysis/hotspots/",
    base_path,
)
Wind = Generation_type(
    "wind",
    ["E-126_7580", "SWT120_3600", "SWT142_3150"],
    "output/wind_power/",
    "plots/analysis/hotspots/",
    base_path,
)

for rel in [True, False]:
    for i, Generation in enumerate([Wind, Solar]):
        for j, scenario in enumerate(Generation.scenarios):
            # load data
            all_power = Generation.open_data(scenario).load()
            # evaluate data
            all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
            delta_power = all_power.max(dim="time") - all_power.min(dim="time")
            if rel:
                delta_power *= 100.0 / all_power.min(dim="time")
            plot_hotspot(Generation, delta_power, all_power, rel, scenario)
            plot_hotspot(
                Generation,
                delta_power,
                all_power,
                rel,
                scenario,
                ensemblemean=True,
                letter_offset=3 * i + j,
            )

###
# CDF analysis
###
plot_CDFs([Wind, Solar], {"solar": [0.1, 0.15, 0.2], "wind": [0.15, 0.25, 0.35]})
