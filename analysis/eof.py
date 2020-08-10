from eofs.xarray import Eof
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import cartopy.crs as ccrs
import numpy as np
import pandas as pd

sys.path.append("../")
from utils import plot_field, Generation_type, add_letters

mpl.rcParams["axes.spines.left"] = True
mpl.rcParams["axes.spines.bottom"] = True


def running(all_power, time_scale_name):
    """
    Calculate the running mean of all_power
    :param all_power:
    :param time_scale_name:
    :return:
    """
    if time_scale_name == "ten":
        all_power = all_power.rolling(time=10, center=True).mean().dropna("time")
    elif time_scale_name == "twenty":
        all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
    return all_power


def perform_EOF(power, N):
    """
    Perform EOF analysis of input data Array
    :param power:
    :param N:
    :return:
    """
    solver = Eof(power)
    eofs = solver.eofs(neofs=N)
    pcs = solver.pcs(npcs=N)
    variance_fraction = solver.varianceFraction()
    return eofs, pcs, variance_fraction


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
    for time_scale_name in ["annual", "ten", "twenty"]:
        for scenario in Generation.scenarios:
            Generation.reset_plot_path(
                "plots/analysis/eof/"
                + time_scale_name
                + "/"
                + Generation.name
                + "_power/"
                + scenario
                + "/"
            )
            # load data
            all_power = Generation.open_data(scenario).load()
            all_power = running(all_power, time_scale_name)
            # evaluate data
            ens_members = list(all_power.number.values)
            ens_members.append("ensemble_mean")
            for number in ens_members:
                # choose ensemble member
                if number == "ensemble_mean":
                    power = all_power.mean("number")
                else:
                    power = all_power.sel({"number": number})
                # EOF analysis
                N_EOFs = 4
                eofs, pcs, variance_fraction = perform_EOF(
                    power[Generation.var], N_EOFs
                )
                # prepare plots
                cmap = mpl.cm.get_cmap("RdBu_r")
                plt.close("all")
                f, ax = plt.subplots(
                    ncols=N_EOFs,
                    nrows=2,
                    subplot_kw={"projection": ccrs.PlateCarree()},
                    gridspec_kw={"height_ratios": [3, 1]},
                    figsize=(12, 6),
                )
                plt.subplots_adjust(
                    left=0.05, right=0.97, bottom=0.15, hspace=0.5, top=0.92
                )
                cbar_ax = f.add_axes([0.2, 0.08, 0.6, 0.02])
                # make plots
                for i in range(N_EOFs):
                    vmax = np.round(
                        np.max([eofs.values.max(), -eofs.values.min()]) * 1.1, 2
                    )
                    if i != 0:
                        plot_field(
                            eofs.sel({"mode": i}, drop=True),
                            ax=ax[0, i],
                            title=str(
                                np.round(
                                    variance_fraction.sel({"mode": i}).values * 100, 1
                                )
                            )
                            + "% variance ",
                            vmin=-vmax,
                            vmax=vmax,
                            add_colorbar=False,
                            cmap=cmap,
                        )
                    else:
                        plot_field(
                            eofs.sel({"mode": i}, drop=True),
                            ax=ax[0, i],
                            title=str(
                                np.round(
                                    variance_fraction.sel({"mode": i}).values * 100, 1
                                )
                            )
                            + "% variance ",
                            vmin=-vmax,
                            vmax=vmax,
                            add_colorbar=True,
                            cbar_ax=cbar_ax,
                            cmap=cmap,
                            cbar_kwargs={
                                "orientation": "horizontal",
                                "label": Generation.name + " power generation [a.u.]",
                            },
                        )
                    ax[1, i] = plt.subplot(
                        2, N_EOFs, N_EOFs + 1 + i
                    )  # override the GeoAxes object
                    pc = pcs.sel({"mode": i}, drop=True)
                    pc = pd.Series(data=pc.values, index=pc.time.values)
                    pc.plot(ax=ax[1, i])
                    pc.rolling(window=5, center=True).mean().plot(
                        ax=ax[1, i], ls="--", color="black", lw=2
                    )
                plt.suptitle(scenario)
                if number == "ensemble_mean":
                    plotname = "/" + Generation.name + "power_eofs_mean.png"
                else:
                    plotname = (
                        "/"
                        + Generation.name
                        + "power_eofs_"
                        + str(int(number))
                        + ".png"
                    )
                add_letters(ax)
                plt.savefig(Generation.plot_path + plotname)

            # add mean timeseries
            f, ax = plt.subplots()
            try:
                all_power[Generation.var].mean(dim=["latitude", "longitude", "number"]).plot(ax=ax)
            except ValueError:
                all_power[Generation.var].mean(dim=["lat", "lon", "number"]).plot(ax=ax)
            ax.set_title(
                scenario + " " + Generation.name + " generation (mean over Europe, " + time_scale_name + " y, ensemble)"
            )
            plt.tight_layout()
            plt.savefig(Generation.plot_path + "/" + scenario + "mean_timeseries.png")
