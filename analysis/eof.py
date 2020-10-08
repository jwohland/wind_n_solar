from eofs.xarray import Eof
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import xarray as xr

sys.path.append("../")
from utils import plot_field, Generation_type, add_letters, add_row_label

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


def perform_EOF_multi(data, N, coords, datanames):
    """
    Perform multivariate EOF analysis of input data Arrays
    :param data: list of data arrays (first dimension must be time)
    :param N: number of EOFs to keep
    :param coords: Grid for eof output (multivariate not defined as xarray operation)
    :param datanames: Names of input data
    :return:
    """
    from eofs.multivariate.standard import MultivariateEof

    solver = MultivariateEof(data)
    eofs = np.asarray(solver.eofs(neofs=N))
    pcs = solver.pcs(npcs=N)
    variance_fraction = solver.varianceFraction()

    # transform eofs into a dataset
    tmp = []
    for i, dataname in enumerate(datanames):
        tmp2 = []
        for j in range(N):
            tmp2.append(
                xr.DataArray(
                    data=eofs[i][j], coords=coords, dims=["lat", "lon"], name=dataname
                )
            )
        tmp.append(xr.concat(tmp2, pd.Index(np.arange(N), name="mode")))
    eofs = xr.Dataset({datanames[0]: tmp[0], datanames[1]: tmp[1]})

    return eofs, pcs, variance_fraction


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
individual = False
joint = True

if individual:
    for Generation in [Solar, Wind]:  # zip?
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
                        plotname = Generation.name + "power_eofs_mean.png"
                    else:
                        plotname = (
                            Generation.name + "power_eofs_" + str(int(number)) + ".png"
                        )
                    add_letters(ax)
                    plt.savefig(Generation.plot_path + plotname)

                # add mean timeseries
                f, ax = plt.subplots()
                try:
                    all_power[Generation.var].mean(
                        dim=["latitude", "longitude", "number"]
                    ).plot(ax=ax)
                except ValueError:
                    all_power[Generation.var].mean(dim=["lat", "lon", "number"]).plot(ax=ax)
                ax.set_title(
                    scenario
                    + " "
                    + Generation.name
                    + " generation (mean over Europe, "
                    + time_scale_name
                    + " y, ensemble)"
                )
                plt.tight_layout()
                plt.savefig(Generation.plot_path + "/" + scenario + "mean_timeseries.png")


"""
multivariate EOF

multi-variate EOF has one set of PCs and different EOFs for each input field
"""
if joint:
    # load solar (almost no variability in solar, so stick to one member)
    solar_scen = "both_constant"
    solar_power = Solar.open_data(solar_scen).mean("number").load()
    solar_power = running(solar_power, "twenty")
    solar_power["time"] += pd.Timedelta(90, unit="m")  # Solar power is 1:30h shifted back
    # normalize with long-term mean
    solar_power /= solar_power.mean("time")
    for wind_scen in Wind.scenarios:
        # load wind
        wind_power = Wind.open_data(wind_scen).mean("number").load()
        wind_power = running(wind_power, "twenty")
        wind_power /= wind_power.mean("time")
        # MEOF analysis
        N_EOFs = 4
        data = [
            data.values for data in xr.align(solar_power[Solar.var], wind_power[Wind.var])
        ]
        eofs, pcs, variance_fraction = perform_EOF_multi(
            data, N_EOFs, solar_power.drop("time").coords, [Solar.var, Wind.var]
        )
        # prepare plots
        cmap = mpl.cm.get_cmap("RdBu_r")
        plt.close("all")
        f = plt.figure(figsize=(12, 8))
        gs = plt.GridSpec(nrows=5, ncols=40, figure=f, height_ratios=[3, 3, 0.3, .7, 2.5], hspace=.1)
        axes, ax_list = {}, []
        for i_axis, name_axis in enumerate(["wind", "solar"]):
            axes[name_axis] = {}
            for j_axis in range(4):
                axes[name_axis][str(j_axis)] = plt.subplot(gs[i_axis, (j_axis)*10:(j_axis+1)*10], projection= ccrs.PlateCarree(), frameon=False)
                ax_list.append(axes[name_axis][str(j_axis)])  # needed to add letter later

        axes["colorbar"] = plt.subplot(gs[2, 9:31])
        axes["PC"] = {}
        for j_axis in range(4):
            if j_axis == 0:
                axes["PC"][str(j_axis)] = plt.subplot(gs[4, 1:9])
            else:
                axes["PC"][str(j_axis)] = plt.subplot(gs[4, (j_axis)*10+1:(j_axis+1)*10-1])
            ax_list.append(axes["PC"][str(j_axis)])  # needed to add letter later

        # make plots
        vmax = .08

        for i in range(N_EOFs):
            if i != 0:
                plot_field(
                    eofs[Solar.var].sel({"mode": i}, drop=True),
                    ax=axes["solar"][str(i)],
                    vmin=-vmax,
                    vmax=vmax,
                    extend='both',
                    add_colorbar=False,
                    cmap=cmap,
                )
            else:
                plot_field(
                    eofs[Solar.var].sel({"mode": i}, drop=True),
                    ax=axes["solar"][str(i)],
                    vmin=-vmax,
                    vmax=vmax,
                    add_colorbar=True,
                    cbar_ax=axes["colorbar"],
                    cmap=cmap,
                    extend='both',
                    cbar_kwargs={
                        "orientation": "horizontal",
                        "label": r"Multivariate EOF of relative wind and solar power generation $\frac{G_{20y}}{G_{\mathrm{mean}}}$ [a.u.]"
                    },
                )
            plot_field(
                eofs[Wind.var].sel({"mode": i}, drop=True),
                ax=axes["wind"][str(i)],
                title=str(np.round(variance_fraction[i] * 100, 1)) + "% joint variance ",
                vmin=-vmax,
                vmax=vmax,
                add_colorbar=False,
                cmap=cmap,
                extend='both',
            )
            pc = pcs[:, i]
            pc = pd.Series(data=pc, index=wind_power.time.values)
            pc.plot(ax=axes["PC"][str(i)], color="grey", lw=1.2)
            axes["PC"][str(i)].set_ylim(ymin=-1.5, ymax=1.5)
            axes["PC"][str(i)].set_yticks([-1, 0, 1])
            if i != 0:
                axes["PC"][str(i)].set_yticklabels([])
        add_row_label(axes["solar"][str(0)], "Solar", 12, x=-0.18)
        add_row_label(axes["wind"][str(0)], "Wind", 12, x=-.18)
        axes["PC"]["0"].set_ylabel("PC timeseries [a.u.]", fontsize=12)
        plotname = wind_scen + "_power_meofs_mean.png"
        add_letters(np.array(ax_list[:-4]), fs=12, x=-0.08, y=1.003)
        add_letters(np.array(ax_list[-4:]), fs=12, letter_offset=8, x=-0.11)

        plt.subplots_adjust(0.03, 0.05, 0.99, 0.97)
        plt.savefig("/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/plots/analysis/eof/twenty/"+ plotname)
