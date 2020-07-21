from eofs.xarray import Eof
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys

sys.path.append("../../")
from utils import plot_field, load_annual_solar_ensemble
import cartopy.crs as ccrs
import numpy as np
import pandas as pd

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
data_path = "output/solar_power/"
panel_names = ["default_panel"]

for time_scale in ["annual", "ten", "twenty"]:
    plot_path = "plots/analysis/eof/" + time_scale + "/solar_power/"

    for panel_name in panel_names:
        filelist = sorted(
            glob.glob(base_path + data_path + panel_name + "/annual/*.nc")
        )
        all_power = load_annual_solar_ensemble(base_path, data_path + panel_name + "/")
        all_power.load()
        if time_scale == "ten":
            all_power = all_power.rolling(time=10, center=True).mean().dropna("time")
        elif time_scale == "twenty":
            all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
        for number in all_power.number.values:
            power = all_power.sel({"number": number})
            solver = Eof(power["PV"])
            N = 5
            eofs = solver.eofs(neofs=N)
            pcs = solver.pcs(npcs=N)
            variance_fraction = solver.varianceFraction()

            cmap = mpl.cm.get_cmap("RdBu_r")
            plt.close("all")
            f, ax = plt.subplots(
                ncols=N,
                nrows=2,
                subplot_kw={"projection": ccrs.PlateCarree()},
                figsize=(12, 5),
            )
            cbar_ax = f.add_axes([0.3, 0.1, 0.4, 0.02])
            for i in range(5):
                vmax = np.round(
                    np.max([eofs.values.max(), -eofs.values.min()]) * 1.1, 2
                )
                if i != 0:
                    plot_field(
                        eofs.sel({"mode": i}, drop=True),
                        ax=ax[0, i],
                        title=str(
                            np.round(variance_fraction.sel({"mode": i}).values * 100, 1)
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
                            np.round(variance_fraction.sel({"mode": i}).values * 100, 1)
                        )
                        + "% variance ",
                        vmin=-vmax,
                        vmax=vmax,
                        add_colorbar=True,
                        cbar_ax=cbar_ax,
                        cmap=cmap,
                        cbar_kwargs={"orientation": "horizontal"},
                    )
                ax[1, i] = plt.subplot(2, N, N + 1 + i)  # override the GeoAxes object
                pc = pcs.sel({"mode": i}, drop=True)
                pc = pd.Series(data=pc.values, index=pc.time.values)
                pc.plot(ax=ax[1, i])
                pc.rolling(window=5, center=True).mean().plot(
                    ax=ax[1, i], ls="--", color="black", lw=2
                )
            plt.subplots_adjust(left=0.05, right=0.92, bottom=0.25)
            plt.suptitle(panel_name)
            plt.savefig(
                base_path
                + plot_path
                + panel_name
                + "/solarpower_eofs_"
                + str(int(number))
                + ".png"
            )
            # add mean timeseries
            mpl.rcParams["axes.spines.left"] = True
            mpl.rcParams["axes.spines.bottom"] = True
            f, ax = plt.subplots()
            all_power.mean(dim=["lat", "lon", "number"]).PV.plot(ax=ax)
            ax.set_title(
                "PV Generation (mean over Europe, " + time_scale + " y, ensemble)"
            )
            plt.tight_layout()
            plt.savefig(base_path + plot_path + panel_name + "/mean_timeseries.png")
