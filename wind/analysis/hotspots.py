"""
Identification of hotspots of multidecacal variabilty
"""
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys

sys.path.append("../../")
from utils import plot_field
import cartopy.crs as ccrs
import numpy as np

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
data_path = "output/wind_power/"
turbine_names = ["E-126_7580", "SWT120_3600", "SWT142_3150"]
plot_path = "plots/analysis/hotspots/"

for turbine_name in turbine_names:
    for rel in [True, False]:
        # load data
        filelist = sorted(
            glob.glob(base_path + data_path + turbine_name + "/annual/*.nc")
        )
        all_power = xr.open_mfdataset(filelist, combine="by_coords")
        all_power.load()
        # prepare plotting
        cmap = mpl.cm.get_cmap("Greens")
        plt.close("all")
        f, ax = plt.subplots(
            ncols=5,
            nrows=2,
            subplot_kw={"projection": ccrs.PlateCarree()},
            figsize=(12, 5),
        )
        cbar_ax = f.add_axes([0.3, 0.1, 0.4, 0.03])
        # evaluate data
        all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
        delta_power = all_power.max(dim="time") - all_power.min(dim="time")
        label_name = "Max - Min 20y wind power"
        if rel:
            delta_power *= 100.0 / all_power.min(dim="time")  # conversion to percent
            # vmax = 20
            label_name = "(Max - Min)/Min 20y wind power [%]"
        vmax = np.round(delta_power["wind_power"].values.max() * 1.1, 2)
        for i, number in enumerate(all_power.number.values):
            if i == 0:
                plot_field(
                    delta_power.sel({"number": number}, drop=True)["wind_power"],
                    ax=ax.flatten()[i],
                    title=str(number),
                    cmap=cmap,
                    add_colorbar=True,
                    cbar_ax=cbar_ax,
                    cbar_kwargs={"orientation": "horizontal", "label": label_name},
                    vmin=0,
                    vmax=vmax,
                )
            else:
                plot_field(
                    delta_power.sel({"number": number}, drop=True)["wind_power"],
                    ax=ax.flatten()[i],
                    title=str(number),
                    cmap=cmap,
                    add_colorbar=False,
                    vmin=0,
                    vmax=vmax,
                )
            # shade low generation areas
            mean_CF = all_power.sel({"number": number}).mean("time").wind_power
            mean_CF.plot.contourf(
                ax=ax.flatten()[i],
                levels=[0, 0.15],
                hatches=["....", ""],
                alpha=0,
                add_colorbar=False,
                add_labels=False,
            )

        plt.subplots_adjust(left=0.05, right=0.92, bottom=0.18)
        plt.suptitle(turbine_name)
        plotname = "_hotspots"
        if rel:
            plotname += "_rel"
        plt.savefig(base_path + plot_path + turbine_name + plotname + ".png")


###
# CDF analysis
###
def plot_CDF(ds, ax, title):
    """
    plots CDF of dataset ds with variable wind_power on axis ax
    :param ds:
    :param ax:
    :param title
    :return:
    """
    values = ds.wind_power.values.flatten()
    ax.hist(
        values[np.isfinite(values)],  # drop nan values that occur in masked calculation
        cumulative=True,
        bins=3000,
        density=True,
        histtype="step",
        label=turbine_name,
    )
    ax.set_title(title)


def relative_change(ds):
    """
    Calculate relative change in per cent between maximum and minimum value
    :return:
    """
    return (ds.max(dim="time") - ds.min(dim="time")) * 100.0 / ds.min(dim="time")


mpl.rcParams["axes.spines.left"] = True
mpl.rcParams["axes.spines.bottom"] = True
f, ax = plt.subplots(ncols=3, figsize=(15, 5))
for i, turbine_name in enumerate(turbine_names):
    filelist = sorted(glob.glob(base_path + data_path + turbine_name + "/annual/*.nc"))
    all_power = xr.open_mfdataset(filelist, combine="by_coords")
    all_power.load()
    all_power = all_power.rolling(time=20, center=True).mean().dropna("time")
    # all sites
    delta_power = relative_change(all_power)
    plot_CDF(delta_power, ax[0], "all locations")
    # mask low-generation sites
    for j, cf_min in enumerate([0.3, 0.2]):
        masked_wind_power = all_power.where(
            (all_power.wind_power.mean(dim=["time", "number"]) > cf_min)
            & (all_power.wind_power.mean(dim=["time", "number"]) < cf_min + 0.1)
        )
        delta_power = relative_change(masked_wind_power)
        plot_CDF(
            delta_power,
            ax[1 + j],
            str(np.round(cf_min + 0.1, 2)) + " > CF > " + str(cf_min),
        )
    # plot 80th percentile
    # for CF_thr in np.arange(0, 0.5, 0.01):
    #    masked_wind_power = all_power.where(
    #        all_power.wind_power.mean(dim=["time", "number"]) > CF_thr
    #    )
    #    delta_power = relative_change(masked_wind_power).wind_power
    #    np.percentile(delta_power[np.isfinite(masked_wind_power)], 80)

ax[0].set_ylabel("Cumulative density function")
for i in range(3):
    ax[i].set_xlabel(r"$\frac{P_{max} - P_{min}}{P_{min}}$ [%]")  #
    ax[i].set_xlim(xmin=2, xmax=15)
    ax[i].grid(True)
ax[0].legend(loc=2)
plt.tight_layout()
plt.savefig(base_path + plot_path + "CDF_rel.png")

# convergence of 80th percentile with higher CF_thresholds?
