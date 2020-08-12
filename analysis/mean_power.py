import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys

sys.path.append("../")
from utils import plot_field, Generation_type

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
Solar = Generation_type(
    "solar",
    ["both_constant", "tilt_constant", "neither_constant"],
    "output/solar_power/",
    "plots/analysis/mean_CF/",
    base_path,
)  # todo update with scenario names once three solar options ready
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
                    "label": "Mean " + Generation.name +" capacity factor",
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
