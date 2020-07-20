import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys

mpl.use("Agg")

sys.path.append("../../")
from utils import plot_field
import cartopy.crs as ccrs

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
data_path = "output/wind_power/"
turbine_names = ["E-126_7580", "SWT120_3600", "SWT142_3150"]
plot_path = "plots/wind_power/"

f, ax = plt.subplots(
    ncols=3, figsize=(12, 5), subplot_kw={"projection": ccrs.PlateCarree()}
)
cbar_ax = f.add_axes([0.3, 0.15, 0.4, 0.03])


for i, turbine_name in enumerate(turbine_names):
    filelist = sorted(glob.glob(base_path + data_path + turbine_name + "/annual/*.nc"))
    all_power = xr.open_mfdataset(filelist, combine="by_coords")
    mean_CF = all_power.wind_power.mean(dim=["time", "number"])
    if i == 1:
        plot_field(
            mean_CF,
            ax[i],
            turbine_name,
            vmin=0.1,
            vmax=0.5,
            add_colorbar= True,
            cbar_ax=cbar_ax,
            cbar_kwargs={"orientation": "horizontal", "label": 'Mean capacity factor'},
            levels=9,
        )
    else:
        plot_field(
            mean_CF,
            ax[i],
            turbine_name,
            vmin=0.1,
            vmax=0.5,
            add_colorbar=False,
            levels=9,
        )
plt.subplots_adjust(0.05, 0.15, 0.95, 0.98)
plt.savefig(base_path + plot_path + "mean_CF.pdf")
