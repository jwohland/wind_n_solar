import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import sys
import numpy as np

mpl.use("Agg")

sys.path.append("../../")
from utils import plot_field
import cartopy.crs as ccrs

base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
data_path = "output/solar_power/"
plot_path = "plots/solar_power/"

f, ax = plt.subplots(
    ncols=1, figsize=(5, 5), subplot_kw={"projection": ccrs.PlateCarree()}
)
cbar_ax = f.add_axes([0.2, 0.1, 0.6, 0.03])


filelist = sorted(glob.glob(base_path + data_path + "/annual/*.nc"))
all_power = xr.open_mfdataset(filelist, combine="by_coords")
mean_CF = all_power.PV.mean(
    dim=["time"]
)  #  todo number processing will be needed once full data available
mean_CF.compute()

plot_field(
    mean_CF,
    ax,
    "default GSEE panel",
    add_colorbar=True,
    cbar_ax=cbar_ax,
    cbar_kwargs={"orientation": "horizontal", "label": "Mean solar capacity factor"},
    levels=np.arange(0.1, 0.24, 0.02),
)

plt.subplots_adjust(0.05, 0.15, 0.9, 0.98)
plt.savefig(base_path + plot_path + "mean_CF.pdf")
