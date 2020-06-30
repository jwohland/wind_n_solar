import xarray as xr
import matplotlib.pyplot as plt
import glob
import matplotlib as mpl
import cartopy.crs as ccrs
import sys
sys.path.insert(1, '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/code')
from utils import plot_field


mpl.use('Agg')


def identify_representative(ds):
    """
    Identification of most representative ensemble member
    Based on a comparison of spatio-temporal averages, one of the ensemble members next to the median are chosen
    :param ds: input dataset
    :return:
    """
    tmp = ds['s100_low'].mean(dim=['time', 'longitude', 'latitude']).load()
    pseudo_median = sorted(tmp.values)[int(tmp.size / 2)]  # not strictly the median because a median is not contained
    # in the input data if the sample size is even
    return tmp.where(tmp == pseudo_median, drop=True).number.values[0]


data_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/output/'
plot_path = '/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/plots/ensemble_spread/'
colors = ['peru', 'olive']

# Identification of representative ensemble members
noaa_data = xr.open_mfdataset(glob.glob(data_path + '20CR*_latindex_[15 25 35]*.nc'), combine='by_coords')
cera_data = xr.open_mfdataset(glob.glob(data_path + 'CERA*_latindex_[15 25 35]*.nc'), combine='by_coords')
rep_noaa, rep_cera = identify_representative(noaa_data), identify_representative(cera_data)

# Investigate spread in full ensemble for three latitude bands
for lat in [15, 25, 35]:
    noaa_data = xr.open_mfdataset(glob.glob(data_path + '20CR*_latindex_' + str(lat) + '*.nc'), combine='by_coords')
    cera_data = xr.open_mfdataset(glob.glob(data_path + 'CERA*_latindex_' + str(lat) + '*.nc'), combine='by_coords')
    for lon in noaa_data.longitude.values:
        # plot absolute timeseries
        print(lon)
        f, ax = plt.subplots()
        for i, data in enumerate([noaa_data, cera_data]):
            name = '20CRv3 ' if i == 0 else 'CERA20C '
            tmp_data = data.sel({'longitude': lon, 'latitude': data.latitude.values[0]})
            for number in data.number.values:
                # grey color for representative member
                color = 'black' if ((i == 0 and number == rep_noaa) or (i == 1 and number == rep_cera)) else colors[i]
                tmp_data['s100_low'].sel({'number': number}).plot(ax=ax,
                                                                  label=name if number == 1 else '',
                                                                  color=color,
                                                                  alpha=.8)
        ax.set_title('Lat =' + str(data.latitude.values[0]) + ', Lon = ' + str(lon))
        ax.legend(loc=1)
        plt.tight_layout()
        ax.set_ylabel('Lanczos-filtered 100m wind speed [m/s]')
        plt.savefig(plot_path + str(lat) + '/Abs_lat_' + str(data.latitude.values[0]) + '_lon_' + str(lon) + '.jpeg',
                    dpi=300)
        plt.close('all')

    for lon in noaa_data.longitude.values:
        # plot differences
        f, ax = plt.subplots()
        tmp_noaa, tmp_cera = noaa_data.sel({'longitude': lon, 'latitude': noaa_data.latitude.values[0]}), \
                             cera_data.sel({'longitude': lon, 'latitude': cera_data.latitude.values[0]})
        for noaa_num in tmp_noaa.number.values:
            for cera_num in tmp_cera.number.values:
                diff = tmp_noaa['s100_low'].sel({'number': noaa_num}) - tmp_cera['s100_low'].sel({'number': cera_num})
                if ((noaa_num == rep_noaa) and (cera_num == rep_cera)):
                    diff.plot(ax=ax, color='black', alpha=.8, lw=3)
                else:
                    diff.plot(ax=ax, color='indigo', alpha=.2, lw=1)
        ax.set_title('Lat =' + str(noaa_data.latitude.values[0]) + ', Lon = ' + str(lon))
        ax.set_ylabel('Difference in filtered 100m wind speed [m/s]')
        ax.set_ylim(ymin=-1, ymax=3)
        ax.axhline(y=0, ls='--', color='grey')
        plt.tight_layout()
        plt.savefig(
            plot_path + str(lat) + '/Diff_lat_' + str(noaa_data.latitude.values[0]) + '_lon_' + str(lon) + '.jpeg')
        plt.close('all')

# Maps using representative members (CERA 7, NOAA 5)
cera_data = xr.open_mfdataset(glob.glob(data_path + 'CERA_7_latindex_*.nc'), combine='by_coords')
noaa_data = xr.open_mfdataset(glob.glob(data_path + '20CRv3_5_latindex_*.nc'), combine='by_coords')

def plot_diff_mean(ds_noaa, ds_cera, year_start, year_duration):
    cmap = mpl.cm.get_cmap('RdBu_r')
    tmp_noaa = ds_noaa.sel({'time':slice(str(year_start), str(year_start + year_duration))})
    tmp_cera = ds_cera.sel({'time':slice(str(year_start), str(year_start + year_duration))})
    f, ax = plt.subplots(ncols=2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    plot_field(tmp_cera['s100_low'].mean(dim='time'), ax=ax[0], title='CERA', vmin=1.5, vmax=12)
    diff = (tmp_noaa.drop('number') - tmp_cera.drop('number')).squeeze()
    plot_field(diff['s100_low'].mean(dim='time'), ax=ax[1], title='NOAA - CERA', vmin=-3.5, vmax=3.4, cmap=cmap)
    for i_lat in [-15, -25, -35]:  # open_mfdataset orders latitude opposite
        ax[1].plot([-15, 35], [noaa_data.latitude[i_lat].values,
                               noaa_data.latitude[i_lat].values], color='grey', alpha=.5)
    plt.suptitle(str(year_start) + ' - ' + str(year_start + year_duration))
    plt.subplots_adjust(left=0.05, right=0.93)
    plt.savefig(plot_path + 'diffmean_' + str(year_start) + '_' + str(year_start + year_duration) + '.jpeg', dpi=300)
    plt.close('all')


for year_start in range(1901, 2001, 10):
    plot_diff_mean(noaa_data, cera_data, year_start, 10)

plot_diff_mean(noaa_data, cera_data, 1901, 108)
