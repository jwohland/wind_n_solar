"""
Correction of the long-term evolution of 100m wind speeds in CERA20C.

The approach consists of

    a) calculation of 3 hourly s_100 from u_100 and v_100 in CERA20C and 20CRv3
    b) regridding of 20CRv3 to the CERA20C grid
    c) calculation of low-pass filtered timeseries (Lanczos)
    d) Computing s_100_CERA_corrected = s_100_CERA + <s_100_20CRv3> - <s_100_CERA20C>

This script does d).
"""
from itertools import permutations
import glob
import sys
import xarray as xr


def build_filedic(data_path, lanczos_path):
    """
    builds a directory with paths to raw data ['CERA'] and the lanczos-filtered versions ['lanczos(X)']
    :param data_path: path to raw data
    :param lanczos_path: path to lanczos filtered data
    :return:
    """
    filedic = {}
    filedic['CERA'] = glob.glob(data_path + 'CERA20C/' + '*instant.nc')
    filedic['lanczos(CERA)'] = glob.glob(lanczos_path + '*CERA20C*.nc')
    filedic['lanczos(20CR)'] = glob.glob(lanczos_path + '*20CR*.nc')
    return filedic


def determine_combination(file_dic, file_index):
    """
    Maps scalar file_index to a combination of a CERA and 20CR ensemble member
    :param file_index:
    :return:
    """
    perms = list(permutations(range(len(file_dic['lanczos(CERA'])), 2))
    if file_index > len(perms):
        raise ValueError('file_index larger than number of permutations')
    cera_index, noaa_index = perms[file_index]
    lanczos_cera_path = file_dic['lanczos(CERA)'][cera_index]
    cera_stream = lanczos_cera_path.split('_')[1]  # todo this assumes that the stream number is given, not the index
    lanczos_20cr_path = file_dic['lanczos(20CR)'][noaa_index]
    return cera_stream, lanczos_cera_path, lanczos_20cr_path


def correct_longterm(s_raw, s_lanczos_cera, s_lanczos_noaa):
    """
    :param s_raw: 3h raw CERA20C wind speed data
    :param s_lanczos_cera: lanczos filtered CERA20C wind speed data
    :param s_lanczos_noaa: lanczos filtered 20CR wind speed data
    :return:
    """
    return s_raw - s_lanczos_cera + s_lanczos_noaa


def open_files(cera_stream, filedic, lanczos_cera_path, lanczos_noaa_path):
    s_raw = xr.open_mfdataset(filedic['CERA'],
                              chunks={'latitude': 5, 'longitude': 5},
                              combine='by_coords')
    s_raw.sel({'number': cera_stream}, drop=True)
    s_lanczos_cera = xr.open_mfdataset(lanczos_cera_path,
                                       chunks={'latitude': 5, 'longitude': 5},
                                       combine='by_coords')
    s_lanczos_noaa = xr.open_mfdataset(lanczos_noaa_path,
                                       chunks={'latitude': 5, 'longitude': 5},
                                       combine='by_coords')
    return s_raw, s_lanczos_cera, s_lanczos_noaa


def main(file_index):
    data_path = '/cluster/work/apatt/wojan/data/'  # todo update path
    lanczos_path = '/cluster/work/apatt/wojan/data/'  # todo update path
    out_path = ''  # todo update path

    filedic = build_filedic(data_path, lanczos_path)
    cera_stream, lanczos_cera_path, lanczos_noaa_path = determine_combination(filedic, file_index)
    s_raw, s_lanczos_cera, s_lanczos_noaa = open_files(cera_stream, filedic, lanczos_cera_path, lanczos_noaa_path)
    s_corrected = correct_longterm(s_raw, s_lanczos_cera, s_lanczos_noaa)
    noaa_member = '3'  # todo infer member number

    name = 'corrected_s100_cera_' + cera_stream + '_noaa_' + noaa_member + '.nc'
    s_corrected.to_netcdf(out_path + name)


file_index = int(sys.argv[1])
main(file_index)