import pandas as pd
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
mpl.use('Agg')


def prepare_CF(CF, Turbine_data):
    """
    Prepares data for subsequent analysis:
    1) # drop small turbines
    2) Drop turbines without power curve
    3) divide power curves by nominal capacity
    :param CF: DataFrame with turbine power curves
    :param Turbine_data: DataFrame with nominal capacities
    :return:
    """
    for turbine_name, values in CF.iterrows():
        if Turbine_data.loc[turbine_name, :]['nominal_power'] < 2.5 * 10**6:
            CF = CF.drop(turbine_name)
        elif Turbine_data.loc[turbine_name, :]['has_power_curve'] == False:
            CF = CF.drop(turbine_name)
        else:
            nominal_capacity = Turbine_data.loc[turbine_name, :]['nominal_power']  # in W
            CF.loc[turbine_name, :] /= nominal_capacity
    return CF


def plot_all_power_curves(CF, min_turbine, max_turbine, median_turbine, plot_path):
    """
    plot all power curves larger than 2.5 MW
    :param CF:
    :param min_turbine:
    :param max_turbine:
    :param median_turbine:
    :param plot_path:
    :return:
    """
    f, ax = plt.subplots(figsize=(12, 10))
    colors = [cm.Paired(x) for x in np.linspace(0, 1, 80)]
    i = 0
    for turbine_name, values in CF.iterrows():
        CF_values = values[np.isfinite(values)]
        if turbine_name in [min_turbine, max_turbine, median_turbine]:
            ax.plot(CF_values.index.astype(float), CF_values.values, label=turbine_name, lw=3, color='black')
        else:
            ax.plot(CF_values.index.astype(float), CF_values.values, alpha=.6, label=turbine_name, color=colors[i])
        i += 1
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.085), fancybox=True, shadow=True, ncol=7)
    ax.set_ylabel('Normalized power curve')
    ax.set_xlabel('Wind speed [m/s]')
    plt.subplots_adjust(top=.97, bottom=0.28)
    plt.savefig(plot_path + 'power_curves.jpeg', dpi=300)


def plot_representative_power_curves(CF, min_turbine, max_turbine, median_turbine, plot_path):
    """
    plots three representative power curves
    :param CF:
    :param min_turbine:
    :param max_turbine:
    :param median_turbine:
    :param plot_path:
    :return:
    """
    f, ax = plt.subplots(figsize=(12, 10))
    for turbine_name, values in CF.iterrows():
        CF_values = values[np.isfinite(values)]
        if turbine_name in [min_turbine, max_turbine, median_turbine]:
            ax.plot(CF_values.index.astype(float), CF_values.values, label=turbine_name, lw=3)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.085), fancybox=True, shadow=True, ncol=3)
    ax.set_ylabel('Normalized power curve')
    ax.set_xlabel('Wind speed [m/s]')
    plt.subplots_adjust(top=.97, bottom=0.15)
    plt.savefig(plot_path + 'power_curves_representative.jpeg', dpi=300)


def determine_representative_turbines(CF):
    """
    determines three representative turbines. They are the one with the highest, median and lowest power value
    at a wind speed of 7 m/s
    :return:
    """
    max_turbine, min_turbine = CF['7.0'].idxmax(), CF['7.0'].idxmin()
    pseudo_median = sorted(CF['7.0'].values)[int(CF['7.0'].size / 2)]
    median_turbine = CF['7.0'][CF['7.0'] == pseudo_median].index[0]
    return max_turbine, median_turbine, min_turbine


base_path = '/home/janwohland/Documents/Europmult/renewable_generation/wind_n_solar/'
windlib_path = base_path + 'data/wind-python-windpowerlib-v0.2.0/windpowerlib/oedb/'
plot_path = base_path + 'plots/wind_power_conversion/'

CF = pd.read_csv(windlib_path + 'power_curves.csv', index_col=0, header=0)
Turbine_data = pd.read_csv(windlib_path + 'turbine_data.csv', index_col=0,
                           usecols=['turbine_type', 'nominal_power', 'has_power_curve'])
CF = prepare_CF(CF, Turbine_data)
max_turbine, median_turbine, min_turbine = determine_representative_turbines(CF)
plot_all_power_curves(CF, min_turbine, max_turbine, median_turbine, plot_path)
plot_representative_power_curves(CF, min_turbine, max_turbine, median_turbine, plot_path)
