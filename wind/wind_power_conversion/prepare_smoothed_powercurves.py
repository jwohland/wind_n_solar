from windpowerlib import wind_turbine as wt
import windpowerlib.power_curves as pc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def interpol(power):
    """
    Cubic spline interpolation of power curve to a resoltion of 0.1 m/s
    :param power: DataFrame with index wind speed and value power generation
    :return:
    """
    start, end = power.first_valid_index(), power.last_valid_index()
    power = power.reindex(np.arange(start, end+0.01, 0.01), tolerance=10**(-7), method='nearest')
    power = power.interpolate(method='spline', order=3, s=0)
    power.index = np.round(power.index, decimals=2)
    return power


base_path = '/home/janwohland/Documents/Europmult/renewable_generation/wind_n_solar/'
plot_path = base_path + 'plots/wind_power_conversion/'
out_path = '/home/janwohland/Documents/Europmult/renewable_generation/wind_n_solar/output/'
powercurve_data = base_path + '/data/wind-python-windpowerlib-v0.2.0/windpowerlib/oedb/power_curves.csv'

rep_turbines = ['E-126/7580', 'SWT120/3600', 'SWT142/3150']
f, ax = plt.subplots(ncols=3, sharey=True, figsize=(14, 5))


for i, turbine_name in enumerate(rep_turbines):
    power_curve = wt.get_turbine_data_from_file(turbine_name, powercurve_data)
    power_curve.value /= power_curve.value.max() # normalization
    ax[i].plot(power_curve.set_index('wind_speed'), label='manufacturer' if i == 1 else '')
    for method in ['Staffell_Pfenninger', 'turbulence_intensity']:
        power_curve_smoothed = \
            pc.smooth_power_curve(power_curve.wind_speed,
                                  power_curve.value,
                                  standard_deviation_method=method,
                                  turbulence_intensity=0.1246 if method == 'turbulence_intensity' else ''
                                  )
        # turbulence intensity taken from eq. 3-13 in Knorr (2016)
        power_curve_smoothed = power_curve_smoothed.set_index('wind_speed')
        power_curve_smoothed = interpol(power_curve_smoothed)
        ax[i].plot(power_curve_smoothed, label=method if i == 1 else '')
        if method == 'turbulence_intensity':
            #save
            df_final = power_curve_smoothed.rename(columns={'value': turbine_name})
            df_final.to_pickle(out_path + 'final_power_curve_' + str(i) + '.p')
    ax[i].set_title(turbine_name)
    ax[i].set_xlabel('Wind speed [m/s]')
ax[1].legend()
ax[0].set_ylabel('Normalized power curve [W/W]')
plt.subplots_adjust(left=0.05, right=0.98, wspace=0.1)
plt.savefig(plot_path + 'smoothing.jpeg', dpi=300)

