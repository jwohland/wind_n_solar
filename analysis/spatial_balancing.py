import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl

sys.path.append("../")
from utils import add_letters

mpl.rcParams["axes.spines.left"] = True
mpl.rcParams["axes.spines.bottom"] = True


def plot_histograms(wind_shares, solar_shares, plotname):
    assert (
        99.8 < np.sum(wind_shares) + np.sum(solar_shares) < 100.2
    )  # allow some rounding errors

    base_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/"
    plot_path = base_path + "plots/analysis/country_assessment/"

    with open(base_path + "output/country_generation.pickle", "rb") as handle:
        country_generation = pickle.load(handle)

    amplitude_dict, results = {}, {}

    for data_timescale in ["seasonal", "multidecadal"]:
        for i, country_group in enumerate(country_generation[data_timescale].keys()):
            G_c = country_generation[data_timescale][country_group]
            if i == 0:
                G_tot = solar_shares[i] * G_c["solar"] + wind_shares[i] * G_c["wind"]
            else:
                G_tot += solar_shares[i] * G_c["solar"] + wind_shares[i] * G_c["wind"]

        amplitude_dict["All"] = {}
        amplitude_dict["All"]["others"] = np.round(
            float((G_tot.max() - G_tot.min()) / G_tot.min() * 100), 1
        )
        for j, exclude in enumerate(country_generation[data_timescale].keys()):
            for i, country_group in enumerate(
                country_generation[data_timescale].keys()
            ):
                if j == 0:
                    amplitude_dict[country_group] = {}
                G_c = country_generation[data_timescale][country_group]
                if exclude == country_group:
                    G_exclude = (
                        solar_shares[i] * G_c["solar"] + wind_shares[i] * G_c["wind"]
                    )
                    amplitude_dict[exclude]["alone"] = np.round(
                        float(
                            (G_exclude.max() - G_exclude.min()) / G_exclude.min() * 100
                        ),
                        1,
                    )
                else:
                    if (i == 0) or (i == 1 and j == 0):
                        G_tot = (
                            solar_shares[i] * G_c["solar"]
                            + wind_shares[i] * G_c["wind"]
                        )
                    else:
                        G_tot += (
                            solar_shares[i] * G_c["solar"]
                            + wind_shares[i] * G_c["wind"]
                        )
                    amplitude_dict[exclude]["others"] = np.round(
                        float((G_tot.max() - G_tot.min()) / G_tot.min() * 100), 1
                    )
        results[data_timescale] = pd.DataFrame.from_dict(
            data=amplitude_dict, orient="columns"
        )

    # plotting
    fs = 12
    f, axs = plt.subplots(ncols=2, figsize=(12, 8))
    for i, data_timescale in enumerate(results.keys()):
        ax = axs[i]
        ax.set_xticks(range(9))
        ticklabels = list(
            results[data_timescale].keys()[1:-1]
        )  # for plotting asthaetics, last group with linebreak
        ticklabels.append(
            "Greece, Bulgaria, Serbia, Croatia \n Bosnia and Herzogovina, Romania, Albania"
        )
        ax.set_xticklabels(
            ticklabels, rotation=45, ha="right", va="center", rotation_mode="anchor",
        )
        ax.axhline(
            results[data_timescale]["All"]["others"],
            color="red",
            ls="--",
            label="Cooperation",
        )
        ax.set_ylabel(
            r"Amplitude of variability   $\frac{G_\mathrm{max} -G_\mathrm{min}}{G_\mathrm{min}}$ [%]",
            fontsize=fs,
        )
        ax.set_title(data_timescale, fontsize=fs + 1)
        weighted_mean_amp = 0
        for x_off, countries in enumerate(results[data_timescale].keys()[1:]):
            w = 0.3
            ax.bar(
                x_off - w / 2,
                results[data_timescale][countries]["alone"],
                w,
                color="grey",
                alpha=0.6,
                label="Isolated" if (i == 1 and x_off == 0) else "",
            )
            ax.bar(
                x_off + w / 2,
                results[data_timescale][countries]["others"],
                w,
                color="blue",
                label="Others" if (i == 1 and x_off == 0) else "",
            )
            weighted_mean_amp += (
                (solar_shares[x_off] + wind_shares[x_off])
                / 100
                * results[data_timescale][countries]["alone"]
            )
        ax.axhline(weighted_mean_amp, color="grey", ls="--", label="Mean Isolated")
    plt.legend(fontsize=fs)
    add_letters(axs)
    plt.subplots_adjust(bottom=0.3, left=0.075, top=0.95, right=0.98)
    plt.savefig(plot_path + "balancing_" + plotname + ".jpeg", dpi=300)


"""
Different scenarios of installed capacities
"""
# current shares (sum(solar_shares) + sum(wind_shares)=1)
wind_shares = [
    8.6,
    9.4,
    7.5,
    20.4,
    4.3,
    3.5,
    2.0,
    1.0,
    2.5,
]
solar_shares = [
    4.1,
    3.6,
    6.7,
    15.3,
    7.7,
    0.2,
    1.6,
    0.1,
    1.6,
]

plot_histograms(wind_shares, solar_shares, "current")
plot_histograms(
    [x / np.sum(wind_shares) * 100 for x in wind_shares], [0] * 9, "current_windonly"
)
plot_histograms(
    [0] * 9, [x / np.sum(solar_shares) * 100 for x in solar_shares], "current_solaronly"
)

# even distribution of wind or solar while keeping total solar to total wind ratio untouched
plot_histograms([np.mean(wind_shares)] * 9, solar_shares, "even_wind")
plot_histograms(wind_shares, [np.mean(solar_shares)]*9, "even_solar")
plot_histograms([np.mean(wind_shares)] * 9, [np.mean(solar_shares)] * 9, "even_both")

# keeping inter-country distribution of wind and solar untouched while changing the continental wind/solar share
def rescale(shares, target_share):
    return [x/np.sum(shares)*target_share for x in shares]

plot_histograms(rescale(wind_shares, 40), rescale(solar_shares, 60), "windshare_40")
plot_histograms(rescale(wind_shares, 80), rescale(solar_shares, 20), "windshare_80")
