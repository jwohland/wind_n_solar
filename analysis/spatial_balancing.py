import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
from mpl_toolkits.axisartist.parasite_axes import SubplotHost
import matplotlib.ticker as ticker

sys.path.append("../")
from utils import add_letters

mpl.rcParams["axes.spines.left"] = True
mpl.rcParams["axes.spines.bottom"] = True

import os

BASE_PATH = os.path.join(os.path.dirname(__file__), "../")


def calc_histograms(wind_shares, solar_shares):
    assert (
        99.8 < np.sum(wind_shares) + np.sum(solar_shares) < 100.2
    )  # allow some rounding errors

    with open(BASE_PATH + "output/country_generation.pickle", "rb") as handle:
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
    return results


def plot_histograms(results, plotname, solar_shares, wind_shares):
    # plotting
    fs = 12
    f, axs = plt.subplots(ncols=2, figsize=(12, 8))
    plot_path = BASE_PATH + "plots/"
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
            "% Amplitude of variability", fontsize=fs,
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
    # plt.savefig(plot_path + "balancing_" + plotname + ".jpeg", dpi=300)
    results["multidecadal"]["Mean_isolated"] = weighted_mean_amp
    return results


COLORS = {
    "isolated": "chocolate",
    "cooperation": "firebrick",
}


def plot_overview_hists(scenario_results, solar_shares, wind_shares):
    fs = 12
    f = plt.figure(figsize=(12, 8))
    axs = np.array(
        [plt.subplot(121), SubplotHost(f, 122)]
    )  # SubplotHost instance needed for double x axis
    for ax in axs:
        f.add_subplot(ax)
    plot_path = BASE_PATH + "plots/"
    # plot hist of currently installed capacities
    results = scenario_results["Both"]
    data_timescale = "multidecadal"
    axs[0].set_xticks(range(9))
    ticklabels = list(
        results[data_timescale].keys()[1:-2]
    )  # for plotting asthaetics, last group with linebreak
    ticklabels.append("GR, BG, RS, HR, BA, RO, AL")
    axs[0].set_xticklabels(
        ticklabels, rotation=45, ha="right", va="center", rotation_mode="anchor",
    )
    axs[0].axhline(
        results[data_timescale]["All"]["others"],
        color=COLORS["cooperation"],
        ls="--",
        label="Cooperation",
    )
    axs[0].set_ylabel(
        "% Amplitude of variability", fontsize=fs,
    )
    weighted_mean_amp = 0
    for x_off, countries in enumerate(results[data_timescale].keys()[1:-1]):
        w = 0.3
        axs[0].bar(
            x_off - w / 2,
            results[data_timescale][countries]["alone"],
            w,
            color=COLORS["isolated"],
            alpha=0.8,
            label="Isolated" if x_off == 0 else "",
        )
        axs[0].bar(
            x_off + w / 2,
            results[data_timescale][countries]["others"],
            w,
            color="blue",
            label="Others" if x_off == 0 else "",
        )
        weighted_mean_amp += (
            (solar_shares[x_off] + wind_shares[x_off])
            / 100
            * results[data_timescale][countries]["alone"]
        )
    axs[0].axhline(
        weighted_mean_amp, color="grey", ls="--", label="Mean Isolated", alpha=0.8
    )
    axs[0].arrow(
        x=x_off + 0.55,
        y=weighted_mean_amp,
        dx=0,
        dy=-1 * (weighted_mean_amp - results[data_timescale]["All"]["others"]),
        color=COLORS["cooperation"],
        length_includes_head=True,
        width=0.2,
        head_width=0.5,
        alpha=0.8,
        head_length=0.35,
    )
    axs[0].set_xlim(xmax=8.85)
    axs[0].legend(fontsize=fs, ncol=2, frameon=False)

    # overview histogram with cooperation and autarky penalty
    axs[1].set_ylabel(
        "% Amplitude of variability", fontsize=fs,
    )
    for x_off, scenario in enumerate(scenario_results.keys()):
        results = scenario_results[scenario]["multidecadal"]
        axs[1].bar(
            x_off,
            results["All"],
            color="firebrick",
            alpha=0.8,
            label="Cooperation" if x_off == 0 else "",
        )
        axs[1].bar(
            x_off,
            results["Mean_isolated"] - results["All"],
            bottom=results["All"],
            color=COLORS["isolated"],
            alpha=0.8,
            label="Autarky penalty" if x_off == 0 else "",
        )

    axs[1].legend(fontsize=fs, frameon=False)
    axs[1].set_xticks(np.arange(8))
    axs[1].set_xticklabels([x for x in scenario_results.keys()],)
    plt.setp(
        axs[1].axis["bottom"].major_ticklabels,
        #     rotation=45,
        #     ha="center",
        va="top",
        #     rotation_mode="anchor",
    )

    # Add a second line of xlabels
    ax2 = axs[1].twiny()
    offset = 0, -30  # Position of the second axis
    new_axisline = ax2.get_grid_helper().new_fixed_axis
    ax2.axis["bottom"] = new_axisline(loc="bottom", axes=ax2, offset=offset)
    ax2.axis["top"].set_visible(False)
    axs[1].axis["right"].set_visible(
        False
    )  # axis cloning activates right hand side y axis

    ax2.set_xticks([-0.5, 2.5, 5.2, 7.5])
    ax2.xaxis.set_major_formatter(ticker.NullFormatter())
    ax2.xaxis.set_minor_locator(ticker.FixedLocator([1, 3.8, 6.4]))
    ax2.axis["bottom"].minor_ticks.set_ticksize(0)
    ax2.xaxis.set_minor_formatter(
        ticker.FixedFormatter(
            [
                "Current capacity\nand distribution",
                "Current capacity\nwith even\ndistribution of plants",
                "Current distribution\nof plants with\nadjusted share of wind",
            ]
        )
    )
    # Add lines for orientation between both x axes
    def add_line(f, x1):
        x1 /= 3600.0
        x0 = x1
        # x0 = x1 + 220.0 / 3600
        y0 = 1 - 1705.0 / 2400
        y1 = 1 - 1925.0 / 2400
        line = plt.Line2D(
            [x0, x1],
            [y0, y1],
            color="grey",
            linestyle="--",
            alpha=0.4,
            transform=f.transFigure,
        )
        f.add_artist(line)

    for x1 in [2610, 3130]:
        add_line(f, x1)

    add_letters(axs, fs=12)
    for i in range(2):
        axs[i].set_ylim(ymax=6.8)
    plt.subplots_adjust(bottom=0.3, left=0.075, top=0.95, right=0.98)
    # plt.savefig(plot_path + "balancing_overview.jpeg", dpi=300)

    return f, axs


def load_scenario_results():
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

    scenario_results = {}

    scenario_results["Both"] = plot_histograms(
        calc_histograms(wind_shares, solar_shares), "current", solar_shares, wind_shares
    )
    scenario_results["Wind"] = plot_histograms(
        calc_histograms([x / np.sum(wind_shares) * 100 for x in wind_shares], [0] * 9),
        "current_windonly",
        solar_shares,
        wind_shares,
    )
    scenario_results["Solar"] = plot_histograms(
        calc_histograms(
            [0] * 9, [x / np.sum(solar_shares) * 100 for x in solar_shares]
        ),
        "current_solaronly",
        solar_shares,
        wind_shares,
    )

    # even distribution of wind or solar while keeping total solar to total wind ratio untouched
    scenario_results["Both\neven"] = plot_histograms(
        calc_histograms([np.mean(wind_shares)] * 9, [np.mean(solar_shares)] * 9),
        "even_both",
        solar_shares,
        wind_shares,
    )
    scenario_results["Wind\neven"] = plot_histograms(
        calc_histograms([np.mean(wind_shares)] * 9, solar_shares),
        "even_wind",
        solar_shares,
        wind_shares,
    )
    scenario_results["Solar\neven"] = plot_histograms(
        calc_histograms(wind_shares, [np.mean(solar_shares)] * 9),
        "even_solar",
        solar_shares,
        wind_shares,
    )

    # keeping inter-country distribution of wind and solar untouched while changing the continental wind/solar share
    def rescale(shares, target_share):
        return [x / np.sum(shares) * target_share for x in shares]

    scenario_results["40%\nWind"] = plot_histograms(
        calc_histograms(rescale(wind_shares, 40), rescale(solar_shares, 60)),
        "windshare_40",
        solar_shares,
        wind_shares,
    )
    scenario_results["80%\nWind"] = plot_histograms(
        calc_histograms(rescale(wind_shares, 80), rescale(solar_shares, 20)),
        "windshare_80",
        solar_shares,
        wind_shares,
    )

    return scenario_results, solar_shares, wind_shares


if __name__ == "__main__":
    # make paper overview Figure
    scenario_results, solar_shares, wind_shares = load_scenario_results()
    plot_overview_hists(scenario_results, solar_shares, wind_shares)
