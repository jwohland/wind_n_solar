import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.path.append("../")
from utils import add_letters
import numpy as np

import spatial_balancing


sns.set_context("notebook", font_scale=1.2)
sns.set_style("ticks", {"axes.linewidth": 1.0, "font": "LiberationSans"})
plt.rcParams.update(
    {
        # Use sans math font
        "text.usetex": False,
        "mathtext.fontset": "stixsans",
        "svg.fonttype": "none",  # Make text editable
    }
)


def get_and_clean_data():
    scenario_results, _, _, = spatial_balancing.load_scenario_results()

    # Subplot A

    SHORT_NAMES = {
        "United Kingdom, Ireland": "UK, IR",
        "Portugal, Spain": "PT, ES",
        "France, Belgium, Netherlands": "FR, BE, NL",
        "Germany, Denmark": "DE, DK",
        "Italy, Austria, Switzerland, Slovenia": "IT, AU, CH, SI",
        "Sweden, Norway": "SE, NO",
        "Poland, Czech Republic, Slovakia, Hungary": "PL, CZ, SK, HU",
        "Lithuania, Latvia, Estonia, Finland": "LT, LV, EE, FI",
        "Greece, Bulgaria, Serbia, Croatia, Bosnia and Herzogovina, Romania, Albania": "GR, BG, RS, HR, BA, RO, AL",
    }

    for k, v in SHORT_NAMES.items():
        SHORT_NAMES[k] = v.replace(", ", "\n")

    df_countries = scenario_results["Both"]["multidecadal"].T

    mean_coop = df_countries.loc["All", "others"]
    mean_isolated = df_countries.loc["Mean_isolated", "others"]

    df_countries = df_countries.rename(
        columns={"others": "Others", "alone": "Isolated"}, index=SHORT_NAMES
    )

    df_countries = df_countries.drop(["All", "Mean_isolated"], axis=0)

    # Subplot B

    df_europe = pd.DataFrame(
        {k: scenario_results[k]["multidecadal"]["All"] for k in scenario_results.keys()}
    ).T
    df_europe_isolated = pd.DataFrame(
        {
            k: scenario_results[k]["multidecadal"]["Mean_isolated"]
            for k in scenario_results.keys()
        }
    ).T
    df_europe["Autarky penalty"] = df_europe_isolated["others"] - df_europe["others"]
    df_europe = df_europe.rename(columns={"others": "Cooperation"}).drop(
        "alone", axis=1
    )

    return df_countries, df_europe, mean_coop, mean_isolated


def draw_plot(df_countries, df_europe, mean_coop, mean_isolated):
    COLORS = {
        "others": "#5489bd",
        "isolated": "#fdc272",
        "cooperation": "#af4d80",
    }

    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

    # SUBPLOT A

    ax = ax0

    df_countries.plot.bar(
        ax=ax,
        rot=0,
        width=0.7,
        legend=False,
        color=[COLORS["others"], COLORS["isolated"]],
    )

    ax.set_ylabel("Amplitude of variability [%]")

    axhline_kwargs = dict(xmax=1.2, clip_on=False)
    ax.axhline(
        mean_coop,
        color=COLORS["cooperation"],
        label="Mean with cooperation",
        **axhline_kwargs
    )
    ax.axhline(
        mean_isolated,
        color=COLORS["isolated"],
        label="Mean with isolation",
        **axhline_kwargs
    )

    # Arrow
    ax.arrow(
        x=9.5,
        y=mean_isolated,
        dx=0,
        dy=-1 * (mean_isolated - mean_coop),
        color=COLORS["cooperation"],
        length_includes_head=True,
        width=0.2,
        head_width=0.5,
        head_length=0.35,
        clip_on=False,
    )

    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[-2], handles[-1], handles[1], handles[0]]
    labels = [labels[-2], labels[-1], labels[1], labels[0]]
    ax.legend(handles, labels, frameon=False, ncol=2, fontsize=12)

    ax.set_title(r"Country groups for current capacities", loc="left")

    # SUBPLOT B

    ax = ax1

    df_europe.plot.bar(
        ax=ax,
        stacked=True,
        rot=0,
        legend=False,
        color=[COLORS["cooperation"], COLORS["isolated"]],
    )

    axhline_kwargs = dict(xmax=0.09)
    ax.axhline(mean_coop, color=COLORS["cooperation"], **axhline_kwargs)
    ax.axhline(mean_isolated, color=COLORS["isolated"], **axhline_kwargs)

    axvline_kwargs = dict(ymin=-0.315, ymax=1.0, clip_on=False, color="grey", alpha=0.2)
    for x in [-0.5, 2.5, 5.5, 7.5]:
        ax.axvline(x, **axvline_kwargs)
    ax.axhline(-2.3, clip_on=False, color="grey", alpha=0.2)

    arrow_marker_kwargs = dict(color=COLORS["cooperation"], alpha=0.4)
    ax.vlines(
        -1.325, ymin=mean_isolated, ymax=6.45, clip_on=False, **arrow_marker_kwargs
    )
    ax.hlines(6.45, xmin=-1.325, xmax=0.3, clip_on=False, **arrow_marker_kwargs)
    ax.text(
        0.4,
        6.6,
        "Benefit of\ncooperation",
        horizontalalignment="left",
        verticalalignment="center",
        fontsize=12,
        bbox=dict(facecolor="none", edgecolor=COLORS["cooperation"], alpha=0.4),
    )

    text_kwargs = dict(
        horizontalalignment="center", verticalalignment="top", fontsize=10
    )
    ax.text(1, -1, "Current capacity\nand distribution", **text_kwargs)
    ax.text(4, -1, "Current capacity\nand even\nplant distribution", **text_kwargs)
    ax.text(
        6.5, -1, "Current plant\ndistribution\nand adjusted\nwind share", **text_kwargs
    )

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], frameon=False, fontsize=12)

    ax.set_title(r"Europe-wide mean for all scenarios", loc="left")

    # CLEAN UP PLOTS

    for ax in [ax0, ax1]:
        ax.set_ylim(0, 7.3)

    add_letters(np.array([ax0, ax1]), fs=12, y=1.03)

    return fig


if __name__ == "__main__":
    df_countries, df_europe, mean_coop, mean_isolated = get_and_clean_data()
    fig = draw_plot(df_countries, df_europe, mean_coop, mean_isolated)
    plot_path = "/cluster/work/apatt/wojan/renewable_generation/wind_n_solar/plots/analysis/country_assessment/"
    plt.savefig(
        plot_path + "balancing_overview.png",
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=300,
    )
