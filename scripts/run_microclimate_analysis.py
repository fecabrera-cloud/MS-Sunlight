#!/usr/bin/env python3
"""Build final anomaly-filtered microclimate figures and reports."""

from itertools import combinations
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.genmod.cov_struct import Autoregressive


REPO_ROOT = Path(__file__).resolve().parents[1]
HERE = REPO_ROOT / "results" / "reproduced" / "microclimate"
HERE.mkdir(parents=True, exist_ok=True)
DATA_ROOT = REPO_ROOT / "data" / "microclimate"
INPUT = DATA_ROOT / "field_only_anomaly_filtered.csv"
TREATMENTS = ["S+", "UV-", "S-"]
COLORS = {"S+": "#E3AC18", "UV-": "#77BCE8", "S-": "#555555"}
LABELS = {
    "S+": "Full sunlight",
    "UV-": "Reduced UV",
    "S-": "Reduced UV and blue–green",
}
LEGEND_TITLE = "Sunlight exposure"


def display_contrast(contrast):
    """Translate internal treatment codes for figure annotations."""
    first, second = contrast.split(" vs ")
    return f"{LABELS[first]} vs {LABELS[second]}"


def save_figure(figure, stem):
    figure.savefig(HERE / f"{stem}.png", dpi=300)
    figure.savefig(HERE / f"{stem}.pdf")
    plt.close(figure)


def format_dates(axis):
    locator = mdates.AutoDateLocator()
    axis.xaxis.set_major_locator(locator)
    axis.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    axis.grid(alpha=0.2)


def daily_replicate(data):
    frame = data.copy()
    frame["date"] = frame.datetime.dt.floor("D")
    daily = frame.groupby(["date", "treatment", "replicate"], as_index=False).agg(
        temperature_min=("temperature_c", "min"),
        temperature_mean=("temperature_c", "mean"),
        temperature_max=("temperature_c", "max"),
        humidity_min=("humidity_percent", "min"),
        humidity_mean=("humidity_percent", "mean"),
        humidity_max=("humidity_percent", "max"),
    )
    daily["unit"] = daily.treatment + "_rep" + daily.replicate.astype(str)
    daily.to_csv(HERE / "daily_replicate_microclimate.csv", index=False)
    return daily


def treatment_time_plot(data, variable, ylabel, stem):
    summary = data.groupby(["period", "treatment", "reading_index"], as_index=False).agg(
        datetime=("datetime", "mean"), mean=(variable, "mean"), sd=(variable, "std"), n=(variable, "count")
    )
    summary["se"] = summary.sd / np.sqrt(summary.n)
    figure, axis = plt.subplots(figsize=(24, 6), constrained_layout=True)
    for treatment in TREATMENTS:
        treatment_data = summary[summary.treatment == treatment]
        for index, (_, period) in enumerate(treatment_data.groupby("period", sort=False)):
            x, mean, se = period.datetime.to_numpy(), period["mean"].to_numpy(), period.se.to_numpy()
            axis.plot(x, mean, color=COLORS[treatment], linewidth=0.9, label=LABELS[treatment] if index == 0 else None)
            axis.fill_between(x, mean - se, mean + se, color=COLORS[treatment], alpha=0.18, linewidth=0)
    axis.set(xlabel="Date", ylabel=ylabel)
    axis.legend(title=LEGEND_TITLE, ncol=3, frameon=False)
    format_dates(axis)
    save_figure(figure, stem)


def treatment_daily(daily):
    return daily.groupby(["date", "treatment"], as_index=False).agg(
        temperature_min=("temperature_min", "mean"), temperature_mean=("temperature_mean", "mean"), temperature_max=("temperature_max", "mean"),
        humidity_min=("humidity_min", "mean"), humidity_mean=("humidity_mean", "mean"), humidity_max=("humidity_max", "mean"),
    )


def rolling_table(daily_treatment):
    rows = []
    metrics = ["temperature_min", "temperature_mean", "temperature_max", "humidity_min", "humidity_mean", "humidity_max"]
    dates = pd.date_range(daily_treatment.date.min(), daily_treatment.date.max(), freq="D")
    for treatment in TREATMENTS:
        frame = daily_treatment[daily_treatment.treatment == treatment].set_index("date").reindex(dates)
        for metric in metrics:
            values = frame[metric].rolling(4, min_periods=4).mean()
            rows.extend({"date": date, "treatment": treatment, "metric": metric, "rolling_4day_mean": value} for date, value in values.dropna().items())
    rolling = pd.DataFrame(rows)
    rolling.to_csv(HERE / "four_day_rolling_values.csv", index=False)
    return rolling


def rolling_lines(rolling):
    figure, axes = plt.subplots(2, 1, figsize=(24, 10), sharex=True, constrained_layout=True)
    for axis, metric, ylabel in zip(axes, ["temperature_mean", "humidity_mean"], ["Four-day mean temperature (°C)", "Four-day mean humidity (%)"]):
        for treatment in TREATMENTS:
            frame = rolling[(rolling.metric == metric) & (rolling.treatment == treatment)]
            axis.plot(frame.date, frame.rolling_4day_mean, color=COLORS[treatment], linewidth=2, label=LABELS[treatment])
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.2)
    axes[0].legend(title=LEGEND_TITLE, ncol=3, frameon=False)
    axes[-1].set_xlabel("Date")
    format_dates(axes[-1])
    save_figure(figure, "four_day_rolling_temperature_humidity")


def daily_extrema_plot(daily_treatment, variable):
    metrics = [f"{variable}_min", f"{variable}_mean", f"{variable}_max"]
    labels = ["Daily minimum", "Daily mean", "Daily maximum"]
    unit = "°C" if variable == "temperature" else "%"
    figure, axes = plt.subplots(3, 1, figsize=(24, 12), sharex=True, constrained_layout=True)
    for axis, metric, label in zip(axes, metrics, labels):
        for treatment in TREATMENTS:
            frame = daily_treatment[daily_treatment.treatment == treatment]
            axis.plot(frame.date, frame[metric], color=COLORS[treatment], linewidth=1, label=LABELS[treatment])
        axis.set_ylabel(f"{label} ({unit})")
        axis.grid(alpha=0.2)
    axes[0].legend(title=LEGEND_TITLE, ncol=3, frameon=False)
    axes[-1].set_xlabel("Date")
    format_dates(axes[-1])
    save_figure(figure, f"daily_{variable}_min_mean_max")


def rolling_boxplots(rolling):
    metrics = [["temperature_min", "temperature_mean", "temperature_max"], ["humidity_min", "humidity_mean", "humidity_max"]]
    figure, axes = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
    rng = np.random.default_rng(42)
    for row, row_metrics in enumerate(metrics):
        for column, metric in enumerate(row_metrics):
            axis = axes[row, column]
            frame = rolling[rolling.metric == metric]
            values = [frame.loc[frame.treatment == treatment, "rolling_4day_mean"].to_numpy() for treatment in TREATMENTS]
            box = axis.boxplot(values, tick_labels=[LABELS[x] for x in TREATMENTS], patch_artist=True, showfliers=False)
            for patch, treatment in zip(box["boxes"], TREATMENTS):
                patch.set_facecolor(COLORS[treatment]); patch.set_alpha(0.65)
            for position, values_one in enumerate(values, 1):
                axis.scatter(rng.normal(position, 0.04, len(values_one)), values_one, s=6, alpha=0.14, color="#222222")
            axis.set_title(metric.replace("_", " ").title())
            axis.grid(axis="y", alpha=0.2)
    axes[0, 0].set_ylabel("Temperature (°C)")
    axes[1, 0].set_ylabel("Relative humidity (%)")
    save_figure(figure, "four_day_rolling_boxplots")


def day_night_data(data):
    frame = data.copy()
    frame["date"] = frame.datetime.dt.floor("D")
    frame["day_period"] = np.where((frame.datetime.dt.hour >= 8) & (frame.datetime.dt.hour < 20), "Day (08:00–20:00)", "Night (20:00–08:00)")
    summary = frame.groupby(["date", "treatment", "replicate", "day_period"], as_index=False).agg(
        temperature_c=("temperature_c", "mean"), humidity_percent=("humidity_percent", "mean")
    )
    summary["unit"] = summary.treatment + "_rep" + summary.replicate.astype(str)
    summary.to_csv(HERE / "daily_day_night_temperature_by_replicate.csv", index=False)
    return summary


def matched_absolute_differences(frame, value_column, grouping_columns):
    rows = []
    grouped = frame.groupby(grouping_columns) if grouping_columns else [((), frame)]
    for group_values, group in grouped:
        group_values = group_values if isinstance(group_values, tuple) else (group_values,)
        wide = group.pivot_table(index="date", columns="treatment", values=value_column, aggfunc="mean")
        for first, second in combinations(TREATMENTS, 2):
            matched = wide[[first, second]].dropna()
            rows.append(
                {
                    **dict(zip(grouping_columns, group_values)),
                    "contrast": f"{first} vs {second}",
                    "n_matched_dates": len(matched),
                    "mean_absolute_difference": (matched[first] - matched[second]).abs().mean(),
                    "mean_signed_difference": (matched[first] - matched[second]).mean(),
                }
            )
    return pd.DataFrame(rows)


def day_night_plots(summary, statistics, absolute_differences):
    treatment_summary = summary.groupby(["date", "treatment", "day_period"], as_index=False).temperature_c.agg(["mean", "std", "count"]).reset_index()
    treatment_summary["se"] = treatment_summary["std"] / np.sqrt(treatment_summary["count"])
    periods = ["Day (08:00–20:00)", "Night (20:00–08:00)"]
    figure, axes = plt.subplots(2, 1, figsize=(24, 10), sharex=True, constrained_layout=True)
    for axis, period in zip(axes, periods):
        for treatment in TREATMENTS:
            frame = treatment_summary[(treatment_summary.day_period == period) & (treatment_summary.treatment == treatment)]
            axis.plot(frame.date, frame["mean"], color=COLORS[treatment], linewidth=1.2, label=LABELS[treatment])
            axis.fill_between(frame.date, frame["mean"] - frame.se, frame["mean"] + frame.se, color=COLORS[treatment], alpha=0.16)
        axis.set_ylabel(f"{period}\nmean ± SE (°C)")
        axis.grid(alpha=0.2)
    axes[0].legend(title=LEGEND_TITLE, ncol=3, frameon=False)
    axes[-1].set_xlabel("Date")
    format_dates(axes[-1])
    save_figure(figure, "temperature_day_night_timeseries")

    figure, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
    rng = np.random.default_rng(42)
    for axis, period in zip(axes, periods):
        frame = summary[summary.day_period == period]
        values = [frame.loc[frame.treatment == treatment, "temperature_c"].to_numpy() for treatment in TREATMENTS]
        box = axis.boxplot(values, tick_labels=[LABELS[x] for x in TREATMENTS], patch_artist=True, showfliers=False)
        for patch, treatment in zip(box["boxes"], TREATMENTS):
            patch.set_facecolor(COLORS[treatment]); patch.set_alpha(0.65)
        for position, treatment_values in enumerate(values, 1):
            axis.scatter(rng.normal(position, 0.055, len(treatment_values)), treatment_values, s=8, alpha=0.17, color="#222222", linewidths=0)
            median = np.nanmedian(treatment_values)
            axis.text(position, median, f"  {median:.2f} °C", ha="left", va="center", fontsize=9, fontweight="bold", bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 1.5})
        axis.set(title=period, ylabel="Daily replicate mean temperature (°C)")
        axis.grid(axis="y", alpha=0.2)
        period_stats = statistics[statistics.day_period == period]
        period_differences = absolute_differences[absolute_differences.day_period == period].set_index("contrast")
        lines = [f"GEE AR(1) global p = {period_stats.global_p.iloc[0]:.3f}"]
        for _, row in period_stats.iterrows():
            absolute = period_differences.loc[row.contrast, "mean_absolute_difference"]
            lines.append(f"{display_contrast(row.contrast)}: adj. p = {row.p_adjusted_BH:.3f}; mean |Δ| = {absolute:.2f} °C")
        axis.text(0.02, 0.98, "\n".join(lines), transform=axis.transAxes, va="top", ha="left", fontsize=9, bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "#777777"})
    save_figure(figure, "temperature_day_night_boxplots")


def humidity_day_night_plots(summary, statistics, absolute_differences):
    treatment_summary = summary.groupby(["date", "treatment", "day_period"], as_index=False).humidity_percent.agg(["mean", "std", "count"]).reset_index()
    treatment_summary["se"] = treatment_summary["std"] / np.sqrt(treatment_summary["count"])
    periods = ["Day (08:00–20:00)", "Night (20:00–08:00)"]
    figure, axes = plt.subplots(2, 1, figsize=(24, 10), sharex=True, constrained_layout=True)
    for axis, period in zip(axes, periods):
        for treatment in TREATMENTS:
            frame = treatment_summary[(treatment_summary.day_period == period) & (treatment_summary.treatment == treatment)]
            axis.plot(frame.date, frame["mean"], color=COLORS[treatment], linewidth=1.2, label=LABELS[treatment])
            axis.fill_between(frame.date, frame["mean"] - frame.se, frame["mean"] + frame.se, color=COLORS[treatment], alpha=0.16)
        axis.set_ylabel(f"{period}\nmean ± SE (%)")
        axis.grid(alpha=0.2)
    axes[0].legend(title=LEGEND_TITLE, ncol=3, frameon=False)
    axes[-1].set_xlabel("Date")
    format_dates(axes[-1])
    save_figure(figure, "humidity_day_night_timeseries")

    figure, axes = plt.subplots(1, 2, figsize=(17, 7), constrained_layout=True)
    rng = np.random.default_rng(42)
    for axis, period in zip(axes, periods):
        frame = summary[summary.day_period == period]
        values = [frame.loc[frame.treatment == treatment, "humidity_percent"].dropna().to_numpy() for treatment in TREATMENTS]
        box = axis.boxplot(values, tick_labels=[LABELS[x] for x in TREATMENTS], patch_artist=True, showfliers=False)
        for patch, treatment in zip(box["boxes"], TREATMENTS):
            patch.set_facecolor(COLORS[treatment]); patch.set_alpha(0.65)
        for position, treatment_values in enumerate(values, 1):
            axis.scatter(rng.normal(position, 0.055, len(treatment_values)), treatment_values, s=8, alpha=0.17, color="#222222", linewidths=0)
            median = np.nanmedian(treatment_values)
            axis.text(position, median, f"  {median:.2f}%", ha="left", va="center", fontsize=9, fontweight="bold", bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 1.5})
        axis.set(title=period, ylabel="Daily replicate mean relative humidity (%)")
        axis.grid(axis="y", alpha=0.2)
        period_stats = statistics[statistics.day_period == period]
        period_differences = absolute_differences[absolute_differences.day_period == period].set_index("contrast")
        lines = [f"GEE AR(1) global p = {period_stats.global_p.iloc[0]:.3f}"]
        for _, row in period_stats.iterrows():
            absolute = period_differences.loc[row.contrast, "mean_absolute_difference"]
            lines.append(f"{display_contrast(row.contrast)}: adj. p = {row.p_adjusted_BH:.3f}; mean |Δ| = {absolute:.2f} pp")
        axis.text(0.02, 0.98, "\n".join(lines), transform=axis.transAxes, va="top", ha="left", fontsize=9, bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "#777777"})
    save_figure(figure, "humidity_day_night_boxplots")


def bh(values):
    values = np.asarray(values); order = np.argsort(values); ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted); result[order] = np.minimum(adjusted, 1); return result


def gee_day_night(summary, value_column, output_stem):
    rows = []
    for period, frame in summary.groupby("day_period"):
        frame = frame.sort_values(["unit", "date"]).copy()
        frame["day_index"] = (frame.date - frame.date.min()).dt.days
        frame["centered"] = frame[value_column] - frame.groupby("date")[value_column].transform("mean")
        design = pd.get_dummies(frame.treatment, dtype=float).reindex(columns=TREATMENTS, fill_value=0).drop(columns="S+")
        design = sm.add_constant(design)
        model = sm.GEE(frame.centered.to_numpy(), design.to_numpy(), groups=frame.unit.to_numpy(), time=frame.day_index.to_numpy(), family=sm.families.Gaussian(), cov_struct=Autoregressive(grid=True)).fit()
        names = list(design.columns); params = pd.Series(model.params, index=names); covariance = pd.DataFrame(model.cov_params(), index=names, columns=names)
        restriction = np.zeros((2, 3)); restriction[0, 1] = 1; restriction[1, 2] = 1
        global_test = model.wald_test(restriction, scalar=True)
        pair_rows = []
        for first, second in combinations(TREATMENTS, 2):
            vector = pd.Series(0.0, index=names)
            if first != "S+": vector[first] += 1
            if second != "S+": vector[second] -= 1
            estimate = float(vector @ params); se = float(np.sqrt(vector @ covariance @ vector)); z = estimate / se
            pair_rows.append({"day_period": period, "contrast": f"{first} vs {second}", "estimated_difference": estimate, "robust_se": se, "z": z, "p_raw": 2 * norm.sf(abs(z)), "global_wald_chi2": float(global_test.statistic), "global_p": float(global_test.pvalue), "ar1_alpha": float(model.cov_struct.dep_params), "n_units": frame.unit.nunique(), "n_observations": len(frame)})
        adjusted = bh([row["p_raw"] for row in pair_rows])
        for row, value in zip(pair_rows, adjusted): row["p_adjusted_BH"] = value; rows.append(row)
    results = pd.DataFrame(rows)
    results.to_csv(HERE / f"{output_stem}_day_night_gee_ar1_results.csv", index=False)
    return results


def copy_statistics():
    for name in ["anomaly_filtering_report.md", "excluded_anomalous_values.csv", "repeated_measures_gee_ar1_results.csv", "repeated_measures_gee_ar1_report.md", "four_day_rolling_statistical_comparisons.csv"]:
        source = DATA_ROOT / name
        if source.exists():
            (HERE / name).write_bytes(source.read_bytes())


def main():
    data = pd.read_csv(INPUT, parse_dates=["datetime"])
    daily = daily_replicate(data)
    daily_treatment = treatment_daily(daily)
    daily_treatment.to_csv(HERE / "daily_treatment_microclimate.csv", index=False)
    rolling = rolling_table(daily_treatment)
    treatment_time_plot(data, "temperature_c", "Mean temperature ± SE (°C)", "temperature_mean_se_wide")
    treatment_time_plot(data, "humidity_percent", "Mean relative humidity ± SE (%)", "humidity_mean_se_wide")
    rolling_lines(rolling)
    daily_extrema_plot(daily_treatment, "temperature")
    daily_extrema_plot(daily_treatment, "humidity")
    rolling_boxplots(rolling)
    day_night = day_night_data(data)
    day_night_statistics = gee_day_night(day_night.dropna(subset=["temperature_c"]), "temperature_c", "temperature")
    day_night_differences = matched_absolute_differences(day_night, "temperature_c", ["day_period"])
    day_night_differences.to_csv(HERE / "temperature_day_night_mean_absolute_differences.csv", index=False)
    overall_daily_differences = matched_absolute_differences(daily_treatment, "temperature_mean", [])
    overall_daily_differences.to_csv(HERE / "daily_temperature_mean_absolute_differences.csv", index=False)
    day_night_plots(day_night, day_night_statistics, day_night_differences)
    humidity_statistics = gee_day_night(day_night.dropna(subset=["humidity_percent"]), "humidity_percent", "humidity")
    humidity_differences = matched_absolute_differences(day_night, "humidity_percent", ["day_period"])
    humidity_differences.to_csv(HERE / "humidity_day_night_mean_absolute_differences.csv", index=False)
    humidity_day_night_plots(day_night, humidity_statistics, humidity_differences)
    copy_statistics()


if __name__ == "__main__":
    main()
