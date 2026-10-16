# Anomaly-filtered microclimate analysis

## Method

The previously defined field-only deployment windows were retained. Quality control was variable-specific: temperatures outside the DS1923 stated measurement range of −40 to 85 °C were set to missing for temperature summaries, and relative-humidity values outside the physically meaningful 0–100% interval were set to missing for humidity summaries. Valid measurements in the same rows were retained. No distribution-based filtering, winsorization, interpolation, or removal of high but plausible solar-heating peaks was performed.

- Temperature values excluded: 1 of 19344
- Humidity values excluded: 2761 of 19344

Daily minimum, mean, and maximum values were calculated by treatment. Tests used the 140 dates containing all three treatments. Global linear models tested treatment while blocking by date (`metric ~ treatment + date`). Pairwise comparisons were paired by date and Benjamini–Hochberg corrected across the three contrasts within each metric. Temporal autocorrelation is not explicitly modeled, so inference remains exploratory.

Four-day rolling analyses use overlapping four-day means of each daily minimum, mean, and maximum metric. Their detailed global and pairwise results are provided in `four_day_rolling_statistical_comparisons.csv`. Because adjacent rolling windows share three days, their p-values are especially sensitive to temporal autocorrelation and should be treated as descriptive/exploratory.

The repeated-measures GEE AR(1) analysis in `repeated_measures_gee_ar1_results.csv` uses replicate-level summaries. It removes shared daily conditions by centering observations on the mean across available experimental units for each date, then estimates treatment effects with treatment-replicate as the repeated cluster and AR(1) within-cluster correlation. This is the preferred inferential analysis.

## Temperature: global tests

| metric | n_complete_dates | F_treatment | p_treatment |
| --- | --- | --- | --- |
| daily_min_c | 140 | 84.468 | 2.1782e-29 |
| daily_mean_c | 140 | 3.5967 | 0.028696 |
| daily_max_c | 140 | 84.815 | 1.7554e-29 |

## Temperature: pairwise tests

| metric | contrast | n_complete_dates | mean_difference | t | p | p_adjusted_BH_within_metric |
| --- | --- | --- | --- | --- | --- | --- |
| daily_min_c | S+ vs UV- | 140 | 1.3156 | 11.136 | 5.511e-21 | 1.6533e-20 |
| daily_min_c | S+ vs S- | 140 | -0.16958 | -2.0506 | 0.042185 | 0.042185 |
| daily_min_c | UV- vs S- | 140 | -1.4852 | -9.1845 | 5.2535e-16 | 7.8802e-16 |
| daily_mean_c | S+ vs UV- | 140 | -0.12231 | -1.3701 | 0.17285 | 0.25343 |
| daily_mean_c | S+ vs S- | 140 | 0.071098 | 1.1468 | 0.25343 | 0.25343 |
| daily_mean_c | UV- vs S- | 140 | 0.19341 | 3.0014 | 0.0031867 | 0.0095601 |
| daily_max_c | S+ vs UV- | 140 | 1.9337 | 5.3815 | 3.0568e-07 | 3.0568e-07 |
| daily_max_c | S+ vs S- | 140 | 4.2249 | 14.595 | 7.7845e-30 | 2.3354e-29 |
| daily_max_c | UV- vs S- | 140 | 2.2912 | 7.1211 | 5.2512e-11 | 7.8768e-11 |

## Humidity: global tests

| metric | n_complete_dates | F_treatment | p_treatment |
| --- | --- | --- | --- |
| daily_min_percent | 140 | 24.643 | 1.4037e-10 |
| daily_mean_percent | 140 | 21.036 | 3.11e-09 |
| daily_max_percent | 140 | 14.771 | 8.0108e-07 |

## Humidity: pairwise tests

| metric | contrast | n_complete_dates | mean_difference | t | p | p_adjusted_BH_within_metric |
| --- | --- | --- | --- | --- | --- | --- |
| daily_min_percent | S+ vs UV- | 140 | 4.0267 | 3.2913 | 0.0012641 | 0.0012641 |
| daily_min_percent | S+ vs S- | 140 | -4.2263 | -5.4471 | 2.2587e-07 | 3.388e-07 |
| daily_min_percent | UV- vs S- | 140 | -8.253 | -5.7671 | 5.0073e-08 | 1.5022e-07 |
| daily_mean_percent | S+ vs UV- | 140 | 2.2823 | 4.7266 | 5.5323e-06 | 8.2984e-06 |
| daily_mean_percent | S+ vs S- | 140 | -0.58685 | -1.4925 | 0.13782 | 0.13782 |
| daily_mean_percent | UV- vs S- | 140 | -2.8691 | -5.5468 | 1.4201e-07 | 4.2602e-07 |
| daily_max_percent | S+ vs UV- | 140 | 0.11189 | 0.37188 | 0.71055 | 0.71055 |
| daily_max_percent | S+ vs S- | 140 | 1.6468 | 4.7768 | 4.468e-06 | 1.3404e-05 |
| daily_max_percent | UV- vs S- | 140 | 1.5349 | 4.184 | 5.0495e-05 | 7.5742e-05 |

## Four-day rolling comparisons

| variable | metric | contrast | treatment_1 | treatment_2 | n_complete_dates | treatment_1_mean | treatment_2_mean | mean_difference | paired_t | p_raw | global_F | global_p | p_adjusted_BH_within_metric |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| humidity | daily_max_percent | S+ vs UV- | S+ | UV- | 131 | 94.984 | 94.895 | 0.089468 | 0.41873 | 0.67611 | 23.874 | 3.0285e-10 | 0.67611 |
| humidity | daily_max_percent | S+ vs S- | S+ | S- | 131 | 94.984 | 93.289 | 1.6948 | 5.6998 | 7.6609e-08 | 23.874 | 3.0285e-10 | 2.2983e-07 |
| humidity | daily_max_percent | UV- vs S- | UV- | S- | 131 | 94.895 | 93.289 | 1.6053 | 5.2214 | 6.8525e-07 | 23.874 | 3.0285e-10 | 1.0279e-06 |
| humidity | daily_mean_percent | S+ vs UV- | S+ | UV- | 131 | 65.078 | 62.849 | 2.2292 | 5.5922 | 1.2659e-07 | 27.724 | 1.219e-11 | 1.8988e-07 |
| humidity | daily_mean_percent | S+ vs S- | S+ | S- | 131 | 65.078 | 65.607 | -0.52861 | -1.6785 | 0.095658 | 27.724 | 1.219e-11 | 0.095658 |
| humidity | daily_mean_percent | UV- vs S- | UV- | S- | 131 | 62.849 | 65.607 | -2.7578 | -6.0825 | 1.2324e-08 | 27.724 | 1.219e-11 | 3.6971e-08 |
| humidity | daily_min_percent | S+ vs UV- | S+ | UV- | 131 | 27.146 | 23.513 | 3.6331 | 4.478 | 1.6308e-05 | 41.335 | 2.5861e-16 | 1.6308e-05 |
| humidity | daily_min_percent | S+ vs S- | S+ | S- | 131 | 27.146 | 31.19 | -4.0444 | -7.0824 | 7.9807e-11 | 41.335 | 2.5861e-16 | 1.1971e-10 |
| humidity | daily_min_percent | UV- vs S- | UV- | S- | 131 | 23.513 | 31.19 | -7.6774 | -7.1384 | 5.9587e-11 | 41.335 | 2.5861e-16 | 1.1971e-10 |
| temperature | daily_max_c | S+ vs UV- | S+ | UV- | 131 | 49.909 | 47.976 | 1.9335 | 6.4416 | 2.102e-09 | 119.63 | 1.4586e-37 | 2.102e-09 |
| temperature | daily_max_c | S+ vs S- | S+ | S- | 131 | 49.909 | 45.664 | 4.2446 | 17.105 | 4.4176e-35 | 119.63 | 1.4586e-37 | 1.3253e-34 |
| temperature | daily_max_c | UV- vs S- | UV- | S- | 131 | 47.976 | 45.664 | 2.3111 | 8.4489 | 5.1089e-14 | 119.63 | 1.4586e-37 | 7.6634e-14 |
| temperature | daily_mean_c | S+ vs UV- | S+ | UV- | 131 | 27.39 | 27.537 | -0.14706 | -1.7069 | 0.090231 | 4.1815 | 0.016315 | 0.13535 |
| temperature | daily_mean_c | S+ vs S- | S+ | S- | 131 | 27.39 | 27.346 | 0.044218 | 0.73792 | 0.46189 | 4.1815 | 0.016315 | 0.46189 |
| temperature | daily_mean_c | UV- vs S- | UV- | S- | 131 | 27.537 | 27.346 | 0.19127 | 3.2923 | 0.0012801 | 4.1815 | 0.016315 | 0.0038403 |
| temperature | daily_min_c | S+ vs UV- | S+ | UV- | 131 | 16.51 | 15.173 | 1.3366 | 12.246 | 2.0922e-23 | 98.72 | 1.2681e-32 | 6.2766e-23 |
| temperature | daily_min_c | S+ vs S- | S+ | S- | 131 | 16.51 | 16.697 | -0.18683 | -2.4406 | 0.016009 | 98.72 | 1.2681e-32 | 0.016009 |
| temperature | daily_min_c | UV- vs S- | UV- | S- | 131 | 15.173 | 16.697 | -1.5234 | -9.7955 | 2.6364e-17 | 98.72 | 1.2681e-32 | 3.9546e-17 |
