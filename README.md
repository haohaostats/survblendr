
# survblendr: Adaptive Spline-Weighted Blending for Survival Extrapolation

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Last commit](https://img.shields.io/github/last-commit/haohaostats/survblendr)

**survblendr** is an R package for survival extrapolation using an adaptive spline-weighted blending method on the cumulative hazard scale. It smoothly combines a semi-parametric model fitted to the observed data with a parametric external model representing the long-term trend.

> The core of this method is a weighted blend on the cumulative hazard scale, combining a piecewise-exponential model (PEM) fitted via INLA with a Gompertz tail model calibrated by a user-defined anchor point.

---

## Key Features

* **Short-term Model**: Fits a flexible baseline hazard using a piecewise-exponential model with a first-order random walk (RW1) prior via `INLA`.
* **Long-term Model**: Utilizes a Gompertz model as the external tail, calibrated by a user-specified anchor point for long-term survival.
* **Adaptive Weighting**: Innovatively uses a monotone P-spline (via the `scam` package) to model the log-cumulative-hazard difference between the two models. The resulting weight is automatically calibrated to ensure the blended hazard is non-negative.
* **Integrated Workflow**: Provides a complete pipeline from data preparation, model fitting, and blending to visualization and result summarization.

## Installation

You can install the development version of `survblendr` from GitHub using the `remotes` package.

```R
# Install remotes if you haven't already
# install.packages("remotes")
# (INLA is hosted outside CRAN; the extra repo line is a safe fallback)
options(repos = c(getOption("repos"),
                  INLA = "https://inla.r-inla-download.org/R/stable"))
# for the stable version:
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

remotes::install_github("haohaostats/survblendr")
```

---

## Quick Start Example

The following example demonstrates the main workflow using the `survblendr_demo` dataset included with the package.

#### 1. Load the Package and Data

First, load the `survblendr` package and the example dataset.

```R
library(survblendr)

# survblendr_demo is a simulated dataset where follow-up is
# administratively censored at time t=10.
data(survblendr_demo)
head(survblendr_demo)
```

#### 2. Run the Extrapolation

Use the core function `survblendr_extrapolate()` to fit the models and perform the blending. We specify that the observation period ends at `t_obs = 10` and we want to extrapolate the survival curve up to `t_max = 25`. The `seed` argument ensures the results are reproducible.

```R
fit <- survblendr_extrapolate(
  survblendr_demo, t_obs = 10, t_max = 25, interval = 1,
  anchor_t = 25, anchor_mean_Sa = 0.035,
  nsim_inla = 2000, nsim_ext = 2000,
  inla_threads = 1,
  seed = 20240901
)
```

#### 3. Visualize the Results

The `plot_curves()` function provides an easy way to visualize the outcome. The plot clearly shows the short-term model (Cox PEM), the long-term external model (Gompertz), and the final smooth blended curve (Spline).

```R
p <- plot_curves(fit)
print(p)
```

![survblendr Plot Example](https://raw.githubusercontent.com/haohaostats/survblendr/main/man/figures/example_plot.svg)

#### 4. Summarize the Results

Finally, use `survblendr_summary_table()` to generate a table of key metrics, such as the predicted survival probability and the Restricted Mean Survival Time (RMST) calculated from time 0 up to each specified year.

```R
# Get an annual summary from year 10 to 25
summary_df <- survblendr_summary_table(fit, years = 10:25)
print(summary_df)
knitr::kable(summary_df, digits = 3, format = "pipe")
```

| time| S_obs| S_ext| S_blend| RMST_obs| RMST_ext| RMST_blend|
|----:|-----:|-----:|-------:|--------:|--------:|----------:|
|   10| 0.274| 0.535|   0.419|    4.963|    6.626|      5.900|
|   11| 0.218| 0.490|   0.364|    4.924|    6.912|      6.016|
|   12| 0.172| 0.445|   0.314|    4.867|    7.134|      6.082|
|   13| 0.136| 0.401|   0.270|    4.808|    7.294|      6.111|
|   14| 0.110| 0.359|   0.233|    4.762|    7.396|      6.119|
|   15| 0.091| 0.317|   0.201|    4.725|    7.442|      6.109|
|   16| 0.076| 0.277|   0.172|    4.696|    7.438|      6.084|
|   17| 0.064| 0.238|   0.148|    4.675|    7.389|      6.047|
|   18| 0.056| 0.202|   0.126|    4.660|    7.301|      6.001|
|   19| 0.049| 0.168|   0.107|    4.650|    7.183|      5.945|
|   20| 0.043| 0.137|   0.090|    4.643|    7.042|      5.882|
|   21| 0.038| 0.110|   0.075|    4.638|    6.888|      5.813|
|   22| 0.035| 0.086|   0.062|    4.635|    6.731|      5.741|
|   23| 0.032| 0.065|   0.050|    4.634|    6.582|      5.670|
|   24| 0.029| 0.048|   0.040|    4.633|    6.447|      5.601|
|   25| 0.026| 0.035|   0.032|    4.632|    6.332|      5.538|

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Authors

Hao Chen.
