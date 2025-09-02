
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

Use the core function `survblendr_extrapolate()` to fit the models and perform the blending. In this example, our trial data (`survblendr_demo`) was observed for 10 years, and we want to project future survival up to a 25-year horizon.

We configure the function with the following key parameters:

* `df = survblendr_demo`: The input dataset containing `time` and `status` columns.
* `t_obs = 10`: The time point where the trial's observation period ends. The model will use data up to this point to learn the short-term survival trend.
* `t_max = 25`: The maximum time horizon for the extrapolation. The survival curve will be generated up to this point.
* `anchor_t = 25`: This is the "anchor point" for the long-term external model. We are providing an external assumption about survival at this future time.
* `anchor_mean_Sa = 0.035`: This is our expert opinion or external evidence for the anchor point. We are assuming that the mean survival probability at 25 years is 3.5%.
* `nsim_inla = 2000` & `nsim_ext = 2000`: The number of posterior simulations to draw for the short-term (INLA) and long-term (Gompertz) models, respectively. Higher numbers lead to more stable estimates.
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

Finally, use `survblendr_summary_table()` to generate a table of key metrics, such as the predicted survival probability and the Restricted Mean Survival Time (RMST) at target times.
The extrapolation-period RMST is 
$$
\mathrm{RMST}(t)=\int_{t_{\mathrm{obs}}}^{t} S(u)\,du.
$$

```R
# Annual summary from year 10 to 25 (extrapolation-only RMST over [t_obs, t])
summary_df <- survblendr_summary_table(fit, years = 10:25)
knitr::kable(summary_df, digits = 3, format = "pipe")
```

| time| S_obs| S_ext| S_blend| RMST_obs| RMST_ext| RMST_blend|
|----:|-----:|-----:|-------:|--------:|--------:|----------:|
|   10| 0.274| 0.535|   0.419|    0.000|    0.000|      0.000|
|   11| 0.218| 0.490|   0.364|    0.246|    0.512|      0.392|
|   12| 0.172| 0.445|   0.314|    0.441|    0.980|      0.731|
|   13| 0.136| 0.401|   0.270|    0.595|    1.403|      1.023|
|   14| 0.110| 0.359|   0.233|    0.718|    1.783|      1.275|
|   15| 0.091| 0.317|   0.201|    0.818|    2.121|      1.492|
|   16| 0.076| 0.277|   0.172|    0.902|    2.418|      1.678|
|   17| 0.064| 0.238|   0.148|    0.972|    2.675|      1.838|
|   18| 0.056| 0.202|   0.126|    1.031|    2.896|      1.975|
|   19| 0.049| 0.168|   0.107|    1.084|    3.081|      2.092|
|   20| 0.043| 0.137|   0.090|    1.129|    3.234|      2.190|
|   21| 0.038| 0.110|   0.075|    1.170|    3.357|      2.272|
|   22| 0.035| 0.086|   0.062|    1.207|    3.455|      2.341|
|   23| 0.032| 0.065|   0.050|    1.240|    3.530|      2.396|
|   24| 0.029| 0.048|   0.040|    1.270|    3.587|      2.441|
|   25| 0.026| 0.035|   0.032|    1.298|    3.629|      2.477|

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Authors

hao hao ([@haohaostats](https://github.com/haohaostats))