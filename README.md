
# survblendr: Adaptive Spline-Weighted Blending for Survival Extrapolation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
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
|   10| 0.274| 0.537|   0.422|    4.958|    6.645|      5.917|
|   11| 0.217| 0.492|   0.366|    4.916|    6.934|      6.034|
|   12| 0.170| 0.448|   0.316|    4.852|    7.159|      6.098|
|   13| 0.133| 0.404|   0.271|    4.783|    7.322|      6.122|
|   14| 0.106| 0.361|   0.233|    4.730|    7.427|      6.126|
|   15| 0.087| 0.320|   0.200|    4.689|    7.477|      6.113|
|   16| 0.072| 0.280|   0.172|    4.661|    7.476|      6.089|
|   17| 0.061| 0.241|   0.147|    4.635|    7.430|      6.049|
|   18| 0.052| 0.205|   0.126|    4.619|    7.345|      6.002|
|   19| 0.046| 0.171|   0.107|    4.608|    7.230|      5.947|
|   20| 0.040| 0.140|   0.090|    4.600|    7.093|      5.886|
|   21| 0.036| 0.113|   0.075|    4.595|    6.943|      5.819|
|   22| 0.033| 0.089|   0.062|    4.594|    6.790|      5.752|
|   23| 0.030| 0.068|   0.051|    4.591|    6.643|      5.682|
|   24| 0.027| 0.052|   0.041|    4.588|    6.511|      5.616|
|   25| 0.025| 0.038|   0.033|    4.586|    6.397|      5.555|

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Authors

Hao Chen.
