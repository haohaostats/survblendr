
# aswb: Adaptive Spline-Weighted Blending for Survival Extrapolation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
**aswb** is an R package for survival extrapolation using an adaptive spline-weighted blending method on the cumulative hazard scale. It smoothly combines a semi-parametric model fitted to the observed data with a parametric external model representing the long-term trend.

> The core of this method is a weighted blend on the cumulative hazard scale, combining a piecewise-exponential model (PEM) fitted via INLA with a Gompertz tail model calibrated by a user-defined anchor point.

---
  
  ## Key Features
  
  * **Short-term Model**: Fits a flexible baseline hazard using a piecewise-exponential model with a first-order random walk (RW1) prior via `INLA`.
* **Long-term Model**: Utilizes a Gompertz model as the external tail, calibrated by a user-specified anchor point for long-term survival.
* **Adaptive Weighting**: Innovatively uses a monotone P-spline (via the `scam` package) to model the log-cumulative-hazard difference between the two models. The resulting weight is automatically calibrated to ensure the blended hazard is non-negative.
* **Integrated Workflow**: Provides a complete pipeline from data preparation, model fitting, and blending to visualization and result summarization.

## Installation

You can install the development version of `aswb` from GitHub using the `remotes` package.

```R
# Install remotes if you haven't already
# install.packages("remotes")

remotes::install_github("haohaostats/aswb")
```

---
  
  ## Quick Start Example
  
  The following example demonstrates the main workflow using the `aswb_demo` dataset included with the package.

#### 1. Load the Package and Data

First, load the `aswb` package and the example dataset.

```R
library(aswb)

# aswb_demo is a simulated dataset where follow-up is 
# administratively censored at time t=10.
data(aswb_demo)
head(aswb_demo)
```

#### 2. Run the Extrapolation

Use the core function `aswb_extrapolate()` to fit the models and perform the blending. Here, we specify that the observation period ends at `t_obs = 10` and we want to extrapolate the survival curve up to `t_max = 25`.

```R
fit <- aswb::aswb_extrapolate(
  aswb_demo, t_obs = 10, t_max = 25, interval = 1,
  anchor_t = 25, anchor_mean_Sa = 0.035,
  nsim_inla = 200, nsim_ext = 200, inla_threads = 1
)
```

#### 3. Visualize the Results

The `plot_curves()` function provides an easy way to visualize the outcome. The plot clearly shows the short-term model (Cox PEM), the long-term external model (Gompertz), and the final smooth blended curve (Spline).

```R
p <- plot_curves(fit)
print(p)
```

![ASWB Plot Example](https://raw.githubusercontent.com/haohaostats/aswb/main/man/figures/README-example-1.png)


#### 4. Summarize the Results

Finally, use `aswb_summary_table()` to generate a table of key metrics, such as the predicted survival probability and the Restricted Mean Survival Time (RMST) calculated from time 0 up to each specified year.

```R
# Get an annual summary from year 10 to 20
summary_df <- aswb_summary_table(fit, years = 10:20)
print(summary_df)
```

---
  
  ## License
  
  This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Authors

Hao Chen and Clara Grazian.