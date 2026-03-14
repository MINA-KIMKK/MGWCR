# MGWCR: Multiscale Geographically Weighted Cox Regression

This repository contains the data and code to reproduce all findings reported in:

**"Multiscale Geographically Weighted Cox Regression with Covariate-Specific Bandwidth"**

## Repository Structure

```
MGWCR/
├── README.md
├── R/                              # Core functions
│   ├── bw.sel.r
│   ├── gw.weight.r
│   ├── gw.dist.r
│   └── mgwcr_func.R
├── src/                            # C++ source files for parallel computing
│   ├── cox_derivatives_parallel.cpp
│   ├── gw_reg.cpp
│   └── MGWmodel.cpp
├── Application/                    # Empirical study (Section 4)
│   ├── final_data.csv
│   ├── 1__GWCox_Application.R
│   ├── 2__MGWCR_Application.R
│   └── 3__Application_figure.R
└── Simulation/                     # Simulation study (Section 3)
    ├── 1__GLBCox_Simulation_Scenario1.R
    ├── 2__GWCox_Simulation_Scenario1.R
    ├── 3__MGWCR_Simulation_Scenario1.R
    ├── 4__GLBCox_Simulation_Scenario2.R
    ├── 5__GWCox_Simulation_Scenario2.R
    ├── 6__MGWCR_Simulation_Scenario2.R
    └── 7__Simulation_figure.R
```

## Software Requirements

- R (>= 4.0)
- A C++ compiler compatible with Rcpp
- Required R packages:

```r
install.packages(c("MASS", "foreach", "doParallel", "cvTools", "dplyr",
                   "survival", "sp", "parallel", "Rcpp", "RcppParallel",
                   "lubridate", "ggplot2", "sf", "tigris"))
```

## Important Notes

- All scripts use **relative paths**. Before running, set your working directory to the root of this repository:
```r
setwd("/path/to/MGWCR")
```
- The `num_threads` parameter in each script controls the number of CPU cores used for parallel computing. Adjust this value based on your system (e.g., `parallel::detectCores() - 1`).
- Each simulation runs 1,000 iterations and may take several hours depending on hardware.
- Scripts within each section (Application or Simulation) should be run **sequentially** in the numbered order, as later scripts depend on outputs from earlier ones.

## Reproduction Instructions

### Simulation Study (Section 3)

The simulation study evaluates the performance of MGWCR under two scenarios.

**Scenario 1: Same Variation Scale** (two covariates with uniform spatial scale)

| Step | Script | Output |
|------|--------|--------|
| 1 | `Simulation/1__GLBCox_Simulation_Scenario1.R` | Global Cox model estimates |
| 2 | `Simulation/2__GWCox_Simulation_Scenario1.R` | GWCR bandwidth and coefficient estimates |
| 3 | `Simulation/3__MGWCR_Simulation_Scenario1.R` | MGWCR covariate-specific bandwidth and coefficient estimates |

**Scenario 2: Varied Spatial Heterogeneity** (three covariates: constant, medium, high variation)

| Step | Script | Output |
|------|--------|--------|
| 4 | `Simulation/4__GLBCox_Simulation_Scenario2.R` | Global Cox model estimates |
| 5 | `Simulation/5__GWCox_Simulation_Scenario2.R` | GWCR bandwidth and coefficient estimates |
| 6 | `Simulation/6__MGWCR_Simulation_Scenario2.R` | MGWCR covariate-specific bandwidth and coefficient estimates |

**Table 1 & Figure 2:**

| Step | Script | Output |
|------|--------|--------|
| 7 | `Simulation/7__Simulation_figure.R` | Reproduces Table 1 (performance metrics) and Figure 2 (estimated coefficient surfaces) |

### Empirical Study (Section 4)

The empirical study applies MGWCR to analyze the spatially varying risk factors of depression using Nurses' Health Study (NHS) data in New York State.

**Table 2 & Figures 3-6:**

| Step | Script | Output |
|------|--------|--------|
| 1 | `Application/1__GWCox_Application.R` | GWCR single bandwidth and coefficient estimates (saves `.RData` files) |
| 2 | `Application/2__MGWCR_Application.R` | MGWCR covariate-specific bandwidths (Table 2), spatially varying coefficients, and p-values |
| 3 | `Application/3__Application_figure.R` | Reproduces Figures 3-6 (coefficient maps and p-value maps) |

