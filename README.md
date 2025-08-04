# EcospaceCal: Automated Calibration of Ecospace Models in R

**EcospaceCal** is an **R-based framework** that automates the calibration of **Ecopath with Ecosim (EwE)** Ecospace models using powerful modern optimization techniques. It provides a **systematic, transparent, and repeatable** process for aligning spatial-temporal models with observational data.

---

## Overview

Calibrating complex ecosystem models is essential for credible projections but often involves **manual, time-consuming**, and **subjective tuning**—particularly in EwE models where **predator-prey vulnerability** settings require careful adjustment.

**EcospaceCal** tackles this by offering an **open-source**, fully-scripted interface that connects the EwE command-line executable with cutting-edge optimization methods. It’s designed to:

- Simplify and streamline the calibration process.
- Enhance reproducibility and transparency.
- Save researchers valuable time and effort.

---

## Key Features

- **Multiple Optimization Algorithms**

  - **Genetic Algorithm (GA):** Population-based, flexible, and widely used.
  - **CMA-ES:** Covariance Matrix Adaptation Evolution Strategy for efficient real-valued optimization.
  - **Bayesian Optimization (BO):** Perfect for minimizing computational load with fewer model runs.

- **Automated Workflow**

  - Manages parameter selection, model runs, and results aggregation in a single loop.

- **Parallel Processing**

  - Speeds up computation by distributing tasks across CPU cores.

- **Smart Caching**

  - Avoids redundant simulations by reusing previously evaluated parameter sets.

- **Robust Error Handling**

  - Includes timeouts and graceful recovery from failed model runs.

---

## System Requirements

- **R:** Version 4.0 or later
- **Ecopath with Ecosim (EwE):** Installed command-line version
- **R Packages:**\
  Install all required packages with:
  ```r
  install.packages(c(
    "GA", "cmaes", "rBayesianOptimization",
    "dplyr", "digest", "readxl",
    "doParallel", "R.utils"
  ))
  ```

---

## Installation & Setup

### 1. Clone the Repository

```bash
git clone https://github.com/mattwoodstock/EcospaceCal.git
cd EcospaceCal
```

### 2. Configure `setup.R`

- Open `setup.R` and edit the necessarry file locations and variables

---

## Quick Start

### 1. Open the Run Script

In RStudio or your preferred IDE, open:

```
R/CalProject.Rproj
analysis/example.R
```

### 2. Choose Your Optimization Method

Uncomment the relevant `run_calibration()` call for:

- `"GA"` – Genetic Algorithm
- `"CMAES"` – CMA-ES
- `"BO"` – Bayesian Optimization

### 3. (Optional) Customize Configuration

Override defaults by passing a config list:

```r
ga_config <- list(
  popSize = 20,
  run = 5
)

run_calibration(
  method = "GA",
  sensitivity_file = "./data/Gulf_of_Mexico/Master Vulnerability Table.xlsx",
  sheet_number = 2,
  output_base = "GA_Calibration_Run",
  n_cores = 8,
  config = ga_config
)
```

### 4. Run the Analysis

Execute the script to begin optimization:

```r
source("analysis/example.R")
```

> Results are stored in a timestamped directory like:\
> `GA_Calibration_Run/GA_Run_20250804_130900/`

---

## Project Structure

```
EcospaceCal/
├── R/                  # Core framework scripts
│   ├── calibration_framework.R
│   └── setup_template.R
│
├── data/               # Input data for case studies
│   ├── Gulf_of_Mexico/
│   └── **To Add More**
│
├── analysis/           # Script to launch analyses
│   └── example.R
│
└── results/            # Output folder (ignored by Git)
```

---

## Citation

If you use **EcospaceCal** in your research, please cite this repo

> Will add full citation once published.
---

## License

This project is released under the **MIT License**. See `LICENSE.txt` for details.

---

## Contact

For questions, issues, or feature requests:\
👉 [Open an Issue](https://github.com/mattwoodstock/EcospaceCal/issues)

