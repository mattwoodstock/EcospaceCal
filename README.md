# 🌍 EcospaceCal: Automated Calibration of Ecospace Models in R

**EcospaceCal** is an **R-based framework** that automates the calibration of **Ecopath with Ecosim (EwE)** Ecospace models using powerful modern optimization techniques. It provides a **systematic, transparent, and repeatable** process for aligning spatial-temporal models with observational data.

---

## 📌 Overview

Calibrating complex ecosystem models is essential for credible projections but often involves **manual, time-consuming**, and **subjective tuning**—particularly in EwE models where **predator-prey vulnerability** settings require careful adjustment.

**EcospaceCal** tackles this by offering an **open-source**, fully-scripted interface that connects the EwE command-line executable with cutting-edge optimization methods. It’s designed to:

- Simplify and streamline the calibration process.
- Enhance reproducibility and transparency.
- Save researchers valuable time and effort.

---

## 🚀 Key Features

- 🔧 **Multiple Optimization Algorithms**

  - **Genetic Algorithm (GA):** Population-based, flexible, and widely used.
  - **CMA-ES:** Covariance Matrix Adaptation Evolution Strategy for efficient real-valued optimization.
  - **Bayesian Optimization (BO):** Perfect for minimizing computational load with fewer model runs.

- ⚙️ **Automated Workflow**

  - Manages parameter selection, model runs, and results aggregation in a single loop.

- ⚡ **Parallel Processing**

  - Speeds up computation by distributing tasks across CPU cores.

- 🧠 **Smart Caching**

  - Avoids redundant simulations by reusing previously evaluated parameter sets.

- 🛡️ **Robust Error Handling**

  - Includes timeouts and graceful recovery from failed model runs.

---

## 💻 System Requirements

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

## 🔧 Installation & Setup

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/AutoCal-Eco.git
cd AutoCal-Eco
```

### 2. Configure `setup.R`

This step links your EwE executable and sets up the calibration command:

- Navigate to `R/setup_template.R`
- Copy and rename it as `setup.R` in the same directory:

```bash
cp R/setup_template.R R/setup.R
```

- Open `setup.R` and edit the following variables:
  - `cmd_base` — Base command to run EwE
  - `fn.runEwE` — Function to launch the model

> ⚠️ Your `setup.R` is `.gitignore`d to protect local paths and credentials.

---

## ⚡ Quick Start

### 1. Open the Run Script

In RStudio or your preferred IDE, open:

```
analysis/run_analysis.R
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
  config = ga_config
)
```

### 4. Run the Analysis

Execute the script to begin optimization:

```r
source("analysis/run_analysis.R")
```

> ⏳ Results are stored in a timestamped directory like:\
> `GA_Calibration_Run/GA_Run_20250804_130900/`

---

## 📁 Project Structure

```
AutoCal-Eco/
├── R/                  # Core framework scripts
│   ├── calibration_framework.R
│   └── setup_template.R
│
├── data/               # Input data for case studies
│   ├── Gulf_of_Mexico/
│   └── West_Florida_Shelf/
│
├── analysis/           # Script to launch analyses
│   └── run_analysis.R
│
└── results/            # Output folder (ignored by Git)
```

---

## 📖 Citation

If you use **AutoCal-Eco** in your research, please cite:

> **Author A, Author B, Author C**. (*Year*). *AutoCal-Eco: An R Framework for Automated Calibration of Ecospace Models Using Modern Optimization Techniques*. **SoftwareX**, Volume, Pages. [https://doi.org/your-doi-here](https://doi.org/your-doi-here)

---

## 📜 License

This project is released under the **MIT License**. See `LICENSE.txt` for details.

---

## 🤝 Contact

For questions, issues, or feature requests:\
👉 [Open an Issue](https://github.com/your-username/AutoCal-Eco/issues)

