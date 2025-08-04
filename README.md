AutoCal-Eco: An R Framework for Automated Calibration of Ecospace Models

 <!--- Replace with your actual DOI once published --->

An R-based framework for the automated calibration of Ecopath with Ecosim (EwE) models using modern optimization techniques. This software provides a systematic, transparent, and repeatable procedure for fitting Ecospace models to temporal and spatial data.

Overview
The calibration of complex ecosystem models is a critical step for ensuring the credibility of their projections, yet it is often a labor-intensive and subjective process. This is particularly true for the widely used Ecopath with Ecosim (EwE) software, where the manual tuning of predator-prey vulnerability parameters presents a significant bottleneck.

AutoCal-Eco addresses this challenge by providing an open-source framework that automates the calibration of spatial-temporal Ecospace models. The software provides a flexible interface between the EwE command-line executable and a suite of powerful optimization algorithms, allowing researchers to efficiently find optimal parameter sets that align model outputs with observed data.

Features
Multiple Optimization Methods: Easily switch between three powerful optimization algorithms:

Genetic Algorithm (GA): A robust, population-based search method.

Covariance Matrix Adaptation Evolution Strategy (CMA-ES): A highly efficient evolutionary algorithm for real-valued problems.

Bayesian Optimization (BO): Ideal for very computationally expensive models, as it intelligently minimizes the number of required model runs.

Automated Workflow: Manages the entire calibration process, from parameter selection to model execution and results aggregation.

Parallel Processing: Significantly reduces runtime by distributing concurrent model evaluations across multiple CPU cores.

Efficient Caching: Avoids re-running simulations with identical parameters, saving valuable computation time.

Robust Error Handling: Includes timeouts and error catching to ensure that a single failed model run does not halt a long optimization process.

System Requirements
R: Version 4.0 or newer.

Ecopath with Ecosim: A working command-line version of the EwE executable.

R Packages: The following R packages are required and can be installed from CRAN:

GA

cmaes

rBayesianOptimization

dplyr

digest

readxl

doParallel

R.utils

You can install all required packages by running the following command in your R console:

install.packages(c("GA", "cmaes", "rBayesianOptimization", "dplyr", "digest", "readxl", "doParallel", "R.utils"))

Installation & Setup
Clone the Repository:
Clone this repository to your local machine using Git:

git clone https://github.com/your-username/AutoCal-Eco.git
cd AutoCal-Eco

Configure the setup.R File:
This is the most important setup step. The framework needs to know where your EwE executable is located and what the base commands for your model are.

In the root directory of the project, find the file R/setup_template.R.

Copy this file and rename the copy to setup.R in the root directory of the project.

Open the new setup.R and edit the variables cmd_base and fn.runEwE according to the instructions within that file.

Note: The setup.R file is listed in .gitignore and will not be tracked by Git, keeping your local file paths private.

Quick Start / Usage
The primary way to use the framework is through the run_analysis.R script located in the analysis/ directory. This script sources the core functions and provides a simple interface to start a calibration run.

Open analysis/run_analysis.R in RStudio or your preferred R editor.

Choose a Method:
Uncomment the run_calibration() call for the method you wish to use (GA, CMAES, or BO).

Customize (Optional):
You can override the default settings for any method by passing a config list. For example, to run the Genetic Algorithm with a smaller population size:

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

Run the Script:
Execute the run_analysis.R script from your R console. The framework will start the optimization process, printing progress to the console. Results will be saved in a new, time-stamped directory (e.g., GA_Calibration_Run/GA_Run_20250804_130900/).

File Structure
AutoCal-Eco/
│
├── R/                  # Core source code for the framework
│   ├── calibration_framework.R
│   └── setup_template.R
│
├── data/               # Input data for case studies
│   ├── Gulf_of_Mexico/
│   └── West_Florida_Shelf/
│
├── analysis/           # Scripts for running analyses
│   └── run_analysis.R
│
└── results/            # Output directory (ignored by Git)

Citation
If you use AutoCal-Eco in your research, please cite our publication:

{Author A, Author B, Author C. (Year). AutoCal-Eco: An R Framework for Automated Calibration of Ecospace Models Using Modern Optimization Techniques. SoftwareX, Volume, Pages. DOI:https://doi.org/your-doi-here}

License
This project is licensed under the MIT License - see the LICENSE.txt file for details.

Contact
For questions, bug reports, or suggestions, please open an issue on the GitHub repository.
