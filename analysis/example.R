# ===================================================================
# Example Script for Running an Automated Calibration
# ===================================================================

rm(list = ls())
# --- Step 1: Source the necessary setup and framework files ---

# This file must define cmd_base, fn.runEwE, and other required variables
source("./setup.R")

# This file contains the main run_calibration function
source("./calibration_pointer.R")

run_calibration(
  method = "GA",
  predprey_pairs = predprey,
  output_base = "GA_Calibration_Run")



