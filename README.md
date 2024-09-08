# WOLCAN

Weighted Overfitted Latent Class Analysis for Non-probability samples (WOLCAN)

This repo contains the following folders and files:

-   `Simulation_Code` folder: Code to generate simulated data

    - `simulate_data.R`: All data simulation scenarios and parameters.
    -   `simulate_data_functions.R`: Helper functions for data simulation.
    - `Old_Files` sub-folder: Old files to create toy simulated datasets.

-   `Model_Code` folder: Code to run the WOLCAN model.
  - `model_sims.R`: Pipeline for running the model using simulated datasets.
  - `model_functions.R`: Helper functions for running the model.
  - `R_batch_submit.sh`: Bash script for running the model simulations on the FASRC cluster.
  - `model_outcome.R`: Code for running outcome model sanity check simulations.
  - `wtd_logistic.stan`: Outcome Bayesian weighted logistic regression model in Stan.
  - `Old_Files` sub-folder: Old files for running the model on toy examples and for troubleshooting.

-   `Summary_Code`: Code to summarize the simulation study results.
  - `summary.R`: Pipeline for summarizing simulation study model results and creating graphics.
  - `summary_functions.R`: Helper functions for summarizing model results.
  - `weights_sim.R`: Code for running and summarizing simulations comparing pseudo-weight estimation models.

-   `Summary_Results`: Stored summary results for the various simulation scenarios.

-   `Application`: Data and code for running the data application with the PROSPECT and PRCS data.
