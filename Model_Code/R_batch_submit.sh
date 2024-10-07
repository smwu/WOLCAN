#!/bin/bash
#SBATCH -J WOLCAN_scen4_        # Job name for the array
#SBATCH -o WOLCAN_scen4_%A.out  # Shared standard output with job ID
#SBATCH -e WOLCAN_scen4_%A.err  # Shared standard error with job ID
#SBATCH -p shared      # Partition to submit to
#SBATCH -c 8	       # Number of cores (for parallelization)
#SBATCH -N 1           # Number of nodes
#SBATCH -t 0-6:00:00  # Runtime (D-HH:MM:SS)
#SBATCH --mem=12000     # Memory request
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load R/4.2.2-fasrc01 gcc/10.2.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
scenario=4
Rscript /n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/Model_Code/model_sims.R ${scenario} ${SLURM_ARRAY_TASK_ID}
