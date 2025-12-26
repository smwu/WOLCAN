#!/bin/bash
#SBATCH -J WOLCAN_wts_scen_24_model11_        # Job name for the array
#SBATCH -o WOLCAN_wts_scen_24_model11_%A.out  # Shared standard output with job ID
#SBATCH -e WOLCAN_wts_scen_24_model11_%A.err  # Shared standard error with job ID
#SBATCH -p shared      # Partition to submit to
#SBATCH -c 8	       # Number of cores (for parallelization)
#SBATCH -N 1           # Number of nodes
#SBATCH -t 0-2:00:00  # Runtime (D-HH:MM:SS)
#SBATCH --mem=10000     # Memory request
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load R/4.2.2-fasrc01 gcc/10.2.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
data_scen=24
model_scen_num=11
Rscript /n/netscratch/stephenson_lab/Lab/stephwu18/WOLCAN/Summary_Code/weights_sims_par_samps.R ${data_scen} ${model_scen_num} ${SLURM_ARRAY_TASK_ID} 
