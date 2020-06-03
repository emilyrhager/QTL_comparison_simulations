#!/bin/bash
#SBATCH -n 1			# Number of cores
#SBATCH -N 1			# All cores on one machine
#SBATCH -t 0-06:00			# Time in D-HH:MM
#SBATCH -p hoekstra,shared	# partition
#SBATCH --mem=8000		# memory
#SBATCH -o /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/job_files/output/step2_maxvert_%j.out
#SBATCH -e /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/job_files/output/step2_maxvert_%j.err

module load Anaconda3/5.0.1-fasrc01
source activate rqtl_env

Rscript simulate_LOD_distrib_given_effect_sizes.R --cross_object /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/data/EPK_20191218_genoprob.rds --phenotype maxvert --folder /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/sims_ERH_phenotype_max_vert.resid.sacrum_chr_15 --effect_ratio 1 --geno_reversed --dom_model additive  

source deactivate
