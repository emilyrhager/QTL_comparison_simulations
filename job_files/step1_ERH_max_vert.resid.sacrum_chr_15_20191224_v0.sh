#!/bin/bash
#SBATCH -n 1			# Number of cores
#SBATCH -N 1			# All cores on one machine
#SBATCH -t 1-08:00    		# Time in D-HH:MM
#SBATCH -p hoekstra,shared	# partition
#SBATCH --mem=8000		# memory
#SBATCH -o /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/job_files/output/step1_max_vert.resid.sacrum_%j.out
#SBATCH -e /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/job_files/output/step1_max_vert.resid.sacrum_%j.err

module load Anaconda3/5.0.1-fasrc01
source activate rqtl_env

mkdir /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/sims_ERH_phenotype_max_vert.resid.sacrum_chr_15

Rscript simulate_effect_size_distrib_given_LOD.R --cross_object /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/data/SWxBK_rqtl_all_object_20191223.RDS --chromosome 15 --phenotype max_vert.resid.sacrum --folder /n/holylfs03/LABS/hoekstra_lab/Users/hager/QTL/compare_QTL_sims/sims_ERH_phenotype_max_vert.resid.sacrum_chr_15 --nsim 5000 --neffect 25 --min_effect 0 --max_effect 2.5 --dom_model additive --lod_window 0.1 --method ehk --model normal  

source deactivate
