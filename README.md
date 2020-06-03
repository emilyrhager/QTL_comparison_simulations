---
title: "Simulation-based approach for comparing QTL mapping results"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#### Purpose

QTL mapping experiments identify regions of the genome where ancestry is associated with variation in some particular trait of interest, usually in a lab-generated hybrid population. 

One challenge for comparing QTL maps (e.g. maps for the same trait in two populations, or for two traits in the same population) is that for typical effect sizes and sample sizes, causative loci can often be missed due to sampling error, and in addition, effect sizes are likely to be overestimated for loci that do pass genome-wide significance thresholds. 

This means that if we detect a QTL for some trait in one experiment, but not in another, it is not trivial to assess whether that result is likely because a true, shared associated locus was missed in one of the experiments, or instead whether the causal locus is in fact not shared.

These scripts implement simulations to provide a quantitative estimate of the probability of obtaining given QTL mapping results if the causal locus were shared with similar effect size in two experiments. Each step can be run from the command line, given two cross objects in the format output by Karl Broman's rqtl package. 

#### Approach

The simulations are broken down into two steps, as follows.

For each QTL detected in one experiment, we would ultimately like to know: If the QTL did represent a shared allele with the same effect size in both experiments, what is the probability that we would have detected no linkage stronger than what we actually observed in the other experiment? We can write that this way:

$$
P_{results|shared} = P( l_2 \le L_2 \mid l_1 = L_1 )
$$
where $l_1$ represents the peak LOD score in experiment 1, $l_2$ is the highest LOD score in some region of the chromosome in experiment 2, and $L_1$ and $L_2$ are the actual experimentally determined values of $l_1$ and $l_2$. In other words, we define $P_{results|shared}$ as the probability that, given that we detected a QTL in experiment 1 with LOD score $L_1$, we would detect no LOD score greater than $L_2$ in experiment 2. 

To find $P_{results|shared}$ we decompose it into two parts: the probability that the true effect size $E$ is some value $e$, given that $l_1 = L_1$, and the probability that $l_2 \le L_2$ given that the true effect size $E$ is $e$:

$$
P( l_2 \le L_2 \mid l_1 = L_1 ) = \sum_e P( l_2 \le L_2 \mid E = e)  P(E = e \mid l_1 = L_1 )
$$
The two distributions $P( l_2 \le L_2 \mid E = e)$ and $P(E = e \mid l_1 = L_1 )$ we find by simulation. 

*Aside:* We need to decide which LOD score from experiment 2 to use for comparison. In practice, the code saves four options: the highest LOD score on the entire chromosome; the highest LOD score located within the Bayes' credible interval from experiment 1; and the LOD score at the two markers adjacent to the peak marker in experiment 1 (which can ultimately be used to calculate the LOD score at the closest marker). For typical cases the Bayes' credible interval definition is recommended. 


#### Running the simulations
The first step in running the simulations is to approximate $P(E = e \mid l_1 = L_1 )$ in experiment 1 (the experiment in which the QTL of interest was detected). This is done using the file **simulate_effect_size_distrib_given_LOD.R**. An example job submission script and a template for submission scripts for similar jobs are located in the job_files subfolder. 

The second step is to approximate $P( l_2 \le L_2 \mid E = e)$ for each effect size that could have generated the LOD score in experiment 1. This is done using the file **simulate_LOD_distrib_given_effect_size.R**. An example job submission script and a template for submission scripts for similar jobs are located in the job_files subfolder. 

The third step is to combine these distributions; because the file names and desired output may depend somewhat on the particular details of the experiment, this is done interactively; an example is shown in the file **gather_QTL_comparison_output_example.R**. 

#### Required packages
To run these scripts, you will need the following packages installed:

--argparse  
--qtl  
--dplyr  
--tidyr

#### Inputs

The scripts are designed to be run such that the simulations for each QTL peak reside in a unique folder. This makes the output of simulations from step 1 unambiguously available as input for step 2. 

**For step 1**, required inputs are:  
-- the cross object for experiment 1 (must have genotype probabilities already calculated, via qtl::calc.genoprob or some other method)  
-- the chromosome where the detected QTL is located (note: X chromosome not supported; assumes just one QTL per chromosome at this point.)  
-- the name of the phenotype (i.e. name of the column in cross_object$pheno)  
-- the folder to save output (should be a separate folder for each QTL peak)  

The default options run 100 simulations per effect size, for 5 effect sizes between 0 and 2.5x the effect size detected in the original experiment, and count LOD scores that fall within +/- 1 of the experimentally determined LOD score. These defaults are designed for ease of testing locally; in practice I use 5000 simulations per effect size for 25 effect sizes in the same range, and count LOD scores that fall within +/- 0.1 of the original LOD score. The default assumes the locus is purely additive; this can be changed to either take the dominance ratio from the experiment or to set d/a to a specific value (e.g. 1 for a purely dominant locus).

Other defaults that may be adjusted include:  
-setting the per-genotype variance to the total variance in the phenotype in the hybrid population; this can be changed so that the per-genotype variance varies with effect size to maintain total variance in the hybrids instead;  
-the method used for QTL mapping; default is extended Haley-Knott ('ehk')

**For step 2**, the required inputs are:  
-- the cross object for experiment 2 (must have genotype probabilities already calculated, via qtl::calc.genoprob or some other method)  
-- the name of the phenotype in experiment 2 (i.e. name of the column in cross_object$pheno)  
-- the folder to save output (should be the output folder used in step 1)  

The default options for the most part draw the parameters for the simulations to match the output from step 1. One exception is that to identify the correct input file, the default for step 2 is to assume the per-genotype variance was constant; this parameter must be changed for simulations where total variance was held constant instead. Similarly, the dominance ratio can be set if needed to distinguish among different runs of step_1 when choosing the input file. 

The default setting assumes the effect size in experiment 2 is equal to that in experiment 1; to perform a sensitivity analysis to see how the results vary if this assumption is relaxed, \--effect_ratio can be set to any value (e.g. setting effect_ratio to 0.5 will perform the simulations with the effect of the locus in experiment 2 equal to half that from experiment 1).

Additional options can be specified; use Rscript <filename.R> --help to see the documentation for these. 

#### Outputs

Each step outputs two csv files. The first (suffix _all_sims.csv) contains the results from each individual simulation run (e.g. the peak marker and LOD scores); to save space, the \--lose_sims flag can be set to prevent saving this file. 

The second (suffix _prob_distrib.csv for step 1 and _summary_sims.csv for step 2) contains summary information about the settings that were used to run the simulations and the output for each effect size e. 

The final step is to combine the information from these two files to calculate the final output probabilities, as shown in gather_QTL_comparison_output_example.R.

