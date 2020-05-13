# QTL comparison simulations

Scripts to run simulations to test the likelihood that detected QTL are NOT shared with similar effect size in two experiments. 

This approach can be applied either for QTL maps generated for two comparable traits in the same mapping population, or for the same trait in two different mapping populations. 

It relies on QTL mapping data represented using the cross object from Karl Broman's rqtl package (TODO: LINK). 

# Overview

## The question
QTL mapping studies often generate association maps for multiple, related traits (e.g. the lengths of several whiskers) in one population, or for the same trait (e.g. total tail length) in two, parallel populations. An interesting question then becomes whether each detected locus is likely independent (i.e. has an effect on one trait/in one population, and a much smaller or no effect in the other) or whether it may be shared. 

Just using a binary comparison (detected above the significance threshold or not; confidence intervals overlap or not) is not sufficient for this purpose, because we know that QTL mapping experiments will likely not detect all causative alleles present in the cross founders, and will likely over-estimate effect sizes for loci that are detected (e.g. Beavis, 1994; Slate, 2013). With this in mind, we would like to have a quantitative estimate, given experimental data, for how likely it is, if a certain detected locus is shared (between two traits in one experiment or between the same trait in two experiments), that we would obtain a pattern of association as or more extreme than our actual results. Here we use a simulation-based approach to address this question.

## Main ideas
We approach this by breaking the question down into two steps: first, we estimate the distribution of true effect sizes that may have generated the detected LOD score for the experiment in which the QTL of interest was detected (experiment 1). Second, we estimate the probability of obtaining LOD scores less than or equal to the values found in experiment 2, if we assume the locus is shared with each possible effect size. We sum these estimates over the distribution of possible effect sizes to get the total probability that the locus is shared with similar effect size given the data. 

# Approach

# To use

# Author
Emily R. Hager, 2020
Hoekstra lab
Harvard University
