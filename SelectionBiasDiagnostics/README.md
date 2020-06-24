# A simulation study of diagnostics for selection bias

### Current Suggested Citation

Boonstra, Philip S, Little, Roderick JA, West, Brady T, Andridge, Rebecca R, and Alvarado-Leiton, Fernanda, "A simulation study of diagnostics for bias in non-probability samples" (2019). *The University of Michigan Department of Biostatistics Working Paper Series*. Working Paper 125.
https://biostats.bepress.com/umichbiostat/paper125

See also:

Roderick J A Little, Brady T West, Philip S Boonstra, Jingwei Hu, 
"Measures of the Degree of Departure from Ignorable Sample Selection" (2019), 
*Journal of Survey Statistics and Methodology*. 
https://doi.org/10.1093/jssam/smz023

Andridge, Rebecca R, West, Brady T, Little, Roderick JA, Boonstra, Philip S, and Alvarado‐Leiton, Fernanda, "Indices of non‐ignorable selection bias for proportions estimated from non‐probability samples" (2020) *Journal of the Royal Statistical Society Series C, Royal Statistical Society*, 68(5), 1465-1483.

## Executive Summary

This manuscript is one element of a multipart project studying estimators of non-ignorable
selection bias. By definition, bias due to non-ignorable selection is not estimable
without making additional assumptions. Little, et al. (2019) outline one such set 
of assumptions based on normal pattern-mixture models that yield estimates of
selection bias. These estimators are the 'standardized measures of unadjusted
bias (SMUB)' and 'standardized measures of adjusted bias (SMAB)'. 

In this manuscript, we describe the results of a comprehensive simulation study 
of the SMUB and SMAB estimators as well as other potential diagnostics for selection 
bias in non-probability sampling mechanisms. Assorted selection mechanisms of both the
ignorable and non-ignorable type are considered

## Further details

To recreate the simulation study in Boonstra et al., run the script <samp>run_sims.R</samp> 837 times, with the variable `array_id` ranging from 
1 to 837. Then, run the script <samp>make_plots.R</samp> to create the
tables and figures. 

The specific R scripts and other files are as follows. 

### <samp>R</samp> files

<samp>example.R</samp> demonstrates how to run a simulation for a single 
scenario

<samp>run_sims.R</samp> is the script we used to do the entire simulation study. 
Specifically, by setting `my.work.computer = F` on line 7 of this script, this 
can be called multiple times by a SLURM batch script (<samp>run_sims.txt</samp>; see below). The results are written to a csv file in a subfolder called 'out'. Setting
`my.work.computer = T` allows for running a single simulation on your local 
machine.

<samp>make_populations.R</samp> contains a function of the same name to simulate
an arbitrary number of populations (and non-probability samples from each 
population). 

<samp>construct_statistics.R</samp> contains a function of the same name, which
takes as input the output from a call to <samp>make_populations</samp> and returns
all of the diagnostics and error measures considered in the manuscript. 

<samp>FMI_v.1.3.R</samp> was provided by Rafael Nishimura and computes the
fraction of missing information diagnostic. The function 
<samp>construct_statistics.R</samp> depends on this function. 

<samp>msb_functions.R</samp>, <samp>nisb_functions.R</samp> are contained in the
parent directory IndicesOfNISB

<samp>make_plots.R</samp> creates the tables and figures for the manuscript. 

### <samp>txt</samp> files

<samp>run_sims.txt</samp> is a SLURM batch script. To rerun the simulation study
reported in Boonstra, et al., ensure that line 11 of the script is set to
<samp>#SBATCH --array=1-837</samp>. Then, from the command line, run 
`sbatch run_sims.txt`, which will call <samp>run_sims.R</samp> 837 
times, with the variable `array_id` ranging from 1 to 837. Each of the 279 unique scenarios is repeated three times with different random seeds per `array_id` and 667 iterations per `array_id`, for a total of 2001 iterations per scenario. 

## Acknowledgments 

This work was supported by the National Institutes of Health (1R21HD090366-01A1)