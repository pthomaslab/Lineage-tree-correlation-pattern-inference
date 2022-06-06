README.txt

This repository contains code in three parts:
1. matlab_correlation_pipeline
2. bayesian_sampling.jl
3. X mathematica

1. Calculates correlations between related pairs of cells
2. Uses MCMC with adaptive Metropolis-within-Gibbs sampling to fit parameters of the two-dimensional inheritance matrix model.
3. Mathematica pipeline

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
matlab_correlation_pipeline

Contains three MATLAB files
* P1_pair_cells.m - pairs cell IDs according to family relation
* P2_pair_colums.m - using paired cell IDs, pairs the chosen measurement
* bootstrapping_single.m - calculated correlation coefficients and 95% boostrapped CI

-------
USAGE GUIDE

P1_pair_cells.m
* INPUT file must be a .txt file with column 1 giving the cell ID and column 2 the parent ID. Columns 3 onwards can be any desired measurement of the cell.
* This pairs the cell IDs together according to family relation.
* OUTPUT file is [celltype]_paired_cell_ids.mat which is used in the next MATLAB file, P2_pair_columns.m

P2_pair_columns.m
* Loads [celltype]_paired_cell_ids.mat and pairs a chosen measurement (specified by the column of the original .txt file).
* Note: when prompted 'What is it that you are pairing?', enter something that relates to the chosen measurement (eg. 'idt' for interdivision time).
* OUTPUT file is [celltype]_[choice]_pairs.mat whch is used in the next MATLAB file, bootstrapping_single.m

bootstrapping_single.m
* Loads [celltype]_[choice]_pairs.mat; computes family correlation coefficients and 95% bootstrapped confidence intervals.
* Note: when prompted 'What did you choose to pair?' this must match exactly the choice you specified in P2_pair_columns.m
* OUTPUT file is [celltype]_[choice]_corrs.txt which can be used in the Julia code (bayesian_sampling.jl).

-------
FORMAT OF [celltype]_[choice]_corrs.txt
* Format must be like this for use in bayesian_sampling.jl, however, only rows 1-6 are required.

column 1: mean value
column 2: distance from mean to lower bound of 95% CI
column 3: distance from mean to upper bound of 95% CI
column 4: standard deviation

row 1: chosen measurement
row 2: variance of chosen measurement
row 3: mother-daughter correlation
row 4: sister-sister correlation
row 5: cousin-cousin correlation
row 6: grandmother-granddaughter correlation
row 7: great-grandmother great-granddaughter correlation
row 8: greatx2-grandmother greatx2-granddaughter correlation
row 9: greatx3-grandmother greatx3-granddaughter correlation
row 10: greatx4-grandmother greatx4-granddaughter correlation
row 11: greatx5-grandmother greatx5-granddaughter correlation
row 12: greatx6-grandmother greatx6-granddaughter correlation

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bayesian_sampling.jl

Takes .txt file [celltype]_[choice]_corrs.txt which is either generated with the MATLAB pipeline or from other sources and fits the two-dimensional inheritance matrix model to the correlation coefficients using MCMC sampling with adaptive Metropolis-Gibb.

-------
USAGE GUIDE

* Only input file required is [celltype]_[choice]_corrs.txt. It should be in the same folder as bayesian_sampling.jl
* Other input parameters are
  * Number of samples desired:
	Burn in added on will be set at 10% of this value.
  * Number of samples desired for thinned output.
	Code produces a thinned output file which will have this many samples.
* When asked 'What did you pair?' the input must be the same as used in the MATLAB pipeline in order to read the .txt file correctly.

This code gives a collection of output files:

* ..._sampling_params.dat - parameter sampling file.

* ..._inf.txt - produces inferred correlation coefficients and 95% CI. Format is as follows:
  column 1 - median value
  column 2: distance from mean to lower bound of 95% CI
  column 3: distance from mean to upper bound of 95% CI
  row 1: mother-daughter
  row 2: grandmother-granddaughter
  row 3: great-grandmother great-granddaughter
  row 4: greatx2-grandmother greatx2-granddaughter
  row 5: sister-sister
  row 6: cousin-cousin

* ..._behaviour_dist.txt - % of samples that fall into each behaviour
 in order: oscillator, alternator, aperiodic

* ..._mlp.txt - maximum likelihood parameter set; in order:
theta11, theta12, theta21, theta22, lambda1, lambda2, gamma12, delta11, delta12, delta22, log likelihood
  theta__ parameters give the parameters of the inheritance matrix
  lambda_ values are the variances of the noise terms
  gamma12 is the correlation between the two noise terms
  delta__ gives the correlations between the noise terms between sisters

* ..._mle.txt - correlations calculated from maximum likelihood parameter set, in order:
mother-daughter, grandmother-granddaughter, great-grandmother great-granddaughter, greatx2-grandmother greatx2-granddaughter, sister-sister, cousin-cousin

* ..._output.txt - full sampling output file, has all relevant outputs
* ..._output_thin.txt - thinned output file
  output files have the following format of columns:

## OUTPUT MATRIX FORMAT
# 1 - 11 | sparas
    # 1 - theta11
    # 2 - theta12
    # 3 - theta21
    # 4 - theta22
    # 5 - lambda1
    # 6 - lambda2
    # 7 - gamma12
    # 8 - delta11
    # 9 - delta12
    # 10 - delta22
    # 11 - log likelihood
# 12 - 17 | correlations
    # 12 - mother | 13 - grandmother | 14 - great-grandmother | 15 - great-great-grandmother
    # 16 - sister | 17 - cousin
# 18 | period of oscillator
# 19 | damping of oscillator
# 20 - 22 | behaviours
    # 20 - oscillator - complex eigenvalues
    # 21 - alternator - at least one negative eigenvalue, period 2
    # 22 - aperiodic - both positive eigenvalues, period inf
# 23 - 32 | inferred periods
# 33 real part eig 1
# 34 complex part eig 1
# 35 real part eig 2
# 36 complex part eig 2
# 37 correlation between 2 factors
# 38 m1d1 mother-daughter factor 1 to factor 1 correlation
# 39 m1d2 mother-daughter factor 1 to factor 2 correlation
# 40 m2d1 mother-daughter factor 2 to factor 1 correlation
# 41 m2d2 mother-daughter factor 2 to factor 2 correlation
