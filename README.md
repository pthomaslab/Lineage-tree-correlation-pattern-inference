#Lineage tree correlations
This repository contains the code required to compute correlations between pairs of cells on a lineage tree, fit a two-dimensional inheritance matrix model to these correlations, and analyse the model fit and output.

This has three parts:
1. matlab_correlation_pipeline
2. bayesian_sampling.jl
3. MATHEMATICA PART

These should be used in order for the best experience.
###General usage
To edit; i.e. what to download, how to download and where to put files etc

------------

### matlab_correlation_pipeline
Contains three MATLAB files
* `P1_pair_cells.m` - pairs cell IDs according to family relation
* `P2_pair_colums.m` - using paired cell IDs, pairs the chosen measurement
* `bootstrapping_single.m` - calculates correlation coefficients and 95% boostrapped confidence intervals

####Usage

`P1_pair_cells.m`
* **INPUT** file of lineage tree data must be a .txt file with column 1 giving the cell ID and column 2 the parent ID. Columns 3 onwards can be any desired measurement of the cell.
* This script pairs the cell IDs together according to family relation.
* **OUTPUT** file is `[celltype]_paired_cell_ids.mat` which is used in the next *MATLAB* file, `P2_pair_columns.m`

`P2_pair_columns.m`
* Loads `[celltype]_paired_cell_ids.mat` and pairs a chosen measurement (specified by the column of the original .txt file).
* **Note: ** when prompted *'What is it that you are pairing?'*, enter something that relates to the chosen measurement (eg. 'idt' for interdivision time).
* **OUTPUT** file is `[celltype]_[choice]_pairs.mat` whch is used in the next *MATLAB* file, `bootstrapping_single.m`

`bootstrapping_single.m`
* Loads `[celltype]_[choice]_pairs.mat`; computes family correlation coefficients and 95% bootstrapped confidence intervals.
* Note: when prompted *'What did you choose to pair?'* this must match exactly the choice you specified previously in `P2_pair_columns.m`
* **OUTPUT** file is `[celltype]_[choice]_corrs.txt` which can be used in the *Julia* code `bayesian_sampling.jl`.

The format of this file `[celltype]_[choice]_corrs.txt` is shown below:

||Mean value|From mean to lower bound|From mean to upper bound|Standard deviation|
| ------------ | ------------ | ------------ | ------------ | ------------ |
|Chosen measure|   |   |   |   |
|Variance of chosen measure|   |   |   |   |
|Mother-daughter correlation|   |   |   |   |
|Sister-sister correlation|   |   |   |   |
|Cousin-cousin correlation|   |   |   |   |
|Grandmother-granddaughter correlation|   |   |   |   |
|Great-grandmother great-granddaughter correlation|   |   |   |   |
|Greatx2-grandmother greatx2-granddaughter correlation|   |   |   |   |
|Greatx3-grandmother greatx3-granddaughter correlation|   |   |   |   |
|Greatx4-grandmother greatx4-granddaughter correlation|   |   |   |   |
|Greatx5-grandmother greatx4-granddaughter correlation|   |   |   |   |
|Greatx6-grandmother greatx4-granddaughter correlation|   |   |   |   |


------------
###bayesian_sampling.jl
Julia script, takes .txt file `[celltype]_[choice]_corrs.txt` which is either generated with the *MATLAB* pipeline `matlab_correlation_pipeline`, or from other sources, and fits the two-dimensional inheritance matrix model to the correlation coefficients using MCMC sampling with adaptive Metropolis-Gibbs.

####Usage
* **INPUT** file required is `[celltype]_[choice]_corrs.txt`. It should be in the same folder as `bayesian_sampling.jl`. This MUST follow the same format as specified in the table above, though only rows 1-6 are needed here.
* Other input parameters are
  * Number of samples desired: *the burn in added on will be set at 10% of this value.*
  * Number of samples desired for thinned output: *the script produces a thinned output file which will have this many samples.*
* When asked *'What did you pair?' * the input **must be the same ** `[choice]` as used in the *MATLAB* pipeline in order to read the .txt file correctly.

This code gives a collection of output files:
`..._sampling_params.dat`- parameter sampling file.
`..._inf.txt `- produces inferred correlation coefficients and 95% CI. Format is as follows:

||Median value|Median to lower bound|Median to Upper Bound|
|------------| ------------ | ------------ | ------------ |
|Mother-daughter correlation|   |   |   |
|Grandmother-granddaughter correlation|   |   |   |
|Great-grandmother great-granddaughter correlation|   |   |   |
|Greatx2-grandmother greatx2-granddaughter correlation|   |   |   |
|Sister-sister correlation|   |   |   |
|Cousin-cousin correlation|   |   |   |

`..._behaviour_dist.txt` - percentage of samples that fall into each behaviour, in order: *oscillator, alternator, aperiodic*
`..._mlp.txt` - maximum likelihood parameter set; in order:
`theta11, theta12, theta21, theta22, lambda1, lambda2, gamma12, delta11, delta12, delta22, log likelihood`
 * `theta__` parameters give the parameters of the inheritance matrix
 * `lambda_ `values are the variances of the noise terms
 *  `gamma12` is the correlation between the two noise terms
 * `delta__` gives the correlations between the noise terms between sisters
`..._mle.txt` - correlations calculated from maximum likelihood parameter set, in order:
mother-daughter, grandmother-granddaughter, great-grandmother great-granddaughter, greatx2-grandmother greatx2-granddaughter, sister-sister, cousin-cousin
`..._output.txt` - full sampling output file, has all relevant outputs
`..._output_thin.txt` - thinned output file

The `output` files have the following format of columns:

## OUTPUT MATRIX FORMAT
##### 1 - 11: Parameters
1\. theta11
2\. theta12
3\. theta21
4\. theta22
5\. lambda1
6\. lambda2
7\. gamma12
8\. delta11
9\. delta12
10\. delta22
11\. log likelihood
##### 12 - 17: Correlations
12\. mother
13\. grandmother
14\. great-grandmother
15\. great-great-grandmother
16\. sister
17\. cousin
#####18 - 41 : Model output
18\. period of oscillator
19\. damping of oscillator
###### 20 - 22 : Behaviours
20\. Oscillator - complex eigenvalues
21\. Alternator - at least one negative eigenvalue, period 2
22\. Aperiodic - both positive eigenvalues, period infinite
###### 23 - 32: Underlying periods
33\. real part eig 1
34\. complex part eig 1
35\. real part eig 2
36\. complex part eig 2
37\. correlation between 2 factors
38\. m1d1 mother-daughter factor 1 to factor 1 correlation
39\. m1d2 mother-daughter factor 1 to factor 2 correlation
40\. m2d1 mother-daughter factor 2 to factor 1 correlation
41\. m2d2 mother-daughter factor 2 to factor 2 correlation
