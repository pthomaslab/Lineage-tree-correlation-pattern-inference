# Lineage tree correlation pattern inference

Scripts written and used as part of the pre-print *[link to be added].*

This repository contains the scripts required to simulate lineage tree data, compute correlations between pairs of cells on a lineage tree, fit a two-dimensional inheritance matrix model to these correlations, and analyse the model fit and output.

This repository has four parts:

1. `sim_using_MLP`
2. `matlab_correlation_pipeline`
3. `bayesian_sampling`
4. `imm_output_analysis`

`sim_using_MLP` can be used to simulate lineage tree data. The following three parts can be used to analyse either real or simulated data. These should be used in order for the best experience.

### General usage

To use these files, download the .zip file and extract. The scripts were written in *MATLAB R2021a*,  *Julia 1.6.2* and *Wolfram Mathematica 12*. They have not been tested for compatability with other versions. For `bayesian_sampling` and `imm_output_analysis`it is best to have your input files in the same folder as the script.

#### Note on parameters:

In `sim_using_MLP`, `bayesian_sampling` and `imm_output_analysis` the individual model parameter names correspond to the inheritance matrix $\theta$ and variance covariance matrices $\boldsymbol{S}_1$ and $\boldsymbol{S}_2$ as follows:

![params_image.png](https://raw.githubusercontent.com/fernhughes/Lineage-tree-correlation-pattern-inference/main/params_image.png)

---

### sim_using_MLP

Contains two files, `generate_ic.m` and `simulate_data2D_MLP.m` which simulate the two dimensional inheritance matrix model for a choice of model parameters.

#### Usage

**Before you begin**

* Download [Random trees - File Exchange - MATLAB Central](https://uk.mathworks.com/matlabcentral/fileexchange/2516-random-trees) and put the three files `branch.m`,`trimtreelayout.m` and `trimtreeplot.m` in the `sim_using_MLP` folder.

* **INPUT** required for these scripts is a .txt file that contains the chosen model parameters in seperate rows in the following order: theta11, theta12, theta21, theta22, lambda1, lambda2, gamma12, delta11, delta12, delta22.

###### Script specific usage

`generate_ic.m`

Generates initial conditions to start simulated trees from. **Must be run first.**

* You will first be prompted to input your *'mlp .txt file name'* which is the .txt file that contains your chosen parameters. This must contain your chosen model parameters in seperate rows in the following order: theta11, theta12, theta21, theta22, lambda1, lambda2, gamma12, delta11, delta12, delta22.

* *'Simulation name:'* choose a name for your simulation.

* *'Number of cells to sim'*- choose total number of cells you want to simulate.

* *'Mean'* is the mean of your desired simulated measurement which you specify here (usually mean interdivision time).

* **OUTPUT** is a .mat file `[simname]_sim_initial_conditions.mat` containing initial conditions which `simulate_data2D_MLP.m` will sample from when producing trees.

`simulate_data2D_MLP.m`

Simulates data using the two-dimensional inheritance matrix model.

* *'Simulation name:'* must match the simulation name `[simname]` you chose before.

* **OUTPUT** is `[simname]_sim_final_data.txt`which has column 1 cell ID, column 2 mother ID and column 3 measurement (usually interdivision time) in the correct format for correlation calculations and model inference. 

* Copy this file `[simname]_sim_final_data.txt` into `matlab_correlation_pipeline` to calculate correlation coefficients.

---

### matlab_correlation_pipeline

Contains three MATLAB files

* `P1_pair_cells.m` - pairs cell IDs according to family relation
* `P2_pair_colums.m` - using paired cell IDs, pairs the chosen measurement
* `bootstrapping_single.m` - calculates correlation coefficients and 95% boostrapped confidence intervals

#### Usage

`P1_pair_cells.m`

* **INPUT** file of lineage tree data must be a .txt file with column 1 giving the cell ID and column 2 the parent ID. Columns 3 onwards can be any desired measurement of the cell.
* When prompted *'Type the file name'* enter the name of the .txt file **without** the .txt extension.
* When prompted *'Type the cell type'* choose a name for your analysis. This will be `[celltype]`.
* This script pairs the cell IDs together according to family relation.
* **OUTPUT** file is `[celltype]_paired_cell_ids.mat` which is used in the next *MATLAB* file, `P2_pair_columns.m`

`P2_pair_columns.m`

* Loads `[celltype]_paired_cell_ids.mat` and pairs a chosen measurement (specified by the column of the original .txt file).

* When prompted *'Type cell name'* this must be `[celltype]`as specified previously.

* **Note:** when prompted *'What is it that you are pairing?'*, enter something that relates to the chosen measurement (eg. 'idt' for interdivision time). This will be `[choice]`.

* **Note:** when prompted *'Which column do you want to pair?'*, usually this is column **3** but you can choose a different column if you have it.

* **OUTPUT** file is `[celltype]_[choice]_pairs.mat` whch is used in the next *MATLAB* file, `bootstrapping_single.m`

`bootstrapping_single.m`

* Loads `[celltype]_[choice]_pairs.mat`; computes family correlation coefficients and 95% bootstrapped confidence intervals.
* Note: when prompted *'What did you choose to pair?'* this must match exactly `[choice]` you specified previously in `P2_pair_columns.m`
* **OUTPUT** file is `[celltype]_[choice]_corrs.txt` which can be used in the *Julia* code `bayesian_sampling.jl`.
* Copy this file `[celltype]_[choice]_corrs.txt` into the folder `bayesian_sampling`.

The format of this file `[celltype]_[choice]_corrs.txt` is shown below:

|                                                         | Mean value | From mean to lower bound | From mean to upper bound | Standard deviation |
| ------------------------------------------------------- | ---------- | ------------------------ | ------------------------ | ------------------ |
| Chosen measure                                          |            |                          |                          |                    |
| Variance of chosen measure                              |            |                          |                          |                    |
| Mother-daughter correlation                             |            |                          |                          |                    |
| Sister-sister correlation                               |            |                          |                          |                    |
| Cousin-cousin correlation                               |            |                          |                          |                    |
| Grandmother-granddaughter correlation                   |            |                          |                          |                    |
| Great-grandmother great-granddaughter correlation       |            |                          |                          |                    |
| Great x2-grandmother great x2-granddaughter correlation |            |                          |                          |                    |
| Great x3-grandmother great x3-granddaughter correlation |            |                          |                          |                    |
| Great x4-grandmother great x4-granddaughter correlation |            |                          |                          |                    |
| Great x5-grandmother great x5-granddaughter correlation |            |                          |                          |                    |
| Great x6-grandmother great x6-granddaughter correlation |            |                          |                          |                    |

------------

### bayesian_sampling.jl

Julia script, takes .txt file `[celltype]_[choice]_corrs.txt` which is either generated with the *MATLAB* pipeline `matlab_correlation_pipeline`, or from other sources, and fits the two-dimensional inheritance matrix model to the correlation coefficients using MCMC sampling with adaptive Metropolis-Gibbs.

#### Usage

* Load the project folder `bayesian_sampling` in Atom or your chosen Julia editor.

* **INPUT** file required is `[celltype]_[choice]_corrs.txt`. It should be in the same folder as `bayesian_sampling.jl`. This **MUST** follow the same format as specified in the table above, though only rows 1-6 are needed here.

* Other input parameters are
  
  * **Number of samples desired**: *the burn in added on will be set at 10% of this value.*
  * **Number of samples desired for thinned output**: *the script produces a thinned output file which will have this many samples.*

* When asked *'What did you pair?'* the input **must be the same** `[choice]` as used in the *MATLAB* pipeline in order to read the .txt file correctly.

* Move the files `..._output_thin.txt``..._behaviour_dist.txt` `..._mlp.txt` and `[celltype]_[choice]_corrs.txt` to the folder `imm_output_analysis`.

This sampling script gives a collection of output files:
`..._sampling_params.dat`- parameter sampling file.
`..._inf.txt `- produces inferred correlation coefficients and 95% CI. Format is as follows:

|                                                       | Median value | Median to lower bound | Median to Upper Bound |
| ----------------------------------------------------- | ------------ | --------------------- | --------------------- |
| Mother-daughter correlation                           |              |                       |                       |
| Grandmother-granddaughter correlation                 |              |                       |                       |
| Great-grandmother great-granddaughter correlation     |              |                       |                       |
| Greatx2-grandmother greatx2-granddaughter correlation |              |                       |                       |
| Sister-sister correlation                             |              |                       |                       |
| Cousin-cousin correlation                             |              |                       |                       |

`..._behaviour_dist.txt` - percentage of samples that fall into each behaviour, in order: *oscillator, alternator, aperiodic*
`..._mlp.txt` - maximum likelihood parameter set; in order:
`theta11, theta12, theta21, theta22, lambda1, lambda2, gamma12, delta11, delta12, delta22, log likelihood`

* `theta__` parameters give the parameters of the inheritance matrix
* `lambda_ `values are the variances of the noise terms
* `gamma12` is the correlation between the two noise terms
* `delta__` gives the correlations between the noise terms between sisters
  `..._mle.txt` - correlations calculated from maximum likelihood parameter set, in order:
  mother-daughter, grandmother-granddaughter, great-grandmother great-granddaughter, greatx2-grandmother greatx2-granddaughter, sister-sister, cousin-cousin
  `..._output.txt` - full sampling output file, has all relevant outputs
  `..._output_thin.txt` - thinned output file

The `output` files have the following format of columns:

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

##### 18 - 41 : Model output

18\. period of oscillator
19\. damping of oscillator

###### 20 - 22 : Behaviours

20\. Oscillator - complex eigenvalues
21\. Alternator - at least one negative eigenvalue, period 2
22\. Aperiodic - both positive eigenvalues, period infinite

###### 23 - 32: Underlying periods

33\. real part eigenvalue1
34\. complex part eigenvalue 1
35\. real part eigenvalue 2
36\. complex part eigenvalue 2
37\. correlation between 2 factors
38\. m1d1: mother-daughter factor 1 to factor 1 correlation
39\. m1d2: mother-daughter factor 1 to factor 2 correlation
40\. m2d1: mother-daughter factor 2 to factor 1 correlation
41\. m2d2: mother-daughter factor 2 to factor 2 correlation

---

### imm_output_analysis

Contains two files, `inheritance_matrix_model_setup.mx` which contains equations for all required correlations in the two-dimensional model. **You do not need to open this.** `imm_output_analysis.nb` loads this file and the output of the MCMC inference and plots the output as seen in the manuscript.

#### Usage

`imm_output_analysis.nb`

Requires the following input files (which should be obtained through `bayesian_sampling.jl`:

* `..._output_thin.txt`

* `..._behaviour_dist.txt`

* `..._mlp.txt`_

and

* `[celltype]_[choice]_corrs.txt`which is obtained through the *MATLAB* pipeline.

You will be asked **four prompts**:

1. *"Name of dataset"* - this is the `[celltype]`and must be exactly as in the input file name.

2. *"What did you choose to measure"* - this is `[choice]` and again must be exactly as in the input file name.

3. *"Number of steps"* - the number of steps specified in the sampling in `bayesian_sampling.jl`.

4. *"Dataset title"* - free form, this will be the title on your plots.

Run the whole notebook to obtain behaviour distribution, correlation visualisations, generalised tree correlation function and a histogram of underlying oscillatory periods (if oscillator behaviour is given.)
