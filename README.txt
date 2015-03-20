COMPLEX TRAIT ARCHITECTURE: THE PLEIOTROPIC MODEL REVISITED
NORTH & BEAUMONT 2015
PUBLISHED IN SCIENTIFIC REPORTS:

North, T.-L. & Beaumont, M.A. Complex trait architecture: the pleiotropic model revisited. Sci. Rep. 5, 9351;(2015)

RELATED REFERENCES:
Eyre-Walker A. Genetic architecture of a complex trait and its implications for fitness and genome-wide 
association studies. Proceedings of the National Academy of Sciences 107, 1752-1756 (2010).

Galassi M., et al. GNU Scientific Library Reference Manual, 3rd edn. (Network Theory Ltd, 2009).

INFORMATION ON SCRIPTS

1/ complex_tld.c

This is the main c code to run the simulation of evolution. The script uses functions from the 
GNU Scientific Library (GSL), which you need to have installed.
 
It requires a parameter text file ("complex_params_tld.txt") with the following information (one value per line):
Psize 
Stoptime 
Gentime 
Murate 
Betashape 
MeanS
Ne
Tau
Trait_is_fit_switch
Ep_mu
Ep_std
Purge_interval

Explanations of each parameter are provided in the script comments. 
Note that Psize should always equal Ne, and that Trait_is_fit_switch should always equal 0


Notes/
-complex_tlc.c and complex_tld_new.c both use the error function
printerr, as defined in myutil.c
-complex_tlc.c and complex_tld_new.c were run 100 times for each model in North and Beaumont (2015)
-complex_tlc.c and complex_tld_new.c output various text files along the course and at the end of the simulation 

2/ complex_tld_new.c

This is an updated version of complex_tld.c to include the new formula for the relationship
between trait and fitness effects

It requires a parameter file with the following information:

Psize
Stoptime
Gentime
Murate
Betashape
MeanS
Ne
A
B
Ep_mu
Ep_std
Purge_interval

Explanations of each parameter are provided in the script comments. 
Note that Psize should always equal Ne


3/ complex_var_plots.R
Aside from some QC checks, complex_var_plots.R analyses the data from the end of the simulation
to output a text file var_plot_x.txt (where x=run) which catalogues the variance contribution of
each mutational frequency class in this run

It requires 8 arguments:
trait_effects.txt
row count of trait_effects.txt
column count of trait_effects.txt
mut_stats_12345.txt
row count of mut_stats_12345.txt
column count of mut_stats_12345.txt
N
Run number

4/ complex_var_plots_iif.R
This script creates the plots shown in the supplement of North and Beaumont (2015)
It is interactive and requires copying and pasting different sections with updates for
input and output filenames
At the beginning you should be in the directory containing all of the var_plot_x.txt files

5/ complex_kurtosis_d.R
This script creates the kurtosis plot in North and Beaumont (2015)
It begins by reading in some of the text files written by complex_var_plots_iif.R.
Like complex_var_plots_iif.R, this script needs updating with folder names etc.
Note that models 1-5 were models a-e, and models 7-9 were models f-h
The script assumes that you have created a folder containing the subfolders
model1, model2, model3, model4, model5, model7, model8 and model9
Each subfolder contains 100 text files 1_run.txt to 100_run.txt which have the
10,000 trait values for each individual at the end of the simulation




