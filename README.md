# CliMatch

Computer code underlying the following publication:
Duarte S. Viana & Jonathan M. Chase. Ecological traits underlying interspecific variation in climate matching of birds. Global Ecology and Biogeography, https://doi.org/10.1111/geb.13480

Author: Duarte S. Viana

The files contained in this repository provide the R code for simulating and analysing the data corresponding to the mentioned paper. Below we provide a brief description of the workflow and respective files. Note that most of the code for data downloading and management was taken from Harris et al. (2018; PeerJ 6, e4278).

"Data_download.R": download the bird and climate data (it uses the R functions in the files "Data_management.R" and "get_prism_data.R").

"Data_preparation.R": prepare and format data for the data analysis.

"Model_fitting_species-climate.R": fit species-climate models (BRT and GAM) with the functions in "BRT&GAM.R" to estimate climate matching.

"Trait_models.R": fit Bayesian multilevel models to estimate the effects of species traits on climate matching.

"Simulations_count_data.R": simulate the effect of mean abundance on the fitting performance of models of count data.
