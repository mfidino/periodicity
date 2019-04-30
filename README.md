# Using Fourier series to predict periodic patterns in dynamic occupancy models

This is model code from:

Fidino, M., and S. B. Magle. (2017). Using Fourier series to predict periodic patterns in dynamic occupancy models. Ecosphere

It can be used to estimate periodic signals within the linear predictor of a dynamic occupancy model. 

All of the files in Data S1 must be within your working directory for the analysis to work.  Furthermore, the analyses have been set up to be done in parallel and uses all but two of the cores in a computer. Therefore, if you have two or less cores on your computer you will need to adjust any function in Fidino_2017_utility_functions.R that uses the detectCores function. This includes zest_posterior, make_pred_PTM_pulse, make_pred_PTM_boom and fit_models. Furthermore the n.chains arguments in Fidino_2017_total_variabiltity.R associated to the run.jags functions would need to be changed.

Once all of the files have been copied, the models can be fit through Fidino_2017_analysis_script.R

# Analysis scripts

There are 4 scripts used for analysis and model summary.

- **Fidino_2017_utility_functions.R:** This script contains utility functions to fit homogeneous time, stochastic time, and periodic time dynamics occupancy models. This script needs to be sourced before using any of the other scripts. Functions in this file are explained within the script itself. 
- **Fidino_2017_analysis_script.R:** This script uses the functions in Fidino_2017_periodic_utility_functions.R to fit periodic time, stochastic time, and homogeneous time dynamic occupancy models to the 9 seasons of Chicago camera trap data for coyote, red fox, striped skunk, raccoon, and Virginia opossum data.
- **Fidino_2017_summarize_and_plot.R:** Uses the saved model outputs from the analysis script and the functions in Fidino_2017_utility_functions.R to calculate summary statistics for plotting (i.e., this script summarizes the MCMC outputs to generate Figure 1 in the manuscript).
- **Fidino_2017_total_temporal_variability.R:** This script uses the functions in Fidino_2017_periodic_utility_functions.R to fit a periodic time model to each species that also includes a random temporal component. The primary variable of interest from these models is the standard deviation associated to the random temporal component, as it can be compared to the standard deviation from a model that only contains a random temporal component (i.e., no periodic element).

# Models

There are 6 JAGS dynamic occupancy models used throughout this analysis.
- **Fidino_2017_periodic_time_model_pulse.R:** The periodic time model used for coyote, striped skunk, and red fox. Used in Fidino_2017_analysis_script.R
- **Fidino_2017_periodic_time_model_boom_bust.R:** The periodic time model used for raccoon and Virginia opossum. Used in Fidino_2017_analysis_script.R
- **Fidino_2017_stochastic_time_model.R:** Used on all species and includes a random temporal effect on colonization rates. Used in Fidino_2017_analysis_script.R
- **Fidino_2017_homogeneous_time_model.R:** Used on all species and does not have any temporally varying parameters on colonization rates. Used in Fidino_2017_analysis_script.R
- **Fidino_2017_periodic_time_full_model_pulse.R:** The periodic time model used for coyote, striped skunk, and red fox that also includes a random temporal component. This requires δ to be supplied as data. Used in Fidino_2017_total_temporal_variability.R
- **Fidino_2017_periodic_time_full_model_boom_bust.R:** The periodic time model used for raccoon and Virginia opossum that also includes a random temporal component. This requires δ to be supplied as data. Used in Fidino_2017_total_temporal_variability.R

# Data
The data used in this analysis comes from a large-scale long-term camera trapping survey in Chicago, Illinois.
- **Fidino_2017_community_incidence_matrix_sp11_sp13.txt:** This includes data on whether or not coyote, red fox, striped skunk, raccoon, and Virginia opossum were observed at the 95 camera trapping sites each season between spring 2010 to spring 2013. If they were detected the cell takes a value of 1, if they were not detected it takes a 0, and if the site was not sampled it takes an NA. It can be converted to a 3-dimensional array using the df_2_array function in Fidino_2017_periodic_utility_functions.R and then used to calculate beta priors and supplied as initial values to model. This is a species by site by season array after being converted. The species order is coyote, red fox, striped skunk, raccoon, and Virginia opossum.
- **Fidino_2017_chicago_detection_data_sp11_sp13.txt:** This includes data on the number of days each species was detected and is supplied as data to the JAGS model so that each species detection probability can be calculated. It can be converted to a species by site by season array using the df_2_array function in Fidino_2017_perioidic_utility_functions.R. If a site was not sampled at a particular season an NA is reported.
- **Fidino_2017_Chicago_days_camera_active_sp11_sp13.txt:** This includes data on the number of days a camera trap was active each site and season. Used with the detection data to calculate detection probabilities, and is a site by season matrix. If a site was not sampled a zero is reported.
- **Fidino_2017_URB_covariate.txt:** The URB covariate for each site, in the same order as all of the other data.
