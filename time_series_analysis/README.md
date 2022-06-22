## *Pycnopodia helianthoides* time series analysis

This folder contains the codeand figures for the population trend analysis for *Pycnopodia helianthoides* in Canada for our upcoming COSEWIC report. Data that we are able to share are included in the root folder of this repo, but it is incomplete as several sources require permission for use - please contact me if you would like to be put in contact with the data providers. Note that model output files (i.e. the `.rds` files) are too large to upload even with Github's large file storage, but I am happy to email them upon request!

*Code*

The code folder contains files labelled with the dataset they pertain to, as well as the model type (if there are multiple models). All `.R` files with the prefix `cleaning_` contain code to process the raw data from each source for use in the analyses. All `.R` files with the prefix `model_` contain the code to execute, validate, and plot the output of the Stan models. The `custom_functions.R` script contains custom functions designed for use with these analyses. These include those created specficially for handling ordinal response output by Dan Greenberg, and my own functions for validating Stan models in R - this must be sourced at the beginning of each R script to be able to access these functions. All `.stan` files contain the hierarchical Bayesian state-space models. All code is heavily annotated to allow for greater readibility, but if you have any questions or concerns, please feel free to get in touch.

*Figures*

The figures folder contains all of the figures pertaining to each analysis. For the time series figures, both the observed annual trends (in blue) and the underlying state (in green) are shown. The underlying state is estimated by partitioning inter-year variation into process noise (i.e., actual changes in the underlying population size) and observation error (i.e., error inherent in the sampling methods), as described by Dennis et al, 2006 (https://doi.org/10.1890/0012-9615(2006)76[323:EDDPNA]2.0.CO;2). The number of surveys conducted each year are shown along the x-axis.
