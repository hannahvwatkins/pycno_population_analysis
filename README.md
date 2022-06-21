## *Pycnopodia helianthoides* population trend analysis

This repository contains the code, model outputs and figures for the population trend analysis for *Pycnopodia helianthoides* in Canada for our upcoming COSEWIC report. Data are stored outside the repository, as most sources require permission for use - please contact me if you would like to be put in contact with the data providers. 

*Data*

Data sources include those used in the IUCN report for *P. helianthoides* as well as from Fisheries and Oceans Canada and Dr. Isabelle Côté.

*Code*

The code folder contains files labelled with the dataset they pertain to, as well as the model type (if there are multiple models). All **.R** files contain cleaning code to format the data in an appropriate way for analysis, and code to execute the Stan models, validate these models, and plot their output. The only exception is **custom_functions.R**, which is a script with custom functions designed for use with these analyses. These include those created specficially for the REEF analyses by Dan Greenberg, and my own functions for validating Stan models in R - this must be sourced at the beginning of each R script to be able to access these functions. All **.stan** files contain the hierarchical Bayesian state-space models. All code is heavily annotated to allow for greater readibility, but if you have any questions or concerns, please feel free to get in touch. Some notes:

-All abundance-related models are based on an index of abundance (count per standardized dive time) rather than a specific density, as most surveys in British Columbia did not explicitly measure the area surveyed. This means these results are still useful as a measure of population decline, but are less useful in estimating total population size.
-REEF survey data have been analyzed separately from other abundance data, as these data are recorded as scores rather than actual abundances. They have been run using an ordered logistic regression model, and presented using a conservative method of converting scores to abundances (i.e., assuming low scores in each category). The relative population index does not change when using less conservative methods, and the estimated population declines stay relatively constant.
-Only data sources that spanned time periods both before and after the crash were included in the sighting frequency model. The exclusion of sources that only surveyed before or after the crash reduced the size of the dataset by <5%, which allowed for better model convergence without sacrificing accuracy. The data sources excluded from this component of the analysis will still be discussed in the report, as a few point to areas with a relatively high potential for recovery.

*Model Output*

The model output folder contains the **.rds** files for each model output, so there is no need to rerun Stan each time you need to extract information from it. The outputs are called in each of the relevant R scripts. The output for the sighting frequency model is currently too large to upload even with Github large file storage, so please email me if you'd like it!

*Figures*

The figures folder contains all of the figures pertaining to each analysis. For the time series figures, both the observed annual trends (in blue) and the underlying state (in green) are shown. The underlying state is estimated by partitioning inter-year variation into process noise (i.e., actual changes in the underlying population size) and observation error (i.e., error inherent in the sampling methods), as described by Dennis et al, 2006 (https://doi.org/10.1890/0012-9615(2006)76[323:EDDPNA]2.0.CO;2). The number of surveys conducted each year are shown along the x-axis.
