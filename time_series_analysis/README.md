## *Pycnopodia helianthoides* time series analysis

This folder contains the code and figures for the population trend analysis for *Pycnopodia helianthoides* in Canada for our upcoming COSEWIC report. Data that we are able to share are included in the root folder of this repo, but it is incomplete as several sources require permission for use - please contact me if you would like to be put in contact with the data providers. Note that model output files (i.e. the `.rds` files) are too large to upload even with Github's large file storage, but I am happy to email them upon request!

*Code*

The code folder contains files labelled with the dataset they pertain to, as well as the model type (if there are multiple models). All `.R` files with the prefix `cleaning_` contain code to process the raw data from each source for use in the analyses. All `.R` files with the prefix `model_` contain the code to execute, validate, and plot the output of the Stan models. The `custom_functions.R` script contains custom functions designed for use with these analyses. These include those created specficially for handling ordinal response output by Dan Greenberg, and functions for validating Stan models in R by myslef and Helen Yan - this must be sourced at the beginning of each R script to be able to access these functions. All `.stan` files contain the hierarchical Bayesian state-space models. All code is heavily annotated to allow for greater readibility, but if you have any questions or concerns, please feel free to get in touch.

*Figures*

The figures folder contains all of the figures pertaining to each analysis. For the time series figures, both the observed annual trends (in blue) and the underlying state (in green) are shown. The underlying state is estimated by partitioning inter-year variation into process noise (i.e., actual changes in the underlying population size) and observation error (i.e., error inherent in the sampling methods), as described by Dennis et al, 2006 (https://doi.org/10.1890/0012-9615(2006)76[323:EDDPNA]2.0.CO;2). The number of surveys conducted each year are shown along the x-axis.


Session Info
```
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252   
[3] LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
[5] LC_TIME=English_Canada.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.1.1      DHARMa_0.4.5         rstan_2.21.5        
 [4] StanHeaders_2.21.0-7 cmdstanr_0.5.2       forcats_0.5.1       
 [7] stringr_1.4.0        dplyr_1.0.7          purrr_0.3.4         
[10] readr_1.4.0          tidyr_1.1.3          tibble_3.1.2        
[13] ggplot2_3.3.5        tidyverse_1.3.1     

loaded via a namespace (and not attached):
 [1] httr_1.4.2           jsonlite_1.7.2       splines_4.1.0       
 [4] modelr_0.1.8         RcppParallel_5.1.4   assertthat_0.2.1    
 [7] posterior_1.1.0      distributional_0.2.2 stats4_4.1.0        
[10] tensorA_0.36.2       cellranger_1.1.0     pillar_1.6.2        
[13] backports_1.2.1      lattice_0.20-44      glue_1.6.2          
[16] checkmate_2.0.0      minqa_1.2.4          rvest_1.0.0         
[19] colorspace_2.0-2     Matrix_1.3-3         pkgconfig_2.0.3     
[22] broom_0.7.8          haven_2.4.1          scales_1.1.1        
[25] processx_3.5.2       lme4_1.1-27.1        generics_0.1.0      
[28] farver_2.1.0         ellipsis_0.3.2       withr_2.5.0         
[31] cli_3.3.0            magrittr_2.0.1       crayon_1.4.1        
[34] readxl_1.3.1         ps_1.6.0             fs_1.5.0            
[37] fansi_0.5.0          nlme_3.1-152         MASS_7.3-54         
[40] xml2_1.3.2           pkgbuild_1.2.0       tools_4.1.0         
[43] loo_2.4.1            prettyunits_1.1.1    hms_1.1.0           
[46] lifecycle_1.0.1      matrixStats_0.59.0   munsell_0.5.0       
[49] reprex_2.0.0         callr_3.7.0          compiler_4.1.0      
[52] rlang_1.0.4          nloptr_1.2.2.2       grid_4.1.0          
[55] rstudioapi_0.13      boot_1.3-28          gtable_0.3.0        
[58] codetools_0.2-18     inline_0.3.19        abind_1.4-5         
[61] DBI_1.1.1            R6_2.5.1             gridExtra_2.3       
[64] lubridate_1.7.10     knitr_1.33           utf8_1.2.1          
[67] stringi_1.6.1        parallel_4.1.0       Rcpp_1.0.7          
[70] vctrs_0.4.1          dbplyr_2.1.1         tidyselect_1.1.1    
[73] xfun_0.24           
```
