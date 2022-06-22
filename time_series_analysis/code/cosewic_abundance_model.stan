//this is a custom prior for the cutpoints in the ordered logistic regression
//derived from Michael Betancourt's blog:
//https://betanalpha.github.io/assets/case_studies/ordinal_regression.html
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
  
}
//upload the data as defined in the list in pycno_abundance.R
data{
  int<lower=1> N_reef;//number of observations (REEF surveys)
  int<lower=1> N_ow;
  int<lower=1> N_iucn;
  int<lower=1> N_aba; 
  int<lower=1> N_multi;
  int y_reef[N_reef]; //abundance category for each survey (1-5)
  int y_ow[N_ow]; //abundance value for each survey (1-7)
  int y_iucn[N_iucn]; //counts 
  int y_aba[N_aba]; 
  int y_multi[N_multi]; 
  int<lower=0> N_site_reef; //number of sites
  int<lower=1,upper=N_site_reef> site_reef[N_reef]; // vector of site identities
  int<lower=0> N_site_ow; 
  int<lower=1,upper=N_site_ow> site_ow[N_ow]; 
  int<lower=0> N_site_iucn; 
  int<lower=1,upper=N_site_iucn> site_iucn[N_iucn]; 
  int<lower=0> N_site_aba; 
  int<lower=1,upper=N_site_aba> site_aba[N_aba]; 
  int<lower=0> N_site_multi; 
  int<lower=1,upper=N_site_multi> site_multi[N_multi]; 
  int<lower=0> N_source_iucn; //number of data sources
  int<lower=1,upper=N_source_iucn> source_iucn[N_iucn]; //vector of source id
  int<lower=0> N_dv_ow; //number of divers
  int<lower=1,upper=N_dv_ow> diver_ow[N_ow]; // vector of diver identities
  int Z_reef; // columns in the covariate matrix
  int Z_ow; 
  int Z_iucn;
  int Z_aba; 
  int Z_multi;
  matrix[N_reef,Z_reef] X_reef; // design matrix X
  matrix[N_ow,Z_ow] X_ow; 
  matrix[N_iucn,Z_iucn] X_iucn;
  matrix[N_aba,Z_aba] X_aba; 
  matrix[N_multi,Z_multi] X_multi;
  int K_reef; //ordinal levels
  int K_ow;
  int TT; // total timespan across all datasets
  int<lower=0> N_yr; //number of years total for REEF 
  int yr_index[N_yr]; //index of years
  int<lower=0> N_yr_ow; //number of years with surveys
  int<lower=0> N_yr_iucn;
  int<lower=0> N_yr_aba; 
  int<lower=0> N_yr_multi;
  int<lower=1,upper=N_yr> year_id_reef[N_reef]; //vector of years
  int<lower=1,upper=N_yr> year_id_ow[N_ow]; 
  int<lower=1,upper=N_yr> year_id_iucn[N_iucn];
  int<lower=1,upper=N_yr_aba> year_id_aba[N_aba]; 
  int<lower=1,upper=N_yr_multi> year_id_multi[N_multi];
  int wasting_index[N_yr];
  int tt_convert_ow[N_yr_ow];
  int tt_convert_iucn[N_yr_iucn];
  int tt_convert_aba[N_yr_aba];
  int tt_convert_multi[N_yr_multi];  
}
parameters {
  ordered[K_reef-1] c_reef; //cutpoints
  ordered[K_ow-1] c_ow; //cutpoints
  real x0; //initial popn size
  real a_ow; //scalars for time-series
  real a_iucn;
  real a_aba;
  real a_multi; 

  real<lower=0> phi_iucn; //dispersion parameters for negative binomials
  real<lower=0> phi_aba;
  real<lower=0> phi_multi;
  real theta_multi; //theta term for zero inflation

  vector[Z_reef] beta_reef; //survey-specific effort coefficients
  vector[Z_ow] beta_ow; 
  vector[Z_iucn] beta_iucn;
  vector[Z_aba] beta_aba; 
  vector[Z_multi] beta_multi; 
  vector[N_site_reef] a_site_reef; //deviation between sites
  vector[N_site_ow] a_site_ow; 
  vector[N_site_iucn] a_site_iucn; 
  vector[N_site_aba] a_site_aba; 
  vector[N_site_multi] a_site_multi; 
  vector[N_source_iucn] a_source_iucn;//deviation between sources within dataset
  vector[N_dv_ow] a_dv_ow; //deviation between divers OW
 
  //variance on the deviance components
  real<lower = 0> sd_site_reef;
  real<lower = 0> sd_site_ow;
  real<lower = 0> sd_site_iucn;
  real<lower = 0> sd_site_aba;
  real<lower = 0> sd_site_multi;
  real<lower = 0> sd_dv_ow;
  real<lower = 0> sd_source_iucn;
  real<lower = 0> sd_r_reef; //variance on observation deviations 
  real<lower = 0> sd_r_ow; 
  real<lower = 0> sd_r_iucn; 
  real<lower = 0> sd_r_aba; 
  real<lower = 0> sd_r_multi;
  real<lower = 0> sd_q; //variance on process deviance for underlying state
  real<lower = 0> sd_wasting; //variance on process in catastrophic years
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev_reef; //observation deviations 
  vector[N_yr_ow] obs_dev_ow; 
  vector[N_yr_iucn] obs_dev_iucn;
  vector[N_yr_aba] obs_dev_aba; 
  vector[N_yr_multi] obs_dev_multi;

}

transformed parameters{
  vector[TT] x; //underlying state (i.e., estimated true population)
  vector[N_yr] a_yr_reef; //year by year estimates
  vector[N_yr_ow] a_yr_ow; 
  vector[N_yr_iucn] a_yr_iucn;
  vector[N_yr_aba] a_yr_aba; 
  vector[N_yr_multi] a_yr_multi;
  
  //this section forms the state-space component of this model - a_yr are the
  //estimated year by year abundances derived from the ordered logistic or
  //negative binomial regressions in the model{} section, and x is the 
  //underlying true population state

  //the underlying state for each year (x[t]), other than the first is based on 
  //the previous year's state (x[t-1]) plus some process deviation
  //this process deviation accounts for all "true" changes that are driven by
  //demographic processes
  //note that the odd notation on the standard deviation is called non-centered
  //parameterization: http://econ.upf.edu/~omiros/papers/val7.pdf
  //https://mc-stan.org/docs/2_28/stan-users-guide/reparameterization.html
  //when pro_dev ~ std_normal(), pro_dev*sd_q is the same as just pro_dev when
  //pro_dev ~ normal(0, sd_q), this notation is just more computationally 
  //efficient
  //the first year follows the same equation, but based on some unknown previous
  //year, x0, which we need to estimate
  x[1] = x0 + pro_dev[1]*sd_q;
  
  //we've also defined years as "normal" (i.e., wasting_index == 0) or 
  //"catastrophic" (i.e., wasting_index == 1) based on whether a given year was
  //during the main part of the SSWD outbreak in Canada (i.e., 2014-2015), and
  //allowed the SD for those to categories of years to vary. This allows process
  //error to be more extreme in the crash years (sd_wasting) without influencing 
  //the estimated process error in normal years (sd_q), meaning we're less 
  //likely to underestimate observation errors - see Appendix in report for more
  for(t in 2:TT){
    if(wasting_index[t] == 0)
      x[t] = x[t-1] + pro_dev[t]*sd_q;
    else
      x[t] = x[t-1] + pro_dev[t]*sd_wasting;
  }

  //the next five equations link the separate year by year estimates to the same
  //underlying state, x. The model looks for the value of x that is the best
  //fit given all the datasets, and estimates the corresponding observation 
  //error for each. To see the links between these equations and the one above 
  //(i.e., how the model can partition fluctuations between both observation 
  //error and process changes, even if there were only one dataset), see Dennis 
  //et al., 2006 https://doi.org/10.1890/0012-9615(2006)76[323:EDDPNA]2.0.CO;2
   
  //for each year, the observed estimate is just the true state plus some 
  //observation error (with the same non-centered parameterization for 
  //efficiency)
  for(i in 1:N_yr){
    a_yr_reef[i] = x[yr_index[i]] + obs_dev_reef[i]*sd_r_reef; 
  }
  //the one key difference here is that we need to scale all of the datasets to 
  //a baseline dataset to account for consistent differences between the five 
  //(otherwise the observation error will likely always be positive for some 
  //datasets and negative for others). Here, we use an additive scalar (i.e., 
  //"+ a") to scale up any datasets that have consistently lower counts than
  //our baseline (in this case, the REEF surveys) or scale down any datasets 
  //that have higher counts
  
  for(i in 1:N_yr_ow){
    a_yr_ow[i] = a_ow + x[tt_convert_ow[i]] + obs_dev_ow[i]*sd_r_ow; 
  }
  for(i in 1:N_yr_iucn){
    a_yr_iucn[i] = a_iucn + x[tt_convert_iucn[i]] + obs_dev_iucn[i]*sd_r_iucn; 
  }
  for(i in 1:N_yr_aba){
    a_yr_aba[i] = a_aba + x[tt_convert_aba[i]] + obs_dev_aba[i]*sd_r_aba; 
  }
  for(i in 1:N_yr_multi){
    a_yr_multi[i] = a_multi + x[tt_convert_multi[i]] + obs_dev_multi[i]*sd_r_multi; 
  }
}  

model{
  //priors on sds for process and observation error
  //gamma priors on variance terms provide reasonable constraints, and the model 
  //does not appear to be strongly influenced by the choice of prior here, as I 
  //have tried a range of options with the same result. Can check the shape of 
  //these priors in R with bayesAB::plotGamma
  sd_q ~ gamma(2,4);    
  sd_wasting ~ gamma(2,2); //wider prior since we want to give the model a lot 
                          //of flexibility to vary in the year wasting hit
  sd_r_reef ~ gamma(2,4);  
  sd_r_ow ~ gamma(2,4);    
  sd_r_iucn ~ gamma(2,4);
  sd_r_aba ~ gamma(2,4);  
  sd_r_multi ~ gamma(2,4);
  x0 ~ normal(0,5); //initial state
  a_ow ~ normal(0,5); //scalars
  a_iucn ~ normal(0,5);
  a_aba ~ normal(0,5);
  a_multi ~ normal(0,5);
  //varying intercepts of observation and process errors
  pro_dev ~ std_normal(); 
  obs_dev_reef ~ std_normal();
  obs_dev_ow ~ std_normal();
  obs_dev_iucn ~ std_normal();
  obs_dev_aba ~ std_normal();
  obs_dev_multi ~ std_normal();
  
  //priors on fixed effects and cutpoints
  c_reef ~ induced_dirichlet(rep_vector(1, K_reef), 0); //custom prior function
  c_ow ~ induced_dirichlet(rep_vector(1, K_ow), 0); 
  beta_reef ~ normal(0,2); 
  beta_ow ~ normal(0,2); 
  beta_iucn ~normal(0,2);
  beta_aba ~ normal(0,2); 
  beta_multi ~ normal(0,2);
  //dispersion parameters for negative binomial models - inv_gamma does a good
  //job of constraining on both the lower bound and to realistic upper values
  //and changing the choice of shape and scale does not change the model fit
  phi_iucn ~ inv_gamma(5,5);
  phi_aba ~ inv_gamma(5,5);
  phi_multi ~ inv_gamma(5,5);
  //zero inflation component of multi dataset
  theta_multi ~ normal(0,3);
  
  //sds on random effects
  sd_site_reef ~ gamma(2,3); 
  sd_site_ow ~ gamma(2,3); 
  sd_site_iucn ~ gamma(2,3);
  sd_site_aba ~ gamma(2,3); 
  sd_site_multi ~ gamma(2,3); 
  sd_dv_ow ~ gamma(2,3);   
  sd_source_iucn ~ gamma(2,3);
  
  //varying intercepts of random effects
  a_site_reef ~ std_normal();
  a_site_ow ~ std_normal();
  a_site_iucn ~ std_normal();
  a_site_aba ~ std_normal();
  a_site_multi ~ std_normal();
  a_dv_ow ~ std_normal();
  a_source_iucn ~ std_normal();
  //note that there's no a_yr prior, since it's in the transformed parameter 
  //sections

  //main ordered logistic/negative binomial regression models for each dataset
  //since there is no intercept in either model and all continuous variables are
  //centered around 0 (as well as the random effects of site and diver), the 
  //a_yr term tells us the expected abundance for each year at the mean value
  //of all the other variables
  y_reef ~ ordered_logistic(a_yr_reef[year_id_reef]+a_site_reef[site_reef]*sd_site_reef+X_reef*beta_reef,c_reef);
  y_ow ~ ordered_logistic(a_yr_ow[year_id_ow]+a_site_ow[site_ow]*sd_site_ow+a_dv_ow[diver_ow]*sd_dv_ow+X_ow*beta_ow,c_ow);
  y_iucn ~ neg_binomial_2_log(a_yr_iucn[year_id_iucn]+a_site_iucn[site_iucn]*sd_site_iucn+a_source_iucn[source_iucn]*sd_source_iucn+X_iucn*beta_iucn, phi_iucn);
  y_aba ~ neg_binomial_2_log(a_yr_aba[year_id_aba]+a_site_aba[site_aba]*sd_site_aba+X_aba*beta_aba, phi_aba);
  
  //based on posterior predictive checks of a negative binomial model, the multi
  //dataset appears to be quite strongly zero inflated (i.e., all iterations 
  //underpredict the number of zeros), so we'll add a simple zero inflation 
  //component to this part of the model to help with that
  for (i in 1:N_multi){
    if(y_multi[i] == 0) {
      target += log_sum_exp(bernoulli_logit_lpmf(1 | theta_multi),
                            bernoulli_logit_lpmf(0 | theta_multi)
                              + neg_binomial_2_log_lpmf(y_multi[i] | a_yr_multi[year_id_multi[i]]+a_site_multi[site_multi[i]]*sd_site_multi+X_multi[i,]*beta_multi, phi_multi));
    }
    else {
       target += bernoulli_logit_lpmf(0 | theta_multi)
                  + neg_binomial_2_log_lpmf(y_multi[i] | a_yr_multi[year_id_multi[i]]+a_site_multi[site_multi[i]]*sd_site_multi+X_multi[i,]*beta_multi, phi_multi);
    }
  }
}
  generated quantities{
  vector[N_reef] y_predict_reef; //generate model predictions to assess fit
  vector[N_ow] y_predict_ow; 
  vector[N_iucn] y_predict_iucn; 
  vector[N_aba] y_predict_aba; 
  vector[N_multi] y_predict_multi; 
  vector[N_reef+N_ow+N_iucn+N_aba+N_multi] log_lik; //generate log likelihoods 

  //make separate y_predicts for each dataset
  for (i in 1:N_reef) y_predict_reef[i] = ordered_logistic_rng(a_yr_reef[year_id_reef[i]]+a_site_reef[site_reef[i]]*sd_site_reef+X_reef[i,]*beta_reef,c_reef); 
  for (z in 1:N_ow) y_predict_ow[z] = ordered_logistic_rng(a_yr_ow[year_id_ow[z]]+a_site_ow[site_ow[z]]*sd_site_ow+a_dv_ow[diver_ow[z]]*sd_dv_ow+X_ow[z,]*beta_ow,c_ow); 
  for (j in 1:N_iucn) y_predict_iucn[j] = neg_binomial_2_log_rng(a_yr_iucn[year_id_iucn[j]]+a_site_iucn[site_iucn[j]]*sd_site_iucn+a_source_iucn[source_iucn[j]]*sd_source_iucn+X_iucn[j,]*beta_iucn, phi_iucn);
  for (k in 1:N_aba) y_predict_aba[k] = neg_binomial_2_log_rng(a_yr_aba[year_id_aba[k]]+a_site_aba[site_aba[k]]*sd_site_aba+X_aba[k,]*beta_aba, phi_aba);
  for (m in 1:N_multi) {
    if(bernoulli_logit_rng(theta_multi) == 1)
    y_predict_multi[m] = 0;
    else
    y_predict_multi[m] = neg_binomial_2_log_rng(a_yr_multi[year_id_multi[m]]+a_site_multi[site_multi[m]]*sd_site_multi+X_multi[m,]*beta_multi, phi_multi);
  }
  //and a single log likelihood for the whole model
  for (i in 1:N_reef) log_lik[i] = ordered_logistic_lpmf(y_reef[i]|a_yr_reef[year_id_reef[i]]+a_site_reef[site_reef[i]]*sd_site_reef+X_reef[i,]*beta_reef, c_reef);
  for (z in 1:N_ow) log_lik[N_reef+z] = ordered_logistic_lpmf(y_ow[z]|a_yr_ow[year_id_ow[z]]+a_site_ow[site_ow[z]]*sd_site_ow+a_dv_ow[diver_ow[z]]*sd_dv_ow+X_ow[z,]*beta_ow, c_ow);
  for (j in 1:N_iucn) log_lik[N_reef+N_ow+j] = neg_binomial_2_log_lpmf(y_iucn[j]|a_yr_iucn[year_id_iucn[j]]+a_site_iucn[site_iucn[j]]*sd_site_iucn+a_source_iucn[source_iucn[j]]*sd_source_iucn+X_iucn[j,]*beta_iucn, phi_iucn);
  for (k in 1:N_aba) log_lik[N_reef+N_ow+N_iucn+k] = neg_binomial_2_log_lpmf(y_aba[k]|a_yr_aba[year_id_aba[k]]+a_site_aba[site_aba[k]]*sd_site_aba+X_aba[k,]*beta_aba, phi_aba);
    for (m in 1:N_multi){
    if(y_multi[m] == 0) {
      log_lik[N_reef+N_ow+N_iucn+N_aba+m] = log_sum_exp(bernoulli_logit_lpmf(1 | theta_multi),
                                                        bernoulli_logit_lpmf(0 | theta_multi)
                                                            + neg_binomial_2_log_lpmf(y_multi[m] | a_yr_multi[year_id_multi[m]]+a_site_multi[site_multi[m]]*sd_site_multi+X_multi[m,]*beta_multi, phi_multi));
    }
    else {
      log_lik[N_reef+N_ow+N_iucn+N_aba+m] = bernoulli_logit_lpmf(0 | theta_multi)
                                                + neg_binomial_2_log_lpmf(y_multi[m] | a_yr_multi[year_id_multi[m]]+a_site_multi[site_multi[m]]*sd_site_multi+X_multi[m,]*beta_multi, phi_multi);
    }
  }
}
//leave an empty line at the end
