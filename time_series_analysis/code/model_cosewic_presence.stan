//upload the data as defined in the list in pycno_presence.R
data{
  int<lower=1> N1;//number of observations 
  int<lower=0> N_site1; //number of sites 
  int<lower=1,upper=N_site1> site1[N1]; // vector of site identities
  int<lower=0> TT; //total timespan
  int<lower=0> yr_index[TT]; //index of years across all datasets
  int<lower=0> N_yr; //number of years with data
  int<lower=1,upper=N_yr> year_id1[N1]; // vector of year
  int presence1[N1]; //was a pycno observed on a survey
  int wasting_index[TT]; //index defining "normal" and "catastrophic" years
  int<lower=0> N_source1;
  int<lower=1,upper=N_source1> source1[N1];
  vector[N1] time; //time spent on survey (scaled)
  vector[N1] area; //area surveyed (scaled)
  vector[N1] effort_type; //whether time (0) or area (1) was used as effort metric
}

parameters {
  real x0; //initial popn size
  vector[N_site1] a_site1; //deviation between sites
  vector[N_source1] a_source; //deviation between sources
  real<lower = 0> sd_source; //deviation between sources
  real beta_time; //effort coefficient for time
  real beta_area; //effort coefficient for area
  real a_area; //scalar for area vs time
 
  //variance on the deviance components
  real<lower = 0> sd_site1; //sites
  real<lower = 0> sd_r1; //observation error
  real<lower = 0> sd_q; //process error - normal years
  real<lower = 0> sd_wasting; //process error - catastrophic
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr] obs_dev1; //observation deviations 
}
transformed parameters{
  vector[TT] x; 
  vector[N_yr] a_yr1;

  //this section forms the state-space component of this model - a_yr is the
  //estimated year by year abundances derived from the logistic regression in 
  //the model{} section, and x is the underlying true population state

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
  for(t in 2:TT){
    if(wasting_index[t] == 0)
      x[t] = x[t-1] + pro_dev[t]*sd_q;
    else
      x[t] = x[t-1] + pro_dev[t]*sd_wasting;
  }
  
  //the next equation links the year by year estimates to the 
  //underlying state, x. The model uses the temporal autocorrelation to 
  //estimate the corresponding observation vs process error for each year. 
  //To see the links between this equation and the one above 
  //(i.e., how the model can partition fluctuations between both observation 
  //error and process changes, even if there were only one dataset), see Dennis 
  //et al., 2006 https://doi.org/10.1890/0012-9615(2006)76[323:EDDPNA]2.0.CO;2
   
  //for each year, the observed estimate is just the true state plus some 
  //observation error (with the same non-centered parameterization for 
  //efficiency) 

  for(i in 1:N_yr){
    a_yr1[i] = x[yr_index[i]] + obs_dev1[i]*sd_r1; 
  }
}  

model{
  //priors
  x0 ~ normal(0,5); 
  beta_time ~ normal(0,2); 
  beta_area ~ normal(0,2); 
  a_area ~ normal(0,2);
  //variance terms - gamma priors seem to constrain the parameters to realistic
  //values well and the model does not appear to be sensitive to the choice of
  //alpha and beta
  sd_q ~ gamma(2,4);   
  sd_wasting ~ gamma(2,4);
  sd_r1 ~ gamma(2,4);  
  sd_site1 ~ gamma(2,3);
  sd_source ~ gamma(2,3);

  //varying intercepts
  a_site1 ~ std_normal();
  a_source ~ std_normal();
  pro_dev ~ std_normal();
  obs_dev1 ~ std_normal(); 

//main logistic regression - note that the same year, site, and source terms are
//being estimated from the full dataset, it's just the time/area predictors and
//the scalar between the two that are split
  for(i in 1:N1){
    if(effort_type[i] == 0){
      presence1[i] ~ bernoulli_logit(a_yr1[year_id1[i]]+a_site1[site1[i]]*sd_site1+a_source[source1[i]]*sd_source+time[i]*beta_time);
    }
    else{
      presence1[i] ~ bernoulli_logit(a_yr1[year_id1[i]]+a_site1[site1[i]]*sd_site1+a_source[source1[i]]*sd_source+area[i]*beta_area+a_area);
    }
  }
}


generated quantities{
  vector[N1] y_predict; //generate model predictions to assess fit
  vector[N1] log_lik; //generate log likelihoods
  
  for(i in 1:N1){
    if(effort_type[i] == 0){
        y_predict[i] = bernoulli_logit_rng(a_yr1[year_id1[i]]+a_site1[site1[i]]*sd_site1+a_source[source1[i]]*sd_source+time[i]*beta_time);
        log_lik[i] = bernoulli_logit_lpmf(presence1[i] | a_yr1[year_id1[i]]+a_site1[site1[i]]*sd_site1+a_source[source1[i]]*sd_source+time[i]*beta_time);
    }
    else{
        y_predict[i] = bernoulli_logit_rng(a_yr1[year_id1[i]]+a_site1[site1[i]]*sd_site1+a_source[source1[i]]*sd_source+area[i]*beta_area+a_area);
        log_lik[i] = bernoulli_logit_lpmf(presence1[i] | a_yr1[year_id1[i]]+a_site1[site1[i]]*sd_site1+a_source[source1[i]]*sd_source+area[i]*beta_area+a_area);
    }
  }
}

//leave empty line at end
