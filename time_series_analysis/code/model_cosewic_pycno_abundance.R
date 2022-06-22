library(tidyverse) #for data wrangling
library(cmdstanr) #for running model
library(rstan) #for extracting and processing model output
library(DHARMa) #for examining residuals
library(patchwork) #for plotting multiple panels
source("code/custom_functions.R")

#load in data------------------------------------------------------------------
reef <- read_csv("../data/cleaned_reef_data.csv")
ow <- read_csv("../data/cleaned_oceanwise_data.csv",
               col_types = list(comments = col_character()))
iucn_pres_com <- read_csv("../data/cleaned_iucn_data.csv") %>% 
  #remove the surveys without an area measurement since we need to be able to 
  #control for survey effort (this just removes some of the Hakai surveys and 
  #the Cote dataset)
  filter(!is.na(area)) %>% 
  filter(!is.na(depth)) %>% 
  mutate(year = event_year) %>%
  #removing pre-2009 data since there are only 7 surveys per year before that
  filter(year > 2008)
aba <- read_csv("../data/cleaned_dfo_abalone.csv")
multi <- read_csv("../data/cleaned_dfo_multispecies.csv") %>% 
  unite("site", c(stat_area, sub_area), remove = FALSE)

#prep data for Stan------------------------------------------------------------
#X1 and X2 are the fixed effects matrices for the survey-level effects
X_reef <- matrix(data=c(scale(as.numeric(reef$btime)),
                        scale(as.numeric(reef$maxdepth)), 
                        scale(as.numeric(reef$visibility)),
                        scale(as.numeric(reef$current)),
                        reef$exp_binary,
                        reef$multiple_surveys),
                 ncol=6,nrow=nrow(reef))
X_ow <- matrix(data=c(scale(as.numeric(ow$bottom_time)),
                      scale(as.numeric(ow$max_depth)),
                      ow$multiple_surveys),
               ncol=3,nrow=nrow(ow))
X_iucn <- matrix(data=c(scale(as.numeric(iucn_pres_com$depth)),
                        scale(as.numeric(iucn_pres_com$area))),
                 ncol=2,
                 nrow=nrow(iucn_pres_com))
X_aba <- matrix(data=c(scale(as.numeric(aba$num_quadrats))),
                ncol=1,nrow=nrow(aba))
X_multi <- matrix(data=c(scale(as.numeric(multi$max_of_cor_depth_m)),
                         scale(as.numeric(multi$num_quadrats))),
                  ncol=2, nrow=nrow(multi))
#create a year-level vector for whether or not a given year was pre or post 
#wasting
wasting_index <- reef %>% 
  select(year, wasting) %>% 
  distinct() %>% 
  arrange(year) %>% 
  select(wasting)
yr_index <- tibble(x = seq(from=1, to=22))

#make an index linking the years in each dfo dataset to the full timespan
tt_ow <- ow %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

tt_iucn <- iucn_pres_com %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

tt_aba <- aba %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

tt_multi <- multi %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

#define all the data for the Stan model
data_pycno_abun = list(TT=22,
                       y_reef = reef$abundance_recode,
                       N_reef = nrow(reef),
                       site_reef=as.numeric(factor(reef$geogr)),
                       N_site_reef=length(unique(reef$geogr)),
                       K_reef=length(unique(reef$abundance)),
                       X_reef=X_reef,
                       Z_reef=ncol(X_reef),
                       N_yr=length(unique(reef$year)),
                       year_id_reef=as.numeric(factor(reef$year)),
                       y_ow = ow$abundance_recode,
                       N_ow = nrow(ow),
                       site_ow=as.numeric(factor(ow$place_name)),
                       N_site_ow=length(unique(ow$place_name)),
                       diver_ow=as.numeric(factor(ow$observer_name)),
                       N_dv_ow=length(unique(ow$observer_name)),
                       K_ow=length(unique(ow$abundance_recode)),
                       X_ow=X_ow,
                       Z_ow=ncol(X_ow),
                       N_yr_ow=length(unique(ow$year)),
                       year_id_ow=as.numeric(factor(ow$year)),
                       y_iucn = iucn_pres_com$totalpycnos,
                       N_iucn = nrow(iucn_pres_com),
                       site_iucn=as.numeric(factor(iucn_pres_com$location_code_numeric)),
                       N_site_iucn=length(unique(iucn_pres_com$location_code_numeric)),
                       X_iucn=X_iucn,
                       Z_iucn=ncol(X_iucn),
                       N_yr_iucn=length(unique(iucn_pres_com$event_year)),
                       year_id_iucn=as.numeric(factor(iucn_pres_com$event_year)),
                       source_iucn=as.numeric(factor(iucn_pres_com$source)),
                       N_source_iucn=length(unique(iucn_pres_com$source)),
                       yr_index=sort(unique(as.numeric(factor(reef$year)))),
                       yr_index_ow=sort(unique(as.numeric(factor(ow$year)))),
                       yr_index_iucn=sort(unique(as.numeric(factor(
                         iucn_pres_com$event_year)))),
                       wasting_index = wasting_index$wasting,
                       tt_convert_ow = tt_ow$tt,
                       tt_convert_iucn = tt_iucn$tt,
                       tt_convert_aba = tt_aba$tt,
                       tt_convert_multi = tt_multi$tt,
                       X_aba = X_aba,
                       X_multi = X_multi,
                       Z_aba = ncol(X_aba),
                       Z_multi = ncol(X_multi),
                       year_id_aba = as.numeric(factor(aba$year)),
                       year_id_multi = as.numeric(factor(multi$year)),
                       N_yr_aba = length(unique(aba$year)),
                       N_yr_multi = length(unique(multi$year)),
                       N_aba = nrow(aba),
                       N_multi = nrow(multi),
                       y_aba = aba$num_pycno,
                       y_multi = multi$num_pycno,
                       site_aba=as.numeric(factor(aba$stat_area)),
                       N_site_aba=length(unique(aba$stat_area)),
                       site_multi=as.numeric(factor(multi$site)),
                       N_site_multi=length(unique(multi$site)))

#note that we're using the cmdstanr package for this model because it uses the
#latest version of Stan and there is a bug in the older version that rstan uses.
#the cutpoints on the ordered logistic regressions are super highly constrained, 
#and sometimes at the beginning of warmup Stan tries out values of the cutpoints
#that don't match those constraints - this is usually no problem but for some
#reason the older version of Stan causes a chain to abort if it tests out these
#inappropriate values (but only when there is an rng function in the generated
#quantities section for some reason). If you try running the model in cmdstanr
#you'll still get the same warnings, but Stan quickly figures out the 
#appropriate values to test out and runs fine - note that these warnings are 
#only cause for concern if they happen after the warmup!
abun_mod_code <- cmdstan_model("code/cosewic_abundance_model.stan")
abun_mod <- abun_mod_code$sample(data = data_pycno_abun,
                                 #set seed for reproducibility - originally just
                                 #used a random seed from cmdstanr and then 
                                 #extracted with the $metadata() function so the
                                 #exact results can be extracted again
                                 seed = 910420862, 
                                 refresh = 100, #how often to print model progress
                                 chains = 4,
                                 parallel_chains = 4, 
                                 iter_warmup = 1000,
                                 iter_sampling = 3000
                                 )
#need to use save_object instead of saveRDS for cmdstanr objects
#and we also need to actually run any of the things we care about if we want 
#them to be saved which is wild but whatever
md <- abun_mod$metadata()
#check for convergence issues
di <- abun_mod$cmdstan_diagnose()
#abun_mod$save_object("model_output/pycno_abundance_model_cmdstanr.rds")

#and we'll convert and save as an rstan object too, just to be safe
abun_mod_rstan <- rstan::read_stan_csv(abun_mod$output_files())
#saveRDS(abun_mod_rstan, "model_output/pycno_abundance_model_rstan.rds")

#run model checks--------------------------------------------------------------
#read in model object if not already run
abun_mod <- readRDS("model_output/pycno_abundance_model_cmdstanr.rds")
abun_mod_rstan <- readRDS("model_output/pycno_abundance_model_rstan.rds")

#check model metadata
metadata <- abun_mod$metadata()
#check for convergence issues
abun_mod$cmdstan_diagnose()
#check r-hat values (should all be 1)
abundance_summary <- summary(abun_mod_rstan)$summary
#look good from a quick check but we can double check to make sure all are 
#exactly 1 when rounded to two decimal places
abundance_summary %>% 
  as_tibble() %>% 
  mutate(rhat = round(Rhat,2)) %>% 
  filter(rhat != 1)
#and there isn't a single one that's not 1 - yay! 

#posterior predictive check to compare model predictions with observed data - 
#need to use custom function to work with rstan object - this function works
#the same way as the ppc and pp_check functions available in other packages
#see custom_functions.R for details
#combine reef and ow for plotting with y_predict
all_pycno <- reef %>% 
  transmute(abundance_recode=abundance_recode, source="REEF") %>% 
  bind_rows(transmute(ow, abundance_recode=abundance_recode, source = "OW")) %>% 
  bind_rows(transmute(iucn_pres_com, abundance_recode = totalpycnos, 
                      source = "IUCN")) %>% 
  bind_rows(transmute(aba, abundance_recode = num_pycno, source="ABA")) %>% 
  bind_rows(transmute(multi, abundance_recode = num_pycno, source="MULTI"))

#to get it to work on all five datasets without needing to copy and paste the 
#code a bunch, we'll make some lists to differentiate the datasets and then 
#run a for loop to store the plots
n_source <- 5
sources <- c("REEF","OW","IUCN","ABA","MULTI")
y_predicts <- c("y_predict_reef", "y_predict_ow", "y_predict_iucn",
                "y_predict_aba", "y_predict_multi")
response <- list(reef$abundance_recode, ow$abundance_recode, 
                 iucn_pres_com$totalpycnos, aba$num_pycno, multi$num_pycno)
ppc_plots <- list()
iter <- 12000

for(i in 1:n_source){
  ppc_plots[[i]] <- ppc_stanfit(model = abun_mod_rstan,
                                iter = iter,
                                response = response[[i]],
                                plot_type = "density",
                                y_predict = y_predicts[[i]])
}
#the posterior predictive checks actually look really good!

#check residuals
for(i in 1:n_source){
  dharma_resids(model = abun_mod_rstan,
                        iter = iter,
                        response = response[[i]],
                        y_predict = y_predicts[[i]])
}
#and the residuals aren't perfect for all five datasets but given the 
#complexity of the model and the massive variation across datasets, they look
#not too shabby
#note also that all of the DHARMa tests are based on significance thresholds, so 
#when the sample size is large, there are nearly always going to be significant 
#p-values! 

#examine pre/post crash changes-------------------------------------------------
abun_mod <- readRDS("model_output/pycno_abundance_model_cmdstanr.rds")
abun_mod_rstan <- readRDS("model_output/pycno_abundance_model_rstan.rds")
params1 <- rstan::extract(abun_mod_rstan)
iter <- 12000

#use ord_to_n to compare popn estimates from different years to look at percent
#declines
pre_count <- data.frame(matrix(NA, nrow = iter, ncol =14))
for(i in 1:14){
  pre_count[,i] <- ord_to_n(params1$x[,i],params1$c_reef)
}
mean_pre <- rowMeans(pre_count)

#we won't include 2014 in either pre or post, since we don't
#know the exact date that SSWD hit
post_count <- data.frame(matrix(NA, nrow = iter, ncol = 7))
for(i in 16:22){
  post_count[,i-15] <- ord_to_n(params1$x[,i],params1$c_reef)
}
mean_post <- rowMeans(post_count)

diff <- as_tibble(cbind(mean_pre, mean_post)) %>% 
  mutate(diff = mean_post/mean_pre)
median(diff$diff - 1)
quantile(diff$diff,0.025)-1
quantile(diff$diff,0.975)-1

#plot model output--------------------------------------------------------------
#pycno_mod_abundance <- 
#  readRDS("model_output/pycno_abundance_model_all_data_additive.rds")
#the way we convert model estimates (i.e., in link space) to real space will be
#different for each dataset (e.g., if we want to plot them on the same scale, we
#need to use the same transformation but with the dataset specific scalars, or 
#if we want to plot them on their own scales, we need to use different 
#transformations for the ordered logistic regressions vs. the negative binomial
#regressions), so we're going to extract the relevant info for each and then 
#make a for loop to help reduce repetition in the code
df_list <- list(reef, ow, iucn_pres_com, aba, multi)
n_dfs <- length(df_list)
iter <- 12000
TT <- data_pycno_abun$TT
params1 <- rstan::extract(abun_mod_rstan)
#need to specify the name of the year by year estimates associated with each df
a_yr <- list(params1$a_yr_reef, params1$a_yr_ow, params1$a_yr_iucn,
             params1$a_yr_aba, params1$a_yr_multi)
a <- list(rep(0, times = iter), params1$a_ow, params1$a_iucn, 
          params1$a_aba, params1$a_multi)

#create arrays to populate
est_x <- array(NA, c(iter, TT))
x_mat <- array(NA, c(TT, 3))

#extract the estimates for the underlying state for each year in the full 
#timespan
for(i in 1:TT){
  x_coef <- data.frame(p_0.x=NA,
                       p_1.x=NA,p_2.x=NA,p_11.x=NA,
                       p_101.x=NA,
                       lambda.x=NA,
                       iter = seq(1,iter))
  x_coef[,1]<- plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,2]<- plogis(params1$c_reef[,2]-params1$x[,i])-plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,3]<- plogis(params1$c_reef[,3]-params1$x[,i])-plogis(params1$c_reef[,2]-params1$x[,i])
  x_coef[,4]<- plogis(params1$c_reef[,4]-params1$x[,i])-plogis(params1$c_reef[,3]-params1$x[,i])
  x_coef[,5]<- 1-plogis(params1$c_reef[,4]-params1$x[,i])
  x_coef[,6]<- apply(x_coef[,1:5],1,abund_tranfs)
  
  est_x[,i] <- x_coef[,6]
  
  x_mat[i,1]=median(est_x[,i])
  x_mat[i,2]=quantile(est_x[,i],0.025)
  x_mat[i,3]=quantile(est_x[,i],0.975)
}

#now make a dataframe with the actual years that we can connect the temporally-
#patchy observation estimates to
dat <- as.data.frame(x_mat) %>% 
  transmute(state_prob = V1,
            l.95.x = V2,
            u.95.x = V3,
            year = seq(from = 2000, to = 2021, by = 1))

#for each dataset, we're going to make some empty arrays to populate for each
#of the years with actual survey data
for(j in 1:n_dfs){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  #and then for each of those years, we'll backtransform based on the REEF scale
  #(i.e., after applying the relevant "a" term)
  for(i in 1:n_yr){
    a_coef <- data.frame(p_0=NA,
                         p_1=NA,p_2=NA,p_11=NA,
                         p_101=NA,
                         lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- plogis(params1$c_reef[,1]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,2]<- plogis(params1$c_reef[,2]-(a_yr[[j]][,i]-a[[j]]))-plogis(params1$c_reef[,1]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,3]<- plogis(params1$c_reef[,3]-(a_yr[[j]][,i]-a[[j]]))-plogis(params1$c_reef[,2]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,4]<- plogis(params1$c_reef[,4]-(a_yr[[j]][,i]-a[[j]]))-plogis(params1$c_reef[,3]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,5]<- 1-plogis(params1$c_reef[,4]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,6]<- apply(a_coef[,1:5],1,abund_tranfs)
    
    est_a[,i] <- a_coef[,6]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  #then connect the correct sampling year with the observation estimate
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  #then add it to the main plotting dataframe
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob = `1`,
              l.95.o = `2`,
              u.95.o = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#and we can extract the info on the number of surveys conducted each year across
#all datasets to add to the plot
num_surveys <- reef %>% 
  bind_rows(ow) %>% 
  bind_rows(transmute(iucn_pres_com, year = event_year)) %>% 
  bind_rows(aba) %>% 
  bind_rows(multi) %>% 
  group_by(year) %>% 
  summarize(n = n())

#now we can plot all the datasets together on top of the estimated underlying
#state
ggplot(data = dat) +
  geom_ribbon(aes(x = year, y = state_prob, ymin = l.95.x, ymax = u.95.x),
              fill='darkcyan', alpha = 0.2) +
  geom_line(aes(year, state_prob, colour = "Estimated\nstate"), lty=5,lwd=0.8) +
  geom_line(aes(year, obs_prob.x, colour = "REEF"), lwd=0.8) +
  geom_point(aes(year, obs_prob.x, fill = "REEF"), col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob.y, colour = "OW"),lwd=0.8) +
  geom_point(aes(year, obs_prob.y, fill = "OW"), col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob.x.x, colour = "Small datasets"),lwd=0.8) +
  geom_point(aes(year, obs_prob.x.x, fill = "Small datasets"), col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob.y.y, colour = "DFO-Abalone"), lwd=0.8) +
  geom_point(aes(year, obs_prob.y.y, fill = "DFO-Abalone"),col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob, colour = "DFO-Multispecies"),lwd=0.8) +
  geom_point(aes(year, obs_prob, fill = "DFO-Multispecies"),col='white',pch=21,cex=3) +
  geom_text(data = num_surveys, aes(y=rep(-0.2,22), x=year, label = n), 
            size = 2.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_colour_manual(name = "", values = c("Estimated\nstate" = 'darkcyan', 
                                            "REEF" = "dodgerblue4", 
                                            "OW" ='dodgerblue',
                                            "DFO-Abalone" = "darkblue",
                                            "DFO-Multispecies" = "mediumpurple",
                                            "Small datasets" = "dodgerblue3"))+
  scale_fill_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                          "OW" ='dodgerblue',
                                          "DFO-Abalone" = "darkblue",
                                          "DFO-Multispecies" = "mediumpurple",
                                          "Small datasets" = "dodgerblue3"),
                    guide = "none")+
  labs(x = "Year", y = "Annual index of abundance\n(Mean count per survey)") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12))
#ggsave("figures/combined_abundance_model_five.png", width = 8, height = 6)  

#examine effort----------------------------------------------------------------
#each dataset has a different metric of "survey effort" (either area or bottom
#time), and since each dataset is scaled independently, what is a "mean effort" 
#in one dataset may be very different than another. The scalar term in the model
#should account for those mean differences (e.g., if all the surveys in one 
#dataset are very small relative to another, the scalar term 'a' for that 
#dataset will likely be large). Here, we'll determine what these mean efforts 
#are.

reef_effort <- scale(reef$btime) %>% 
  attributes()
ow_effort <- scale(ow$bottom_time) %>% 
  attributes()
iucn_effort <- scale(iucn_pres_com$area) %>% 
  attributes()
aba_effort <- scale(aba$num_quadrats) %>%
  attributes()
multi_effort <- scale(multi$num_quadrats) %>% 
  attributes()

effort_info <- tibble(df = c("reef","ow","iucn","aba","multi"),
                      metric = c("time","time","area","area","area"),
                      mean = NULL,
                      sd = NULL)  
efforts <- list(reef_effort,ow_effort,iucn_effort,aba_effort,multi_effort)

for(i in 1:5){
  effort_info[i,3] <- efforts[[i]][2]
  effort_info[i,4] <- efforts[[i]][3]
}

#plot individual trends with error---------------------------------------------
df_list <- list(reef, ow, iucn_pres_com, aba, multi)
n_dfs <- length(df_list)
iter <- 12000
TT <- data_pycno_abun$TT
params1 <- rstan::extract(abun_mod_rstan)
#need to specify the name of the year by year estimates associated with each df
a_yr <- list(params1$a_yr_reef, params1$a_yr_ow, params1$a_yr_iucn,
             params1$a_yr_aba, params1$a_yr_multi)
a <- list(rep(0, times = iter), params1$a_ow, params1$a_iucn, 
          params1$a_aba, params1$a_multi)

est_x <- array(NA, c(iter, TT))
x_mat <- array(NA, c(TT, 3))

num_surveys_plotting <- list()
for(i in 1:n_dfs){
  num_surveys_plotting[[i]] <- df_list[[i]] %>% 
    group_by(year) %>% 
    summarize(n = n())
}

#for each year in the full timespan
for(i in 1:TT){
  x_coef <- data.frame(p_0.x=NA,
                       p_1.x=NA,p_2.x=NA,p_11.x=NA,
                       p_101.x=NA,
                       lambda.x=NA,
                       iter = seq(1,iter))
  x_coef[,1]<- plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,2]<- plogis(params1$c_reef[,2]-params1$x[,i])-plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,3]<- plogis(params1$c_reef[,3]-params1$x[,i])-plogis(params1$c_reef[,2]-params1$x[,i])
  x_coef[,4]<- plogis(params1$c_reef[,4]-params1$x[,i])-plogis(params1$c_reef[,3]-params1$x[,i])
  x_coef[,5]<- 1-plogis(params1$c_reef[,4]-params1$x[,i])
  x_coef[,6]<- apply(x_coef[,1:5],1,abund_tranfs)
  
  est_x[,i] <- x_coef[,6]
  
  x_mat[i,1]=median(est_x[,i])
  x_mat[i,2]=quantile(est_x[,i],0.025)
  x_mat[i,3]=quantile(est_x[,i],0.975)
}

dat <- as.data.frame(x_mat) %>% 
  transmute(state_prob = V1,
            l.95.x = V2,
            u.95.x = V3,
            year = seq(from = 2000, to = 2021, by = 1))

#these don't need to be for loops but I got lazy and it was easier to cut and
#paste this way
#REEF only
for(j in 1:1){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(p_0=NA,
                         p_1=NA,p_2=NA,p_11=NA,
                         p_101=NA,
                         lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- plogis(params1$c_reef[,1]-a_yr[[j]][,i])
    a_coef[,2]<- plogis(params1$c_reef[,2]-a_yr[[j]][,i])-plogis(params1$c_reef[,1]-a_yr[[j]][,i])
    a_coef[,3]<- plogis(params1$c_reef[,3]-a_yr[[j]][,i])-plogis(params1$c_reef[,2]-a_yr[[j]][,i])
    a_coef[,4]<- plogis(params1$c_reef[,4]-a_yr[[j]][,i])-plogis(params1$c_reef[,3]-a_yr[[j]][,i])
    a_coef[,5]<- 1-plogis(params1$c_reef[,4]-a_yr[[j]][,i])
    a_coef[,6]<- apply(a_coef[,1:5],1,abund_tranfs)
    
    est_a[,i] <- a_coef[,6]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_r = `1`,
              l.95.r = `2`,
              u.95.r = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#OW only
for(j in 2:2){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(p_0=NA,
                         p_1=NA,p_11=NA,p_26=NA,
                         p_51=NA,p_101=NA,p_1001=NA,
                         lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- plogis(params1$c_ow[,1]-a_yr[[j]][,i])
    a_coef[,2]<- plogis(params1$c_ow[,2]-a_yr[[j]][,i])-plogis(params1$c_ow[,1]-a_yr[[j]][,i])
    a_coef[,3]<- plogis(params1$c_ow[,3]-a_yr[[j]][,i])-plogis(params1$c_ow[,2]-a_yr[[j]][,i])
    a_coef[,4]<- plogis(params1$c_ow[,4]-a_yr[[j]][,i])-plogis(params1$c_ow[,3]-a_yr[[j]][,i])
    a_coef[,5]<- plogis(params1$c_ow[,5]-a_yr[[j]][,i])-plogis(params1$c_ow[,4]-a_yr[[j]][,i])
    a_coef[,6]<- plogis(params1$c_ow[,6]-a_yr[[j]][,i])-plogis(params1$c_ow[,5]-a_yr[[j]][,i])
    a_coef[,7]<- 1-plogis(params1$c_ow[,6]-a_yr[[j]][,i])
    a_coef[,8]<- apply(a_coef[,1:7],1,abund_tranfs_ow)
    
    est_a[,i] <- a_coef[,8]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_o = `1`,
              l.95.o = `2`,
              u.95.o = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#IUCN only
for(j in 3:3){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- exp(a_yr[[j]][,i])
    
    est_a[,i] <- a_coef[,1]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_i = `1`,
              l.95.i = `2`,
              u.95.i = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}


#abalone only
for(j in 4:4){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- exp(a_yr[[j]][,i])
    
    est_a[,i] <- a_coef[,1]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_a = `1`,
              l.95.a = `2`,
              u.95.a = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#multi only
for(j in 5:5){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- exp(a_yr[[j]][,i])
    
    est_a[,i] <- a_coef[,1]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_m = `1`,
              l.95.m = `2`,
              u.95.m = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}


#REEF ONLY
reef_plot <- ggplot(data = dat) +
  #geom_ribbon(aes(x = year, y = state_prob, ymin = l.95.x, ymax = u.95.x),
  # fill='darkcyan', alpha = 0.2) +
  #geom_line(aes(year, state_prob, colour = "Estimated\nstate"), lty=5,lwd=0.8) +
  geom_line(aes(year, obs_prob_r, colour = "REEF"), lwd=0.8) +
  geom_point(aes(year, obs_prob_r, fill = "REEF"), col='white',pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_r, ymin = l.95.r, ymax = u.95.r,
                  fill = 'REEF'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[1]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[1]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                            "OW" ='dodgerblue',
                                            "IUCN" = "dodgerblue3",
                                            "DFO-Abalone" = "darkblue",
                                            "DFO-Multispecies" = "lightblue"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                          "OW" ='dodgerblue',
                                          "IUCN" = "dodgerblue3",
                                          "DFO-Abalone" = "darkblue",
                                          "DFO-Multispecies" = "lightblue"),
                    guide = "none")+
  labs(x = "", y = "",
       title = "REEF dive surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

#OW ONLY
ow_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_o, colour = "OW"), lwd=0.8) +
  geom_point(aes(year, obs_prob_o, fill = "OW"), col='white',pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_o, ymin = l.95.o, ymax = u.95.o,
                  fill = 'OW'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[2]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[2]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                            "OW" ='dodgerblue',
                                            "IUCN" = "dodgerblue3",
                                            "DFO-Abalone" = "darkblue",
                                            "DFO-Multispecies" = "lightblue"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                          "OW" ='dodgerblue',
                                          "IUCN" = "dodgerblue3",
                                          "DFO-Abalone" = "darkblue",
                                          "DFO-Multispecies" = "lightblue"),
                    guide = "none")+
  labs(x = "", y = "", title = "OW Pacific Marine\nLife surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

#IUCN ONLY
iucn_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_i, colour = "IUCN"), lwd=0.8) +
  geom_point(aes(year, obs_prob_i, fill = "IUCN"), col='white',
             pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_i, ymin = l.95.i, ymax = u.95.i,
                  fill = 'IUCN'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[3]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[3]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                            "OW" ='dodgerblue',
                                            "IUCN" = "dodgerblue3",
                                            "DFO-Abalone" = "darkblue",
                                            "DFO-Multispecies" = "lightblue"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                          "OW" ='dodgerblue',
                                          "IUCN" = "dodgerblue3",
                                          "DFO-Abalone" = "darkblue",
                                          "DFO-Multispecies" = "lightblue"),
                    guide = "none")+
  labs(x = "Year", y = "", title="Combined small datasets") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)


#abalone ONLY
#the geom_ribbon only works if there are consecuctive years in the model - since
#there were no surveys in 2020, the 2021 point is on it's own and we need to 
#add in a CI manually
aba_subset <- dat %>% 
  select(year, obs_prob_a:u.95.a) %>% 
  filter(year == 2021)


aba_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_a, colour = "DFO-Abalone"), lwd=0.8) +
  geom_point(aes(year, obs_prob_a, fill = "DFO-Abalone"), 
             col='white',pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_a, ymin = l.95.a, ymax = u.95.a,
                  fill = 'DFO-Abalone'),
              alpha = 0.2) +
  geom_linerange(data=aba_subset, aes(x=year,  
                                      ymin = l.95.a, ymax = u.95.a,  
                                      colour = 'DFO-Abalone'), 
                 alpha = 0.2, size = 3) +
  geom_text(data = num_surveys_plotting[[4]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[4]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                            "OW" ='dodgerblue',
                                            "IUCN" = "dodgerblue3",
                                            "DFO-Abalone" = "darkblue",
                                            "DFO-Multispecies" = "lightblue"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                          "OW" ='dodgerblue',
                                          "IUCN" = "dodgerblue3",
                                          "DFO-Abalone" = "darkblue",
                                          "DFO-Multispecies" = "lightblue"),
                    guide = "none")+
  labs(x = "Year", y = "",
       title="DFO-Abalone dive surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

#multi ONLY
multi_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_m, colour = "DFO-Multispecies"), lwd=0.8) +
  geom_point(aes(year, obs_prob_m, fill = "DFO-Multispecies"), col='white',
             pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_m, ymin = l.95.m, ymax = u.95.m,
                  fill = 'DFO-Multispecies'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[5]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[5]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                            "OW" ='dodgerblue',
                                            "IUCN" = "dodgerblue3",
                                            "DFO-Abalone" = "darkblue",
                                            "DFO-Multispecies" = "mediumpurple"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("REEF" = "dodgerblue4", 
                                          "OW" ='dodgerblue',
                                          "IUCN" = "dodgerblue3",
                                          "DFO-Abalone" = "darkblue",
                                          "DFO-Multispecies" = "mediumpurple"),
                    guide = "none")+
  labs(x = "", y = "", 
       title = "DFO-Multispecies\ndive surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

common_y_axis <- ggplot(data.frame(l = "Annual index of abundance (mean count per survey)", 
                                   x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  #theme() +
  coord_cartesian(clip = "off")

common_y_axis+(reef_plot+ow_plot)/
  (multi_plot+aba_plot)/
  (iucn_plot + plot_spacer()) + 
  #now make the label v small relative to the plots
  plot_layout(widths = c(1,25))
#ggsave("figures/individual_dataset_abundance_trends.png", width = 7, height = 8)  

