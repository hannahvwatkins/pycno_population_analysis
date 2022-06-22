library(tidyverse) #for wrangling data
library(rstan) #for running model
library(DHARMa) #for examining residuals
source("code/custom_functions.R") #for processing ordinal data and stan output
rstan_options(auto_write = TRUE)
options(mc.cores = 4) #to parallelize stan model

#load in data------------------------------------------------------------------
#select and rename only the variables that are consistent across all datasets
#since we'll be combining them here
reef <- read_csv("../data/cleaned_reef_data.csv") %>% 
  transmute(year = year, 
            presence = presence,
            site = as.character(geogr),
            source = "REEF",
            time = btime)
ow <- read_csv("../data/cleaned_oceanwise_data.csv",
               col_types = list(comments = col_character())) %>% 
  transmute(year = year, 
            presence = pres_abs,
            site = as.character(location_code),
            source = "OW",
            time = bottom_time)
iucn <- read_csv("../data/cleaned_iucn_data.csv") %>% 
  transmute(year = event_year, 
            presence = pres_abs,
            site = as.character(location_code_numeric),
            source = source,
            time = b_time,
            area = area,
            source = "small")
dfo <- read_csv("../data/cleaned_dfo_dives.csv", guess_max = 3000) %>% 
  unite("location_code", c(lat,lon)) %>% 
  transmute(year = year, 
            presence = presence,
            site = case_when(!is.na(stat_area) ~ as.character(stat_area),
                             TRUE ~ location_code),
            source = source,
            area = num_quadrats)

#prep data for stan---------------------------------------------------------
df <- bind_rows(reef,ow,iucn,dfo) %>% 
  mutate(scaled_time = scale(time),
         scaled_area = scale(area),
         scaled_effort = case_when(!is.na(area) ~ scaled_area,
                                   TRUE ~ scaled_time),
         effort_type = case_when(!is.na(area) ~ 0,
                                 TRUE ~ 1)) %>% 
  filter(!is.na(scaled_effort)) %>% 
  #need to replace NAs with 0s so Stan doesn't get mad at us - we will be using
  #an ifelse statement to make sure none of the NAs are used anyways so it
  #doesn't matter what we replace them with as long as it's numeric
  mutate(time_stan = case_when(is.na(time) ~ 0,
                               TRUE ~ scaled_time),
         area_stan = case_when(is.na(area) ~ 0,
                               TRUE ~ scaled_area))

wasting_index <- read_csv("../data/cleaned_reef_data.csv") %>% 
  select(year, wasting) %>% 
  distinct() %>% 
  arrange(year) %>% 
  select(wasting)
yr_index <- tibble(x = seq(from=1, to=22))


data_presence = list(presence1 = df$presence,
                          N1 = nrow(df),
                          time=as.numeric(df$time_stan),
                          area=as.numeric(df$area_stan),
                          effort_type=df$effort_type,
                          site1=as.numeric(factor(df$site)),
                          N_site1=length(unique(df$site)),
                          yr_index=yr_index$x,
                          TT=22,
                          N_yr=length(unique(df$year)),
                          year_id1=as.numeric(factor(df$year)),
                          wasting_index = wasting_index$wasting,
                          source1=as.numeric(factor(df$source)),
                          N_source1=length(unique(df$source))
)

presence_mod <- rstan::stan(file = "code/cosewic_presence_model.stan", 
                            data = data_presence,
                            seed = 559570178,
                            pars = c('x0',
                                     'sd_site1',
                                     'sd_r1','sd_q', 
                                     'sd_wasting',
                                     'x',
                                     'a_yr1',
                                     'a_site1', 
                                     'pro_dev', 'obs_dev1',
                                     'sd_source', 'a_source',
                                     'beta_time', 'beta_area', 'a_area',
                                     'log_lik', 'y_predict'
                                ),
                            control = list(adapt_delta = 0.83),
                            iter = 4000)
#saveRDS(presence_mod, "model_output/pycno_presence_model.rds")
presence_mod <- readRDS("model_output/pycno_presence_model.rds")
summ <- summary(presence_mod)$summary

#test model fit----------------------------------------------------------------
#check r-hat values (should all be 1)
presence_summary <- summary(presence_mod)$summary
#look good from a quick check but we can double check to make sure all are 
#exactly 1 when rounded to two decimal places
presence_summary %>% 
  as_tibble() %>% 
  mutate(rhat = round(Rhat,2)) %>% 
  filter(rhat != 1)
#and there isn't a single one that's not 1 - yay! 

#check posterior predictions
iter <- 8000
all_dat <- df %>% 
  select(presence)
#custom function
ppc_stanfit(model = presence_mod,
            iter = iter,
            response = all_dat$presence,
            plot_type = "density")
#the posterior predictions looks great

#custom function for residual analysis
dharma_resids(model = presence_mod,
              iter = iter,
              response = all_dat$presence)
#and the residuals aren't perfect but given the massive variation across 
#datasets, they look not too shabby
#note also that all of the DHARMa tests are based on significance thresholds, so 
#when the sample size is large, there are nearly always going to be significant 
#p-values! 

#examine pre/post crash changes-------------------------------------------------
params1  <- rstan::extract(presence_mod)
iter <- 8000

pre_count <- data.frame(matrix(NA, nrow = iter, ncol = 14))
for(i in 1:14){
  pre_count[,i] <-plogis(params1$x[,i])
}
mean_pre <- rowMeans(pre_count)

#we won't include 2014 in either pre or post, since we don't
#know the exact date that SSWD hit
post_count <- data.frame(matrix(NA, nrow = iter, ncol = 7))
for(i in 16:22){
  post_count[,i-15] <- plogis(params1$x[,i])
}
mean_post <- rowMeans(post_count)

diff <- as_tibble(cbind(mean_pre, mean_post)) %>% 
  mutate(diff = mean_post/mean_pre)
median(diff$diff - 1)
quantile(diff$diff,0.025)-1
quantile(diff$diff,0.975)-1

#plot model output-------------------------------------------------------------
params1 <- rstan::extract(presence_mod)
TT <- 22
iter <- 8000

#convert model output to real estimates
est_probs<- list()
#split these calculations into 3 separate for loops since not all datasets
#appear in all years and need NAs for years they don't appear
for(i in 1:22){
  presence_coef <- data.frame(obs_prob_r=NA, 
                          state_prob=NA, 
                          iter=seq(1,iter))
  presence_coef[,1] <- plogis(params1$a_yr1[,i]) #convert from link space to real
  presence_coef[,2] <- plogis(params1$x[,i])
  
  est_probs[[i]] = presence_coef
}

#extract the median and 95% CI for each year in each dataset
x_mat<- data.frame(median.iucn.x=NA,l.95.x=NA,u.95.x=NA)
for(i in 1:TT){
  x_mat[i,1]=median(est_probs[[i]]$state_prob)
  x_mat[i,2]=quantile(est_probs[[i]]$state_prob,0.025)
  x_mat[i,3]=quantile(est_probs[[i]]$state_prob,0.975)
}

r_mat<- data.frame(median.iucn.r=NA,l.95.r=NA,u.95.r=NA)
for(i in 1:TT){
  r_mat[i,1]=median(est_probs[[i]]$obs_prob_r)
  r_mat[i,2]=quantile(est_probs[[i]]$obs_prob_r,0.025, na.rm=TRUE)
  r_mat[i,3]=quantile(est_probs[[i]]$obs_prob_r,0.975, na.rm=TRUE)
}

#and combine all 4
dat <-  bind_cols(x_mat, r_mat) %>% 
  mutate(year = seq(from=2000, to = 2021, by = 1))
num_surveys <- df %>% 
  group_by(year) %>% 
  summarize(n = n())

#plot!
ggplot(data = dat) +
  geom_ribbon(aes(x = year, y = median.iucn.x, ymin = l.95.x, ymax = u.95.x),
              fill='darkcyan', alpha = 0.2) +
  geom_line(aes(year, median.iucn.x, colour = "Estimated\nstate"), 
            lty=5,lwd=0.8) +
  geom_line(aes(year, median.iucn.r, colour = "Observed"), lwd=0.8) +
  geom_point(aes(year, median.iucn.r, fill = "Observed"), 
             col='white',pch=21,cex=3) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_colour_manual(name = "", values = c("Estimated\nstate" = 'darkcyan', 
                                            "Observed" = "dodgerblue4"))+
  scale_fill_manual(name = "", values = c("Estimated\nstate" = 'darkcyan', 
                                          "Observed" ='dodgerblue4'),
                    guide = "none") +
  geom_text(data = num_surveys, aes(y=rep(0,22), 
                                    x=year, label = n), 
            size = 2.5) +
  labs(x = "Year", y = "Sighting frequency") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,1)
#ggsave("figures/combined_presence_model.png", width = 8, height = 6)  


#raw trends for each dataset---------------------------------------------------
iucn_sources <- read_csv("../data/cleaned_iucn_data.csv") %>% 
  transmute(year = event_year, 
            presence = pres_abs,
            site = as.character(location_code_numeric),
            source = source,
            time = b_time,
            area = area)
df_sources <- bind_rows(reef,ow,iucn_sources,dfo)

source_summary <- df_sources %>% 
  group_by(source, year) %>% 
  summarize(mean_pres = mean(presence),
            n = n()) %>% 
  ungroup()

source_labs <- c(abalone = "DFO-Abalone", 
                 bhm = "DFO-Benthic Habitat Modeling", 
                 CCIRA_Dive_BC = "CCIRA",
                 cote = "SFU-Cote", 
                 Hakai_Dive_BC = "Hakai", 
                 Lee_Dive_BC = "SFU-Lee", 
                 `MARINe_CitSciDive_AK-BC`="MARINe",
                 multi="DFO-Multispecies",
                 OW="OW", 
                 ParksCanada_Dive_Haida="Gwaii Haanas",
                 REEF="REEF",
                 Salomon_Dive_HaidaGwaii="SFU-Salomon",
                 Schultz_Dive_BC="SFU-Schultz",
                 Watson_Dive_Vanls="VIU")

ggplot(source_summary, aes(year, mean_pres, size = n)) +
  geom_point() +
  facet_wrap(~source,
             labeller = as_labeller(source_labs),
             ncol=3) +
  theme_bw() +
  labs(x = "Year", y = "Mean sighting frequency",
       size = "Number of\nsurveys")
#ggsave("figures/source_frequency_trends.png", width = 7.5, height = 8)  
