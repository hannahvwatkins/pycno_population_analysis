#formating data for spatial analyses
library(tidyverse)
library(sp)
library(janitor)
library(measurements)

locations <- read_csv("../data/REEF_locations.csv") %>% 
  filter(lat != "NULL") %>% 
  #this function is being annoying so this is going to be hella janky
  separate(lat, into = c("lat_degree","lat_minute","lat_decimal")) %>% 
  separate(lon, into = c(NA,"lon_degree","lon_minute","lon_decimal")) %>% 
  unite("lat_minute", lat_minute:lat_decimal, sep = ".") %>% 
  unite("lon_minute", lon_minute:lon_decimal, sep = ".") %>% 
  mutate(lat_min = as.numeric(lat_minute)/60,
         lat = as.numeric(lat_degree) + lat_min,
         lon_min = as.numeric(lon_minute)/60,
         lon = (as.numeric(lon_degree) + lon_min)*-1,
         geogr = geogid) %>% 
  select(lat, lon, geogr)
         
reef <- read_csv("../data/cleaned_reef_data.csv") %>%  
  #since not many sites were surveyed more than once on the same day, including
  #a random effect of site_day is really challenging and the model can't
  #converge properly
  #so we're just going to take the highest score for a given site_day to avoid
  #pseudoreplication, since it's more likely that one diver missed individuals 
  #than that the other falsely counted more
#however, we already did that in the cleaning code, so there's no need to do it
#again
  #group_by(site_ymd) %>% 
  #mutate(num_surveys = n()) %>% 
  #slice_max(abundance) %>% 
  #ungroup() %>% 
  left_join(locations, by = "geogr") %>% 
  #turn into counts using the lowest value in each bin as a conservative 
  #estimate
  mutate(count = case_when(abundance == 0 ~ 0,
                           abundance == 1 ~ 1,
                           abundance == 2 ~ 2,
                           abundance == 3 ~ 11,
                           abundance == 4 ~ 101)) %>% 
  group_by(lat, lon, year) %>% 
  summarize(num_visits = n(),
            mean_count = mean(count, na.rm = TRUE),
            median_count = median(count, na.rm = TRUE),
            presence = max(presence)) %>% 
  ungroup() %>% 
  mutate(source = "REEF") 

ow_2001 <- read_csv("../data/cleaned_oceanwise_data.csv",
               guess_max = 10000) %>% 
  rename(lat = latitude, lon = longitude) %>% 
  #turn into counts using the lowest value in each bin as a conservative 
  #estimate
  mutate(count = case_when(abundance_recode == 1 ~ 0,
                           abundance_recode == 2 ~ 1,
                           abundance_recode == 3 ~ 11,
                           abundance_recode == 4 ~ 26,
                           abundance_recode == 5 ~ 51,
                           abundance_recode == 6 ~ 101,
                           abundance_recode == 7 ~ 1001)) %>% 
  group_by(lat, lon, year) %>% 
  summarize(num_visits = n(),
            mean_count = mean(count, na.rm = TRUE),
            median_count = median(count, na.rm = TRUE),
            presence = max(pres_abs)) %>% 
  ungroup() %>% 
  mutate(source = "OW") 

ow_early <- read_csv("../data/cleaned_oceanwise_data_for_spatial.csv",
                   guess_max = 10000) %>% 
  rename(lat = latitude, lon = longitude) %>% 
  filter(pres_abs == 1) %>% 
  group_by(lat, lon, year) %>% 
  summarize(num_visits = n(),
            presence = 1) %>% 
  ungroup() %>% 
  mutate(source = "OW_pre_2001") %>% 
  filter(year < 2001) 

dfo <- read_csv("../data/cleaned_dfo_dives.csv", guess_max = 3000) %>% 
  filter(source != "abalone") %>% 
  group_by(source, year, lat, lon) %>% 
  summarize(num_visits = n(),
            mean_count = mean(num_pycno),
            median_count = median(num_pycno, na.rm = TRUE),
            presence = max(presence))

iucn <- read_csv("../data/cleaned_iucn_data_for_spatial.csv") %>% 
  group_by(source, event_year, dec_lat, dec_long) %>% 
  summarize(num_visits = n(),
            mean_count = mean(totalpycnos, na.rm = TRUE),
            median_count = median(totalpycnos, na.rm = TRUE),
            presence = max(pres_abs)) %>% 
  arrange(event_year) %>% 
  ungroup() %>% 
  rename(lat = dec_lat,
         lon = dec_long,
         year = event_year) %>% 
  #fix up the NaNs that form instead of NAs
  mutate(mean_count = na_if(mean_count, 'NaN'))

#we'll make a separate dataset with the iNat observations since these just
#presence only data and shouldn't be used in the occupancy models
inat <- read_csv("../data/inat_pycnos.csv") %>% 
  transmute(lat = latitude,
            lon = longitude,
            year = lubridate::year(observed_on),
            source = "iNaturalist",
            totalpycnos = 1,
            presence = 1) %>% 
  group_by(source, year, lat, lon) %>% 
  summarize(num_visits = n(),
            mean_count = mean(totalpycnos, na.rm = TRUE),
            median_count = median(totalpycnos, na.rm = TRUE),
            presence = max(presence))

#and another for the MARINe observations, which are presence/absence
marine <- read_csv("../data/cleaned_marine_global.csv") %>% 
  group_by(source, event_year, dec_lat, dec_long) %>% 
  summarize(num_visits = n(),
            presence = max(pres_abs)) %>% 
  arrange(event_year) %>% 
  ungroup() %>% 
  rename(lat = dec_lat,
         lon = dec_long,
         year = event_year)


latlon <- iucn %>% 
  bind_rows(reef) %>%
  bind_rows(ow_2001) %>% 
  bind_rows(ow_early) %>% 
  bind_rows(dfo) %>% 
  bind_rows(inat) %>% 
  bind_rows(marine) %>% 
  distinct()
#write_csv(latlon, "../data/pycno_dive_locations_full.csv")
sample_size_spatial_full <- latlon %>% group_by(source) %>% summarize(n = n())


latlon_counts <- latlon %>% 
  filter(!is.na(mean_count))
#write_csv(latlon_counts, "../data/pycno_dive_locations_counts.csv")
sample_size_counts <- latlon_counts %>% group_by(source) %>% summarize(n = n())


#bycatch data
#these are observations from bycatch surveys so there won't be zeros included
#here, just observations, locations, and depths
crab <- read_csv("../data/dfo_crab.csv") %>% 
  janitor::clean_names() %>% 
  filter(species == "4XE") %>% 
  #convert lat/lon minutes to degrees
  mutate(lat_deg = start_lat_min/60,
         lon_deg = start_long_min/60,
         lat = start_lat_deg + lat_deg,
         lon = start_long_deg + lon_deg,
         presence = 1,
         source = "crab_bycatch",
         depth = case_when(depth_unit == 'F' ~ max_depth/3.28084,
                                            TRUE ~ max_depth)) %>% 
  dplyr::select(source, year, lat, lon, presence, depth) #depth in m

prawn <- read_csv("../data/dfo_prawn.csv") %>% 
  janitor::clean_names() %>% 
  filter(species == '4XE') %>% 
  #convert lat/lon minutes to degrees
  mutate(lat_deg = start_lat_min/60,
         lon_deg = start_long_min/60,
         lat = start_lat_deg + lat_deg,
         lon = start_long_deg + lon_deg,
         presence = 1,
         source = "prawn_bycatch",
         #convert the depths measures in feet to m
         depth = case_when(depth_unit == 'F' ~ max_depth/3.28084,
                           TRUE ~ max_depth)) %>% 
  dplyr::select(source, year, lat, lon, presence, depth)

groundfish <- read_csv("../data/dfo_groundfish_research_db.csv") %>% 
  janitor::clean_names() %>% 
  filter(species_code == '4XE') %>% 
  transmute(source = "groundfish_bycatch",
            year = year,
            lat = latitude,
            lon = longitude,
            presence = 1,
            depth = depth_m)

longline <- read_csv("../data/dfo_groundfish_longline.csv") %>% 
  janitor::clean_names() %>% 
  filter(species_code == '4XE') %>% 
  filter(catch_count > 0) %>% 
  transmute(source = "groundfish_longline_bycatch",
            year = year,
            lat = latitude,
            lon = longitude,
            presence = 1,
            depth = depth_m)


bycatch <- bind_rows(crab, prawn, groundfish, longline)
#write_csv(bycatch, "../data/full_pycno_bycatch_presence_only.csv")

#the crab and longline data have zeros associated with them so we can also do
#a presence/absence dataframe
crab_full <- read_csv("../data/dfo_crab.csv")
crab_soak <- readxl::read_excel("../data/dfo_crab_soak_times.xlsx") %>% 
  left_join(crab_full, by = "H_key") %>% 
  janitor::clean_names() %>% 
  mutate(presence = case_when(species == "4XE" ~ 1,
                              TRUE ~ 0),
         lat_deg = start_lat_min/60,
         lon_deg = start_long_min/60,
         lat = start_lat_deg + lat_deg,
         lon = start_long_deg + lon_deg,
         depth = case_when(depth_unit == 'F' ~ max_depth/3.28084,
                                  TRUE ~ max_depth)) %>% 
  transmute(source = "crab_bycatch",
            year = year_x,
            lat = lat,
            lon = lon,
            presence = presence,
            depth = depth)
#jk the zeros don't have lat/lons
#just the longline then
longline_full <- read_csv("../data/dfo_groundfish_longline.csv") %>% 
  janitor::clean_names() %>% 
  filter(species_code == "4XE") %>% 
  mutate(presence = case_when(catch_count == 0 ~ 0,
                              TRUE ~ 1)) %>% 
  transmute(source = "groundfish_longline_bycatch",
            year = year,
            lat = latitude,
            lon = longitude,
            presence = presence,
            depth = depth_m)
#write_csv(longline_full, "../data/full_pycno_bycatch_presence_absence.csv")
