###

# R script concerning the spatial analyses for the 
# 2022 COSEWIC assessment report on Pycnopodia helianthoides

# Contact: Steven Brownlee (email: steven.fr.brownlee@gmail.com)
# Date completed: May 18 2023

###

# Note for users: see 'Spatial data dictionary' document for full descriptions
# of the files used and where to access them. 

##

# Optional installation code:

# install.packages(c('terra', 'tidyverse', 'sf', 'sp', 'PNWColors', 'MetBrewer',
#                  'patchwork'))

library(terra)
library(tidyverse)
library(sf)
library(sp)
library(PNWColors)
library(MetBrewer)
library(patchwork)

###

# Script requires an installation of QGIS and SAGA GIS on your computer in order
# to run functions from QGISProcess. For most users installing QGIS by itself
# from the official website (https://www.qgis.org/en/site/forusers/download.html)
# should be sufficient to meet all of the requirements. 

# Optional installation code:

# install.packages("remotes")
# remotes::install_github("r-spatial/qgisprocess")

library(qgisprocess)

# The package should automatically detect where your QGIS installation is on 
# your computer and configure the 'qgisprocess' R package accordingly. 

# If there are issues detecting your QGIS installation please see the package
# Github repo here for troubleshooting:

# https://github.com/r-spatial/qgisprocess

###

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###

# Note for users: This script assumes you have downloaded all of the required
# datasets into one directory in a known location. To run the script, simply
# uncomment the following setwd() command and fill it in with the path to your
# spatial data directory. 

#setwd('')

###

# Step 1: Extent of occurrence and occupancy workflow

###

## Area of occupancy

# Conduct joins to bring together bycatch and dive survey data. This data is 
# identical to that used in the main non-spatial analysis pipeline: see 
# appropriate documentation elsewhere for more information.

bycatch <- read_csv('full_pycno_bycatch_presence_only.csv')

dive_survey_inaturalist <- read_csv('pycno_dive_locations_full_final.csv')

# Remove depth information from bycatch data for cleanliness' sake

bycatch <- bycatch %>% select(-depth)

# Remove number of visits and mean counts from dive survey data

dive_survey_inaturalist <- dive_survey_inaturalist %>% select(-num_visits, 
                                                              -mean_count)
# Add identifying column to signify data source

bycatch <- bycatch %>% add_column(dataset = 'bycatch')

dive_survey_inaturalist <- dive_survey_inaturalist %>% 
  add_column(dataset = 'dive_surveys_and_inaturalist')

# Bind two bycatch and dive survey/iNaturalist datasets together

pycno_presence_dataset <- bind_rows(dive_survey_inaturalist, bycatch)

# Remove rows with missing coordinates

pycno_presence_dataset <- pycno_presence_dataset %>% 
  drop_na(lat)

# Remove rows where no pycnopodia were sighted

pycno_presence_dataset_zer_remov <- pycno_presence_dataset %>% 
  filter(presence != 0)

# Convert combined datasets into spatial objects

pycno_crs <- st_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

pycno_presence_sf <- st_as_sf(pycno_presence_dataset_zer_remov,
                              coords = c('lon', 'lat'),
                              crs = pycno_crs)
pycno_presence_sf_zeroes <- st_as_sf(pycno_presence_dataset,
                              coords = c('lon', 'lat'),
                              crs = pycno_crs)

# Note for users: Initial datasets' GPS coordinates are stored in WGS 1984,
# then projected into NAD 1983 BC Environment Albers (EPSG:3005) to 
# match auxiliary data.

# Read in marine EEZ shapefile

bc_eez_path <- 'bc_eez_bc_albers.shp'

bc_eez_sf <- read_sf(bc_eez_path)

# Read in boundary file used to remove non-BC EEZ from the marine EEZ file. 

# Note for users: This boundary file was created manually in QGIS to crop out
# the portions of the Canada-wide EEZ not relevant to the analysis at hand,
# and consists of a roughly square-shaped polygon centered around the BC coast.

bc_boundary_path <- 'bc_eez_boundary.shp'

bc_boundary <- read_sf(bc_boundary_path)

# Crop Canada-wide EEZ using boundary file

bc_eez_sf <- bc_eez_sf %>% st_crop(bc_boundary)

# Reproject occurrence data into EPSG:3005

bc_env_albs_crs <- '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 
+x_0=1000000 +y_0=0 +ellps=GRS80 +datpum=NAD83 +units=m +no_defs'

pycno_presence_sf <- pycno_presence_sf %>% st_transform(bc_env_albs_crs)

bc_eez_sf <- bc_eez_sf %>% st_transform(bc_env_albs_crs)

pycno_presence_sf_zeroes <- pycno_presence_sf_zeroes %>% st_transform(bc_env_albs_crs)

# Crop our combined occurrence datasets to remove records that have erroneous
# coordinates - specifically ones on land. 

# Note for users: No additional filtering for the correctness of coordinates 
# was conducted other than this simple check.

pycno_presence_sf <- pycno_presence_sf %>% st_intersection(bc_eez_sf)
pycno_presence_sf_zeroes <- pycno_presence_sf_zeroes %>% st_intersection(bc_eez_sf)

pycno_presence_summarized <- pycno_presence_sf %>% 
  group_by(source) %>% summarize(n = n())

pycno_presence_summarized_zeroes <- pycno_presence_sf_zeroes %>% 
  group_by(source) %>% summarize(n = n())


# Keep only necessary data columns, remove extraneous ones carried over from
# EEZ file.

pycno_presence_sf <- pycno_presence_sf %>% select(source,
                                                  year,
                                                  presence,
                                                  dataset)

pycno_presence_sf_zeroes <- pycno_presence_sf_zeroes %>% select(source,
                                                  year,
                                                  presence,
                                                  dataset)

# Create unique ID column for final filtering

pycno_presence_sf <- pycno_presence_sf %>% mutate(id = row_number())

pycno_presence_sf_zeroes <- pycno_presence_sf_zeroes %>% mutate(id = row_number())

###

# Note for users: quick visual inspection, optional. Run if you want to make sure 
# everything is working as it should!

# Simple plot of occurrence data by data source

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'white') +
  geom_sf(data = pycno_presence_sf, aes(colour = source))

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'white') +
  geom_sf(data = pycno_presence_sf_zeroes, aes(colour = source))

# Simple plot of occurrence data by dataset

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = pycno_presence_sf, aes(colour = year))

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = pycno_presence_sf, aes(colour = dataset))

# Simple plot of all survey effort - presence/absence

pycno_presence_sf_zeroes <- pycno_presence_sf_zeroes %>% mutate_at('presence',
                                                                   as.factor)

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = pycno_presence_sf_zeroes, aes(colour = presence))

###

# Remove row 2203 - single dive survey by OceanWise in 2003 to the 
# Bowie Seamount that was never repeated. Would artificially inflate
# pre-wasting area of occupancy.

# Note for users: Double check in QGIS or your preferred GIS program to ensure
# the correct coordinate was removed. Updates to the dataset or changes to the 
# composition of the one used here *will* change the row number you need to remove.

pycno_presence_sf <- pycno_presence_sf %>% filter(id != 2203)

pycno_presence_sf_zeroes <- pycno_presence_sf_zeroes %>% filter(id != 3143)

###

# Note for users: Use of zero-included dataset ends here - carried through to 
# final figure generation. 

# Split occurrence dataset into pre and post wasting subsets.

# 1972 - 2013: pre-wasting
# 2015 - 2021: post-wasting

pre_wasting_pycno <- pycno_presence_sf %>% filter(year < 2013)
post_wasting_pycno <- pycno_presence_sf %>% filter(year >= 2015)

# Note for users: Optional visual inspection step

# Pre-wasting

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = pre_wasting_pycno, aes(colour = year))

# Post-wasting

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = post_wasting_pycno, aes(colour = year))

# Construct convex hulls around pre- and post-wasting datasets - first convert
# from sf to SpatVector formats from Terra package

pre_wasting_pycno <- vect(pre_wasting_pycno)
post_wasting_pycno <- vect(post_wasting_pycno)

pre_wasting_hull <- convHull(pre_wasting_pycno)
post_wasting_hull <- convHull(post_wasting_pycno)

total_vect <- vect(pycno_presence_sf)
total_hull <- convHull(total_vect)

# Adjust convex hulls such that segments of hull outside of Canadian jurisdiction
# are removed - use traced segments drawn in QGIS as clipping templates.

# Note to users, *important*: These adjustments were done manually in QGIS by 
# creating polygons that correspond to the portions of the convex hulls that are 
# outside the Canadian EEZ in the Pacific. If you adapt this code to other datasets
# or filter the existing datasets differently than we have here, you *will* need to
# draw your own adjustment polygons in QGIS, as the convex hulls you generate will
# be different than the ones generated here! Contact the authors if you need 
# assistance with this step. 

pre_wasting_hull <- st_as_sf(pre_wasting_hull)
post_wasting_hull <- st_as_sf(post_wasting_hull)
total_hull <- st_as_sf(total_hull)

pre_wasting_adjustment <- read_sf('pre_wasting_eez_clip_template.shp')

post_wasting_adjustment <- read_sf('post_wasting_eez_clip_template.shp')

total_hull_adjustment <- read_sf('total_hull_eez_template.shp')

pre_wasting_hull_clipped <- qgis_run_algorithm(
  'native:difference',
  INPUT = pre_wasting_hull,
  OVERLAY = pre_wasting_adjustment 
)

post_wasting_hull_clipped <- qgis_run_algorithm(
  'native:difference',
  INPUT = post_wasting_hull,
  OVERLAY = post_wasting_adjustment 
)

total_hull_clipped <- qgis_run_algorithm(
  'native:difference',
  INPUT = total_hull,
  OVERLAY = total_hull_adjustment 
)

pre_wasting_hull_clipped <- sf::read_sf(qgis_output(pre_wasting_hull_clipped, "OUTPUT"))

post_wasting_hull_clipped <- sf::read_sf(qgis_output(post_wasting_hull_clipped, "OUTPUT"))

total_hull_clipped <- sf::read_sf(qgis_output(total_hull_clipped, "OUTPUT"))

# Convert to 'terra' formatting

pre_wasting_hull <- vect(pre_wasting_hull_clipped)

post_wasting_hull <- vect(post_wasting_hull_clipped)

total_hull <- vect(total_hull_clipped)

# Calculate area of pre- and post-wasting Extent of Occupancy (EOO)

pre_wasting_area <- expanse(pre_wasting_hull, unit = 'km', transform = TRUE)

post_wasting_area <- expanse(post_wasting_hull, unit = 'km', transform = TRUE)

total_area <- expanse(total_hull, unit = 'km', transform = TRUE)

###

# Note for users: Optional visual inspection step. 

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = pre_wasting_hull_clipped, aes(colour = 'blue'))

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = post_wasting_hull_clipped, aes(colour = 'blue'))

ggplot() +
  geom_sf(data = bc_eez_sf, fill = 'slategray2') +
  geom_sf(data = post_wasting_hull_clipped, aes(colour = 'blue'))


## 

# Note for users: final outputs of this step of the analysis. If you modified
# the filtering steps above or changed the source datasets, your results will
# vary from those noted here!

pre_wasting_area

# Pre-wasting area:  200547.9 sq. km

post_wasting_area

# Post-wasting area: 212516.7 sq. km

total_area

# Total area: 217731 sq. km 

# Note for user: the variable 'total_area' is the EOO reported in the manuscript.
# Pre- and post-wasting values were generated as part of an exploratory analysis
# that did not make it into final manuscript.

###

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###

## Step 2: Index of area of occupancy

# Read in bathymetric raster for extent template:

bathy_path <- 'bc_eez_100m1.tif'

bathy_eez <- rast(bathy_path)

bc_ext <- ext(bathy_eez)

# Create 2 km sq. grid for BC coast

iao_grid <- st_make_grid(bc_ext, cellsize = 2000, crs = bc_env_albs_crs,
                         square = TRUE)

iao_grid <- st_as_sf(iao_grid)

# Intersect grid cells that contain occurrences pre- and post-wasting

pre_wasting_pycno <- st_as_sf(pre_wasting_pycno)

post_wasting_pycno <- st_as_sf(post_wasting_pycno)

total_pycno <- pycno_presence_sf

pre_wasting_grid_iao <- qgis_run_algorithm(
  'native:extractbylocation',
  INPUT = iao_grid,
  PREDICATE = 'intersect',
  INTERSECT = pre_wasting_pycno
)

post_wasting_grid_iao <- qgis_run_algorithm(
  'native:extractbylocation',
  INPUT = iao_grid,
  PREDICATE = 'intersect',
  INTERSECT = post_wasting_pycno
)

total_pycno_grid_iao <- qgis_run_algorithm(
  'native:extractbylocation',
  INPUT = iao_grid,
  PREDICATE = 'intersect',
  INTERSECT = total_pycno
)

# Note for users: this step is made more convoluted by the necessity of 
# having to convert back and forth between different spatial data formats.

pre_wasting_grid_iao <- read_sf(qgis_output(pre_wasting_grid_iao, 'OUTPUT'))

post_wasting_grid_iao <- read_sf(qgis_output(post_wasting_grid_iao, 'OUTPUT'))

total_pycno_grid_iao <- read_sf(qgis_output(total_pycno_grid_iao, 'OUTPUT'))

# Calculate area of cells pre- and post-wasting

pre_wasting_vect_grid <- vect(pre_wasting_grid_iao)

post_wasting_vect_grid <- vect(post_wasting_grid_iao)

total_pycno_vect_grid <- vect(total_pycno_grid_iao)

#

iao_pre <- expanse(pre_wasting_vect_grid, unit = 'km', transform = TRUE)

iao_pre_sum <- sum(iao_pre)

iao_post <- expanse(post_wasting_vect_grid, unit = 'km', transform = TRUE)

iao_post_sum <- sum(iao_post)

iao_total <- expanse(total_pycno_vect_grid, unit = 'km', transform = TRUE)

iao_total_sum <- sum(iao_total)

#

# IAO pre-wasting: 4739.99

iao_pre_sum

# IAO post-wasting: 2323.999

iao_post_sum

# IAO total pycno dataset: 7651.999

iao_total_sum

# Note for user: IAO referenced in the manuscript is the 'iao_total_sum'
# variable. Pre- and post-wasting values were generated as part of an exploratory 
# analysis that did not make it into final manuscript.

###

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###

# Step 3: Change in 'occupancy' map generation workflow.

# Note for the user: an earlier version of this analysis incorporated two spatial
# scales: a 2 km sq. grid and a 10 km sq. grid, in order to test the sensitivity
# of the results to different units of aggregation. The 2 km sq. analyses are what 
# is reported in the manuscript, as that is the standard resolution for COSEWIC
# spatial analyses. The 10 sq. km results are reported here for completeness' sake
# and for the user's potential interest.

dive_survey_inaturalist_occ <- read_csv('pycno_dive_locations_full_final.csv')

###

# Remove extraneous columns and filter for data with missing coordinates.

dive_survey_inaturalist_occ <- dive_survey_inaturalist_occ %>% 
  select(-num_visits)

dive_survey_inaturalist_occ <- dive_survey_inaturalist_occ %>% 
  drop_na(lat)

dive_survey_inaturalist_occ <- dive_survey_inaturalist_occ %>% 
  mutate(id = row_number())

# Define coordinate system for dataset, specifically WGS 1984. 

pycno_crs <- st_crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

occ_dataset <- st_as_sf(dive_survey_inaturalist_occ,
                        coords = c('lon', 'lat'),
                        crs = pycno_crs)

# Project dataset into BC Environment Albers (EPSG:3005) to match auxiliary data

occ_dataset <- occ_dataset %>% st_transform(bc_env_albs_crs)

occ_dataset <- occ_dataset %>% st_intersection(bc_eez_sf)

# Remove Bowie Seamount occurrence, select only needed columns

occ_dataset <- occ_dataset %>% filter(id != 4173)

occ_dataset <- occ_dataset %>% select(source,
                                      year,
                                      median_count,
                                      presence)

# Filter into pre-wasting and post-wasting datasets.

occ_pre_wasting <- occ_dataset %>% filter(year < 2013)

occ_post_wasting <- occ_dataset %>% filter(year >= 2015)

# Retrieve only needed column, specifically median counts. 

occ_pre_wasting <- occ_pre_wasting %>% select(median_count)

occ_post_wasting <- occ_post_wasting %>% select(median_count)

# Rename columns to allow for the concatenation of the two datasets further
# in the script. 

occ_pre_wasting <- occ_pre_wasting %>% 
  rename(pre_wasting_median_count = median_count)

occ_post_wasting <- occ_post_wasting %>% 
  rename(post_wasting_median_count = median_count)

###

# 2 km occupancy workflow

###

# Create 2 km grid that encompasses BC

delta_occ_grid_2 <- st_make_grid(bc_ext, 
                                 cellsize = 2000, 
                                 crs = bc_env_albs_crs,
                                 square = TRUE)

delta_occ_grid_2 <- st_as_sf(delta_occ_grid_2)

# Transpose median counts from points into grid cells - where more than one count
# is transposed into a cell, an average is calculated and reported for the cell.

pre_wasting_grid_2 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_2,
  JOIN = occ_pre_wasting,
  SUMMARIES = 7
)

post_wasting_grid_2 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_2,
  JOIN = occ_post_wasting,
  SUMMARIES = 7
)

pre_wasting_grid_2 <- read_sf(qgis_output(pre_wasting_grid_2, 'OUTPUT'))
post_wasting_grid_2 <- read_sf(qgis_output(post_wasting_grid_2, 'OUTPUT'))

# Remove extraneous columns generated during the joining process. 

pre_wasting_grid_2 <- pre_wasting_grid_2 %>% 
  select(-fid_median) %>% 
  rename(pre_wasting_median = pre_wasting_median_count_median)

post_wasting_grid_2 <- post_wasting_grid_2 %>% 
  select(-fid_median) %>% 
  rename(post_wasting_median = post_wasting_median_count_median)

# Join pre- and post-wasting grids together.

comb_wasting_grid_2 <- qgis_run_algorithm(
  'native:joinattributesbylocation',
  PREDICATE = 'equal',
  INPUT = pre_wasting_grid_2,
  JOIN = post_wasting_grid_2
)

comb_wasting_grid_2 <- read_sf(qgis_output(comb_wasting_grid_2, 'OUTPUT'))

# Remove extraneous columns generated during the joining process. 

comb_wasting_grid_2 <- comb_wasting_grid_2 %>% 
  select(-fid_2)

# Calculate change in median count

delta_median_grid_2 <- comb_wasting_grid_2 %>% 
  mutate(change_median_count  
         = post_wasting_median - pre_wasting_median)

# Drop NAs where sampling did not take place in both time periods.

delta_median_grid_2 <- delta_median_grid_2 %>% 
  drop_na(change_median_count)

# Split out outcomes into 'gain', 'loss' and 'no change' datasets. 

delta_median_pos_2 <- delta_median_grid_2 %>% 
  filter(change_median_count > 0)

delta_median_neg_2 <- delta_median_grid_2 %>% 
  filter(change_median_count < 0)

delta_median_zero_2 <- delta_median_grid_2 %>% 
  filter(change_median_count == 0)

# Calculate percent change pre- and post-wasting for each set of outcomes. 

delta_median_pos_2 <- delta_median_pos_2 %>% 
  mutate(pct_change_median_count  = 
           (change_median_count/post_wasting_median)*100)

delta_median_neg_2 <- delta_median_neg_2 %>% 
  mutate(pct_change_median_count  = 
           (change_median_count/pre_wasting_median)*100)

delta_median_zero_2 <- delta_median_zero_2 %>% 
  mutate(pct_change_median_count  = 0)

# Bind datasets back together

delta_median_comb_2 <- bind_rows(delta_median_pos_2, 
                               delta_median_neg_2, 
                               delta_median_zero_2)


delta_median_test_2 <- delta_median_comb_2 %>% select(pct_change_median_count)

delta_median_display_2 <- st_as_sf(delta_median_comb_2)

###

# 10 km occupancy

###

delta_occ_grid_10 <- st_make_grid(bc_ext, cellsize = 10000, crs = bc_env_albs_crs,
                                 square = TRUE)

delta_occ_grid_10 <- st_as_sf(delta_occ_grid_10)

# Transpose median counts from points into grid cells - where more than one count
# is transposed into a cell, an average is calculated and reported for the cell.

pre_wasting_grid_10 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_10,
  JOIN = occ_pre_wasting,
  SUMMARIES = 7
)

post_wasting_grid_10 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_10,
  JOIN = occ_post_wasting,
  SUMMARIES = 7
)

pre_wasting_grid_10<- read_sf(qgis_output(pre_wasting_grid_10, 'OUTPUT'))
post_wasting_grid_10 <- read_sf(qgis_output(post_wasting_grid_10, 'OUTPUT'))

# Remove extraneous columns generated during the joining process. 

pre_wasting_grid_10 <- pre_wasting_grid_10 %>% 
  select(-fid_median) %>% 
  rename(pre_wasting_median = pre_wasting_median_count_median)

post_wasting_grid_10 <- post_wasting_grid_10 %>% 
  select(-fid_median) %>% 
  rename(post_wasting_median = post_wasting_median_count_median)

# Join pre- and post-wasting grids together.

comb_wasting_grid_10 <- qgis_run_algorithm(
  'native:joinattributesbylocation',
  PREDICATE = 'equal',
  INPUT = pre_wasting_grid_10,
  JOIN = post_wasting_grid_10
)

comb_wasting_grid_10 <- read_sf(qgis_output(comb_wasting_grid_10, 'OUTPUT'))

# Remove extraneous columns generated during the joining process. 

comb_wasting_grid_10 <- comb_wasting_grid_10 %>% 
  select(-fid_2)

# Calculate change in median count

delta_median_grid_10 <- comb_wasting_grid_10 %>% 
  mutate(change_median_count  
         = post_wasting_median - pre_wasting_median)

# Drop NAs where sampling did not take place in both time periods.

delta_median_grid_10 <- delta_median_grid_10 %>% 
  drop_na(change_median_count)

# Split out outcomes into 'gain', 'loss' and 'no change' datasets. 

delta_median_pos_10 <- delta_median_grid_10 %>% 
  filter(change_median_count > 0)

delta_median_neg_10 <- delta_median_grid_10 %>% 
  filter(change_median_count < 0)

delta_median_zero_10 <- delta_median_grid_10 %>% 
  filter(change_median_count == 0)

# Calculate percent change pre- and post-wasting for each set of outcomes. 

delta_median_pos_10 <- delta_median_pos_10 %>% 
  mutate(pct_change_median_count  = 
           (change_median_count/post_wasting_median)*100)

delta_median_neg_10 <- delta_median_neg_10 %>% 
  mutate(pct_change_median_count  = 
           (change_median_count/pre_wasting_median)*100)

delta_median_zero_10 <- delta_median_zero_10 %>% 
  mutate(pct_change_median_count  = 0)

# Bind datasets back together

delta_median_comb_10 <- bind_rows(delta_median_pos_10, 
                                 delta_median_neg_10, 
                                 delta_median_zero_10)


delta_median_test_10 <- delta_median_comb_10 %>% select(pct_change_median_count)

delta_median_display_10 <- st_as_sf(delta_median_comb_10)

###

# Calculate change in presence-absence only. 

# Note for reader: this analysis incorporates both presence/absence only records
# *and* occurrences (ie counts). Occurrences are treated as presences.

pres_abs_dataset <- pycno_presence_sf_zeroes

# Split dataset into pre- and post-wasting datasets.

pres_pre_wasting <- pres_abs_dataset %>% filter(year < 2013)

pres_post_wasting <- pres_abs_dataset %>% filter(year >= 2015)

###

# 2 km sq. grid

###

# Intersect pre- and post-wasting presence/absence data with the template grid

pres_pre_wasting_grid_2 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_2,
  JOIN = pres_pre_wasting,
  SUMMARIES = 3
)

pres_post_wasting_grid_2 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_2,
  JOIN = pres_post_wasting,
  SUMMARIES = 3
)

# 10 km sq. grid

# Join pre- and post-wasting presence/absence data with the template grid

pres_pre_wasting_grid_10 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_10,
  JOIN = pres_pre_wasting,
  SUMMARIES = 3
)

pres_post_wasting_grid_10 <- qgis_run_algorithm(
  'qgis:joinbylocationsummary',
  PREDICATE = 'intersect',
  INPUT = delta_occ_grid_10,
  JOIN = pres_post_wasting ,
  SUMMARIES = 3
)

pres_pre_wasting_grid_2 <- read_sf(qgis_output(pres_pre_wasting_grid_2, 'OUTPUT'))
pres_post_wasting_grid_2 <- read_sf(qgis_output(pres_post_wasting_grid_2, 'OUTPUT'))

pres_pre_wasting_grid_10 <- read_sf(qgis_output(pres_pre_wasting_grid_10, 'OUTPUT'))
pres_post_wasting_grid_10 <- read_sf(qgis_output(pres_post_wasting_grid_10, 'OUTPUT'))

# Remove extraneous columns generated during joining process 

pres_pre_wasting_grid_2 <- pres_pre_wasting_grid_2 %>% 
  select(-fid_max) %>% 
  rename(pre_wasting_presence = presence_max)

pres_post_wasting_grid_2 <- pres_post_wasting_grid_2 %>% 
  select(-fid_max) %>% 
  rename(post_wasting_presence = presence_max)

pres_pre_wasting_grid_10 <- pres_pre_wasting_grid_10 %>% 
  select(-fid_max) %>% 
  rename(pre_wasting_presence = presence_max)

pres_post_wasting_grid_10 <- pres_post_wasting_grid_10 %>% 
  select(-fid_max) %>% 
  rename(post_wasting_presence = presence_max)

# Join pre- and post-wasting grids together at each spatial scale. 

pres_comb_wasting_grid_2 <- qgis_run_algorithm(
  'native:joinattributesbylocation',
  PREDICATE = 'equal',
  INPUT = pres_pre_wasting_grid_2,
  JOIN = pres_post_wasting_grid_2
)

pres_comb_wasting_grid_10 <- qgis_run_algorithm(
  'native:joinattributesbylocation',
  PREDICATE = 'equal',
  INPUT = pres_pre_wasting_grid_10,
  JOIN = pres_post_wasting_grid_10
)

pres_comb_wasting_grid_2 <- read_sf(qgis_output(pres_comb_wasting_grid_2, 'OUTPUT'))

pres_comb_wasting_grid_10 <- read_sf(qgis_output(pres_comb_wasting_grid_10, 'OUTPUT'))

# Remove extraneous columns generated during joining process.

pres_comb_wasting_grid_2 <- pres_comb_wasting_grid_2 %>% 
  select(-fid_2)

pres_comb_wasting_grid_10 <- pres_comb_wasting_grid_10 %>% 
  select(-fid_2)

# Convert presence data to 'numeric' data type. 

pres_comb_wasting_grid_2$pre_wasting_presence <- as.numeric(pres_comb_wasting_grid_2$pre_wasting_presence)
pres_comb_wasting_grid_2$post_wasting_presence <- as.numeric(pres_comb_wasting_grid_2$post_wasting_presence)

pres_comb_wasting_grid_10$pre_wasting_presence <- as.numeric(pres_comb_wasting_grid_10$pre_wasting_presence)
pres_comb_wasting_grid_10$post_wasting_presence <- as.numeric(pres_comb_wasting_grid_10$post_wasting_presence)

# Calculate change in presence at both spatial scales. 

pres_delta_grid_2 <- pres_comb_wasting_grid_2 %>% 
  mutate(change_presence = post_wasting_presence - pre_wasting_presence)

pres_delta_grid_10 <- pres_comb_wasting_grid_10 %>% 
  mutate(change_presence = post_wasting_presence - pre_wasting_presence)

pres_delta_grid_2 <- pres_delta_grid_2 %>% 
  select(pre_wasting_presence, post_wasting_presence, change_presence) %>% 
  drop_na(change_presence)

# Remove data where sampling did not occur in both time periods. 

pres_delta_grid_10 <- pres_delta_grid_10 %>% 
  select(pre_wasting_presence, post_wasting_presence, change_presence) %>% 
  drop_na(change_presence)

# Gather different outcomes and summarize. 

pres_abs_10_summ <- pres_delta_grid_10 %>% count(change_presence)

pres_abs_10_summ

# 10 sq. km
# n = 226, loss presence: 40, no change: 184, gain presence: 2

pres_abs_2_summ <- pres_delta_grid_2 %>% count(change_presence)

pres_abs_2_summ

# 2 km sq. km
# n = 300, loss presence: 92, no change: 202, gain presence: 6

# Code numeric outcomes into text for display purposes. 

pres_display_10 <- pres_delta_grid_10 %>% 
  mutate(change_presence = case_when(change_presence == -1 ~ 'loss of presence', 
                   change_presence == 0 ~ 'no change', 
                   change_presence == 1 ~ 'gain of presence'))

pres_display_2 <- pres_delta_grid_2 %>% 
  mutate(change_presence = case_when(change_presence == -1 ~ 'loss of presence', 
                                     change_presence == 0 ~ 'no change', 
                                     change_presence == 1 ~ 'gain of presence'))

pres_display_10$change_presence <- as.factor(pres_display_10$change_presence)

pres_display_10$change_presence <- factor(pres_display_10$change_presence,
                                              levels = c('loss of presence', 'no change',
                                              'gain of presence'))


pres_display_2$change_presence <- as.factor(pres_display_2$change_presence)

pres_display_2$change_presence <- factor(pres_display_2$change_presence,
                                             levels = c('loss of presence', 'no change',
                                              'gain of presence'))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###

# Figure generation section.

###

# Bring in background GIS files:

na_lakes_01_p = 'na_10m_lakes_bc_env_albs.shp'

na_lakes_02_p = 'na_10m_lakes_northam_bc_env_albs.shp'

poli_boundaries_p = 'poli_boundaries_bc_env_albs.shp'

bc_eez_template_p = 'bc_eez_bc_albers.shp'

canada_range = 'pycno_canadian_range_455m.shp'


na_lakes_01 <- read_sf(na_lakes_01_p)

na_lakes_02 <- read_sf(na_lakes_02_p)

poli_boundaries <- read_sf(poli_boundaries_p)

bc_eez_template <- read_sf(bc_eez_template_p)

coord_ext <- ext(bc_boundary)

coord_ext

bc_albs <- st_crs(bc_eez_template)

canada_range <- read_sf(canada_range)

###

# Figure 1b: Inset map of Pycno range in BC.

# Note for user: Annotations of placenames done via image editing software 
# post-generation. 

# Figure 1a taken from Gravem paper referenced in the manuscript. 

v <- ggplot() +
  theme_dark() +
  theme(legend.position = c(0.17, 0.2)) +
  geom_sf(data = canada_range,
          fill = 'goldenrod4',
          colour = 'grey21',
          alpha = 0.7) +
  geom_sf(data = bc_eez_template, 
          fill = 'grey89', 
          alpha = 0.2, 
          colour = 'gray15') +
  geom_sf(data = poli_boundaries, 
          fill = 'grey32', 
          colour = 'gray15') +
  geom_sf(data = na_lakes_01, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = na_lakes_02, 
          fill = 'grey50',
          colour = 'gray15') +
  coord_sf(xlim = c(1275058.50591579, 475058.7033865398),
           ylim = c(1206292.83107712, 266292.3782252073), 
           expand = FALSE)

ggsave(filename = 'figure_1b.png',
       plot = v,
       width = 12,
       height = 8,
       units = 'in')

###

# Figure 2 - overview of occurrence/presence dataset

pycno_presence_sf_zeroes_display <- pycno_presence_sf_zeroes 

pycno_presence_sf_zeroes_display$presence <- recode_factor(pycno_presence_sf_zeroes_display$presence, 
                                                  `0` = 'Absence', 
                                                  `1`= 'Presence')
pycno_presence_sf_zeroes_display <- pycno_presence_sf_zeroes_display %>% 
  arrange(desc(presence))

pre_wasting_pycno_zeroes <- pycno_presence_sf_zeroes_display %>% 
  filter(year < 2013)

post_wasting_pycno_zeroes <- pycno_presence_sf_zeroes_display %>% 
  filter(year >= 2015)

pycno_total_occ_pres <- pycno_presence_sf_zeroes_display %>% 
  filter(presence == 'Presence')


# Total dataset presence/absence

h <- ggplot() +
  theme_dark() +
  theme(legend.position = c(0.17, 0.2)) +
  theme(legend.title=element_blank()) +
  geom_sf(data = bc_eez_template, 
          fill = 'grey89', 
          alpha = 0.2, 
          colour = 'gray15') +
  geom_sf(data = poli_boundaries, 
          fill = 'grey32', 
          colour = 'gray15') +
  geom_sf(data = na_lakes_01, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = na_lakes_02, 
          fill = 'grey50',
          colour = 'gray15') + 
  geom_sf(data = pycno_presence_sf_zeroes_display, 
          aes(colour = presence)) +
  scale_colour_manual(values = c('sienna4', 'sienna1')) +
  coord_sf(xlim = c(1275058.50591579, 475058.7033865398),
           ylim = c(1206292.83107712, 266292.3782252073), 
           expand = FALSE) +
  ggtitle('Total presence/absence dataset \n 1972-2013')

ggsave(filename = 'figure_2_total_occurrence_presence.png',
       plot = h,
       width = 12,
       height = 8,
       units = 'in')

# Pre-wasting presence/absence


i <- ggplot() +
  theme_dark() +
  theme(legend.position = c(0.2, 0.2)) +
  theme(legend.title=element_blank()) +
  geom_sf(data = bc_eez_template, 
          fill = 'grey89', 
          alpha = 0.2, 
          colour = 'gray15') +
  geom_sf(data = poli_boundaries, 
          fill = 'grey32', 
          colour = 'gray15') +
  geom_sf(data = na_lakes_01, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = na_lakes_02, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = pre_wasting_pycno_zeroes, aes(colour = presence)) +
  scale_colour_manual(values = c('sienna4', 'sienna1')) +
  coord_sf(xlim = c(1275058.50591579, 475058.7033865398),
           ylim = c(1206292.83107712, 266292.3782252073), 
           expand = FALSE) +
  ggtitle('Pre-SSWD presences and absences \n 1972-2013')

ggsave(filename = 'figure_3a_pre_wasting_pres_abs.png',
       plot = i,
       width = 12,
       height = 8,
       units = 'in')

# Post-wasting presence/absence

j <- ggplot() +
  theme_dark() +
  theme(legend.position = c(0.2, 0.2)) +
  theme(legend.title=element_blank()) +
  geom_sf(data = bc_eez_template, 
          fill = 'grey89', 
          alpha = 0.2, 
          colour = 'gray15') +
  geom_sf(data = poli_boundaries, 
          fill = 'grey32', 
          colour = 'gray15') +
  geom_sf(data = na_lakes_01, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = na_lakes_02, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = post_wasting_pycno_zeroes, aes(colour = presence)) +
  scale_colour_manual(name = "Presence/Absence", values = c('sienna4', 'sienna1')) +
  coord_sf(xlim = c(1275058.50591579, 475058.7033865398),
           ylim = c(1206292.83107712, 266292.3782252073), 
           expand = FALSE) +
  ggtitle('Post-SSWD presences and absences \n 2015-2021')

ggsave(filename = 'figure_3b_post_wasting_pres_abs.png',
       plot = j,
       width = 12,
       height = 8,
       units = 'in')

###

k <- ggplot() +
  theme_dark() +
  theme(legend.position = c(0.2, 0.2)) +
  theme(legend.title=element_blank()) +
  geom_sf(data = bc_eez_template, 
          fill = 'grey89', 
          alpha = 0.2, 
          colour = 'gray15') +
  geom_sf(data = poli_boundaries, 
          fill = 'grey32', 
          colour = 'gray15') +
  geom_sf(data = na_lakes_01, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = na_lakes_02, 
          fill = 'grey50',
          colour = 'gray15') +
  geom_sf(data = total_hull_clipped,
          fill = 'darkgoldenrod4',
          alpha = 0.1) +
  geom_sf(data = pycno_total_occ_pres,
          colour = 'darkgoldenrod2') +
  coord_sf(xlim = c(1275058.50591579, 475058.7033865398),
           ylim = c(1206292.83107712, 266292.3782252073), 
           expand = FALSE) +
  ggtitle('Total occurrence and presence dataset \n 1972-2021')

ggsave(filename = 'figure_4_total_pres_occ_eoo.png',
       plot = k,
       width = 12,
       height = 8,
       units = 'in')

###

# Note to user: These figures did not make it into the final analysis, but 
# display the data cited in the analysis. They are presented here for completeness'
# sake and for the user's interest. 

###

# Figure 4 - comparative declines in median count, 2 sq. km vs 10 sq. km

d_2 <- ggplot(delta_median_display_2, aes(x = pct_change_median_count)) +
  theme_classic() +
  geom_histogram(color="darkblue", fill="lightblue",
                 linetype="solid",
                 bins = 20) +
  xlab('Percent change in median count') +
  ylab('Number of cells') +
  ggtitle('2 sq. km grid cell aggregation (n = 220)')

d_10 <- ggplot(delta_median_display_10, aes(x = pct_change_median_count)) +
  theme_classic() +
  geom_histogram(color="darkorange3", fill="orange",
                 linetype="solid",
                 bins = 20) +
  xlab('Percent change in median count') +
  ylab('Number of cells') +
  ggtitle('10 sq. km grid cell aggregation (n = 117)')

d_complete <- d_2 / d_10

d_complete

#

ggsave(filename = 'figure_x_spatial_comparison.png',
       plot = d_complete,
       width = 12,
       height = 8,
       units = 'in')

###

# Change in presence-absence, 2 sq. km vs 10 sq. km

d_a_2 <- ggplot(pres_display_2, aes(x = change_presence)) +
  theme_classic() +
  geom_bar(color="darkblue", fill="lightblue",
                 linetype="solid") +
  xlab('Change in presence') +
  ylab('Number of cells') +
  ggtitle('2 sq. km grid cell aggregation (n = 300)')

d_a_10 <- ggplot(pres_display_10, aes(x = change_presence)) +
  theme_classic() +
  geom_bar(color="darkorange3", fill="orange",
                 linetype="solid") +
  xlab('Change in presence') +
  ylab('Number of cells') +
  ggtitle('10 sq. km grid cell aggregation (n = 226)')

d_a_complete <- d_a_2 / d_a_10

d_a_complete

ggsave(filename = 'figure_x_spatial_comparison_pres_abs.png',
       plot = d_a_complete,
       width = 12,
       height = 8,
       units = 'in')

###

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###