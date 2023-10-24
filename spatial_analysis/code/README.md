## *Pycnopodia helianthoides* spatial data dictionary

**Author:** Steven Brownlee <br>
**Contact:** steven.fr.brownlee@gmail.com

# Note for the user

Data dictionary describing all of the spatial datasets used in the spatial component 
of the COSEWIC assessment and in generating maps. 

# File list

**full_pycno_bycatch_presence_only.csv**

Non-spatial CSV file consisting of bycatch data generated from a variety of sources
for *Pycnopodia helianthoides*. See 'source' field and non-spatial analysis documentation
for further information. 'lat' and 'lon' fields encoded as GPS coordinates, WGS 1984. 

Note for user: some coordinates are missing or otherwise invalid (ie on land). Analysis
script filters for these cases, but if the user is interested in running their own analyses
they will need to be removed.

**pycno_dive_locations_full_final.csv**

Non-spatial CSV file consisting of dive survey data generated from a variety of sources
for *Pycnopodia helianthoides*.  See 'source' field and non-spatial analysis documentation
for further information. 'lat' and 'lon' fields encoded as GPS coordinates, WGS 1984.

Note for user: some coordinates are missing or otherwise invalid (ie on land). Analysis
script filters for these cases, but if the user is interested in running their own analyses
they will need to be removed.

**bc_eez_bc_albers.shp**

ESRI shapefile containing polygons delineating the Canadian marine EEZ on the 
Pacific coast, as well as the shorelines of the coastline and islands themselves. 
Downloaded from marineregions.org and clipped from the overall Canadian EEZ to 
only include the Pacific coast. CRS: EPSG:3005.					(https://www.marineregions.org/gazetteer.php?p=details&id=8493)

**bc_eez_boundary.shp**

ESRI shapefile consisting of a square shape centered on the coastline of BC.
Generated in QGIS to facilitate the clipping of the Canada-wide EEZ to only
include the BC coastine. CRS: EPSG:3005.

**pre_wasting_eez_clip_template.shp**

ESRI shapefile generated in QGIS to clip the portions of the convex hull 
generated for pre-wasting occurrences/presences that are outside of the Canadian 
EEZ. CRS: EPSG:3005.

**post_wasting_eez_clip_template.shp**

ESRI shapefile generated in QGIS to clip the portions of the convex hull 
generated for post-wasting occurrences/presences that are outside of the Canadian 
EEZ. CRS: EPSG:3005.

**total_hull_eez_template.shp**

ESRI shapefile generated in QGIS to clip the portions of the convex hull 
generated for pre- and post-wasting occurrences/presences that are outside of 
the Canadian EEZ. CRS: EPSG:3005.

**bc_eez_100m1.tif**

File too large to upload to Github, can source from below or contact author for a copy.
.tif raster file containing bathymetry from ETOPO1 clipped to Canada’s marine 
EEZ on the Pacific coast. 100 m spatial resolution. 
Sourced from the BC Marine Conservation Analysis (BCMA) website. CRS: EPSG: 3005				(https://bcmca.ca/data/eco_physical_bathymetry/)

**na_10m_lakes_bc_env_albs.shp**

ESRI shapefile adapted from the 'Natural Earth' dataset of physical vector files. 
Consists of lake polygons and reprojected into shared coordinate system.
CRS: EPSG: 3005
(https://www.naturalearthdata.com/downloads/10m-physical-vectors/)

**na_10m_lakes_northam_bc_env_albs.shp**

ESRI shapefile adapted from the 'Natural Earth' dataset of physical vector files.
Supplement of additional lakes and rivers for North America, reprojected into shared coordinate system. CRS: EPSG: 3005
(https://www.naturalearthdata.com/downloads/10m-physical-vectors/)

**poli_boundaries_bc_env_albs.shp**

ESRI shapefile adapted from the 'Natural Earth' dataset of cultural vector files.
Consists of polygons delineating province-level political units in North and Central
America. Reprojected into shared coordinate system. CRS: EPSG: 3005
(https://www.naturalearthdata.com/downloads/10m-cultural-vectors/)

**pycno_canadian_range_455m.shp**

ESRI shapefile generated in QGIS. Consists of polygon delineating the total range 
of *Pycnopodia helianthoides* in Canada, with a cutoff of all marine habitat under
444 m in depth, the recorded depth limit for the species. Exported in shared
coordinate system. CRS: EPSG: 3005
