# Load libraries
library(tidyverse)
library(sf)
library(data.table)
library(rgdal)
library(rgeos)
library(maps)
library(maptools)
library(raster)

# Pacific centered world basemap ----
# Read world basemap
NE_countries <- readOGR(dsn = "./ne_110m_land", layer = "ne_110m_land")

# Split world map by "split line"

# shift central/prime meridian towards west - positive values only
shift <- 180+20

# create "split line" to split country polygons
WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
split.line = SpatialLines(list(Lines(list(Line(cbind(180-shift,c(-90,90)))), ID="line")), 
                          proj4string=WGS84)

# intersecting line with country polygons
line.gInt <- gIntersection(split.line, NE_countries)

# create a very thin polygon (buffer) out of the intersecting "split line"
bf <- gBuffer(line.gInt, byid=TRUE, width=0.01)  

# split country polygons using intersecting thin polygon (buffer)
NE_countries.split <- gDifference(NE_countries, bf, byid=TRUE)
# plot(NE_countries.split) # check map
class(NE_countries.split) # is a SpatialPolygons object

# Create graticules

# create a bounding box - world extent
b.box <- as(raster::extent(-180, 180, -90, 90), "SpatialPolygons")
# assign CRS to box
proj4string(b.box) <- WGS84
# create graticules/grid lines from box
grid <- gridlines(b.box, 
                  easts  = seq(from=-180, to=180, by=20),
                  norths = seq(from=-90, to=90, by=10))

# create labels for graticules
grid.lbl <- labels(grid, side = 1:4)

# transform labels from SpatialPointsDataFrame to a data table that ggplot can use
grid.lbl.DT <- data.table(grid.lbl@coords, grid.lbl@data)

# prepare labels with regular expression:
# - delete unwanted labels
grid.lbl.DT[, labels := gsub(pattern="180\\*degree|90\\*degree\\*N|90\\*degree\\*S", replacement="", x=labels)]
# - replace pattern "*degree" with "�" (* needs to be escaped with \\)
grid.lbl.DT[, lbl := gsub(pattern="\\*degree", replacement="�", x=labels)]
# - delete any remaining "*"
grid.lbl.DT[, lbl := gsub(pattern="*\\*", replacement="", x=lbl)]

# adjust coordinates of labels so that they fit inside the globe
grid.lbl.DT[, long := ifelse(coords.x1 %in% c(-180,180), coords.x1*175/180, coords.x1)]
grid.lbl.DT[, lat  := ifelse(coords.x2 %in% c(-90,90), coords.x2*82/90, coords.x2)]

# Prepare data for ggplot, shift & project coordinates

# give the PORJ.4 string for WGS84
PROJ <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

# transform graticules from SpatialLines to a data table that ggplot can use
grid.DT <- data.table(map_data(SpatialLinesDataFrame(sl=grid, 
                                                     data=data.frame(1:length(grid)), 
                                                     match.ID = FALSE)))
# project coordinates
# assign matrix of projected coordinates as two columns in data table
grid.DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]

# project coordinates of labels
grid.lbl.DT[, c("X","Y") := data.table(project(cbind(long, lat), proj=PROJ))]

# transform split country polygons in a data table that ggplot can use
Country.DT <- data.table(map_data(as(NE_countries.split, "SpatialPolygonsDataFrame")))
# Shift coordinates
Country.DT[, long.new := long + shift]
Country.DT[, long.new := ifelse(long.new > 180, long.new-360, long.new)]
# project coordinates 
Country.DT[, c("X","Y") := data.table(project(cbind(long.new, lat), proj=PROJ))]

# Project world basemap ----
worldMap <- Country.DT %>% 
  st_as_sf(coords = c('long.new', 'lat'), 
           crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160') %>%
  group_by(group) %>%
  summarise(do_union = FALSE) %>%
  st_cast('POLYGON') %>%
  st_union()
worldMap_robin <- worldMap %>% st_transform(crs = '+proj=robin')
# Read Graticules and project ----
grats <- st_read('./ne_110m_graticules_all/ne_110m_graticules_15.shp', 
                 crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160')
grats_robin <- st_transform(grats, crs = '+proj=robin')
# Read Bounding box and project ----
bbox <- st_read('./ne_110m_graticules_all/ne_110m_wgs84_bounding_box.shp', 
                crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160')
bbox_robin <- st_transform(bbox, '+proj=robin')

# Read global heat flow database (IHFC 2018) ----
global_hf <- st_read('./IHFC_global_dataset/HeatFlowIHFC_global_original_2018.shp', 
                     coords = c(3,2), 
                     crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160')
# Clean data
# Remove useless columns and rows with 0 heat flow
global_hf <- global_hf %>% dplyr::select(-Data_Numbe, -Codes, -No__Temps, -No_Heat_Pr, 
                                  -Heat_Prod_, -No__sites, -Year_of_Pu, -Comments, 
                                  -F21, -F22, -F23, -F24, -F25, -F26, -ben, -europe) %>%
  filter(Heat_Flow > 0)
# Filter dataset for heat flow values between 0 and 200 mW/m^2
global_hf_filter <- filter(global_hf, Heat_Flow <= 200 & Heat_Flow >= 0)
# Robinson projection to get units into meters, rather than lat/long
global_hf_robin <- st_transform(global_hf, crs = '+proj=robin +lon_0=160')
# Filter dataset for heat flow values between 0 and 200 mW/m^2
global_hf_robin_filter <- filter(global_hf_robin, Heat_Flow <= 200 & Heat_Flow >= 0)

# Read in plate boundary data ----
# P. Bird dataset (Bird, 2003, G3)
pb_bound <- st_read('./pbird_2003_boundaries/tectonicplates-master/PB2002_boundaries.shp', 
                    crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160')
# Clean up data
# Rename boundaries to be consistent and remove useless columns
pb_bound$Name <- gsub('\\\\', '/', pb_bound$Name)
pb_bound$Name <- gsub('-', '/', pb_bound$Name)
pb_bound <- arrange(pb_bound, Name) %>%
  dplyr::select(-LAYER) %>%
  rename(Boundary = Name) 
# Robinson projection to get units into meters, rather than lat/long
pb_bound_robin <- st_transform(pb_bound, crs = '+proj=robin +lon_0=160')
plot(worldMap_robin)
plot(pb_bound_robin, col = 'blue', add = TRUE)
pb_bound_robin_sz <- filter(pb_bound_robin, Type == 'subduction') %>%
  mutate('segment_index' = 1:65) %>%
  arrange(segment_index) %>%
  dplyr::select(-Type)
# Draw 100km buffers around plate boundaries
pb_bound_robin_buffer <- st_buffer(pb_bound_robin, 1000000)
pb_bound_robin_buffer_sz <- filter(pb_bound_robin_buffer, Type == 'subduction') %>%
  mutate('segment_index' = 1:65) %>%
  arrange(segment_index) %>%
  dplyr::select(-Type)
# Subset heat flow data that lie within (inner join) the buffer
hf_in_pb_buffer <- st_join(global_hf_robin_filter, left = FALSE, 
                           pb_bound_robin_buffer[c('Boundary','Type', 'PlateA', 'PlateB')])
hf_in_pb_buffer_sz <- st_join(global_hf_robin_filter, left = FALSE,
                              pb_bound_robin_buffer_sz[
                                c('Boundary','PlateA','PlateB','segment_index')
                                ]) %>%
  arrange(segment_index)

# USGS Plate Boundary KMZ (added TUESDAY, AUGUST 27, 2019; ----
# https://www.usgs.gov/media/files/plate-boundaries-kmz-file)
usgs_bound <- st_read('./USGS_boundaries/plate-boundaries.kml',
                      crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160')
# Clean up data
# Relabel descriptions from raw HTML to plain text
usgs_bound$Description[grepl('Transform', usgs_bound$Description)] <- 'transform'
usgs_bound$Description[grepl('Divergent', usgs_bound$Description)] <- 'divergent'
usgs_bound$Description[grepl('Other', usgs_bound$Description)] <- 'other'
usgs_bound$Description[grepl('Convergent', usgs_bound$Description)] <- 'convergent'
usgs_bound <- usgs_bound %>% rename(Type = Description, Boundary = Name)
# Robinson projection to get units into meters, rather than lat/long
usgs_bound_robin <- usgs_bound %>% 
  st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")) %>% 
  st_transform('+proj=robin +lon_0=160')
usgs_bound_robin_sz <- filter(usgs_bound_robin, Type == 'convergent') %>%
  arrange(Boundary) %>%
  mutate('segment_index' = 1:66) %>%
  arrange(segment_index) %>%
  dplyr::select(-Type)
# Draw 100km buffers around plate boundaries
usgs_bound_robin_buffer <- st_buffer(usgs_bound_robin, 1000000)
usgs_bound_robin_buffer_sz <- filter(usgs_bound_robin_buffer, Type == 'convergent') %>%
  arrange(Boundary) %>%
  mutate('segment_index' = 1:66) %>%
  arrange(segment_index) %>%
  dplyr::select(-Type)
# Subset heat flow data that lie within (inner join) the buffer
hf_in_usgs_buffer <- st_join(global_hf_robin_filter, left = FALSE, 
                             usgs_bound_robin_buffer[c('Boundary','Type')])
hf_in_usgs_buffer_sz <- st_join(global_hf_robin_filter, left = FALSE, 
                                usgs_bound_robin_buffer_sz[c('Boundary', 'segment_index')]) %>%
  arrange(segment_index)

# University of Texas Institute for Geophysics "PLATES" dataset ----
utig_bound <- st_read('./UTIG_PLATES_boundaries/UTIG_trench.kml',
                      crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_0=160') %>%
  dplyr::select(-Description) %>%
  rename(Type = Name)
# Robinson projection to get units into meters, rather than lat/long
utig_bound_robin <- st_transform(utig_bound, crs = '+proj=robin +lon_0=160')
utig_bound_robin_sz <- utig_bound_robin %>% 
  filter(Type == 'TR') %>%
  mutate(segment_index = 1:87) %>%
  arrange(segment_index)
# Draw 100km buffers around plate boundaries
utig_bound_robin_buffer <- st_buffer(utig_bound_robin, 1000000)
utig_bound_robin_buffer_sz <- utig_bound_robin_buffer %>%
  filter(Type == 'TR') %>%
  mutate(segment_index = 1:87) %>%
  arrange(segment_index)
# Subset heat flow data that lie within (inner join) the buffer
hf_in_utig_buffer <- st_join(global_hf_robin_filter, left = FALSE, 
                             utig_bound_robin_buffer[c('Type')])
hf_in_utig_buffer_sz <- st_join(global_hf_robin_filter, left = FALSE, 
                                utig_bound_robin_buffer_sz[c('Type', 'segment_index')]) %>%
  arrange(segment_index)

# Read in arc segment data from Syracuse and Abers (2006) ----
files <- list.files('./syracuse_abers_2006/gmts/')
for(i in 1:length(files)){
  fname <- files[i]
  tmp <- st_read(paste0('./syracuse_abers_2006/gmts/', fname), crs = 4326)
  # Robinson projection and set dateline wrapping for lines that cross ± 180 long
  tmp <- tmp %>% st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")) %>%
    st_transform(crs = '+proj=robin +lon_0=160')
  # 1000 km buffer
  tmpbuf <- tmp[[1]][1] %>% st_buffer(dist = 1000000)
  # Subset heat flow data that lie within (inner join) the buffer
  tmphf <- st_join(global_hf_robin_filter, left = FALSE, st_as_sf(tmpbuf))
  # Assign line and buffer objects to segment names
  assign(str_replace(fname, '^*.gmt', ''), tmp)
  assign(str_replace(fname, '^*.gmt', '_buffer'), tmpbuf)
  assign(str_replace(fname, '^*.gmt', '_hf'), tmphf)
}
# Clean up environment
rm(tmp, tmpbuf, i, files, fname)
# Create list of segments
segs <- list(alaska.aleutians.contours, andes.contours, central.america.contours,
             kamchatka.marianas.contours, kyushu.ryukyu.contours,
             lesser.antilles.contours, n.philippines.contours, new.britain.solomon.contours,
             s.philippines.contours, scotia.contours, sumatra.banda.sea.contours,
             tonga.new.zealand.contours, vanuatu.contours)
# Create list of buffers
segs_buffer <- list(alaska.aleutians.contours_buffer, andes.contours_buffer, central.america.contours_buffer,
                    kamchatka.marianas.contours_buffer, kyushu.ryukyu.contours_buffer,
                    lesser.antilles.contours_buffer, n.philippines.contours_buffer, new.britain.solomon.contours_buffer,
                    s.philippines.contours_buffer, scotia.contours_buffer, sumatra.banda.sea.contours_buffer,
                    tonga.new.zealand.contours_buffer, vanuatu.contours_buffer)
# Create list of heat flow within buffers
segs_hf <- list(alaska.aleutians.contours_hf, andes.contours_hf, central.america.contours_hf,
                kamchatka.marianas.contours_hf, kyushu.ryukyu.contours_hf,
                lesser.antilles.contours_hf, n.philippines.contours_hf, new.britain.solomon.contours_hf,
                s.philippines.contours_hf, scotia.contours_hf, sumatra.banda.sea.contours_hf,
                tonga.new.zealand.contours_hf, vanuatu.contours_hf)

# Save Data ----
# Cleanup environment
rm(PROJ, shift, b.box, bf, Country.DT, grid, grid.lbl, grid.DT, grid.lbl.DT, line.gInt, NE_countries, NE_countries.split, split.line, tmphf, WGS84)
save.image('backarc_heatflow_data.RData')
