# Load packages and functions
source('functions.R')

# Projections
# WGS84
proj4.wgs <- '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs'

# Robinson Pacific centered
proj4.robin.pacific <- "+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Print Projections
cat('Defining CRS (projections) ',
    'WGS:',
    proj4.wgs, 'Robinson Pacific Centered:',
    proj4.robin.pacific,
    sep = '\n')

# Country Boundaries
# Create sliver around twenty-five degree long to
# cut country boundaries so they don't project
# as straight lines across the globe
st_polygon(
  x = list(rbind(
    c(25.000001, 90),
    c(25, 90),
    c(25, -90),
    c(25.000001, -90),
    c(25.000001, 90)))) %>%
  st_sfc() %>%
  st_set_crs(proj4.wgs) -> shp.sliv

# Trim the countries boundaries overlapping at zero degree long and reproject
rnaturalearth::ne_countries(returnclass = 'sf') %>%
  st_difference(shp.sliv) %>%
  tibble::as_tibble() %>%
  st_as_sf() %>%
  st_transform(proj4.robin.pacific) -> shp.world.robin.pacific

# Read syracuse et al 2006 Segments
cat('\nLoading Syracuse and Abers (2006) segment boundaries ...')
c('Alaska_Aleutians',
  'Andes',
  'Central_America',
  'Kamchatka_Marianas',
  'Kyushu_Ryukyu',
  'Lesser_Antilles',
  'N_Philippines',
  'New_Britain_Solomon',
  'S_Philippines',
  'Scotia',
  'Sumatra_Banda_Sea',
  'Tonga_New_Zealand',
  'Vanuatu') -> seg.names
files <- list.files('data/sa2006/gmts', full.names = TRUE)
cat('\nSegments:', seg.names, sep = '\n')

# Contours Robinson pacific centered
purrr::map2_df(
  files,
  seg.names,
  ~read_latlong(
  .x,
  .y,
  proj4.robin.pacific)) -> shp.sa.contours.robin.pacific

# Filter for segment boundary
shp.sa.contours.robin.pacific %>%
  group_by(segment) %>%
  filter(row_number() == 1) -> shp.sa.segs.robin.pacific

# Buffer
buf.dist <- 1000000
cat('\nDrawing', buf.dist/1000, 'km buffers around segments ...')
shp.sa.segs.robin.pacific %>%
  st_buffer(buf.dist, endCapStyle = 'ROUND') %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES")) -> shp.sa.segs.robin.pacific.buffer

# Read global heat flow database (NGHF: Lucazeau 2019)
cat('\nReading global heat flow data (NGHF: Lucazeau 2019)',
    'Projecting to:',
    proj4.robin.pacific,
    sep = '\n')

# Read csv
hf <- read_csv('data/nghf/NGHF.csv', col_type = cols())

cat('\nRemoving NA heat flow values:',
    nrow(hf[is.na(hf$`heat-flow (mW/m2)`),]),
    '\nRemoving poor quality data (code6 = D):',
    nrow(hf[hf$code6 == 'D',]),
    '\nNo geographic info:',
    nrow(hf[is.na(hf$latitude) | is.na(hf$longitude),]))

# Filter
hf %>% 
  filter(!is.na(`heat-flow (mW/m2)`)) %>% # remove NA
  filter(code6 != 'D') %>% # remove poor quality data
  filter(!is.na(longitude)) %>% # remove no lat
  filter(!is.na(latitude)) %>% # remove no long
  st_as_sf(coords = c(9,10), # make simple feature object
           crs = proj4.wgs) %>% 
  st_transform(proj4.robin.pacific) -> shp.hf

# Check for duplicate measurements and remove
dup <- sp::zerodist(sf::as_Spatial(shp.hf))
cat('\nParsing duplicate locations:', length(dup),
    '\nIf x is better quality than y, keep x',
    '\notherwise randomly keep either x or y')

# For duplicate measurements, select the best quality data point,
# if the quality is the same, randomly select one measurement
rid <- vector('integer', nrow(dup))

for(i in 1:nrow(dup)) {
  if(shp.hf$code6[dup[i,1]] != shp.hf$code6[dup[i,2]] & shp.hf$code6[dup[i,1]] > shp.hf$code6[dup[i,2]]) {
    rid[i] <- dup[i,1]
  } else if(shp.hf$code6[dup[i,1]] != shp.hf$code6[dup[i,2]] & shp.hf$code6[dup[i,1]] < shp.hf$code6[dup[i,2]]) {
    rid[i] <- dup[i,2]
  } else {
    rid[i] <- dup[i,sample(1:2, 1)]
  }
}

cat('\nParsed', length(rid), 'duplicates')

shp.hf <- shp.hf[-rid,]

cat('\nFiltered heat flow:', nrow(shp.hf))

# Read predicted heat flow from Lucazeau 2019
read_delim('data/nghf/HFgrid14.csv',
	   delim = ';',
	   col_types = c('ddddd')) -> hf.pred 

# Make simple feature
hf.pred %>% 
  st_as_sf(
    coords = c(1,2),
    crs = proj4.wgs) %>% 
  st_transform(proj4.robin.pacific) -> shp.hf.pred

# Extract 0.5˚ x 0.5˚ grid from Lucazeau 2019
shp.grid <- st_geometry(shp.hf.pred)

# Read syracuse et al 2006 volcanoes
# Header
cat('\nReading Syracuse and Abers (2006) volcano data ...')
h.volc <- read_table('data/sa2006/volcanoes.txt', n_max = 1, col_type = cols())

# Read table without header
volc <- read_table('data/sa2006/volcanoes.txt', skip = 2, col_type = cols())

# Combine
names(volc) <- names(h.volc)

# Tidy
volcano <- volc %>%
  mutate(
    'ArcName' = as.factor(ArcName),
    'Age' = as.numeric(Age),
    'Phi/100' = as.numeric(`Phi/100`)) %>%
  filter(!is.na(H))

# Header
h.volc.no.spread <- read_table('data/sa2006/volcanoes_nospread.txt', n_max = 1, col_type = cols())

# Read table without header
volc.no.spread <- read_table('data/sa2006/volcanoes_nospread.txt', skip = 2, col_type = cols())

# Combine
names(volc.no.spread) <- names(h.volc.no.spread)

# Tidy
volcano.no.spread <- volc.no.spread %>%
  mutate(
    'ArcName' = as.factor(ArcName),
    'Age' = as.numeric(Age),
    'H' = as.numeric(H),
    'Dip' = as.numeric(Dip),
    'DesRat' = as.numeric(DesRat),
    'Phi/100' = as.numeric(`Phi/100`)
  ) %>%
  filter(!is.na(H))

# Change to simple features object and project
cat('\nProjecting volcano to:', proj4.robin.pacific, sep = '\n')
shp.volc <- volcano %>%
  st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj4.wgs) %>%
  st_transform(proj4.robin.pacific)
shp.volc.no.spread <- volcano.no.spread %>%
  mutate(
    'ArcName' = as.factor(ArcName),
    'H' = as.numeric(H),
    'Dip' = as.numeric(Dip),
    'DesRat' = as.numeric(DesRat),
    'Phi/100' = as.numeric(`Phi/100`),
  ) %>%
  rename(Vcorr = Vc) %>%
  filter(!is.na(H)) %>%
  st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj4.wgs) %>%
  st_transform(proj4.robin.pacific)

# Clean up environment
cat('\nCleaning up environment ...')
rm(list = lsf.str())
rm(files,
   i,
   h.volc,
   h.volc.no.spread,
   hf,
   shp.sliv,
   volc,
   volcano,
   dup,
   rid,
   volc.no.spread,
   volcano.no.spread,
   buf.dist,
   hf.pred)

# Heat flow table (without geometry column)
d.hf <- shp.hf %>% st_set_geometry(NULL)
d.volc <- shp.volc %>% st_set_geometry(NULL)
d.volc.no.spread <- shp.volc.no.spread %>% st_set_geometry(NULL)

# Save
fname.save <- 'data/hf.RData'
cat('\nSaving data to', fname.save, '\n', sep = ' ')
save.image(fname.save)
