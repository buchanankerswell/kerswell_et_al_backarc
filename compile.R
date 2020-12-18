rm(list = ls())
source('functions.R')

# Projections
# WGS84 Pseudo-Mercator [Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI]
proj <- 3857
proj4 <- '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'

# WGS84
proj.wgs <- 4326
proj4.wgs <- '+proj=longlat +datum=WGS84 +no_defs'

# Robinson
proj4.robin <- '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs'
proj4.robin.pacific <- '+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs'

# Graticules 30 degree
shp.grats30 <- st_read('data/grats30/', quiet = T) %>% st_transform(proj4.robin.pacific)

# Country Boundaries
# Create sliver around zero degree long
shp.zero.sliv <- st_polygon(x = list(rbind(c(-0.0001, 90), c(0, 90), c(0, -90), c(-0.0001, -90), c(-0.0001, 90)))) %>%
  st_sfc() %>%
  st_set_crs(proj.wgs)
# Trim the countries boundaries overlapping at zero degree long and reproject
shp.world.robin.pacific <- rnaturalearth::ne_countries(returnclass = 'sf') %>%
  st_difference(shp.zero.sliv) %>%
  tibble::as_tibble() %>%
  st_as_sf() %>% st_transform(proj4.robin.pacific)
shp.world.robin.pacific10 <- rnaturalearth::ne_countries(returnclass = 'sf', scale = 10) %>%
  st_difference(shp.zero.sliv) %>%
  tibble::as_tibble() %>%
  st_as_sf() %>% st_transform(proj4.robin.pacific)

# Read syracuse et al 2006 Segments
# seg.names <- list.files('./sa2006/gmts/') %>% purrr::map_chr(~ gsub('.contours.gmt', '', .x))
seg.names <- c('Alaska Aleutians', 'Andes', 'Central America', 'Kamchatka Marianas', 'Kyushu Ryukyu', 'Lesser Antilles', 'N. Philippines', 'New Britain Solomon', 'S. Philippines', 'Scotia', 'Sumatra Banda Sea', 'Tonga New Zealand', 'Vanuatu')
files <- list.files('data/sa2006/gmts', full.names = TRUE)

# Contours Robinson pacific centered
shp.sa.countours.robin.pacific <- read_wrap_latlong(files, seg.names, proj4.robin.pacific)
shp.sa.segs.robin.pacific <- shp.sa.countours.robin.pacific %>% group_by(segment) %>% filter(row_number() == 1)
shp.sa.segs.robin.pacific.buffer <- shp.sa.segs.robin.pacific %>% st_buffer(500000, endCapStyle = 'ROUND') %>% st_wrap_dateline(options = c("WRAPDATELINE=YES"))

# Make boxes around segments
box.segs.all <- st_bbox(shp.sa.segs.robin.pacific.buffer) %>% bbox_widen(proj4.robin.pacific, c('left' = 0, 'right' = 0, 'top' = 0, 'bottom' = 0))
box.segs.all.wide <- st_bbox(shp.sa.segs.robin.pacific.buffer) %>% bbox_widen(proj4.robin.pacific, c('left' = 0.02, 'right' = 0.01, 'top' = 0.008, 'bottom' = 0.008))
box.segs.bbox <- shp.sa.segs.robin.pacific.buffer %>% tidyr::nest() %>% purrr::pmap(~st_bbox(..2)) %>% purrr::set_names(unique(shp.sa.segs.robin.pacific.buffer$segment))
box.segs <- box.segs.bbox  %>% purrr::map(~bbox_widen(.x, crs = proj4.robin.pacific, borders = c('left' = 0, 'right' = 0, 'top' = 0, 'bottom' = 0)))
box.segs.wide <- box.segs.bbox %>% purrr::map(~bbox_widen(.x, crs = proj4.robin.pacific, borders = c('left' = 0.05, 'right' = 0.05, 'top' = 0.05, 'bottom' = 0.05)))

# Draw box to stuff pacific labels into
box.lab <- st_bbox(shp.sa.segs.robin.pacific) %>% bbox_widen(proj4.robin.pacific, borders = c('left' = -0.3, 'right' = -0.25, 'top' = -0.12, 'bottom' = -0.40))

# Subsegments
shp.sa.segs.robin.pacific.subseg <- purrr::map_df(seg.names, ~{
  shp.sa.segs.robin.pacific %>% filter(segment == .x) %>%
    splt(cut.length = 1000000, buffer = F, buffer.dist = 500000) %>% mutate(segment = .x, .before = geometry)
})

# Subsegment buffers
shp.sa.segs.robin.pacific.subseg.buffer <- purrr::map_df(seg.names, ~{
  shp.sa.segs.robin.pacific %>% filter(segment == .x) %>%
    splt(cut.length = 1000000, buffer = T, buffer.dist = 500000) %>% mutate(segment = .x, .before = geometry)
})

# Read global heat flow database (IHFC 2010), turn into tibble,
# remove useless columns & filter heat flow between 0 and 200, then crop to buffers
hf <- st_read('data/hf/', coords = c(3,2), crs = proj.wgs, quiet = T) %>%
  st_transform(proj4.robin.pacific) %>%
  tibble::as_tibble() %>% st_as_sf() %>%
  select(-Data_Numbe, -Codes, -Year_of_Pu, -Comments, 
         -F21, -F22, -F23, -F24, -F25, -F26, -ben, -europe) %>%
  filter(Heat_Flow <= 200 & Heat_Flow >= 0) %>%
  mutate('No__Temps' = as.numeric(No__Temps), 'No_Heat_Pr' = as.numeric(No_Heat_Pr), 'Heat_Prod_' = as.numeric(Heat_Prod_))
shp.hf <- hf %>% st_join(shp.sa.segs.robin.pacific.buffer, left = F)

# Read global heat flow database (IHFC 2010), turn into tibble,
# remove useless columns & filter heat flow between 0 and 200, then crop to subsegment buffers
shp.hf.subseg <- hf %>% st_join(shp.sa.segs.robin.pacific.subseg.buffer, left = F)

# Read syracuse et al 2006 volcanoes
# Header
h.volc <- readr::read_table('data/sa2006/volcanoes.txt', n_max = 1)
# Read table without header
volc <- readr::read_table('data/sa2006/volcanoes.txt', skip = 2)
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
h.volc.no.spread <- readr::read_table('data/sa2006/volcanoes_nospread.txt', n_max = 1)
# Read table without header
volc.no.spread <- readr::read_table('data/sa2006/volcanoes_nospread.txt', skip = 2)
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
shp.volc <- volcano %>%
  st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj.wgs) %>%
  st_transform(proj4.robin.pacific)
shp.volc.no.spread <- volcano.no.spread %>%
  mutate(
    'ArcName' = as.factor(ArcName),
    'H' = as.numeric(H),
    'Dip' = as.numeric(Dip),
    'DesRat' = as.numeric(DesRat),
    '`Phi/100`' = as.numeric(`Phi/100`),
  ) %>% 
  filter(!is.na(H)) %>%
  st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj.wgs) %>%
  st_transform(proj4.robin.pacific)

# Clean up environment
rm(list = lsf.str())
rm(files, h.volc, h.volc.no.spread)

# Save
save.image('data/hf.RData')
