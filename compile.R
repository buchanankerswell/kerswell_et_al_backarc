suppressMessages(source('functions.R'))
suppressWarnings({
# Projections
# WGS84 Pseudo-Mercator [Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI]
proj <- 3857
proj4 <- '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'

# WGS84
proj.wgs <- 4326
proj4.wgs <- '+proj=longlat +datum=WGS84 +no_defs'

# Robinson
proj4.robin <- '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
proj4.robin.pacific <- "+proj=robin +lon_0=-155 +lon_wrap=-155 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Print Projections
cat('Defining CRS (projections): ', proj4, proj4.wgs, proj4.robin.pacific, sep = '\n')

# Graticules 30 degree
cat('Loading country boundaries and graticules ...', sep = '\n')
shp.grats30 <- st_read('data/grats30/', quiet = T) %>% st_transform(proj4.robin.pacific)

# Country Boundaries
# Create sliver around zero degree long
shp.sliv <- st_polygon(x = list(rbind(c(25.000001, 90),
                                      c(25, 90),
                                      c(25, -90),
                                      c(25.000001, -90),
                                      c(25.000001, 90)))) %>%
  st_sfc() %>%
  st_set_crs(proj.wgs)

# Trim the countries boundaries overlapping at zero degree long and reproject
suppressMessages(
  shp.world.robin.pacific <- rnaturalearth::ne_countries(returnclass = 'sf') %>%
    st_difference(shp.sliv) %>%
    tibble::as_tibble() %>%
    st_as_sf() %>% st_transform(proj4.robin.pacific)
)

# Read syracuse et al 2006 Segments
cat('Loading Syracuse and Abers (2006) segment boundaries ...', sep = '\n')
seg.names <- c('Alaska Aleutians', 'Andes', 'Central America', 'Kamchatka Marianas', 'Kyushu Ryukyu', 'Lesser Antilles', 'N. Philippines', 'New Britain Solomon', 'S. Philippines', 'Scotia', 'Sumatra Banda Sea', 'Tonga New Zealand', 'Vanuatu')
files <- list.files('data/sa2006/gmts', full.names = TRUE)
seg.names

# Contours Robinson pacific centered
shp.sa.countours.robin.pacific <- read_latlong(files, seg.names, proj4.robin.pacific)

# Filter for segment boundary
shp.sa.segs.robin.pacific <- shp.sa.countours.robin.pacific %>% group_by(segment) %>% filter(row_number() == 1)

# Buffer
buf.dist <- c(500000, 750000, 1000000, 1500000)
cat('Drawing', buf.dist/1000, 'km buffers around segments ...\n', sep = ' ')
shp.sa.segs.robin.pacific.buffer <- purrr::map(buf.dist,
                                               ~shp.sa.segs.robin.pacific %>%
                                                 st_buffer(.x, endCapStyle = 'SQUARE') %>%
                                                 st_wrap_dateline(options = c("WRAPDATELINE=YES"))) %>%
  purrr::set_names(nm = paste0('km.', buf.dist/1000))

# Make boxes around segments
box.segs.all <- st_bbox(shp.sa.segs.robin.pacific.buffer$km.1500) %>% bbox_widen(proj4.robin.pacific, c('left' = 0, 'right' = 0, 'top' = 0, 'bottom' = 0))
box.segs.all.wide <- st_bbox(shp.sa.segs.robin.pacific.buffer$km.1500) %>% bbox_widen(proj4.robin.pacific, c('left' = 0.02, 'right' = 0.01, 'top' = 0.008, 'bottom' = 0.008))
box.segs.bbox <- shp.sa.segs.robin.pacific.buffer$km.1500 %>% tidyr::nest() %>% purrr::pmap(~st_bbox(..2)) %>% purrr::set_names(unique(shp.sa.segs.robin.pacific.buffer$km.1500$segment))
box.segs <- box.segs.bbox  %>% purrr::map(~bbox_widen(.x, crs = proj4.robin.pacific, borders = c('left' = 0, 'right' = 0, 'top' = 0, 'bottom' = 0)))
box.segs.wide <- box.segs.bbox %>% purrr::map(~bbox_widen(.x, crs = proj4.robin.pacific, borders = c('left' = 0.05, 'right' = 0.05, 'top' = 0.05, 'bottom' = 0.05)))

# Draw box to stuff pacific labels into
box.lab <- bbox_widen(st_bbox(box.segs.all), crs = proj4.robin.pacific, borders = c('left' = -0.5, 'right' = -0.1, 'top' = -0.1, 'bottom' = -0.1))

# Read global heat flow database (IHFC 2010), turn into tibble,
# remove useless columns, then crop to buffers
cat('Reading global heat flow data and projecting to:', proj4.robin.pacific, sep = '\n')
hf <- st_read('data/hf/', coords = c(3,2), crs = proj.wgs, quiet = T) %>%
  st_transform(proj4.robin.pacific) %>%
  tibble::as_tibble() %>% st_as_sf() %>%
  select(-Data_Numbe, -Codes, -Year_of_Pu, -Comments,
         -F21, -F22, -F23, -F24, -F25, -F26, -ben, -europe) %>%
  mutate('No__Temps' = as.numeric(No__Temps),
         'No_Heat_Pr' = as.numeric(No_Heat_Pr),
         'Heat_Prod_' = as.numeric(Heat_Prod_)) %>%
  rename(id = Object_ID,
         site = Site_Name,
         ele = Elevation,
         tempNo = No__Temps,
         grad = Gradient,
         condNo = No__Cond_,
         con = Conductivi,
         heatPrNo = No_Heat_Pr,
         heatPr = Heat_Prod_,
         hf = Heat_Flow,
         siteNo = No__sites,
         ref = Reference)

# Remove no heat flow values
cat('Removing zero heat flow and duplicate values ...', sep = '\n')
cat('NAs: ', length(hf[!is.na(hf$hf), ]), '\n', sep = ' ')
cat('Zero heat flow: ', length(hf[hf$hf != 0, ]), '\n', sep = ' ')
hf <- hf[!is.na(hf$hf), ]
hf <- hf[hf$hf != 0, ]

# Check for duplicate measurements and remove
cat('Duplicates: ', length(sp::zerodist(sf::as_Spatial(hf))), '\n', sep = ' ')
if(nrow(sp::zerodist(sf::as_Spatial(hf))) != 0) {
  hf <- hf[-sp::zerodist(sf::as_Spatial(hf))[,1],]
}

# Filter hf
cat('Filtering heat flow:\n grad < 500 & grad > 0 & con < 8 & tempNo < 250', sep = '\n')
hf <- hf %>% filter(grad < 500 & grad > 0 & con < 8 & tempNo < 250 & hf < 500)
shp.hf <- hf

# Read syracuse et al 2006 volcanoes
# Header
cat('Reading Syracuse and Abers (2006) volcano data ...', sep = '\n')
h.volc <- suppressMessages(
  readr::read_table('data/sa2006/volcanoes.txt', n_max = 1)
)

# Read table without header
volc <- suppressMessages(
  readr::read_table('data/sa2006/volcanoes.txt', skip = 2)
)

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
h.volc.no.spread <- suppressMessages(
  readr::read_table('data/sa2006/volcanoes_nospread.txt', n_max = 1)
)

# Read table without header
volc.no.spread <- suppressMessages(
  readr::read_table('data/sa2006/volcanoes_nospread.txt', skip = 2)
)

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
cat('Projecting volcano to:', proj4.robin.pacific, sep = '\n')
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
    'Phi/100' = as.numeric(`Phi/100`),
  ) %>%
  rename(Vcorr = Vc) %>%
  filter(!is.na(H)) %>%
  st_as_sf(coords = c('Lon', 'Lat')) %>%
  st_set_crs(proj.wgs) %>%
  st_transform(proj4.robin.pacific)

# Clean up environment
cat('Cleaning up environment ...', sep = '\n')
rm(list = lsf.str())
rm(files, h.volc, h.volc.no.spread, hf)

# Heat flow table (without geometry column)
d.hf <- shp.hf %>% st_set_geometry(NULL)
d.volc <- shp.volc %>% st_set_geometry(NULL)
d.volc.no.spread <- shp.volc.no.spread %>% st_set_geometry(NULL)

# Save
fname.save <- 'data/hf.RData'
cat('Saving data to', fname.save, '\n', sep = ' ')
save.image(fname.save)
})
