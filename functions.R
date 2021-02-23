# Use cairo for png device
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}
cat('Loading libraries\n')
suppressMessages({
  # Load packages
  library(rgeos)
  library(sp)
  library(raster)
  library(sf)
  library(stars)
  library(gstat)
  library(ggrepel)
  library(ggsflabel)
  library(cowplot)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(geosphere)
  library(metR)
})

# Rotate projection
rotate_proj = function(obj, angle) {
  # get bounding box as spatial points object
  box.center <- bbox_widen(st_bbox(obj),
                      crs = st_crs(obj),
                      borders = c('left' = 0, 'right' = 0, 'top' = 0, 'bottom' = 0)) %>% 
    st_centroid() %>% 
    st_transform(crs = 4326) %>% 
    st_coordinates()
  # construct the proj4 string
  paste0("+proj=omerc +lat_0=", box.center[2], " +lonc=", box.center[1], " +x_0=0 +y_0=0 +alpha=", 
               angle, " +gamma=0 +k_0=1 +ellps=WGS84 +datum=WGS84 +units=m +no_def")
  # return as a CRS:
  # st_crs(prj)
}

# Draw a widened box from a st_bbox object
bbox_widen <- function(bbox, crs, borders = c('left' = 0.5, 'right' = 0.5, 'top' = 0, 'bottom' = 0)) {
  b <- bbox # current bounding box
  xrange <- b$xmax - b$xmin # range of x values
  yrange <- b$ymax - b$ymin # range of y values
  b[1] <- b[1] - (borders['left'] * xrange) # xmin - left
  b[3] <- b[3] + (borders['right'] * xrange) # xmax - right
  b[2] <- b[2] - (borders['bottom'] * yrange) # ymin - bottom
  b[4] <- b[4] + (borders['top'] * yrange) # ymax - top
  box <- st_polygon(list(matrix(c(b$xmin, b$ymax, b$xmin, b$ymin, b$xmax, b$ymin, b$xmax, b$ymax, b$xmin, b$ymax), ncol = 2, byrow = TRUE))) %>%
    st_sfc(crs = crs)
  return(box)
}

# Read gmt files, wrap the dateline to avoid plotting horizontal lines on map,
# make into tibble, add segment names, transform projection, and bind into one sf object
read_latlong <- function(files, filenames = NULL, crs) {
  if(is.null(filenames)) {
    shp <- purrr::map2(files, filenames,
                       ~ st_read(.x, crs = proj.wgs, quiet = TRUE) %>%
                         st_transform(crs) %>%
                         tibble::as_tibble() %>%
                         tibble::add_column('segment' = .y) %>%
                         group_by(segment) %>%
                         st_as_sf()) %>%
      bind_rows()
  } else {
    shp <- purrr::map2(files, filenames,
                       ~ st_read(.x, crs = proj.wgs, quiet = TRUE) %>%
                         st_transform(crs) %>%
                         tibble::as_tibble() %>%
                         tibble::add_column('segment' = .y) %>%
                         group_by(segment) %>%
                         st_as_sf()) %>%
      purrr::set_names(nm = filenames) %>%
      bind_rows()
  }
  return(shp)
}

# Variogram fitting and kriging
model_variogram <- function(
  data,
  lags = 30,
  lag.cutoff = 5,
  param,
  ngrid = 1000,
  grid.method = 'regular',
  grid.rotate = 0,
  grid.shape = c('left' = 0.1, 'right' = 0.1, 'top' = 0.1, 'bottom' = 0.1),
  v.mod = c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
  krige = FALSE,
  interp.power = 2,
  proj = st_crs(data)){
  # Rotated projection
  proj.rot <- rotate_proj(data, angle = grid.rotate)
  # Data
  d <- data %>% st_transform(proj.rot)
  # Box
  b <- bbox_widen(st_bbox(d),
                  crs = st_crs(proj.rot),
                  borders = grid.shape) %>% 
    st_as_sf()
  # Drop simple features
  d.t <- d %>% st_set_geometry(NULL)
  # Variogram cutoff
  cutoff <- max(st_distance(data))/lag.cutoff
  # Sample within box for calculating kriging weights
  s <- b %>%
    st_sample(size = ngrid, type = grid.method)
  s.reg <- b %>% 
    st_sample(size = ngrid, type = 'regular')
  # Experimental variogram
  v <- variogram(get(param)~1,
                 locations = d,
                 cutoff = as.vector(cutoff),
                 width = as.vector(cutoff/lags))
  # Model variogram
  f <- fit.variogram(v, v.mod)
  # Variogram plot
  p <- ggplot() +
    geom_path(data = variogramLine(f, maxdist = as.vector(cutoff)),
              aes(x = dist/1000, y = gamma)) +
    geom_point(data = v, aes(x = dist/1000, y = gamma), size = 0.6) +
    geom_text_repel(data = v, aes(x = dist/1000, y = gamma, label = np), size = 3) +
    labs(x = 'Lag (km)', y = 'Semivariance') +
    theme_classic() +
    theme(plot.background = element_rect(fill = 'transparent', color = NA),
          panel.background = element_rect(fill = 'transparent', color = NA),
          axis.text = element_text(color = 'black'))
  if(krige == TRUE){
    p
    # Kriging
    k <- krige(formula = get(param)~1,
               locations = d,
               newdata = s,
               model = f) %>%
      as_tibble() %>% 
      st_as_sf() %>% 
      mutate('var1.pred' = round(var1.pred)) %>% 
      rename(variance := var1.var,
             !!param := var1.pred)
    i <- idw(formula = get(param)~1,
             k,
             newdata = s.reg,
             idp = interp.power) %>% 
      as_tibble() %>% 
      st_as_sf() %>% 
      select(-var1.var) %>% 
      mutate('var1.pred' = round(var1.pred)) %>% 
      rename(!!param := var1.pred)
    cat('New projection:', proj.rot, sep = '\n')
    return(list(variogram = v %>% 
                  as_tibble() %>% 
                  rename(semivar = gamma) %>% 
                  select(np, dist, semivar),
                model.variogram = f,
                variogram.plot = p,
                krige.results = k,
                interp.results = i))
  } else {
    p
    return(list(variogram = v %>% 
                  as_tibble() %>% 
                  rename(semivar = gamma) %>% 
                  select(np, dist, semivar),
                model.variogram = f,
                variogram.plot = p))
  }
}

# Split segments and buffers into equidistant subsegments
splt <- function(object,
                 cut.prop = 4,
                 buffer = TRUE,
                 buffer.dist = 500000,
                 cap.style = 'ROUND', ...) {
  # Cast to points
  obj.pnts <- suppressWarnings(object %>% st_cast("POINT"))
  # Calculate distances between consecutive points
  dst <- as.vector(c(units::as_units(0, 'm'),
                     st_distance(obj.pnts[-1,],
                                 obj.pnts[-nrow(obj.pnts),],
                                 by_element = TRUE)))
  # Calculate cut length
  cut.length <- cumsum(dst)[length(cumsum(dst))] / cut.prop
  # If segment length is shorter than cut length
  if(cumsum(dst)[length(cumsum(dst))] < cut.length) {
    if (buffer != TRUE) {
      warning('Segment length is shorter than cut length. Returning original segment')
      return(object)
    } else if (buffer == TRUE) {
      warning('Segment length is shorter than cut length. Returning original segment with buffer')
      return(st_buffer(object, buffer.dist, endCapStyle = cap.style))
    }
  } else {
    # Find indices to split points into groups with equal distances
    cut.ind <- c(1, rep(NA, ceiling(cut.prop)))
    save.ind <- 2
    start.ind <- 1
    for(i in 1:length(dst)) {
      cmsm <- cumsum(dst[start.ind:i])
      tot <- cmsm[length(cmsm)]
      if(tot < cut.length) {
        i <- i + 1
      } else {
        cut.ind[save.ind] <- i
        start.ind <- i
        save.ind <- save.ind + 1
      }
    }
    # Last cut should be end of segment line
    cut.ind[length(cut.ind)] <- length(dst)
    # Split the object using the saved cut indices, group by subsegment, and cast back into linestrings
    obj.splt <- suppressWarnings(
      purrr::pmap(list(cut.ind[-length(cut.ind)],
                          c(cut.ind[-1]),
                          1:(length(cut.ind)-1)),
                     ~{obj.pnts[..1:..2,] %>%
                         mutate(subseg = ..3, .before = geometry) %>% 
                         summarise(do_union = F, .groups = 'keep') %>%
                         st_cast("LINESTRING")})
    )
    # Return result
    if(buffer != TRUE) {
      return(obj.splt)
    } else if (buffer == TRUE) {
      return(obj.splt %>% purrr::map(~{
        st_buffer(.x, buffer.dist, endCapStyle = cap.style)
        })
        )
    }
  }
}

# Find point positions relative to the trench
pts_position <- function(sf.pts, sf.line, offset.x, offset.y, direction = c('updown', 'leftright'), arc.direction = c('up', 'down', 'left', 'right')){
  if(direction == 'leftright'){
    # Make polygon to left of segment
    poly.left <- rbind(c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1], st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2] + offset.y),
                    c(st_bbox(sf.line)['xmin'] - offset.x, st_bbox(sf.line)['ymax'] + offset.y),
                    c(st_bbox(sf.line)['xmin'] - offset.x, st_coordinates(sf.line)[1, 2] - offset.y),
                    c(st_coordinates(sf.line)[1, 1], st_coordinates(sf.line)[1, 2] - offset.y),
                    as.data.frame(st_coordinates(sf.line))[,c(1,2)],
                    c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1], st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2] + offset.y)) %>%
      as.matrix() %>% list() %>% st_polygon() %>% st_sfc(crs = st_crs(sf.line)) %>% st_sf()
    # Make polygon to right of segment
    poly.right <- rbind(c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1], st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2] + offset.y),
                    c(st_bbox(sf.line)['xmin'] + offset.x, st_bbox(sf.line)['ymax'] + offset.y),
                    c(st_bbox(sf.line)['xmin'] + offset.x, st_coordinates(sf.line)[1, 2] - offset.y),
                    c(st_coordinates(sf.line)[1, 1], st_coordinates(sf.line)[1, 2] - offset.y),
                    as.data.frame(st_coordinates(sf.line))[,c(1,2)],
                    c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1], st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2] + offset.y)) %>%
      as.matrix() %>% list() %>% st_polygon() %>% st_sfc(crs = st_crs(sf.line)) %>% st_sf()
    # Find points to the left of segment
    pts.left <- st_intersects(sf.pts, poly.left, sparse = F)
    # Determine position of points relative to the segment on the arc-side or outbound of the trench
    if(arc.direction == 'left'){
      return(ifelse(pts.left, 'arc-side', 'outboard') %>% as.vector())
    } else if(arc.direction == 'right'){
      return(ifelse(pts.left, 'outboard', 'arc-side') %>% as.vector())
    }
  } else if(direction == 'updown'){
    # Make polygon above segment
    poly.up <- rbind(c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1] + offset.x, st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2]),
                    c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1] + offset.x, st_bbox(sf.line)['ymax'] + offset.y),
                    c(st_bbox(sf.line)['xmin'] - offset.x, st_coordinates(sf.line)[1, 2] + offset.y),
                    c(st_coordinates(sf.line)[1, 1] - offset.x, st_coordinates(sf.line)[1, 2]),
                    as.data.frame(st_coordinates(sf.line))[,c(1,2)],
                    c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1] + offset.x, st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2])) %>%
      as.matrix() %>% list() %>% st_polygon() %>% st_sfc(crs = st_crs(sf.line)) %>% st_sf()
    # Make polygon below segment
    poly.down <- rbind(c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1] + offset.x, st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2]),
                    c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1] + offset.x, st_bbox(sf.line)['ymax'] - offset.y),
                    c(st_bbox(sf.line)['xmin'] - offset.x, st_coordinates(sf.line)[1, 2] - offset.y),
                    c(st_coordinates(sf.line)[1, 1] - offset.x, st_coordinates(sf.line)[1, 2]),
                    as.data.frame(st_coordinates(sf.line))[,c(1,2)],
                    c(st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 1] + offset.x, st_coordinates(sf.line)[nrow(st_coordinates(sf.line)), 2])) %>%
      as.matrix() %>% list() %>% st_polygon() %>% st_sfc(crs = st_crs(sf.line)) %>% st_sf()
    # Find points above segment
    pts.up <- st_intersects(sf.pts, poly.up, sparse = F)
    if(arc.direction == 'up'){
      return(ifelse(pts.up, 'arc-side', 'outboard') %>% as.vector())
    } else if(arc.direction == 'down'){
      return(ifelse(pts.up, 'outboard', 'arc-side') %>% as.vector())
    }
  }
}

# Calculate point distances relative to trench (positive and negative)
trench_distance <- function(sf.pts, pos, sf.trench){
  dist <- st_distance(sf.pts, sf.trench) %>% as.vector()
  dist[pos == 'outboard'] <- -dist[pos == 'outboard']
  return(dist)
}
