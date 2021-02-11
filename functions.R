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
})

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

# Kriging
krg <- function(data, lags = 100, lag.cutoff = 3, param, buf, seg, contours, volc, crs, ngrid = 3e5, grid.method = 'hexagonal', v.mod = 'Sph', plot = TRUE){
  # Buffer
  b <- buf
  # Filter data inside buffer
  d <- data %>% st_join(b, left = F)
  # Check for duplicate measurements and remove
  if(nrow(sp::zerodist(sf::as_Spatial(d))) != 0) {
    d <- d[-sp::zerodist(sf::as_Spatial(d))[,1],]
  }
  # Drop simple features
  d.t <- d %>% st_set_geometry(NULL)
  # Variogram cutoff
  cutoff <- max(st_distance(data))/lag.cutoff
  if(cutoff <= units::as_units(0, 'm')) {
    cutoff <- units::as_units(1, 'm')
  }
  # Volcanoes
  volc <- volc
  # Sample within buffer for calculating kriging weights
  s <- b %>% st_sample(size = ngrid, grid.method)
  # Experimental variogram
  v <- variogram(get(param)~1, locations = d, cutoff = as.vector(cutoff), width = as.vector(cutoff/lags))
  if(is.null(v)) {
    warning('Only one data point. Cannot calculate variogram: returning nothing')
    return(NULL)
  }
  # Model variogram
  f <- fit.variogram(v, vgm(v.mod))
  # Model variogram line for plotting
  v.m <- variogramLine(f, maxdist = max(v$dist))
  # Box to crop simple features within for plotting
  crp <- bbox_widen(st_bbox(b), crs = crs, c('left' = 0.1, 'right' = 0.1, 'top' = 0.1, 'bottom' = 0.1))
  # Kriging
  k <- krige(formula = get(param)~1, locations = d, newdata = s, model = f) %>%
    select(-var1.var) %>%
    mutate('var1.pred' = round(var1.pred))
  # Plotting
  if(plot == TRUE){
    p.h <- ggplot() +
      geom_histogram(data = d, aes_string(x = param), alpha = 0.8, bins = 30) +
      labs(x = bquote(Heat~Flow~~mWm^-2), y = 'Frequency') +
      theme_classic(base_size = 14) +
      theme(axis.text = element_text(color = 'black'), axis.ticks = element_line(color = 'black'))
    p.v <- ggplot() +
      geom_point(data = v, aes(x = dist/1000, y = gamma), size = 1) +
      geom_path(data = v.m, aes(x = dist/1000, y = gamma)) +
      labs(x = bquote(Lag~~km), y = 'Semivariance') +
      coord_cartesian(ylim = c(0, NA)) +
      theme_classic(base_size = 14) +
      theme(axis.text = element_text(color = 'black'), axis.ticks = element_line(color = 'black'))
    p <- ggplot() +
      geom_sf(data = bbox_widen(st_bbox(b), crs = crs, c('left' = 0.1, 'right' = 0.1, 'top' = 0.1, 'bottom' = 0.1)), fill = 'cornflowerblue', alpha = 0.2, color = NA) +
      geom_stars(data = st_rasterize(k) %>% st_crop(b)) +
      geom_sf(data = b %>% st_crop(crp), fill = NA, color = 'black', size = 0.25, alpha = 0.8) +
      geom_sf(data = contours %>% st_crop(crp), color = 'white', alpha = 0.25, size = 0.25) +
      geom_sf(data = seg, fill = NA, color = 'black', size = 1.25) +
      geom_sf(data = volc %>% st_crop(crp), aes(shape = 'volcano'), color = 'deeppink4', alpha = 0.5) +
      geom_sf(data = d, aes_string(color = param), size = 0.5, show.legend = F) +
      labs(x = NULL, y = NULL) +
      scale_shape_manual(name = NULL, values = c('volcano' = 2)) +
      scale_color_viridis_c(option = 1) +
      scale_fill_viridis_c(name = bquote(mWm^-2), option = 1, na.value = 'transparent') +
      theme_map(font_size = 11) +
      theme(
        axis.text = element_text(),
        panel.border = element_blank(),
        panel.grid = element_line(size = 0.25, color = rgb(0.1, 0.1, 0.1, 0.5)),
        panel.background = element_blank(),
        panel.ontop = TRUE,
        plot.background = element_rect(fill = "transparent", color = NA)
      )
    return(list(shp.d = d, data = d.t, variogram = v, v.model = f, krg = k, hist = p.h, k.plot = p, v.plot = p.v))
  } else {
    return(list(shp.d = d, data = d.t, variogram = v, v.model = f, krg = k))
  }
}

# Split segments and buffers into equidistant subsegments
splt <- function(object, cut.prop = 4, buffer = TRUE, buffer.dist = 500000, cap.style = 'ROUND', ...) {
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
      purrr::pmap_df(list(cut.ind[-length(cut.ind)],
                          c(cut.ind[-1]),
                          1:(length(cut.ind)-1)),
                     ~{obj.pnts[..1:..2,] %>% mutate(subseg = ..3, .before = geometry)}) %>%
        group_by(subseg) %>% summarise(do_union = F, .groups = 'keep') %>% st_cast("LINESTRING")
    )
    # Return result
    if(buffer != TRUE) {
      return(obj.splt)
    } else if (buffer == TRUE) {
      return(st_buffer(obj.splt, buffer.dist, endCapStyle = cap.style, ...))
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
