# Load packages
cat('Loading libraries\n')

# Quiet loading
function(a.package){
	suppressWarnings(
	suppressPackageStartupMessages(
	library(a.package, character.only=TRUE)))
} -> sshhh

# Package list
c('magrittr', 'dplyr', 'ggplot2', 'tidyr', 'readr', 'purrr',
	'rgeos', 'sp', 'raster', 'stars', 'gstat', 'geosphere',
	'ggsflabel', 'metR', 'sf', "devtools", "magrittr", "dplyr",
	"tidyr", "readr", "ggplot2", "purrr", "knitr", "ggrepel", "patchwork",
	"cowplot") -> p.list

# auto-load quietly
sapply(p.list, sshhh)

# Rotate projection
rotate_proj <- function(obj, angle) {
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
      purrr::map2(files,
                  filenames,
                  ~ st_read(.x,
                            crs = '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs',
                            quiet = TRUE) %>%
                    st_transform(crs) %>%
                    tibble::as_tibble() %>%
                    tibble::add_column('segment' = .y) %>%
                    dplyr::group_by(segment) %>%
                    st_as_sf()) %>%
      bind_rows()
  } else {
      purrr::map2(files,
                  filenames,
                  ~ st_read(.x,
                            crs = '+proj=longlat +lon_wrap=180 +ellps=WGS84 +datum=WGS84 +no_defs',
                            quiet = TRUE) %>%
                    st_transform(crs) %>%
                    tibble::as_tibble() %>%
                    tibble::add_column('segment' = .y) %>%
                    dplyr::group_by(segment) %>%
                    st_as_sf()) %>%
      purrr::set_names(nm = filenames) %>%
      dplyr::bind_rows()
  }
}

# Krige optimization algorithmcross-validation used for optimizing variogram model
f.obj <- function(
  data,
  lags = 20,
  lag.cutoff = 3,
  param,
  v.mod = 0,
  v.sill = 1,
  v.range = 1,
	v.nub = 0,
  maxdist = Inf,
  fold = 10
) {
  # Projection
  st_crs(data) -> proj
  # Data
  data -> d
  # Drop simple features
  d %>% st_set_geometry(NULL) -> d.t
  # Variogram cutoff
  max(st_distance(data))/lag.cutoff -> cutoff
  # Experimental variogram
  variogram(get(param)~1,
            locations = d,
            cutoff = as.vector(cutoff),
            width = as.vector(cutoff/lags)) -> v
  # Variogram model discritization formula
  if(v.mod >= 0 && v.mod < 1) {
    v.mod <- 'Sph'
  } else if(v.mod >= 1 && v.mod < 2) {
    v.mod <- 'Exp'
  } else if(v.mod >= 2 && v.mod <= 3) {
    v.mod <- 'Gau'
  }
  # Model variogram
  fit.variogram(v, 
		vgm(psill = v.sill, 
		    model = v.mod, 
		    range = v.range,
				nugget = v.nug), 
		fit.method = 7) -> f
  # Kriging
  krige.cv(formula = get(param)~1,
           locations = d,
           model = f,
           maxdist = maxdist,
           verbose = T) %>%
    as_tibble() %>% 
    st_as_sf() %>% 
    st_set_crs(st_crs(d)) %>% 
    mutate('var1.pred' = round(var1.pred, 1)) %>% 
    rename(variance := var1.var,
           !!param := var1.pred) %>% 
    drop_na()-> k
  # Calculating cost function after Li et al 2018
  # Simultaneously minimizes misfit on variogram model and kriged interpolation errors
  # Weights
  wi <- 0.5 # interpolation error
  wf <- 0.5 # variogram fit error
  # Calculate variogram fit error
  # Root mean weighted (Nj/hj^2) sum of squared errors * (1-wi)/v.sigma
  sqrt((attr(f, "SSErr")^2)/nrow(v)) * ((1-wi)/sd(v$gamma)) -> v.s
  # Calculate interpolation error
  # RMSE * wi/i.sigma
  sqrt(sum((k$residual^2))/nrow(k)) * wi/sd(k$hf) -> i.s
  # Cost
  v.s + i.s
}
    
# Variogram fitting and kriging
krige_interp <- function(
  data,
  lags = 20,
  lag.cutoff = 3,
  param,
  grid = NULL,
  ngrid = 1e3,
  grid.method = 'regular',
  grid.rotate = 0,
  grid.shape = c('left' = 0.1, 'right' = 0.1, 'top' = 0.1, 'bottom' = 0.1),
  v.mod = c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
  krige = FALSE,
  plot = FALSE,
  interp = FALSE,
  interp.power = 2){
  # Rotated projection
  if(grid.rotate != 0) {
    rotate_proj(data, angle = grid.rotate) -> proj.rot
  } else {
    st_crs(data) -> proj.rot
  }
  # Data
  data %>% st_transform(proj.rot) -> d
  # Box
  bbox_widen(st_bbox(d),
             crs = st_crs(proj.rot),
             borders = grid.shape) %>% 
    st_as_sf() -> b
  # Drop simple features
  d %>% st_set_geometry(NULL) -> d.t
  # Variogram cutoff
  max(st_distance(data))/lag.cutoff -> cutoff
  # Sample within box for calculating kriging weights
  if(is.null(grid)) {
    b %>%
      st_sample(size = ngrid,
                type = grid.method) -> s
    b %>% 
      st_sample(size = ngrid,
                type = 'regular') -> s.reg
  } else {
    grid %>%
      st_transform(st_crs(d)) %>%
      st_crop(b) -> s
  }
  # Experimental variogram
  variogram(get(param)~1,
            locations = d,
            cutoff = as.vector(cutoff),
            width = as.vector(cutoff/lags)) -> v
  # Model variogram
  fit.variogram(v, vgm(model = v.mod)) -> f
  # Variogram plot
  if (plot == TRUE) {
    ggplot() +
      geom_path(data = variogramLine(f, maxdist = as.vector(cutoff)),
                aes(x = dist/1000, y = gamma)) +
      geom_point(data = v, aes(x = dist/1000, y = gamma), size = 0.6) +
      geom_text_repel(data = v, aes(x = dist/1000, y = gamma, label = np), size = 3) +
      labs(x = 'Lag (km)', y = 'Semivariance') +
      theme_classic() +
      theme(plot.background = element_rect(fill = 'transparent', color = NA),
            panel.background = element_rect(fill = 'transparent', color = NA),
            axis.text = element_text(color = 'black')) -> p
    p
  } else {
    NULL -> p
  }
  if(krige == TRUE){
    # Kriging
    krige(formula = get(param)~1,
          locations = d,
          newdata = s,
          model = f) %>%
      as_tibble() %>% 
      st_as_sf() %>% 
      mutate('var1.pred' = round(var1.pred, 1)) %>% 
      rename(variance := var1.var,
             !!param := var1.pred) -> k
    if(is.null(grid)) {
      idw(formula = get(param)~1,
          k,
          newdata = s.reg,
          idp = interp.power) %>% 
        as_tibble() %>% 
        st_as_sf() %>% 
        dplyr::select(-var1.var) %>% 
        mutate('var1.pred' = round(var1.pred, 1)) %>% 
        rename(!!param := var1.pred) -> i
    } else {
      NULL -> i
    }
    if(grid.rotate != 0) {
      cat('New projection:', proj.rot, sep = '\n')
    }
    return(list(variogram = v %>% 
                  as_tibble() %>% 
                  rename(semivar = gamma) %>% 
                  dplyr::select(np, dist, semivar),
                model.variogram = f,
                variogram.plot = p,
                krige.results = k,
                interp.results = i))
  } else {
    return(list(variogram = v %>% 
                  as_tibble() %>% 
                  rename(semivar = gamma) %>% 
                  dplyr::select(np, dist, semivar),
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
