if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}
library(tidyr)
library(leaflet)
library(ggsflabel)
library(ggplot2)
library(dplyr)
library(sf)
library(ggmap)
library(cowplot)
library(gstat)
library(sp)
library(stars)
library(patchwork)
library(ggrepel)

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
read_wrap_latlong <- function(files, filenames = NULL, crs) {
  if(is.null(filenames)) {
    shp <- purrr::map2(files, filenames, 
                       ~ st_read(.x, crs = proj.wgs, quiet = TRUE) %>%
                         st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")) %>%
                         st_transform(crs) %>%
                         tibble::as_tibble() %>%
                         tibble::add_column('segment' = .y) %>%
                         group_by(segment) %>%
                         st_as_sf()) %>%
      bind_rows()
  } else {
    shp <- purrr::map2(files, filenames, 
                       ~ st_read(.x, crs = proj.wgs, quiet = TRUE) %>%
                         st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")) %>%
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

# Heat flow map plot
plot_hf <- function(
  grats = NULL,
  world = NULL,
  box = NULL,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = NULL,
  seg.contours = NULL,
  seg.buffer = NULL,
  seg.names = TRUE,
  seg.label.size = 3,
  seg.label.alpha = 0.7,
  heatflow = NULL,
  hf.size = 0.1,
  viridis.pal = 'viridis',
  volcano = NULL,
  volc.size.range = c(1, 6),
  volc.size = 0.5,
  volc.alpha = 1,
  volc.color = 'black',
  box.labels = NULL,
  box.limits = NULL,
  legend.pos = 'right',
  legend.just = 'left',
  legend.size = unit(8, 'lines'),
  plot.labs = NULL) {
  if("patchwork" %in% (.packages())){
    detach("package:patchwork", unload=TRUE) 
  }
  volcs <- volcano %>% bind_rows()
  p <- ggplot()
  if(!is.null(grats)) {p <- p + geom_sf(data = grats, size = 0.15)}
  if(!is.null(world)) {p <- p + geom_sf(data = world, alpha = 0.8, color = NA)}
  if(!is.null(box)) {p <- p + geom_sf(data = box, fill = 'cornflowerblue', alpha = 0.15, color = NA)}
  if(!is.null(seg.contours)) {p <- p + geom_sf(data = seg.contours, size = 0.25, color = 'black')}
  if(!is.null(seg)) {p <- p + geom_sf(data = seg, size = 1, color = 'black')}
  if(!is.null(seg.buffer)) {p <- p + geom_sf(data = seg.buffer, fill = 'black', color = NA, alpha = 0.1)}
  if(!is.null(heatflow)) {p <- p + geom_sf(data = heatflow, aes(color = Heat_Flow), size = hf.size)}
  if(!is.null(volcano)) {p <- p + geom_sf(data = volcs, aes(size = H), shape = 2, color = volc.color, alpha = volc.alpha, stroke = volc.size)}
  if(!is.null(seg) && seg.names) {p <- p + geom_sf_label_repel(data = seg, aes(label = segment), size = seg.label.size, alpha = seg.label.alpha, xlim = c(st_bbox(box.labels)$xmin, st_bbox(box.labels)$xmax), ylim = c(st_bbox(box.labels)$ymin, st_bbox(box.labels)$ymax))}
  if(!is.null(box.limits)) {p <- p + coord_sf(crs = proj4.robin.pacific, xlim = c(st_bbox(box.limits)$xmin, st_bbox(box.limits)$xmax), ylim = c(st_bbox(box.limits)$ymin, st_bbox(box.limits)$ymax))}
  p <- p + scale_color_viridis_c(option = viridis.pal) + scale_size_continuous(range = volc.size.range)
  ifelse(!is.null(plot.labs), {p <- p + labs(title = plot.labs, caption = 'Data: Syracuse & Abers (2006); IHFC (2010)', y = '', x = '')}, {p <- p + labs(caption = 'Data: Syracuse & Abers (2006); IHFC (2010)', y = '', x = '')})
  ifelse(legend.pos %in% c('bottom', 'top'), {p <- p + guides(color = guide_colorbar(title = bquote(mW/m^2), barwidth = legend.size, ticks = F))}, {p <- p + guides(color = guide_colorbar(title = bquote(mW/m^2), barlength = legend.size, ticks = F))})
  p <- p + guides(size = guide_legend(override.aes = list(alpha = 1)))
  p <- p + theme_map(font_size = 11) + 
    theme(
      axis.text = element_text(),
      legend.justification = legend.just,
      legend.position = legend.pos,
      panel.border = element_blank(),
      panel.grid = element_line(size = 0.25, color = rgb(0.1, 0.1, 0.1, 0.5)),
      panel.background = element_blank(),
      panel.ontop = TRUE,
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  if (legend.pos %in% c('bottom', 'top')) {p <- p + theme(legend.box.margin = margin(0, 8, 0, 0), legend.margin = margin(0, 8, 0, 8))}
  if (any(grepl('*Andes.*', plot.labs) == TRUE)) {p <- p + theme(plot.title = element_text(hjust = 1))}
  return(p)
}

# Kriging
krg <- function(data, lags = 50, lag.cutoff = 5, param, krg.shp, seg, contours, crs, ngrid = 3e5, grid.method = 'hexagonal', v.mod = 'Sph', plot = TRUE){
  d <- data
  cutoff <- max(st_distance(data))/lag.cutoff
  b <- krg.shp
  s <- b %>% st_sample(size = ngrid, grid.method)
  v <- variogram(get(param)~1, locations = d, cutoff = as.vector(cutoff), width = as.vector(cutoff/lags))
  f <- fit.variogram(v, vgm(NA, v.mod, NA, NA))
  v.m <- variogramLine(f, maxdist = max(v$dist))
  k <- krige(formula = get(param)~1, locations = d, newdata = s, model = f) %>% select(-var1.var)
  if(plot == TRUE){
    p.h <- ggplot() +
      geom_histogram(data = d, aes_string(x = param, y = '..density..'), alpha = 0.8, color = 'white') +
      labs(x = bquote(atop('Heat Flow', mW/m^2)), y = 'Density') +
      theme_bw()
    p.v <- ggplot() +
      geom_point(data = v, aes(x = dist/1000, y = gamma), size = 0.5) +
      geom_path(data = v.m, aes(x = dist/1000, y = gamma)) +
      labs(x = 'Lag (km)', y = 'Semivariance') +
      coord_cartesian(ylim = c(0, NA)) +
      theme_bw()
    p <- ggplot() +
      geom_sf(data = bbox_widen(st_bbox(b), crs = crs, c('left' = 0.1, 'right' = 0.1, 'top' = 0.1, 'bottom' = 0.1)), fill = 'cornflowerblue', alpha = 0.2, color = NA) +
      geom_stars(data = st_rasterize(k) %>% st_crop(b)) +
      geom_sf(data = d, aes_string(color = param), size = 0.5, show.legend = F) +
      geom_sf(data = contours, color = 'white', alpha = 0.8, size = 0.3) +
      geom_sf(data = seg, fill = NA, color = 'black', size = 1.25) +
      labs(x = NULL, y = NULL) +
      scale_color_viridis_c(option = 1) +
      scale_fill_viridis_c(name = bquote(atop('Heat Flow', mW/m^2)), option = 1, na.value = 'transparent') +
      theme_map(font_size = 11) + 
      theme(
        axis.text = element_text(),
        panel.border = element_blank(),
        panel.grid = element_line(size = 0.25, color = rgb(0.1, 0.1, 0.1, 0.5)),
        panel.background = element_blank(),
        panel.ontop = TRUE,
        plot.background = element_rect(fill = "transparent", color = NA)
      )
    return(list(variogram = v, v.model = f, krg = k, hist = p.h, k.plot = p, v.plot = p.v))
  } else {
    return(list(variogram = v, v.model = f, krg = k))
  }
}

plot.subseg <- function(object.subseg, object.hf, param, krg.shp, seg, contours, crs) {
  d <- 
  b <- 
  s <- 
  
  
}