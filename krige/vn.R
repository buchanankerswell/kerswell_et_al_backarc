suppressMessages(source('../functions.R'))
load('../data/hf.RData')

# Compile and transform data
cat('Compiling and transforming spatial data for Vanuatu:
    draw segment
    draw buffer
    draw contours
    crop heat flow
    draw bounding box
    calculate distances from data points to trench
    get volcano positions
    \n')
# Vanuatu

# Kriging options
lags <- 15
lag.cutoff <- 3
ngrid <- 1e3
grid.method <- 'hexagonal'
rotate.angle <- -20
v.mod <- c('Sph', 'Gau', 'Exp', 'Lin')
idp <- 2

cat('Kriging Vanuatu heat flow...\nKriging options:',
    'lags:', lags,
    'lag cutoff (proportion of max distance):', paste0('1/', lag.cutoff),
    'grid size:', ngrid,
    'sampling method:', grid.method,
    'variogram model:', v.mod,
    'inverse weighting power:', idp,
    sep = '\n')

# Segment
shp.seg.vn <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Vanuatu')

# Flat buffer
shp.buf.vn <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'Vanuatu') %>%
  st_buffer(1000000, endCapStyle = 'FLAT')

# Contours
shp.con.vn <- shp.sa.countours.robin.pacific %>%
  filter(segment == 'Vanuatu')

# Heat flow
shp.hf.vn <- shp.hf %>% st_join(shp.buf.vn, left = F) %>%
  relocate(segment, .before = geometry)

# Volcanoes
suppressWarnings(
  shp.volc.vn <- rbind(
    shp.volc %>% st_intersection(shp.buf.vn),
    shp.volc.no.spread %>% st_intersection(shp.buf.vn)))

# Rotated projection
proj.rot <- rotate_proj(shp.hf.vn, rotate.angle)

# Kriging
# Krige entire segment
k <- model_variogram(
  data = shp.hf.vn,
  lags = lags,
  lag.cutoff = lag.cutoff,
  param = 'hf',
  ngrid = ngrid,
  grid.method = grid.method,
  grid.rotate = rotate.angle,
  grid.shape = c('left' = 0.2, 'right' = 0.1, 'top' = 0.05, 'bottom' = 0.05),
  v.mod = vgm(model = v.mod),
  krige = T,
  interp.power = 2
)
cat('\nSaving variogram to figs/vn/variogram.png\n')
ggsave(filename = '../figs/vn/variogram.png',
       plot = k$variogram.plot + labs(title = 'Vanuatu Variogram'),
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
cat('Drawing interpolation sampling grid\n')
p1 <- k$interp.results %>% 
  ggplot() +
  geom_sf(data = shp.seg.vn) +
  geom_sf(data = shp.hf.vn, aes(color = hf), size = 0.3) +
  geom_sf(data = shp.con.vn, color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.vn, color = 'black') +
  geom_sf(data = shp.volc.vn, color = 'black', shape = 2) +
  geom_sf(size = 0.2, shape = 15) +
  labs(color = bquote(mWm^-2),
       shape = 'Position',
       title = 'Vanuatu Interpolation Grid') +
  coord_sf(crs = proj.rot) +
  scale_color_viridis_c(option = 'magma') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.vn)$xmax/3.2e5 -
                                                      st_bbox(shp.hf.vn)$xmin/3.2e5)/3,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow grid to figs/vn/hf_grid.png\n')
ggsave(filename = paste0('../figs/vn/hf_grid.png'),
       plot = p1,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.vn)$xmax/3.2e5 - st_bbox(shp.hf.vn)$xmin/3.2e5),
       height = abs(st_bbox(shp.hf.vn)$ymax/3.2e5 - st_bbox(shp.hf.vn)$ymin/3.2e5),
       bg = 'transparent')
# Draw interpolated results
cat('Drawing interpolated results\n')
p2 <- k$interp.results %>%
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = k$interp.results$hf) %>%
  ggplot() +
  geom_contour_fill(aes(x = X, y = Y, z = hf)) +
  geom_sf(data = shp.hf.vn %>% st_transform(proj.rot), aes(color = hf), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.vn %>% st_transform(proj.rot), color = 'white', alpha = 0.5) +
  geom_sf(data = shp.seg.vn %>% st_transform(proj.rot), color = 'white') +
  geom_sf(data = shp.volc.vn %>% st_transform(proj.rot), color = 'white', shape = 2) +
  labs(fill = bquote(mWm^-2), title = 'Vanuatu Interpolated Heat Flow') +
  scale_color_viridis_c(option = 'magma') +
  scale_fill_viridis_c(option = 'magma') +
  guides(fill = guide_colorbar(title.vjust = 1,
                               barwidth = unit(abs(st_bbox(shp.hf.vn)$xmax/3.2e5 -
                                                     st_bbox(shp.hf.vn)$xmin/3.2e5)/3,
                                               'in'))) +
  coord_sf(crs = proj.rot) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving interpolated heat flow results to figs/vn/hf_interp.png\n')
ggsave(filename = paste0('../figs/vn/hf_interp.png'),
       plot = p2,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.vn)$xmax/3.2e5 - st_bbox(shp.hf.vn)$xmin/3.2e5),
       height = abs(st_bbox(shp.hf.vn)$ymax/3.2e5 - st_bbox(shp.hf.vn)$ymin/3.2e5),
       bg = 'transparent')
# Heat flow transects
cat('Drawing heat flow transects\n')
t <- k$interp.results %>%
  st_coordinates() %>%
  as_tibble() %>%
  mutate(hf = k$interp.results$hf) %>%
  group_by(Y) %>%
  mutate(group = cur_group_id()) %>% 
  st_as_sf(coords = c(1,2), crs = proj.rot) %>% 
  st_transform(crs = proj4.wgs)
p3 <- t %>% ggplot() +
  geom_sf(aes(color = group), size = 0.2, shape = 15) +
  geom_sf(data = shp.hf.vn %>% st_transform(proj.rot), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.vn %>% st_transform(proj.rot), color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.vn %>% st_transform(proj.rot), color = 'black') +
  geom_sf(data = shp.volc.vn %>% st_transform(proj.rot), color = 'black', shape = 2) +
  coord_sf(crs = proj.rot) +
  labs(color = 'Transect', title = 'Vanuatu Transects') +
  scale_color_viridis_c(option = 'viridis') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.vn)$xmax/3.2e5 -
                                                      st_bbox(shp.hf.vn)$xmin/3.2e5)/3,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow transect map to figs/vn/hf_trans_map.png\n')
ggsave(filename = paste0('../figs/vn/hf_trans_map.png'),
       plot = p3,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.vn)$xmax/3.2e5 - st_bbox(shp.hf.vn)$xmin/3.2e5),
       height = abs(st_bbox(shp.hf.vn)$ymax/3.2e5 - st_bbox(shp.hf.vn)$ymin/3.2e5),
       bg = 'transparent')
p4 <- t %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = t$hf,
         'X' = ifelse(X < 0, X + 360, X),
         group = t$group) %>% 
  group_by(group) %>%
  ggplot() +
  geom_path(aes(x = X, y = hf, group = group, color = group), size = 0.3, show.legend = F) +
  scale_color_viridis_c(option = 'viridis') +
  labs(x = 'Longitude',
       y = bquote(Heat~Flow~mWm^-2),
       color = 'Transect',
       title = 'Vanuatu Transects') +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        title = element_text(face = 'bold'))
cat('Saving heat flow transects to figs/vn/hf_trans.png\n')
ggsave(filename = paste0('../figs/vn/hf_trans.png'),
       plot = p4,
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
# Composition
cat('Drawing composition\n')
p <- (p2 + p3) / ((k$variogram.plot +
                     labs(title = 'Vanuatu Variogram') +
                     theme(title = element_text(face = 'bold'))) + p4) +
  plot_annotation(tag_levels = 'a',
                  theme = theme(plot.margin = margin(),
                                panel.background = element_rect(fill = 'transparent', color = NA),
                                plot.background = element_rect(fill = 'transparent', color = NA),
                                legend.background = element_rect(fill = 'transparent', color = NA))) +
  plot_layout(guides = 'collect',
              heights = c(3, 1)) &
  theme(legend.position = 'bottom',
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA))
cat('Saving composition to figs/vn/comp.png\n')
ggsave(filename = paste0('../figs/vn/comp.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 10,
       height = 10,
       bg = 'transparent')
