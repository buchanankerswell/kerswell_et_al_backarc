suppressMessages(source('../functions.R'))
load('../data/hf.RData')

# Compile and transform data
cat('Compiling and transforming spatial data for N. Philippines:
    draw segment
    draw buffer
    draw contours
    crop heat flow
    draw bounding box
    calculate distances from data points to trench
    get volcano positions
    \n')
# N. Philippines

# Kriging options
lags <- 20
lag.cutoff <- 6
ngrid <- 1e3
grid.method <- 'hexagonal'
rotate.angle <- 8
v.mod <- c('Sph', 'Gau', 'Exp', 'Lin')
idp <- 2

cat('Kriging N. Philippines heat flow...\nKriging options:',
    'lags:', lags,
    'lag cutoff (proportion of max distance):', paste0('1/', lag.cutoff),
    'grid size:', ngrid,
    'sampling method:', grid.method,
    'variogram model:', v.mod,
    'inverse weighting power:', idp,
    sep = '\n')

# Segment
shp.seg.np <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'N. Philippines')

# Flat buffer
shp.buf.np <- shp.sa.segs.robin.pacific %>%
  filter(segment == 'N. Philippines') %>%
  st_buffer(1000000, endCapStyle = 'FLAT')

# Contours
shp.con.np <- shp.sa.countours.robin.pacific %>%
  filter(segment == 'N. Philippines')

# Heat flow
shp.hf.np <- shp.hf %>% st_join(shp.buf.np, left = F) %>%
  relocate(segment, .before = geometry)

# Volcanoes
suppressWarnings(
  shp.volc.np <- rbind(
    shp.volc %>% st_intersection(shp.buf.np),
    shp.volc.no.spread %>% st_intersection(shp.buf.np)))

# Rotated projection
proj.rot <- rotate_proj(shp.hf.np, rotate.angle)

# Kriging
# Krige entire segment
k <- model_variogram(
  data = shp.hf.np,
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
cat('\nSaving variogram to figs/np/variogram.png\n')
ggsave(filename = '../figs/np/variogram.png',
       plot = k$variogram.plot + labs(title = 'N. Philippines Variogram'),
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
cat('Drawing interpolation sampling grid\n')
p1 <- k$interp.results %>% 
  ggplot() +
  geom_sf(data = shp.seg.np) +
  geom_sf(data = shp.hf.np, aes(color = hf), size = 0.3) +
  geom_sf(data = shp.con.np, color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.np, color = 'black') +
  geom_sf(data = shp.volc.np, color = 'black', shape = 2) +
  geom_sf(size = 0.2, shape = 15) +
  labs(color = bquote(mWm^-2),
       shape = 'Position',
       title = 'N. Philippines Interpolation Grid') +
  coord_sf(crs = proj.rot) +
  scale_color_viridis_c(option = 'magma') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.np)$xmax/1.5e5 -
                                                      st_bbox(shp.hf.np)$xmin/1.5e5)/3,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow grid to figs/np/hf_grid.png\n')
ggsave(filename = paste0('../figs/np/hf_grid.png'),
       plot = p1,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.np)$xmax/1.5e5 - st_bbox(shp.hf.np)$xmin/1.5e5),
       height = abs(st_bbox(shp.hf.np)$ymax/1.5e5 - st_bbox(shp.hf.np)$ymin/1.5e5),
       bg = 'transparent')
# Draw interpolated results
cat('Drawing interpolated results\n')
p2 <- k$interp.results %>%
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = k$interp.results$hf) %>%
  ggplot() +
  geom_contour_fill(aes(x = X, y = Y, z = hf)) +
  geom_sf(data = shp.hf.np %>% st_transform(proj.rot), aes(color = hf), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.np %>% st_transform(proj.rot), color = 'white', alpha = 0.5) +
  geom_sf(data = shp.seg.np %>% st_transform(proj.rot), color = 'white') +
  geom_sf(data = shp.volc.np %>% st_transform(proj.rot), color = 'white', shape = 2) +
  labs(fill = bquote(mWm^-2), title = 'N. Philippines Interpolated Heat Flow') +
  scale_color_viridis_c(option = 'magma') +
  scale_fill_viridis_c(option = 'magma') +
  guides(fill = guide_colorbar(title.vjust = 1,
                               barwidth = unit(abs(st_bbox(shp.hf.np)$xmax/1.5e5 -
                                                     st_bbox(shp.hf.np)$xmin/1.5e5)/3,
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
cat('Saving interpolated heat flow results to figs/np/hf_interp.png\n')
ggsave(filename = paste0('../figs/np/hf_interp.png'),
       plot = p2,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.np)$xmax/1.5e5 - st_bbox(shp.hf.np)$xmin/1.5e5),
       height = abs(st_bbox(shp.hf.np)$ymax/1.5e5 - st_bbox(shp.hf.np)$ymin/1.5e5),
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
  geom_sf(data = shp.hf.np %>% st_transform(proj.rot), size = 0.3, show.legend = F) +
  geom_sf(data = shp.con.np %>% st_transform(proj.rot), color = 'black', alpha = 0.5) +
  geom_sf(data = shp.seg.np %>% st_transform(proj.rot), color = 'black') +
  geom_sf(data = shp.volc.np %>% st_transform(proj.rot), color = 'black', shape = 2) +
  coord_sf(crs = proj.rot) +
  labs(color = 'Transect', title = 'N. Philippines Transects') +
  scale_color_viridis_c(option = 'viridis') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.hf.np)$xmax/1.5e5 -
                                                      st_bbox(shp.hf.np)$xmin/1.5e5)/3,
                                                'in'))) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        panel.ontop = T)
cat('Saving heat flow transect map to figs/np/hf_trans_map.png\n')
ggsave(filename = paste0('../figs/np/hf_trans_map.png'),
       plot = p3,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.hf.np)$xmax/1.5e5 - st_bbox(shp.hf.np)$xmin/1.5e5),
       height = abs(st_bbox(shp.hf.np)$ymax/1.5e5 - st_bbox(shp.hf.np)$ymin/1.5e5),
       bg = 'transparent')
p4 <- t %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  mutate(hf = t$hf,
         group = t$group) %>% 
  group_by(group) %>% 
  ggplot() +
  geom_path(aes(x = X, y = hf, group = group, color = group), size = 0.3, show.legend = F) +
  scale_color_viridis_c(option = 'viridis') +
  labs(x = 'Longitude',
       y = bquote(Heat~Flow~mWm^-2),
       color = 'Transect',
       title = 'N. Philippines Transects') +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color = rgb(0.5, 0.5, 0.5, 0.5)),
        panel.background = element_rect(fill = 'transparent', color = NA),
        plot.background = element_rect(fill = 'transparent', color = NA),
        legend.background = element_rect(fill = 'transparent', color = NA),
        title = element_text(face = 'bold'))
cat('Saving heat flow transects to figs/np/hf_trans.png\n')
ggsave(filename = paste0('../figs/np/hf_trans.png'),
       plot = p4,
       device = 'png',
       type = 'cairo',
       width = 7,
       height = 4,
       bg = 'transparent')
# Composition
cat('Drawing composition\n')
p <- (p2 + p3) / ((k$variogram.plot +
                     labs(title = 'N. Philippines Variogram') +
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
cat('Saving composition to figs/np/comp.png\n')
ggsave(filename = paste0('../figs/np/comp.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = 13,
       height = 7,
       bg = 'transparent')
