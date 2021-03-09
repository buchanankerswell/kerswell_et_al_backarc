source('functions.R')
cat('Loading segment and heat flow data\n')
load('data/hf.RData')

# Define bounding boxes
# All segments
shp.pacific.box <-
  bbox_widen(st_bbox(shp.sa.segs.robin.pacific),
             proj4.robin.pacific,
             borders = c('left' = 0.05, 'right' = 0.05, 'top' = 0.05, 'bottom' = 0.05))

# Intividual segments
shp.segs.box <- 
  shp.sa.segs.robin.pacific.buffer$geometry %>% 
  purrr::map(~st_bbox(.x)) %>% 
  purrr::set_names(shp.sa.segs.robin.pacific$segment)

cat('\nPlotting ...', 'Saving figures to figs/base/')
cat('\nPlotting world map with segments...')
# World map
p <- ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.grats30, alpha = 0.1) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 3,
                      alpha = 0.8,
                      max.overlaps = 30) +
  labs(title = 'Subduction Zone Segments', subtitle = 'from Syracuse & Abers (2006)') +
  coord_sf(expand = F) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.grid = element_line(size = 0.2, color = rgb(0,0,0,0.2)),
        panel.ontop = T)
ggsave(filename = paste0('figs/base/world.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 - st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

cat('\nPlotting world heat flow map...')
# All Heat Flow
p <- ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.hf %>% filter(`heat-flow (mW/m2)` < 250), aes(color = `heat-flow (mW/m2)`), size = 0.1) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(option = 'magma') +
  labs(title = 'New Global Heat Flow', subtitle = bquote(HF < 250 ~ mWm^-2 ~~ 'from Lucazeau (2019)')) +
  guides(color = guide_colorbar(title.vjust = 1, barwidth = unit(abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6)/5, 'in'))) +
  coord_sf(expand = F) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'center',
        panel.grid = element_line(size = 0.2, color = rgb(0,0,0,0.2)),
        panel.ontop = T)
ggsave(filename = paste0('figs/base/hf.png'),
       plot = p,
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 - st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

cat('\nPlotting world heat flow map zoomed to Pacific...')
# All Heat Flow Pacific Zoom
p <- ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf(data = shp.hf %>% filter(`heat-flow (mW/m2)` < 250),
          aes(color = `heat-flow (mW/m2)`),
          size = 0.01) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 3,
                      alpha = 0.8,
                      force = 50,
                      max.overlaps = 30) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(option = 'magma') +
  guides(color = guide_colorbar(title.vjust = 1, barwidth = unit(abs(st_bbox(shp.pacific.box)$xmax/2.5e6 - st_bbox(shp.pacific.box)$xmin/2.5e6)/3, 'in'))) +
  coord_sf(expand = F,
           xlim = c(st_bbox(shp.pacific.box)$xmin, st_bbox(shp.pacific.box)$xmax),
           ylim = c(st_bbox(shp.pacific.box)$ymin, st_bbox(shp.pacific.box)$ymax)) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'center',
        panel.grid = element_line(size = 0.2, color = rgb(0,0,0,0.2)),
        panel.ontop = T)
ggsave(filename = paste0('figs/base/hf_zoom.png'),
       device = 'png',
       plot = p,
       type = 'cairo',
       width = abs(st_bbox(shp.pacific.box)$xmax/2.5e6 - st_bbox(shp.pacific.box)$xmin/2.5e6),
       height = abs(st_bbox(shp.pacific.box)$ymax/2.5e6 - st_bbox(shp.pacific.box)$ymin/2.5e6))

cat('\nPlotting Lucazeau (2019) predicted heat flow ...')
# Lucazeau (2019) predicted
p <- shp.hf.pred %>%
  ggplot() +
  geom_sf(aes(color = HF_pred)) +
  geom_sf(data = shp.world.robin.pacific, color = rgb(0, 0, 0, 0.7), fill = NA, size = 0.2) +
  labs(title = 'Predicted Heat Flow', subtitle = 'from Lucazeau (2019)', color = bquote(mWm^-2)) +
  scale_color_viridis_c(option = 'magma', trans = 'log10') +
  guides(color = guide_colorbar(title.vjust = 1, barwidth = unit(abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6)/5, 'in'))) +
  coord_sf(expand = F) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'center',
        panel.grid = element_line(size = 0.2, color = rgb(1,1,1,0.2)),
        panel.ontop = T)
ggsave(filename = paste0('figs/base/hf_pred.png'),
       device = 'png',
       plot = p,
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 - st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

cat('\nPlotting individual segments ...')
# Individual Segments
purrr::walk2(shp.segs.box, shp.sa.segs.robin.pacific$segment, ~{
  cat('\nPlotting', .y)
  hf <- shp.hf %>% st_crop(.x) %>% filter(`heat-flow (mW/m2)` < 250)
  p <- ggplot() +
    geom_sf(data = shp.world.robin.pacific, color = NA) +
    geom_sf(data = shp.sa.segs.robin.pacific %>% filter(segment == .y), size = 1.1) +
    geom_sf(data = shp.sa.countours.robin.pacific %>% filter(segment == .y), size = 0.1) +
    geom_sf(data = hf,
            aes(color = `heat-flow (mW/m2)`),
            size = 1) +
    labs(color = bquote(mWm^-2), title = .y, subtitle = bquote(NGHF<~250~mWm^-2~~n==.(nrow(hf)))) +
    scale_color_viridis_c(option = 'magma') +
    guides(color = guide_colorbar(title.vjust = 1, barwidth = unit(abs(st_bbox(.x)$xmax/1e6 - st_bbox(.x)$xmin/1e6)/1.5, 'in'))) +
    coord_sf(expand = F,
             xlim = c(st_bbox(.x)$xmin, st_bbox(.x)$xmax),
             ylim = c(st_bbox(.x)$ymin, st_bbox(.x)$ymax)) +
    theme_map(font_size = 11) +
    theme(axis.text = element_text(),
          legend.position = 'bottom',
          legend.justification = 'center',
          legend.direction = 'horizontal',
          panel.grid = element_line(size = 0.1, color = rgb(0,0,0,0.2)),
          panel.ontop = T)
  ggsave(filename = paste0('figs/base/', .y, '.png'),
         plot = p,
         device = 'png',
         type = 'cairo',
         width = abs(st_bbox(.x)$xmax/5e5 - st_bbox(.x)$xmin/5e5),
         height = abs(st_bbox(.x)$ymax/5e5 - st_bbox(.x)$ymin/5e5))
})
cat('\nDone')