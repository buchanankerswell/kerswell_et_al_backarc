source('functions.R')
cat('Loading segment and heat flow data\n')
load('data/hf.RData')

# World map
ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.grats30, alpha = 0.1) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 3,
                      alpha = 0.5,
                      max.overlaps = 30) +
  coord_sf(expand = F) +
  theme_map(font_size = 11) +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.background = element_rect(fill = rgb(0.392, 0.584, 0.929, 0.1)))
ggsave(filename = paste0('figs/base/world.png'),
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 - st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

# All Heat Flow
ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.grats30, alpha = 0.1) +
  geom_sf(data = shp.hf, aes(color = hf), size = 0.1) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(option = 'magma', trans = 'log10') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 -
                                                      st_bbox(shp.world.robin.pacific)$xmin/2.5e6)/5,
                                                'in'))) +
  coord_sf(expand = F) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.background = element_rect(fill = rgb(0.392, 0.584, 0.929, 0.1)))
ggsave(filename = paste0('figs/base/all_hf.png'),
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(shp.world.robin.pacific)$xmax/2.5e6 - st_bbox(shp.world.robin.pacific)$xmin/2.5e6),
       height = abs(st_bbox(shp.world.robin.pacific)$ymax/2.5e6 - st_bbox(shp.world.robin.pacific)$ymin/2.5e6))

# All Heat Flow Zoomed
ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.grats30, alpha = 0.1) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf(data = shp.hf, aes(color = hf), size = 0.1) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 3,
                      alpha = 0.8,
                      force = 50,
                      max.overlaps = 30) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(option = 'magma', trans = 'log10') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(box.segs.all.wide)$xmax/2.5e6 -
                                                      st_bbox(box.segs.all.wide)$xmin/2.5e6)/5,
                                                'in'))) +
  coord_sf(expand = F,
           xlim = c(st_bbox(box.segs.all.wide)$xmin, st_bbox(box.segs.all.wide)$xmax),
           ylim = c(st_bbox(box.segs.all.wide)$ymin, st_bbox(box.segs.all.wide)$ymax)) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.background = element_rect(fill = rgb(0.392, 0.584, 0.929, 0.1)))
ggsave(filename = paste0('figs/base/all_hf_zoom.png'),
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(box.segs.all.wide)$xmax/2.5e6 - st_bbox(box.segs.all.wide)$xmin/2.5e6),
       height = abs(st_bbox(box.segs.all.wide)$ymax/2.5e6 - st_bbox(box.segs.all.wide)$ymin/2.5e6))

# Buffered Heat Flow Zoomed
ggplot() +
  geom_sf(data = shp.world.robin.pacific, color = NA) +
  geom_sf(data = shp.grats30, alpha = 0.1) +
  geom_sf(data = shp.sa.segs.robin.pacific.buffer$km.1000, fill = NA, color = 'black', size = 0.1) +
  geom_sf(data = shp.sa.segs.robin.pacific, size = 1.1) +
  geom_sf(data = shp.sa.countours.robin.pacific, size = 0.1) +
  geom_sf(data = shp.hf %>% st_join(shp.sa.segs.robin.pacific.buffer$km.1000, left = F), aes(color = hf), size = 0.1) +
  geom_sf_label_repel(data = shp.sa.segs.robin.pacific,
                      aes(label = segment),
                      size = 3,
                      alpha = 0.8,
                      force = 50,
                      max.overlaps = 30) +
  labs(color = bquote(mWm^-2)) +
  scale_color_viridis_c(option = 'magma', trans = 'log10') +
  guides(color = guide_colorbar(title.vjust = 1,
                                barwidth = unit(abs(st_bbox(box.segs.all.wide)$xmax/2.5e6 -
                                                      st_bbox(box.segs.all.wide)$xmin/2.5e6)/5,
                                                'in'))) +
  coord_sf(expand = F,
           xlim = c(st_bbox(box.segs.all.wide)$xmin, st_bbox(box.segs.all.wide)$xmax),
           ylim = c(st_bbox(box.segs.all.wide)$ymin, st_bbox(box.segs.all.wide)$ymax)) +
  theme_map() +
  theme(axis.text = element_text(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        panel.background = element_rect(fill = rgb(0.392, 0.584, 0.929, 0.1)))
ggsave(filename = paste0('figs/base/hf_buff.png'),
       device = 'png',
       type = 'cairo',
       width = abs(st_bbox(box.segs.all.wide)$xmax/2.5e6 - st_bbox(box.segs.all.wide)$xmin/2.5e6),
       height = abs(st_bbox(box.segs.all.wide)$ymax/2.5e6 - st_bbox(box.segs.all.wide)$ymin/2.5e6))

# Individual Segments
purrr::walk2(box.segs.wide, shp.sa.segs.robin.pacific$segment, ~{
  ggplot() +
    geom_sf(data = shp.world.robin.pacific, color = NA) +
    geom_sf(data = shp.grats30, alpha = 0.1) +
    geom_sf(data = shp.sa.segs.robin.pacific.buffer$km.1000 %>% filter(segment == .y),
            fill = NA, color = 'black', size = 0.1) +
    geom_sf(data = shp.sa.segs.robin.pacific %>% filter(segment == .y), size = 1.1) +
    geom_sf(data = shp.sa.countours.robin.pacific %>% filter(segment == .y), size = 0.1) +
    geom_sf(data = shp.hf %>%
              st_join(shp.sa.segs.robin.pacific.buffer$km.1000 %>% filter(segment == .y), left = F),
            aes(color = hf), size = 1) +
    labs(color = bquote(mWm^-2), title = .y) +
    scale_color_viridis_c(option = 'magma') +
    guides(color = guide_colorbar(title.vjust = 1, barwidth = unit(abs(st_bbox(.x)$xmax/1e6 - st_bbox(.x)$xmin/1e6)/3, 'in'))) +
    coord_sf(expand = F,
             xlim = c(st_bbox(.x)$xmin, st_bbox(.x)$xmax),
             ylim = c(st_bbox(.x)$ymin, st_bbox(.x)$ymax)) +
    theme_map(font_size = 11) +
    theme(axis.text = element_text(),
          legend.position = 'bottom',
          legend.justification = 'right',
          legend.direction = 'horizontal',
          panel.grid = element_blank(),
          panel.background = element_rect(fill = rgb(0.392, 0.584, 0.929, 0.1)),
          panel.ontop = T)
  ggsave(filename = paste0('figs/base/', .y, '.png'),
         device = 'png',
         type = 'cairo',
         width = abs(st_bbox(.x)$xmax/1e6 - st_bbox(.x)$xmin/1e6),
         height = abs(st_bbox(.x)$ymax/1e6 - st_bbox(.x)$ymin/1e6))
})