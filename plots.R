rm(list = ls())
source('functions.R')
load('data/hf.RData')

# Open PDF device
cairo_pdf('figs/maps/hf.pdf', onefile = T, width = 10, height = 7)

# Pacific centered world map
# World plot [all]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.buffer = shp.sa.segs.robin.pacific.buffer,
  seg.label.size = 2,
  seg.names = TRUE,
  heatflow = shp.hf,
  hf.size = 0.005,
  volcano = list(shp.volc, shp.volc.no.spread),
  volc.size.range = c(0.01, 3),
  volc.size = 0.1,
  volc.alpha = 1,
  volc.color = 'deeppink',
  box.labels = box.segs.all.wide,
  box.limits = shp.grats30,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/full.comp/global.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# World plot [segments]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.buffer = shp.sa.segs.robin.pacific.buffer,
  seg.label.size = 2,
  box.labels = box.segs.all.wide,
  box.limits = shp.grats30,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/segments/global_segments.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# World plot [segments + heatflow]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.buffer = shp.sa.segs.robin.pacific.buffer,
  seg.label.size = 2,
  heatflow = shp.hf,
  hf.size = 0.005,
  box.labels = box.segs.all.wide,
  box.limits = shp.grats30,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/heatflow/global_segs_hf.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# World plot [segments + volcanoes]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.label.size = 2,
  volcano = list(shp.volc, shp.volc.no.spread),
  volc.size.range = c(0.01, 3),
  volc.size = 0.3,
  volc.alpha = 1,
  volc.color = 'deeppink',
  box.labels = box.segs.all.wide,
  box.limits = shp.grats30,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/volcano/global_segs_volc.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# Pacific centered regional map
# Pacific Region Map [All]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.buffer = shp.sa.segs.robin.pacific.buffer,
  seg.label.size = 2,
  heatflow = shp.hf,
  hf.size = 0.01,
  volcano = list(shp.volc, shp.volc.no.spread),
  volc.size.range = c(0.25, 3),
  volc.size = 0.2,
  volc.alpha = 1,
  volc.color = 'deeppink',
  box.labels = box.segs.all,
  box.limits = box.segs.all.wide,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/full.comp/pacific.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# Pacific Region Map [Segments]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.label.size = 2,
  box.labels = box.segs.all,
  box.limits = box.segs.all.wide,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/segments/pacific_segments.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# Pacific Region Map [Segments + heatflow]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.2,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.buffer = shp.sa.segs.robin.pacific.buffer,
  seg.label.size = 2,
  heatflow = shp.hf,
  hf.size = 0.01,
  box.labels = box.segs.all,
  box.limits = box.segs.all.wide,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/heatflow/pacific_segs_hf.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# Pacific Region Map [Segments + volcanoes]
plot_hf(
  grats = shp.grats30,
  world = shp.world.robin.pacific,
  box = box.segs.all.wide,
  box.col = 'cornflowerblue',
  box.alpha = 0.15,
  seg = shp.sa.segs.robin.pacific,
  seg.contours = shp.sa.countours.robin.pacific,
  seg.label.size = 2,
  volcano = list(shp.volc, shp.volc.no.spread),
  volc.size.range = c(0.25, 3),
  volc.size = 0.3,
  volc.alpha = 1,
  volc.color = 'deeppink',
  box.labels = box.segs.all,
  box.limits = box.segs.all.wide,
  legend.pos = 'top',
  legend.just = 'right',
  legend.size = unit(8, 'lines')
)
ggsave('figs/maps/volcano/pacific_segs_volc.png', device = 'png', dpi = 320, bg = 'transparent', width = 10, height = 7)

# Individual arc segment closeups

# Split segments into long (portrait) and wide (landscape) groups
# Portrait list
portrait.list.boxes <- box.segs %>% purrr::keep(names(.) %in% c('Andes', 'Kamchatka Marianas', 'Lesser Antilles', 'N. Philippines', 'S. Philippines', 'Tonga New Zealand', 'Vanuatu', 'Kyushu Ryukyu', 'Scotia'))
portrait.list.boxes.wide <- box.segs.wide %>% purrr::keep(names(.) %in% c('Andes', 'Kamchatka Marianas', 'Lesser Antilles', 'N. Philippines', 'S. Philippines', 'Tonga New Zealand', 'Vanuatu', 'Kyushu Ryukyu', 'Scotia'))
portrait.list.names <- names(portrait.list.boxes)
portrait.list.fnames <- portrait.list.names %>% purrr::map_chr(~{.x %>% gsub(' ', '', .)})
# Landscape list
landscape.list.boxes <- box.segs %>% purrr::keep(!(names(.) %in% c('Andes', 'Kamchatka Marianas', 'Lesser Antilles', 'N. Philippines', 'S. Philippines', 'Tonga New Zealand', 'Vanuatu', 'Kyushu Ryukyu', 'Scotia')))
landscape.list.boxes.wide <- box.segs.wide %>% purrr::keep(!(names(.) %in% c('Andes', 'Kamchatka Marianas', 'Lesser Antilles', 'N. Philippines', 'S. Philippines', 'Tonga New Zealand', 'Vanuatu', 'Kyushu Ryukyu', 'Scotia')))
landscape.list.names <- names(landscape.list.boxes)
landscape.list.fnames <- landscape.list.names %>% purrr::map_chr(~{.x %>% gsub(' ', '', .)})

# Plotting
# Portrait segments [all]
purrr::pmap(
  list(
    portrait.list.boxes,
    portrait.list.boxes,
    portrait.list.names,
    portrait.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.buffer = shp.sa.segs.robin.pacific.buffer,
      seg.names = FALSE,
      heatflow = shp.hf,
      hf.size = 1,
      volcano = list(shp.volc, shp.volc.no.spread),
      volc.size.range = c(1, 6),
      volc.size = 0.4,
      volc.alpha = 0.5,
      volc.color = 'deeppink',
      box.labels = box.lab,
      box.limits = ..2,
      legend.pos = 'right',
      legend.just = 'left',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/full.comp/', ..4, '.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 7,
      height = 10
    )
  }
)

# Portrait segments [segments]
purrr::pmap(
  list(
    portrait.list.boxes,
    portrait.list.boxes,
    portrait.list.names,
    portrait.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.names = FALSE,
      box.limits = ..2,
      legend.pos = 'right',
      legend.just = 'left',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/segments/', ..4, '_segment.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 7,
      height = 10
    )
  }
)

# Portrait segments [segments + heatflow]
purrr::pmap(
  list(
    portrait.list.boxes,
    portrait.list.boxes,
    portrait.list.names,
    portrait.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.buffer = shp.sa.segs.robin.pacific.buffer,
      seg.names = FALSE,
      heatflow = shp.hf,
      hf.size = 1,
      box.limits = ..2,
      legend.pos = 'right',
      legend.just = 'left',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/heatflow/', ..4, '_heatflow.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 7,
      height = 10
    )
    }
)

# Portrait segments [segments + volcanoes]
purrr::pmap(
  list(
    portrait.list.boxes,
    portrait.list.boxes,
    portrait.list.names,
    portrait.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.names = FALSE,
      volcano = list(shp.volc, shp.volc.no.spread),
      volc.size.range = c(1, 6),
      volc.size = 0.4,
      volc.alpha = 0.5,
      volc.color = 'deeppink',
      box.limits = ..2,
      legend.pos = 'right',
      legend.just = 'left',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/volcano/', ..4, '_volc.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 7,
      height = 10
    )
    }
)

# Plotting
# Landscape segments [all]
purrr::pmap(
  list(
    landscape.list.boxes,
    landscape.list.boxes,
    landscape.list.names,
    landscape.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.buffer = shp.sa.segs.robin.pacific.buffer,
      seg.names = FALSE,
      heatflow = shp.hf,
      hf.size = 1,
      volcano = list(shp.volc, shp.volc.no.spread),
      volc.size.range = c(1, 6),
      volc.size = 0.4,
      volc.alpha = 1,
      volc.color = 'deeppink',
      box.labels = box.lab,
      box.limits = ..2,
      legend.pos = 'top',
      legend.just = 'right',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/full.comp/', ..4, '.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 10,
      height = 7
    )
  }
)

# Landscape segments [segments]
purrr::pmap(
  list(
    landscape.list.boxes,
    landscape.list.boxes,
    landscape.list.names,
    landscape.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.names = FALSE,
      box.limits = ..2,
      legend.pos = 'top',
      legend.just = 'right',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/segments/', ..4, '_segment.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 10,
      height = 7
    )
  }
)

# Landscape segments [segments + heatflow]
purrr::pmap(
  list(
    landscape.list.boxes,
    landscape.list.boxes,
    landscape.list.names,
    landscape.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.buffer = shp.sa.segs.robin.pacific.buffer,
      seg.names = FALSE,
      heatflow = shp.hf,
      hf.size = 1,
      box.limits = ..2,
      legend.pos = 'top',
      legend.just = 'right',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/heatflow/', ..4, '_heatflow.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 10,
      height = 7
    )
  }
)

# Landscape segments [segments + volcanoes]
purrr::pmap(
  list(
    landscape.list.boxes,
    landscape.list.boxes,
    landscape.list.names,
    landscape.list.fnames
  ),
  ~ {
    plot_hf(
      grats = shp.grats30,
      world = shp.world.robin.pacific10,
      box = ..1,
      box.col = 'cornflowerblue',
      box.alpha = 0.15,
      seg = shp.sa.segs.robin.pacific,
      seg.contours = shp.sa.countours.robin.pacific,
      seg.names = FALSE,
      volcano = list(shp.volc, shp.volc.no.spread),
      volc.size.range = c(1, 6),
      volc.size = 0.4,
      volc.alpha = 1,
      volc.color = 'deeppink',
      box.limits = ..2,
      legend.pos = 'top',
      legend.just = 'right',
      legend.size = unit(8, 'lines'),
      plot.labs = ..3
    )
    ggsave(
      paste0('figs/maps/volcano/', ..4, '_volc.png'),
      device = 'png',
      dpi = 320,
      bg = 'transparent',
      width = 10,
      height = 7
    )
  }
)

# Volcano plots
# Tidy Volcano data
d.volcano <- volcano %>% full_join(volcano.no.spread) %>% group_by(ArcName) %>% mutate('ArcName' = gsub('_', ' ', ArcName))
order.volcano.H <- forcats::fct_reorder(d.volcano$ArcName, d.volcano$H, median)
d.volcano.long <- d.volcano %>% select(-c('Lat', 'Lon', 'Vtot', 'Vcorr', 'AZconv', 'Strike', 'Dip')) %>% tidyr::pivot_longer(cols = where(is.numeric))

# 1D
# Slab depth
ggplot(d.volcano) +
  geom_boxplot(aes(x = forcats::fct_reorder(ArcName, H, median), y = H, group = ArcName)) +
  labs(title = 'Slab Depth beneath Volcanoes', x = 'Arc Segment', y = 'H (km)', caption = 'Data: Syracuse & Abers (2006)') +
  scale_y_reverse() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor.y = element_blank()
  )
ggsave(
  paste0('figs/H_boxplot.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 10,
  height = 5.5
)

# Phi
ggplot(d.volcano) +
  geom_boxplot(aes(x = forcats::fct_reorder(ArcName, H, median), y = `Phi/100`, group = ArcName)) +
  labs(title = 'Slab Thermal Parameter', x = 'Arc Segment', y = bquote(over(Phi, 100)~~~bgroup('(', over(km, 100), ')')), caption = 'Data: Syracuse & Abers (2006)') +
  scale_y_reverse() +
  theme_bw() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor.y = element_blank()
  )
ggsave(
  paste0('figs/Phi_boxplot.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 10,
  height = 5.5
)

# Tdis
ggplot(d.volcano) +
  geom_boxplot(aes(x = forcats::fct_reorder(ArcName, H, median), y = Tdis, group = ArcName)) +
  labs(title = 'Arc-Trench Distance', x = 'Arc Segment', y = 'Distance from Trench (km)', caption = 'Data: Syracuse & Abers (2006)') +
  scale_y_reverse() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor.y = element_blank()
  )
ggsave(
  paste0('figs/Tdis_boxplot.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 10,
  height = 5.5
)

# All important parameters vs. arc segment in 1D
ggplot(d.volcano.long) +
  geom_boxplot(aes(x = factor(ArcName, levels = levels(order.volcano.H)), y = value, group = ArcName)) +
  labs(title = 'Volcano & Slab Parameters by Arc Segment', x = 'Arc Segment', y = NULL, caption = 'Data: Syracuse & Abers (2006)') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.minor.y = element_blank()
  ) +
  facet_wrap(~ factor(name, levels = c('H', 'Phi/100', 'Vc', 'Age', 'DesRat', 'Tdis')), scales = 'free_y', ncol = 2)
ggsave(
  paste0('figs/facet_boxplot.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 11,
  height = 11
)

# 2D
# Slab Depth vs. Trench Distance
ggplot(d.volcano) +
  geom_point(aes(x = Tdis, y = H)) +
  scale_y_reverse() +
  labs(title = 'Slab Depth beneath Volcanoes', x = 'Distance from Trench (km)', y = 'H (km)', caption = 'Data: Syracuse & Abers (2006)') +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent", color = NA))
ggsave(
  paste0('figs/H_Tdis.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 7,
  height = 5.5
)

# Slab Depth vs. Phi
ggplot(d.volcano) +
  geom_point(aes(x = `Phi/100`, y = H)) +
  scale_y_reverse() +
  labs(title = 'Slab Depth beneath Volcanoes', x = bquote(over(Phi, 100)~~~bgroup('(', over(km, 100), ')')), y = 'H (km)', caption = 'Data: Syracuse & Abers (2006)') +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent", color = NA))
ggsave(
  paste0('figs/H_Phi.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 7,
  height = 5.5
)

# Slab depth vs. trench distance faceted by segment
1:5 %>% purrr::map(~ {
  ggplot(d.volcano) +
    geom_point(aes(x = Tdis, y = H)) +
    scale_y_reverse() +
    labs(title = 'Slab Depth beneath Volcanoes', x = 'Distance from Trench (km)', y = 'H (km)', caption = 'Data: Syracuse & Abers (2006)') +
    theme_bw() +
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    ggforce::facet_wrap_paginate(
      ~ forcats::fct_reorder(ArcName, H, median),
      ncol = 3,
      nrow = 3,
      page = .x
    )
  ggsave(
    paste0('figs/H_Tdis_', .x, '.png'),
    device = 'png',
    dpi = 320,
    bg = 'transparent',
    width = 7,
    height = 7
  )
})

# Slab depth vs. Phi faceted by segment
1:5 %>% purrr::map(~ {
  ggplot(d.volcano) +
    geom_point(aes(x = `Phi/100`, y = H)) +
    scale_y_reverse() +
    labs(title = 'Slab Depth beneath Volcanoes', x = bquote(over(Phi, 100)~~~bgroup('(', over(km, 100), ')')), y = 'H (km)', caption = 'Data: Syracuse & Abers (2006)') +
    theme_bw() +
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    ggforce::facet_wrap_paginate(
      ~ forcats::fct_reorder(ArcName, H, median),
      ncol = 3,
      nrow = 3,
      page = .x
    )
  ggsave(
    paste0('figs/H_Phi_', .x, '.png'),
    device = 'png',
    dpi = 320,
    bg = 'transparent',
    width = 7,
    height = 7
  )
})

# Age, Vc, and H
ggplot(d.volcano) +
  geom_point(aes(x = Vc, y = Age, color = H), shape = 15) +
  labs(title = 'Slab Depth beneath Volcanoes', x = 'Convergence Velocity (km/Ma)', y = 'Age (Ma)', caption = 'Data: Syracuse & Abers (2006)') +
  scale_color_viridis_c() +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent", color = NA))
ggsave(
  paste0('figs/Age_Vc_H.png'),
  device = 'png',
  dpi = 320,
  bg = 'transparent',
  width = 7,
  height = 5.5
)

# Kriging
load('data/krig.RData')

# Composition layout
lyt <- '
1111122
1111122
3333333
3333333
3333333
3333333
3333333
'

# Draw and save compositions for segments with np > 30
purrr::pmap(list(purrr::map(k100.w[seg.30np$segment], ~ .x$v.plot), purrr::map(k100.w[seg.30np$segment], ~ .x$hist), purrr::map(k100.w[seg.30np$segment], ~ .x$k.plot), seg.30np$segment),
            ~ {..1 +
                (..2 + ylab(bquote(Heat~Flow~(mW/m^2))) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())) +
                ..3 +
                plot_layout(design = lyt) +
                plot_annotation(tag_levels = 'a', title = ..4) &
                theme(plot.background = element_rect(fill = "transparent", color = NA))
              ggsave(filename = paste0('figs/maps/kriged/k100w.', ..4, '.png'), device = 'png')
            })

# Draw and save compositions for segments with np > 30
purrr::pmap(list(purrr::map(k100.n[seg.30np$segment], ~ .x$v.plot), purrr::map(k100.n[seg.30np$segment], ~ .x$hist), purrr::map(k100.n[seg.30np$segment], ~ .x$k.plot), seg.30np$segment),
            ~ {..1 +
                (..2 + ylab(bquote(Heat~Flow~(mW/m^2))) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())) +
                ..3 +
                plot_layout(design = lyt) +
                plot_annotation(tag_levels = 'a', title = ..4) &
                theme(plot.background = element_rect(fill = "transparent", color = NA))
              ggsave(filename = paste0('figs/maps/kriged/k100n.', ..4, '.png'), device = 'png')
            })

# Close device
dev.off()