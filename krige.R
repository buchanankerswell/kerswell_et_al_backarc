rm(list = ls())
source('functions.R')
load('data/hf.RData')
rm.list <- ls()

# This script uses the `krg` function to fit calculate variograms, fit variogram models, krige,
# and interpolate for a sampled area near arc segments

# Each individual arc segment is kriged individually as the spatial correlation of heat flow data
# varies from segment to segment. Experimental variograms for some segments do not produce good results,
# so these segments are split into equidistant subsegments before using the `krg` function.

# Arc Segments
# - Alaska Aleutians
# - Andes
# - Central America
# - Kamchatka Marianas
# - Kyushu Ryukyu
# - Lesser Antilles
# - N. Philippines
# - New Britain Solomon
# - S. Philippines
# - Scotia
# - Sumatra Banda Sea
# - Tonga New Zealand
# - Vanuatu

# Filter heat flow data
shp.hf.filtered <- shp.hf %>%
  filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500, No__Temps < 250) %>%
  rename('Conductivity' = Conductivi, 'Heat Flow' = Heat_Flow, 'Min Depth' = minD, 'Max Depth' = maxD, 'No. of Cond.' = No__Cond_, 'No. of Temps' = No__Temps)
d.hf.filtered <- shp.hf.filtered %>% st_set_geometry(NULL)

# Summarise global hf data ----
d.hf.filtered %>%
  group_by(segment) %>%
  summarise(n = n())
# Summary plots
# Number of hf data points
p.no <- d.hf.filtered %>%
  group_by(segment) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(x = n, y = forcats::fct_rev(factor(segment)), fill = segment), stat = 'identity', alpha = 0.8, show.legend = F) +
  labs(x = NULL, y = NULL, title = 'Number of Heat Flow Data Points') +
  labs(y = 'Segment', x = NULL, title = bquote(Data~Points)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Print
p.no
# Heat flow
p.hf <- d.hf.filtered %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = `Heat Flow`, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Heat~Flow~~mWm^-2)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Print
p.hf
# Gradient
p.grad <- d.hf.filtered %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = Gradient, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Gradient~~mKm^-1)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Print
p.grad
# Conductivity
p.con <- d.hf.filtered %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = Conductivity, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Conductivity~~Wm^-1~K^-1)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Print
p.con
# Elevation
p.el <- d.hf.filtered %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = Elevation, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Elevation~~m)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Print
p.el
# No. of Temps
p.nt <- d.hf.filtered %>%
  group_by(segment) %>%
  ggplot() +
  geom_boxplot(aes(x = `No. of Temps`, y = forcats::fct_rev(factor(segment)), fill = segment), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = 'Segment', x = NULL, title = bquote(Number~of~Temps)) +
  scale_fill_viridis_d(option = 'D') +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(hjust = 1),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Print
p.nt
# Composite summary plot
(p.no + theme(plot.margin = margin(8, 0, 8, 8), axis.text.y = element_text(margin = margin(0, 4, 0, 0)))) +
  (p.hf + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) +
  (p.grad + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) +
  (p.con + theme(plot.margin = margin(8, 0, 8, 8), axis.text.y = element_text(margin = margin(0, 4, 0, 0)))) +
  (p.el + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) +
  (p.nt + theme(axis.text.y = element_blank(), plot.margin = margin(8, 0, 8, 0))) &
  theme(panel.spacing.x = unit(0, 'pt'), axis.line.y.left = element_line(color = rgb(0, 0, 0, 0.2)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/summary.png', device = 'png', dpi = 330, width = 11, height = 7, bg = 'transparent')
# Summary by variable
d.hf.filtered %>%
  select(c(`Heat Flow`, Gradient, Conductivity, Elevation, `No. of Temps`, segment)) %>%
  pivot_longer(-segment, 'variable') %>%
  mutate(variable = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))) %>%
  group_by(segment, variable) %>%
  summarise(
    min = min(value, na.rm = T),
    max = max(value, na.rm = T),
    median = median(value, na.rm = T),
    mean = mean(value, na.rm = T))
# Plot all (unfiltered) heat flow relative to trench
purrr::map(unique(shp.hf$segment), ~shp.hf %>% filter(segment == .x) %>%
             mutate(pos = pts_position(shp.hf %>% filter(segment == .x), shp.sa.segs.robin.pacific %>% filter(segment == .x), 5e6, 5e6, 'updown', 'up'),
                    t.dist = trench_distance(shp.hf %>% filter(segment == .x), pos, shp.sa.segs.robin.pacific %>% filter(segment == .x))) %>%
             ggplot() + geom_point(aes(x = t.dist/1000, y = Heat_Flow)) + labs(x = bquote(Distance~~km), y = bquote(Heat~Flow~~mWm^-2), title = .x) + theme_classic(base_size = 14))

# Alaska Aleutians ----
# Data
d.hf.filtered %>%
  filter(segment == 'Alaska Aleutians') %>%
  select(c(`Heat Flow`, Gradient, Conductivity, Elevation, `No. of Temps`))
# Summary
d.hf.filtered %>%
  filter(segment == 'Alaska Aleutians') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  mutate(variable = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))) %>%
  group_by(variable) %>%
  summarise(
            min = min(value, na.rm = T),
            max = max(value, na.rm = T),
            median = median(value, na.rm = T),
            mean = mean(value, na.rm = T))
# Summary Plot
d.hf.filtered %>%
  filter(segment == 'Alaska Aleutians') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  group_by(variable) %>%
  ggplot() +
  geom_boxplot(aes(x = value, group = variable, fill = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = NULL, x = NULL, title = paste0('Alaska Aleutians | n = ', d.hf.filtered %>% filter(segment == 'Alaska Aleutians') %>% count())) +
  scale_fill_viridis_d(option = 'D') +
  facet_wrap(~factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps')), ncol = 1, scales = 'free') +
  theme_void() +
  theme(strip.text = element_text(size = 14, margin = margin(0, 0, 8, 0)),
        axis.text.x = element_text(),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
        )
# Save Plot
ggsave('figs/aa-summary.png', device = 'png', dpi = 330, width = 7, height = 7, bg = 'transparent')
# Plot all (unfiltered) heat flow vs. distance to trench
shp.hf %>% filter(segment == 'Alaska Aleutians') %>%
  mutate(pos = pts_position(shp.hf %>% filter(segment == 'Alaska Aleutians'), shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians'), 5e6, 5e6, 'updown', 'up'),
         t.dist = trench_distance(shp.hf %>% filter(segment == 'Alaska Aleutians'), pos, shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians'))) %>%
  ggplot() + geom_point(aes(x = t.dist, y = Heat_Flow))
shp.hf %>% filter(segment == 'Alaska Aleutians') %>%
  mutate(pos = pts_position(shp.hf %>% filter(segment == 'Alaska Aleutians'), shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians'), 5e6, 5e6, 'updown', 'up'),
         t.dist = trench_distance(shp.hf %>% filter(segment == 'Alaska Aleutians'), pos, shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians'))) %>%
  ggplot() + geom_sf(data = shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians')) + geom_sf(aes(color = Heat_Flow, shape = pos))
# Krige entire segment
k.aa <- krg(
  data = shp.hf.filtered %>% filter(segment == 'Alaska Aleutians') %>% rename('Heat_Flow' = `Heat Flow`),
  lags = 100,
  lag.cutoff = 3,
  param = 'Heat_Flow',
  krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == 'Alaska Aleutians'),
  seg = shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians'),
  contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Alaska Aleutians'),
  crs = proj4.robin.pacific,
  ngrid = 3e5,
  grid.method = 'hexagonal',
  v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
  plot = T
)
# Composite summary plot
k.aa$k.plot / (k.aa$v.plot + k.aa$hist) + plot_layout(heights = c(3, 1)) +
  plot_annotation(title = 'Alaska Aleutians') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), panel.spacing.y = unit(0, 'pt'), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/aa.png', device = 'png', dpi = 330, width = 11, height = 7, bg = 'transparent')
# Summarise experimental variogram
k.aa$variogram %>%
  summarise(np.max = max(np), np.min = min(np), np.mean = mean(np), np.median = median(np))
# Variogram model
k.aa$v.model

# Subsegment kriging
# Segment length
l.aa <- st_length(shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians'))
# Split segment
shp.aa.splt.seg <- shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians') %>%
  splt(cut.length = as.vector(l.aa)/4, buffer = F, buffer.dist = 500000, cap.style = 'ROUND') %>% mutate(segment = 'Alaska Aleutians', .before = geometry)
shp.aa.splt.buffer <- shp.sa.segs.robin.pacific %>% filter(segment == 'Alaska Aleutians') %>%
  splt(cut.length = as.vector(l.aa)/4, buffer = T, buffer.dist = 500000, cap.style = 'ROUND') %>% mutate(segment = 'Alaska Aleutians', .before = geometry)
# Split hf data
shp.hf.aa.splt <- shp.hf.filtered %>% st_join(shp.aa.splt.buffer, left = F)
# Drop simple features geometry
hf.aa.splt <- shp.hf.aa.splt %>% st_set_geometry(NULL)
# Summarise
hf.aa.splt %>%
  group_by(subseg) %>%
  summarise(n = n())
# Krige individual segments
k.aa.subseg <- purrr::map(unique(shp.hf.aa.splt$subseg) %>% sort(), ~{
  krg(
    data = shp.hf.aa.splt %>% filter(subseg == .x) %>% rename('Heat_Flow' = `Heat Flow`),
    lags = 100,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = shp.aa.splt.buffer %>% filter(subseg == .x),
    seg = shp.aa.splt.seg %>% filter(subseg == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Alaska Aleutians'),
    crs = proj4.robin.pacific,
    ngrid = 3e5,
    grid.method = 'hexagonal',
    v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
    plot = T
  )
})
# Summarise experimental variograms
purrr::map_df(k.aa.subseg, ~.x$variogram %>% summarise(np.max = max(np), np.min = min(np), np.mean = round(mean(np)), np.median = median(np)), .id = 'subseg')
# Summarise variogram models
purrr::map_df(k.aa.subseg, ~.x$v.mod, .id = 'subseg')
# Composite summary plot (krige results)
k.aa$k.plot / ((k.aa.subseg[[1]]$k.plot + k.aa.subseg[[2]]$k.plot) / (k.aa.subseg[[3]]$k.plot + k.aa.subseg[[4]]$k.plot)) +
  plot_annotation(title = 'Alaska Aleutians') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/aa-subseg-krig.png', device = 'png', dpi = 330, width = 11, height = 10, bg = 'transparent')
# Composite summary plot (variograms)
k.aa$v.plot / ((k.aa.subseg[[1]]$v.plot + k.aa.subseg[[2]]$v.plot + theme(axis.title.y = element_blank())) / (k.aa.subseg[[3]]$v.plot + k.aa.subseg[[4]]$v.plot + theme(axis.title.y = element_blank()))) +
  plot_annotation(title = 'Alaska Aleutians') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/aa-subseg-variogram.png', device = 'png', dpi = 330, width = 11, height = 10, bg = 'transparent')

# Andes ----
# Data
d.hf.filtered %>%
  filter(segment == 'Andes') %>%
  select(c(`Heat Flow`, Gradient, Conductivity, Elevation, `No. of Temps`))
# Summary
d.hf.filtered %>%
  filter(segment == 'Andes') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  mutate(variable = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))) %>%
  group_by(variable) %>%
  summarise(
    min = min(value, na.rm = T),
    max = max(value, na.rm = T),
    median = median(value, na.rm = T),
    mean = mean(value, na.rm = T))
# Summary Plot
d.hf.filtered %>%
  filter(segment == 'Andes') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  group_by(variable) %>%
  ggplot() +
  geom_boxplot(aes(x = value, group = variable, fill = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = NULL, x = NULL, title = paste0('Andes | n = ', d.hf.filtered %>% filter(segment == 'Andes') %>% count())) +
  scale_fill_viridis_d(option = 'D') +
  facet_wrap(~factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps')), ncol = 1, scales = 'free') +
  theme_void() +
  theme(strip.text = element_text(size = 14, margin = margin(0, 0, 8, 0)),
        axis.text.x = element_text(),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Save Plot
ggsave('figs/an-summary.png', device = 'png', dpi = 330, width = 7, height = 7, bg = 'transparent')
# Krige entire segment
k.an <- krg(
  data = shp.hf.filtered %>% filter(segment == 'Andes') %>% rename('Heat_Flow' = `Heat Flow`),
  lags = 100,
  lag.cutoff = 3,
  param = 'Heat_Flow',
  krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == 'Andes'),
  seg = shp.sa.segs.robin.pacific %>% filter(segment == 'Andes'),
  contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Andes'),
  crs = proj4.robin.pacific,
  ngrid = 3e5,
  grid.method = 'hexagonal',
  v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
  plot = T
)
# Composite summary plot
k.an$k.plot + (k.an$v.plot / k.an$hist) +
  plot_annotation(title = 'Andes') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), panel.spacing.y = unit(0, 'pt'), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/an.png', device = 'png', dpi = 330, width = 10, height = 7, bg = 'transparent')
# Summarise experimental variogram
k.an$variogram %>%
  summarise(np.max = max(np), np.min = min(np), np.mean = mean(np), np.median = median(np))
# Variogram model
k.an$v.model

# Subsegment kriging
# Segment length
l.an <- st_length(shp.sa.segs.robin.pacific %>% filter(segment == 'Andes'))
# Split segment
shp.an.splt.seg <- shp.sa.segs.robin.pacific %>% filter(segment == 'Andes') %>%
  splt(cut.length = as.vector(l.an)/9, buffer = F, buffer.dist = 500000, cap.style = 'ROUND') %>% mutate(segment = 'Andes', .before = geometry)
shp.an.splt.buffer <- shp.sa.segs.robin.pacific %>% filter(segment == 'Andes') %>%
  splt(cut.length = as.vector(l.an)/9, buffer = T, buffer.dist = 500000, cap.style = 'ROUND') %>% mutate(segment = 'Andes', .before = geometry)
# Split hf data
shp.hf.an.splt <- shp.hf.filtered %>% st_join(shp.an.splt.buffer, left = F)
# Drop simple features geometry
hf.an.splt <- shp.hf.an.splt %>% st_set_geometry(NULL)
# Summarise
hf.an.splt %>%
  group_by(subseg) %>%
  summarise(n = n())
# Krige individual segments
k.an.subseg <- purrr::map(unique(shp.hf.an.splt$subseg) %>% sort(), ~{
  krg(
    data = shp.hf.an.splt %>% filter(subseg == .x) %>% rename('Heat_Flow' = `Heat Flow`),
    lags = 100,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = shp.an.splt.buffer %>% filter(subseg == .x),
    seg = shp.an.splt.seg %>% filter(subseg == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Andes'),
    crs = proj4.robin.pacific,
    ngrid = 3e5,
    grid.method = 'hexagonal',
    v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
    plot = T
  )
})
# Summarise experimental variograms
purrr::map_df(k.an.subseg, ~.x$variogram %>% summarise(np.max = max(np), np.min = min(np), np.mean = round(mean(np)), np.median = median(np)), .id = 'subseg')
# Summarise variogram models
purrr::map_df(k.an.subseg, ~.x$v.mod, .id = 'subseg')
# Composite summary plot (krige results)
k.an$k.plot + ((k.an.subseg[[1]]$k.plot + k.an.subseg[[2]]$k.plot + k.an.subseg[[3]]$k.plot) / (k.an.subseg[[4]]$k.plot + k.an.subseg[[5]]$k.plot + k.an.subseg[[6]]$k.plot) / (k.an.subseg[[7]]$k.plot + k.an.subseg[[8]]$k.plot + k.an.subseg[[9]]$k.plot)) +
  plot_annotation(title = 'Andes') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/an-subseg-krig.png', device = 'png', dpi = 330, width = 13, height = 10, bg = 'transparent')
# Composite summary plot (variograms)
k.an$v.plot / ((k.an.subseg[[1]]$v.plot + k.an.subseg[[2]]$v.plot + theme(axis.title.y = element_blank()) + k.an.subseg[[3]]$v.plot + theme(axis.title.y = element_blank())) / (k.an.subseg[[4]]$v.plot + k.an.subseg[[5]]$v.plot + theme(axis.title.y = element_blank()) + k.an.subseg[[6]]$v.plot + theme(axis.title.y = element_blank())) / (k.an.subseg[[7]]$v.plot + k.an.subseg[[8]]$v.plot + theme(axis.title.y = element_blank()) + k.an.subseg[[9]]$v.plot + theme(axis.title.y = element_blank()))) +
  plot_annotation(title = 'Andes') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/an-subseg-variogram.png', device = 'png', dpi = 330, width = 13, height = 10, bg = 'transparent')

# Central America ----
# Data
d.hf.filtered %>%
  filter(segment == 'Central America') %>%
  select(c(`Heat Flow`, Gradient, Conductivity, Elevation, `No. of Temps`))
# Summary
d.hf.filtered %>%
  filter(segment == 'Central America') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  mutate(variable = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))) %>%
  group_by(variable) %>%
  summarise(
    min = min(value, na.rm = T),
    max = max(value, na.rm = T),
    median = median(value, na.rm = T),
    mean = mean(value, na.rm = T))
# Summary Plot
d.hf.filtered %>%
  filter(segment == 'Central America') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  group_by(variable) %>%
  ggplot() +
  geom_boxplot(aes(x = value, group = variable, fill = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = NULL, x = NULL, title = paste0('Central America | n = ', d.hf.filtered %>% filter(segment == 'Central America') %>% count())) +
  scale_fill_viridis_d(option = 'D') +
  facet_wrap(~factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps')), ncol = 1, scales = 'free') +
  theme_void() +
  theme(strip.text = element_text(size = 14, margin = margin(0, 0, 8, 0)),
        axis.text.x = element_text(),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Save Plot
ggsave('figs/ca-summary.png', device = 'png', dpi = 330, width = 7, height = 7, bg = 'transparent')
# Krige entire segment
k.ca <- krg(
  data = shp.hf.filtered %>% filter(segment == 'Central America') %>% rename('Heat_Flow' = `Heat Flow`),
  lags = 100,
  lag.cutoff = 3,
  param = 'Heat_Flow',
  krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == 'Central America'),
  seg = shp.sa.segs.robin.pacific %>% filter(segment == 'Central America'),
  contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Central America'),
  crs = proj4.robin.pacific,
  ngrid = 3e5,
  grid.method = 'hexagonal',
  v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
  plot = T
)
# Composite summary plot
k.ca$k.plot / (k.ca$v.plot + k.ca$hist) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(title = 'Central America') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), panel.spacing.y = unit(0, 'pt'), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/ca.png', device = 'png', dpi = 330, width = 10, height = 10, bg = 'transparent')
# Summarise experimental variogram
k.ca$variogram %>%
  summarise(np.max = max(np), np.min = min(np), np.mean = mean(np), np.median = median(np))
# Variogram model
k.ca$v.model

# Kamchatka Marianas ----
# Data
d.hf.filtered %>%
  filter(segment == 'Kamchatka Marianas') %>%
  select(c(`Heat Flow`, Gradient, Conductivity, Elevation, `No. of Temps`))
# Summary
d.hf.filtered %>%
  filter(segment == 'Kamchatka Marianas') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  mutate(variable = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))) %>%
  group_by(variable) %>%
  summarise(
    min = min(value, na.rm = T),
    max = max(value, na.rm = T),
    median = median(value, na.rm = T),
    mean = mean(value, na.rm = T))
# Summary Plot
d.hf.filtered %>%
  filter(segment == 'Kamchatka Marianas') %>%
  select(c(`Heat Flow`, Conductivity, Elevation, Gradient, `No. of Temps`)) %>%
  pivot_longer(everything(), 'variable') %>%
  group_by(variable) %>%
  ggplot() +
  geom_boxplot(aes(x = value, group = variable, fill = factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps'))), color = 'black', alpha = 0.8, outlier.alpha = 1, show.legend = F) +
  labs(y = NULL, x = NULL, title = paste0('Kamchatka Marianas | n = ', d.hf.filtered %>% filter(segment == 'Kamchatka Marianas') %>% count())) +
  scale_fill_viridis_d(option = 'D') +
  facet_wrap(~factor(variable, levels = c('Heat Flow', 'Gradient', 'Conductivity', 'Elevation', 'No. of Temps')), ncol = 1, scales = 'free') +
  theme_void() +
  theme(strip.text = element_text(size = 14, margin = margin(0, 0, 8, 0)),
        axis.text.x = element_text(),
        plot.margin = margin(12, 14, 12, 14),
        panel.spacing.y = unit(14, 'pt'),
        plot.title = element_text(size = 18, margin = margin(0, 0, 14, 0), hjust = 0.5),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.background = element_rect(fill = 'transparent', color = NA)
  )
# Save Plot
ggsave('figs/km-summary.png', device = 'png', dpi = 330, width = 7, height = 7, bg = 'transparent')
# Krige entire segment
k.km <- krg(
  data = shp.hf.filtered %>% filter(segment == 'Kamchatka Marianas') %>% rename('Heat_Flow' = `Heat Flow`),
  lags = 100,
  lag.cutoff = 3,
  param = 'Heat_Flow',
  krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == 'Kamchatka Marianas'),
  seg = shp.sa.segs.robin.pacific %>% filter(segment == 'Kamchatka Marianas'),
  contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Kamchatka Marianas'),
  crs = proj4.robin.pacific,
  ngrid = 3e5,
  grid.method = 'hexagonal',
  v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
  plot = T
)
# Composite summary plot
k.km$k.plot + (k.km$v.plot / k.km$hist) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(title = 'Kamchatka Marianas') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), panel.spacing.y = unit(0, 'pt'), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA), legend.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/km.png', device = 'png', dpi = 330, width = 11, height = 7, bg = 'transparent')
# Summarise experimental variogram
k.km$variogram %>%
  summarise(np.max = max(np), np.min = min(np), np.mean = mean(np), np.median = median(np))
# Variogram model
k.km$v.model

# Subsegment kriging
# Segment length
l.km <- st_length(shp.sa.segs.robin.pacific %>% filter(segment == 'Kamchatka Marianas'))
# Split segment
shp.km.splt.seg <- shp.sa.segs.robin.pacific %>% filter(segment == 'Kamchatka Marianas') %>%
  splt(cut.length = as.vector(l.km)/3, buffer = F, buffer.dist = 500000, cap.style = 'ROUND') %>% mutate(segment = 'Kamchatka Marianas', .before = geometry)
shp.km.splt.buffer <- shp.sa.segs.robin.pacific %>% filter(segment == 'Kamchatka Marianas') %>%
  splt(cut.length = as.vector(l.km)/3, buffer = T, buffer.dist = 500000, cap.style = 'ROUND') %>% mutate(segment = 'Kamchatka Marianas', .before = geometry)
# Split hf data
shp.hf.km.splt <- shp.hf.filtered %>% st_join(shp.km.splt.buffer, left = F)
# Drop simple features geometry
hf.km.splt <- shp.hf.km.splt %>% st_set_geometry(NULL)
# Summarise
hf.km.splt %>%
  group_by(subseg) %>%
  summarise(n = n())
# Krige individual segments
k.km.subseg <- purrr::map(unique(shp.hf.km.splt$subseg) %>% sort(), ~{
  krg(
    data = shp.hf.km.splt %>% filter(subseg == .x) %>% rename('Heat_Flow' = `Heat Flow`),
    lags = 100,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = shp.km.splt.buffer %>% filter(subseg == .x),
    seg = shp.km.splt.seg %>% filter(subseg == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == 'Kamchatka Marianas'),
    crs = proj4.robin.pacific,
    ngrid = 3e5,
    grid.method = 'hexagonal',
    v.mod <- c('Sph', 'Nug', 'Gau', 'Mat', 'Lin', 'Cir', 'Per', 'Wav', 'Log'),
    plot = T
  )
})
# Summarise experimental variograms
purrr::map_df(k.km.subseg, ~.x$variogram %>% summarise(np.max = max(np), np.min = min(np), np.mean = round(mean(np)), np.median = median(np)), .id = 'subseg')
# Summarise variogram models
purrr::map_df(k.km.subseg, ~.x$v.mod, .id = 'subseg')
# Composite summary plot (krige results)
k.km$k.plot + ((k.km.subseg[[1]]$k.plot + k.km.subseg[[2]]$k.plot + k.km.subseg[[3]]$k.plot) / (k.km.subseg[[4]]$k.plot + k.km.subseg[[5]]$k.plot + k.km.subseg[[6]]$k.plot) / (k.km.subseg[[7]]$k.plot + k.km.subseg[[8]]$k.plot + k.km.subseg[[9]]$k.plot)) +
  plot_annotation(title = 'Kamchatka Marianas') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/km-subseg-krig.png', device = 'png', dpi = 330, width = 13, height = 10, bg = 'transparent')
# Composite summary plot (variograms)
k.km$v.plot / ((k.km.subseg[[1]]$v.plot + k.km.subseg[[2]]$v.plot + theme(axis.title.y = element_blank()) + k.km.subseg[[3]]$v.plot + theme(axis.title.y = element_blank())) / (k.km.subseg[[4]]$v.plot + k.km.subseg[[5]]$v.plot + theme(axis.title.y = element_blank()) + k.km.subseg[[6]]$v.plot + theme(axis.title.y = element_blank())) / (k.km.subseg[[7]]$v.plot + k.km.subseg[[8]]$v.plot + theme(axis.title.y = element_blank()) + k.km.subseg[[9]]$v.plot + theme(axis.title.y = element_blank()))) +
  plot_annotation(title = 'Kamchatka Marianas') &
  theme(plot.title = element_text(size = 18, hjust = 0.5, margin = margin(14, 0, 0, 0)), plot.background = element_rect(fill = 'transparent', color = NA), panel.background = element_rect(fill = 'transparent', color = NA))
# Save plot
ggsave('figs/km-subseg-variogram.png', device = 'png', dpi = 330, width = 10, height = 10, bg = 'transparent')

