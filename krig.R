rm(list = ls())
source('functions.R')
load('data/hf.RData')

# Kriging
# Filter data
dat <- shp.hf %>% filter(Heat_Flow > 10 & Heat_Flow < 120 & minD < maxD & Gradient > 0 & Gradient < 500)

# Calculate experimental variogram, fit experimental variogram, krig, and interpolate
k100.w <- purrr::map(
  unique(shp.hf$segment),
  ~ krg(
    data = dat %>% filter(segment == .x),
    lags = 100,
    lag.cutoff = 3,
    param = 'Heat_Flow',
    krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == .x),
    seg = shp.sa.segs.robin.pacific %>% filter(segment == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == .x),
    crs = proj4.robin.pacific,
    ngrid = 3.2e5,
    grid.method = 'hexagonal',
    v.mod = 'Sph',
    plot = T
  )
)
# Save
save(k100.w, file = 'data/k100w.RData')

# Calculate experimental variogram, fit experimental variogram, krig, and interpolate
k100.n <- purrr::map(
  unique(shp.hf$segment),
  ~ krg(
    data = dat %>% filter(segment == .x),
    lags = 100,
    lag.cutoff = 6,
    param = 'Heat_Flow',
    krg.shp = shp.sa.segs.robin.pacific.buffer %>% filter(segment == .x),
    seg = shp.sa.segs.robin.pacific %>% filter(segment == .x),
    contours = shp.sa.countours.robin.pacific %>% filter(segment == .x),
    crs = proj4.robin.pacific,
    ngrid = 3.2e5,
    grid.method = 'hexagonal',
    v.mod = 'Sph',
    plot = T
  )
)
# Save
save(k100.n, file = 'data/k100n.RData')

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

# Draw and save compositions
purrr::pmap(list(
  purrr::map(k100.w, ~ .x$v.plot),
  purrr::map(k100.w, ~ .x$hist),
  purrr::map(k100.w, ~ .x$k.plot),
  1:13
),
~ {
  ..1 + (
    ..2 + ylab(bquote(Heat~Flow~(mW/m^2))) +
      theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  ) +
    ..3 + plot_layout(design = lyt) + plot_annotation(tag_levels = 'a') + theme(plot.background = element_rect(fill = "transparent", color = NA))
  ggsave(filename = paste0('figs/maps/kriged/k100.w', ..4, '.png'), device = 'png')
})

# Draw and save compositions
purrr::pmap(list(
  purrr::map(k100.n, ~ .x$v.plot),
  purrr::map(k100.n, ~ .x$hist),
  purrr::map(k100.n, ~ .x$k.plot),
  1:13
),
~ {
  ..1 + (
    ..2 + ylab(bquote(Heat~Flow~(mW/m^2))) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  ) +
    ..3 + plot_layout(design = lyt) + plot_annotation(tag_levels = 'a') + theme(plot.background = element_rect(fill = "transparent", color = NA))
  ggsave(filename = paste0('figs/maps/kriged/k100.n', ..4, '.png'), device = 'png')
})