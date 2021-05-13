# Load functions and libraries
source('functions.R')

# Load decoded genes
cat('\nLoading decoded chromosomes from data/genes_decoded.RData')
load('data/genes_decoded.RData')

v.grms.ga <- v.grms
v.fits.ga <- v.fits

# Define paths and names
list.files('data/diff', pattern = '.RData', full.names = T) -> fpath
purrr::map_chr(
  list.files('data/diff', pattern = '.RData'),
  ~.x %>%
  stringr::str_replace('.RData', '')) -> fname

# Load data
cat('\nLoading variograms defined in data/diff')
for (i in fpath) load(i)

purrr::map(fname, ~get(.x)$v.grm) %>%
purrr::set_names(fname) -> v.grms

purrr::map(fname, ~get(.x)$v.mod) %>%
purrr::set_names(fname) -> v.fits

purrr::pmap(list(v.grms[names(v.grms) %in% names(v.grms.ga)],
                 v.fits[names(v.fits) %in% names(v.fits.ga)],
                 v.grms.ga, v.fits.ga,
                 names(v.grms[names(v.grms) %in% names(v.grms.ga)])), ~{
  cat('\nPlotting', ..5, 'variograms')
  v.grm <- ..1
  v.fit <- ..2
  v.grm.ga <- ..3
  v.fit.ga <- ..4
  seg.name <- ..5 %>% stringr::str_replace_all('_', ' ')
  ggplot() +
  geom_line(data = variogramLine(v.fit, maxdist = max(v.grm$dist)),
    aes(x = dist/1000, y = gamma)) +
  geom_line(data = variogramLine(v.fit.ga, maxdist = max(v.grm.ga$dist)),
    aes(x = dist/1000, y = gamma), color = 'deeppink') +
  geom_point(data = v.grm, aes(x = dist/1000, y = gamma), size = 0.6) +
  geom_point(data = v.grm.ga, aes(x = dist/1000, y = gamma), size = 0.6, color = 'deeppink') +
  labs(x = NULL, y = 'Semivariance', title = seg.name) +
  theme_classic() +
  theme(
    axis.text = element_text(color = 'black'),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank())}) %>%
wrap_plots() -> p

cat('\nSaving variogram plots to figs/vgrms_compare.png')
ggsave(
  'figs/vgrms_compare.png',
  plot = p,
  device = 'png',
  type = 'cairo',
  width = 7,
  height = 7
)

cat('\nDone')