# Load functions and libraries
source('../functions.R')
load('../data/hf.RData')

# Load genetic algorithm results and decode chromosome into a variogram model
flist <- list.files(pattern = 'opt.RData', full.name = T)
fnames <- list.files(pattern = 'opt.RData') %>%
	purrr::map(~.x %>% stringr::str_replace('_opt.RData', '_vmod'))

# Save environment
init.list <- ls()

# Load ga results
for(i in flist) load(i)

# Decode chromosomes
v.mods <- purrr::map(ls()[!ls() %in% c(init.list, 'i', 'init.list')],
					 ~{
	# Get ga result
	o <- get(.x)@solution
  # Variogram model discritization formula
  if(o[1] >= 0 && o[1] < 1) {
    m <- 'Sph'
  } else if(o[1] >= 1 && o[1] < 2) {
    m <- 'Exp'
  } else if(o[1] >= 2 && o[1] <= 3) {
    m <- 'Gau'
  }
	# Construct variogram model
	vgm(psill = o[2], model = m, range = o[3], nugget = o[4])
}) %>% purrr::set_names(fnames)

# Krige segments and compare with Lucazeau (2019) predicted heat flow
s.names <- purrr::map(names(v.mods), ~.x %>%
											stringr::str_replace('_vmod', '') %>%
											stringr::str_replace('_', ' '))

purrr::walk2(v.mods, s.names, ~Krige_diff(.x, .y))