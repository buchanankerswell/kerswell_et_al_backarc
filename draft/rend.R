args <- commandArgs(trailingOnly = TRUE)
library(rmarkdown)
render(args[1])
system2('open', args = paste0(substr(args[1], 1, nchar(args[1])-4), '.pdf'), wait = FALSE)