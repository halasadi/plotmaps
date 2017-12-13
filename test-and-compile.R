setwd("/Users/halasadi/plotmaps/")
devtools::document()
devtools::load_all()

#source("~/plotmaps/R/plotmaps.R")

path = "~/eems2/test-plotting/test-data/"
mcmcpath=path
longlat=TRUE
is.mrates=TRUE
dimns <- read_dimns(path, longlat)
params <- list(is.mrates=is.mrates, add.pts = TRUE, add.graph = TRUE, add.countries=FALSE,
               plot.median=TRUE, longlat=longlat)
g <- make_map(dimns, params)
g <- add_contours(mcmcpath, dimns, g, params)
ggsave("median-simulation.pdf", width = 10, height = 6)

