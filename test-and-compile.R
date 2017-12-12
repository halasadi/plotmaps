setwd("/Users/halasadi/plotmaps/")
devtools::document()
devtools::load_all()

source("~/plotmaps/R/plotmaps.R")

path = "~/eems2/data/overall/r2/popressard_2_Inf/output/"
mcmcpath=path
longlat=FALSE
is.mrates=TRUE
dimns <- read_dimns(path, longlat)
params <- list(is.mrates=is.mrates, add.pts = TRUE, add.graph = TRUE, add.countries=TRUE)
g <- make_map(dimns, params)
g <- add_contours(mcmcpath, dimns, g, params)
