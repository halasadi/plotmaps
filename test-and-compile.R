setwd("/Users/halasadi/plotmaps/")
devtools::document()
devtools::load_all()

source("~/plotmaps/R/plotmaps.R")

path = "~/eems2/test-plotting/test-data/"
mcmcpath=path
longlat=TRUE
is.mrates=TRUE
dimns <- read_dimns(path, longlat)

plot_contour(mcmcpath, dimns, is.mrates)
