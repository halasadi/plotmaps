library(devtools)

devtools::document("/Users/halasadi/plotmaps/")
devtools::load_all("/Users/halasadi/plotmaps/")

path = "~/eems2/test-plotting/test-data/"
mcmcpath=path
longlat=TRUE
is.mrates=TRUE
dimns <- read_dimns(path, longlat)
params <- list(is.mrates=is.mrates, add.pts = TRUE, add.graph = TRUE, add.countries=FALSE,
               plot.median=FALSE, longlat=longlat, plot.mean = FALSE, plot.sign = TRUE)
g <- make_map(dimns, params)
g <- add_contours(mcmcpath, dimns, g, params)
ggsave("median-simulation.pdf", width = 10, height = 6)

