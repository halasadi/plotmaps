library(devtools)

devtools::document("/Users/halasadi/plotmaps/")
devtools::load_all("/Users/halasadi/plotmaps/")

path = "~/eems2/test-plotting/test-data/"
params <- list(is.mrates=FALSE, add.pts = TRUE, add.graph = TRUE, add.countries=FALSE,
               plot.median=FALSE, longlat=TRUE, plot.mean = FALSE, plot.sign = TRUE, 
               mcmcpath = path, outpath = "garbage/", width = 10, height = 6)

plotmaps(params)
plot_trace(params) 
plot_fit_data(params)
