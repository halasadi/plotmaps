library(devtools)

devtools::document("/Users/halasadi/plotmaps/")
devtools::load_all("/Users/halasadi/plotmaps/")

path = "~/eems2/test-plotting/test-data/"

plot_all(add.pts = TRUE, add.graph = TRUE, add.countries = FALSE,
         plot.median = TRUE, longlat = TRUE, mcmcpath = path, outpath = "garbage/",
         width = 10, height = 6)
