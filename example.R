# install package
#library(devtools)
#install("~/plotmaps")

library(plotmaps)
path = "~/MAPS/test-plotting/test-run"

plot_all(add.pts = TRUE, add.graph = TRUE, add.countries = FALSE,
         plot.median = TRUE, longlat = TRUE, mcmcpath = path, outpath = "garbage/",
         width = 10, height = 6)
