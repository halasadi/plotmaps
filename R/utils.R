#' @import sp

#' @title default_eems_colors
#' @return a vector of the default eems colors
#' @export
default_eems_colors <- function( ) {
  ## To reproduce the default eems colors:
  ## Oranges <- dichromat::colorschemes$BluetoDarkOrange.12[12:7]
  ## Blues <- dichromat::colorschemes$BrowntoBlue.12[7:12]
  ## eems.colors <- c(Oranges, "#FBFBFB", Blues)
  eems.colors <- c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", ## orange sequence
                   "#FBFBFB", ## very slightly off-white
                   "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99") ## blue sequence
  return (eems.colors)
}


#' reads the underlying graph
#' @param path path to mcmc output files
#' @param longlat a binary variable, TRUE if long appears before lat
#' @export
read_graph <- function(path, longlat) {
  eems.output <- NULL
  if (file.exists(paste0(path, '/demes.txt')) &&
      file.exists(paste0(path, '/ipmap.txt')) &&
      file.exists(paste0(path, '/outer.txt'))) {
    eems.output <- TRUE
  } 
  if (is.null(eems.output)) {
    stop(paste0(path, ' is not a mcmcpath.'))
  }
  ## Read the assigned sample coordinates
  ipmap <- scan(paste0(path, '/ipmap.txt'), what = numeric(), quiet = TRUE)
  demes <- scan(paste0(path, '/demes.txt'), what = numeric(), quiet = TRUE)
  outer <- scan(paste0(path, '/outer.txt'), what = numeric(), quiet = TRUE)
  demes <- matrix(demes, ncol = 2, byrow = TRUE)
  outer <- matrix(outer, ncol = 2, byrow = TRUE)
  edges <- read.table(paste0(path, '/edges.txt'), colClasses = numeric())
  edges <- as.matrix(edges)
  if (!longlat) {
    demes <- demes[, c(2, 1)]
    outer <- outer[, c(2, 1)]
  }
  sizes <- table(ipmap)
  alpha <- as.numeric(names(sizes))
  sizes <- as.numeric(sizes)
  return(list(ipmap = ipmap, demes = demes, edges = edges, alpha = alpha, sizes = sizes, outer = outer))
}

#'  checks if files exist
#' @param files a list of files to check if they exist
#' @export
check_files_at_path <- function(files, paths=".", strict=F){
  # checks if ech subfolder `paths` contains all the files in `files`.
  file_paths <- c(outer(paths, files, paste, sep=.Platform$file.sep))
  file_exist <- file.exists(file_paths)
  if( all(file_exist) ) return( paths )
  
  print(paste0(file_paths[!file_exist], " not found", collate="\n") )
  
  if( strict ){
    stop( )
  }
  
  f <- matrix(file_exist, nrow=length(paths))
  return(paths[ rowSums(f) == ncol(f)])
}


geo_distm <- function(coord, longlat) {
  if (!longlat) {
    long <- coord[, 2]; lat <- coord[, 1]
  } else {
    long <- coord[, 1]; lat <- coord[, 2]
  }
  x <- cbind(long, lat)
  Dist <- sp::spDists(x, x, longlat = TRUE)
  Dist <- Dist[upper.tri(Dist, diag = FALSE)]
  return(Dist)
}

#' plots the diagnostic files to assess whether MAPS fits the data "well"
#' @param params list of parameter options
#' @export
plot_fit_data <- function(mcmcpath, outpath, longlat) {
  oDemes <- scan(paste0(mcmcpath, '/rdistoDemes.txt'), quiet = TRUE)
  oDemes <- matrix(oDemes, ncol = 3, byrow = TRUE)
  Sizes <- oDemes[, 3]
  nPops <- nrow(oDemes)
  Demes <- seq(nPops)
  Sobs <- as.matrix(read.table(paste0(mcmcpath, '/rdistJtDobsJ.txt'), header = FALSE))
  Shat <- as.matrix(read.table(paste0(mcmcpath, '/rdistJtDhatJ.txt'), header = FALSE))
  colnames(Sobs) <- Demes
  rownames(Sobs) <- Demes
  colnames(Shat) <- Demes
  rownames(Shat) <- Demes
  Dist = geo_distm(oDemes[, 1:2], longlat)
  df.between <- data.frame(Dist=Dist, Sobs = Sobs[upper.tri(Sobs)], Shat = Shat[upper.tri(Shat)],
        row.names = NULL)
        
  df.within <- data.frame(Sobs.within = diag(Sobs), Shat.within = diag(Shat), row.names=NULL)
  
  plot_pw(df.between)
  ggsave(filename = paste0(outpath, "/observed_vs_fitted-between.pdf"), width = 4, height = 4)
  plot_variogram(df.between)
  ggsave(paste0(outpath, "/semivariogam.pdf"), width = 4, height = 4)
  plot_ww(df.within)
  ggsave(paste0(outpath, "/observed_vs_fitted-within.pdf"), width = 4, height = 4)

}


plot_pw <- function(df){
  P <- ggplot(df) + geom_point(aes(y=Sobs, x=Shat), alpha=0.6)  
  P <- P + theme_classic() + 
    geom_abline(intercept=0) +
    theme(legend.position=0) +
    ylab("genetic (observed) similarity between demes") +
    xlab("fitted similarity between demes")  
  P
}

plot_ww <- function(df){
  P <- ggplot(df) + geom_point(aes(y=Sobs.within, x=Shat.within), alpha=0.6)  
  P <- P + theme_classic() + 
    geom_abline(intercept=0) +
    theme(legend.position=0) +
    ylab("genetic (observed) similarity within demes") +
    xlab("fitted similarity within demes")  
  P
}

plot_variogram <- function(df){
  P <- ggplot(df) + geom_point(aes(y=Sobs, x=Dist), alpha=0.6)
  P + theme_classic() + geom_abline(intercept=0) +
    theme(legend.position=0) + ylab("genetic similarity") + 
    xlab("geographic distance")
  P
}

#' plots the log-likelihood and prior as a function of MCMC iteration (thinned)
#' @param params list of parameter options
#' @export
plot_trace <- function(mcmcpath, outpath){
  pilogl <- scan(paste0(mcmcpath, '/mcmcpilogl.txt'), quiet = TRUE)
  pilogl <- data.frame(matrix(pilogl, ncol = 2, byrow = TRUE))
  colnames(pilogl) <- c("prior", "logll")
  
  ggplot(pilogl, aes(y = logll, x = seq(1, length(logll)))) + geom_line() + xlab("iteration (thinned)")
  ggsave(paste0(outpath, "/logll.pdf"), width = 5, height = 3)

  ggplot(pilogl, aes(y = prior, x = seq(1, length(prior)))) + geom_line() + xlab("iteration (thinned)")
  ggsave(paste0(outpath, "/prior.pdf"), width = 5, height = 3)
  
}

#' computes scaling factors so results are interpretable in physical quanitites.
#' @param mcmcpath path to mcmc files
#' @export
compute_scaling <- function(mcmcpath){
  demes = read.table(paste0(mcmcpath, "/demes.txt"))
  odemes = nrow(read.table(paste0(mcmcpath, "/rdistoDemes.txt")))
  ndemes = nrow(demes)
  
  dx = min(dist(demes)) * 110.57
  m.scalingfactor = (dx)^2
  
  total.area = Polygon(read.table(paste0(mcmcpath, "/outer.txt")))@area * (110.57)^2
  N.scalingfactor = total.area / ndemes
  
  return(list(m.scalingfactor = m.scalingfactor, N.scalingfactor = N.scalingfactor))
}
