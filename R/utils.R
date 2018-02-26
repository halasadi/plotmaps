#' @import sp
#' @importFrom gridExtra arrangeGrob

#' @title default_eems_colors
#' @return a vector of the default eems colors
#' @export
default_eems_colors <- function( ) {
  eems.colors <- c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", ## orange sequence
                   "#FBFBFB", ## very slightly off-white
                   "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99") ## blue sequence
  return (eems.colors)
}

#' @title inferno_colors
#' @return a vector of colors
#' @export
inferno_colors <- function(){
  inferno.colors <- c("#000004FF", "#140B35FF", "#3A0963FF", "#60136EFF", "#85216BFF", "#A92E5EFF", 
                      "#FBFBFB", "#FBFBFB", "#FBFBFB", ## very slightly off-white
                      "#CB4149FF", "#E65D2FFF", "#F78311FF", "#FCAD12FF", "#F5DB4BFF", "#FCFFA4FF")
  return(inferno.colors)
  
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

#' checks if files exist
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
  oDemes <- scan(paste0(mcmcpath[1], '/rdistoDemes.txt'), quiet = TRUE)
  oDemes <- matrix(oDemes, ncol = 3, byrow = TRUE)
  sizes <- as.matrix(oDemes[, 3])
  nPops <- nrow(oDemes)
  Demes <- seq(nPops)
  Sobs <- as.matrix(read.table(paste0(mcmcpath[1], '/rdistJtDobsJ.txt'), header = FALSE))
  
  Shat = matrix(nrow = nrow(Sobs), ncol = ncol(Sobs), 0)
  for (path in mcmcpath){
    Shat <- Shat + as.matrix(read.table(paste0(path, '/rdistJtDhatJ.txt'), header = FALSE))
  }
  Shat <- Shat / length(mcmcpath)
  
  colnames(Sobs) <- Demes
  rownames(Sobs) <- Demes
  colnames(Shat) <- Demes
  rownames(Shat) <- Demes
  Dist = geo_distm(oDemes[, 1:2], longlat)
  Sizes <- sizes %*% t(sizes)
  diag(Sizes) <- (sizes * (sizes-1))/2
  df.between <- data.frame(Dist=Dist, Sobs = Sobs[upper.tri(Sobs)], 
                           Shat = Shat[upper.tri(Shat)], Sizes = Sizes[upper.tri(Sizes)],
                          row.names = NULL)
        
  df.within <- data.frame(Sobs.within = diag(Sobs), Shat.within = diag(Shat), 
                          Sizes = diag(Sizes), row.names=NULL)
  
  plot_pw(df.between)
  ggsave(filename = paste0(outpath, "/observed_vs_fitted-between.pdf"), width = 4, height = 4)
  plot_variogram(df.between)
  ggsave(paste0(outpath, "/semivariogam.pdf"), width = 4, height = 4)
  plot_ww(df.within)
  ggsave(paste0(outpath, "/observed_vs_fitted-within.pdf"), width = 4, height = 4)

}


plot_pw <- function(df){
  P <- ggplot(df) + geom_point(aes(y=Sobs, x=Shat, size = Sizes), alpha=0.6)  
  P <- P + scale_size_continuous(range = c(1, 10), name = "# pairs") + theme_classic() + 
    geom_abline(intercept=0) +
    ylab("genetic (observed) similarity between demes") +
    xlab("fitted similarity between demes")  
  P
}

plot_ww <- function(df){
  P <- ggplot(df) + geom_point(aes(y=Sobs.within, x=Shat.within, size = Sizes), alpha=0.6)  
  P <- P  + scale_size_continuous(range = c(1, 10), name = "# pairs") + theme_classic() + 
    geom_abline(intercept=0) +
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

#' plots the log-likelihood and log-posterior as a function of MCMC iteration (thinned)
#' @param mcmcpath path to output mcmc files
#' @param outpath where to save the file
#' @export
plot_trace <- function(mcmcpath, outpath){
  logll <- c()
  for (path in mcmcpath){
    pilogl <- scan(paste0(path, '/mcmcpilogl.txt'), quiet = TRUE)
    pilogl <- data.frame(matrix(pilogl, ncol = 2, byrow = TRUE))
    logll      <- cbind(logll, pilogl[,2])
  }
  
  colnames(logll)   <- paste0("run", 1:length(mcmcpath))
  itrs    <- 1:dim(logll)[1]
  logll   <- as.data.frame(cbind(logll, itrs)) %>% gather(variable, value, -itrs)
  
  
  g1 <- ggplot(logll, aes(x=itrs, y=value, color=variable)) + geom_line(size=0.5, alpha = 0.5) +
    xlab("iteration (thinned)") + ylab("logll")
  g2 <- ggplot(logll, aes(x=itrs, y=value, color=variable)) + geom_line() + facet_wrap(~variable) + 
    xlab("iteration (thinned)") + ylab("logll")
  ggsave(paste0(outpath, "/logll_trace.pdf"), arrangeGrob(g1, g2))
  
}


## not used currently
compute_e_summaries <- function(mcmcpath, outpath, type = "m"){
  mcmchyper <- read.table(paste0(mcmcpath, '/mcmc', type, 'hyper.txt'))
  colnames(mcmchyper) <- c("mu", "omega")
  rates <- scan(paste(mcmcpath,'/mcmc', type, 'rates.txt',sep=''),what=numeric(),quiet=TRUE)
  tiles <- read.table(paste0(mcmcpath, '/mcmc', type, 'tiles.txt'))$V1
  niter <- length(tiles)
  
  e_mean <- rep(0, niter)
  e_sd   <- rep(0, niter)
  count <- 0
  for (i in 1:niter) {
    now.tiles <- tiles[i]
    now.rates <- rates[(count+1):(count+now.tiles)]
    e.rates <- (log10(now.rates) - mcmchyper$mu[i]) / 10^(mcmchyper$omega[i])
    e_mean[i] <- mean(e.rates)
    e_sd[i]   <- sd(e.rates)
  }
  
  return(list(e_mean = e_mean, e_sd = e_sd))
  
}



#' plots the parameters of MAPS as function of iteration (thinned)
#' @param mcmcpath path to output mcmc files
#' @param outpath where to save the file
#' @export
plot_parameter_trace <- function(mcmcpath, outpath) {
  
  logll <- c()
  for (path in mcmcpath){
    pilogl <- scan(paste0(path, '/mcmcpilogl.txt'), quiet = TRUE)
    pilogl <- data.frame(matrix(pilogl, ncol = 2, byrow = TRUE))
    logll      <- cbind(logll, pilogl[,2])
  }
  
  colnames(logll)   <- paste0("run", 1:length(mcmcpath))
  itrs    <- 1:dim(logll)[1]
  logll   <- as.data.frame(cbind(logll, itrs)) %>% gather(variable, value, -itrs)
  
  mu_m <- c()
  omega_m <- c()
  mu_q <- c()
  omega_q <- c()
  qtiles <- c()
  mtiles <- c()
  for (path in mcmcpath){
    mcmcmhyper <- read.table(paste0(path, '/mcmcmhyper.txt'))
    mcmcqhyper <- read.table(paste0(path, '/mcmcqhyper.txt'))
    mcmcqtiles <- read.table(paste0(path, '/mcmcqtiles.txt'))
    mcmcmtiles <- read.table(paste0(path, '/mcmcmtiles.txt'))
    
    mu_m       <- cbind(mu_m, mcmcmhyper[,1])
    omega_m    <- cbind(omega_m, mcmcmhyper[,2])
    mu_q      <- cbind(mu_q, mcmcqhyper[,1])
    omega_q   <- cbind(omega_q, mcmcqhyper[,2])
    mtiles     <- cbind(mtiles, mcmcmtiles[,1])
    qtiles     <- cbind(qtiles, mcmcqtiles[,1])
    
  }
  
  colnames(omega_m) <- paste0("run", 1:length(mcmcpath))
  colnames(omega_q) <- paste0("run", 1:length(mcmcpath))
  colnames(mu_m)    <- paste0("run", 1:length(mcmcpath))
  colnames(mu_q)    <- paste0("run", 1:length(mcmcpath))
  colnames(mtiles)  <- paste0("run", 1:length(mcmcpath))
  colnames(qtiles)  <- paste0("run", 1:length(mcmcpath))
  itrs <- 1:dim(mcmcmtiles)[1]
  
  omega_m <- as.data.frame(cbind(omega_m, itrs)) %>% gather(variable, value, -itrs)
  omega_q <- as.data.frame(cbind(omega_q, itrs)) %>% gather(variable, value, -itrs)
  mu_m    <- as.data.frame(cbind(mu_m, itrs)) %>% gather(variable, value, -itrs)
  mu_q    <- as.data.frame(cbind(mu_q, itrs)) %>% gather(variable, value, -itrs)
  logll   <- as.data.frame(cbind(logll, itrs)) %>% gather(variable, value, -itrs)
  mtiles  <- as.data.frame(cbind(mtiles, itrs)) %>% gather(variable, value, -itrs)
  qtiles  <- as.data.frame(cbind(qtiles, itrs)) %>% gather(variable, value, -itrs)
  
  
  g1 <- ggplot(omega_m, aes(x=value, color=variable)) + geom_freqpoly(bins=30) + xlab("log10(omega_m)") + theme(legend.position="none")
  g2 <- ggplot(omega_q, aes(x=value, color=variable)) + geom_freqpoly(bins=30) + xlab("log10(omega_q)") + theme(legend.position="none")
  g3 <- ggplot(mu_m, aes(x=value, color=variable)) + geom_freqpoly(bins=30)    + xlab("mu_m") + theme(legend.position="none")
  g4 <- ggplot(mu_q, aes(x=value, color=variable)) + geom_freqpoly(bins=30)    + xlab("mu_q") + theme(legend.position="none")
  g5 <- ggplot(mtiles, aes(x=value, color=variable)) + geom_freqpoly(bins=30)  + xlab("mtiles") + theme(legend.position="none")
  g6 <- ggplot(qtiles, aes(x=value, color=variable)) + geom_freqpoly(bins=30)  + xlab("qtiles") + theme(legend.position="none")
  
  ggsave(paste0(outpath, "/parameters_trace.pdf"), arrangeGrob(g1, g2, g3, g4, g5, g6, ncol=2))
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
  N.scalingfactor = ndemes / total.area
  
  return(list(m.scalingfactor = m.scalingfactor, N.scalingfactor = N.scalingfactor))
}
