read_dimns <- function(path, longlat) {
  eems.output <- NULL
  nxmrks = NULL
  nymrks = NULL
  datapath.outer <- paste0(path, '.outer')
  mcmcpath.outer <- paste0(path, '/outer.txt')
  if (file.exists(mcmcpath.outer)) {
    eems.output <- TRUE
  } else if (file.exists(datapath.outer)) {
    eems.output <- FALSE
  }
  if (is.null(eems.output)) {
    stop(paste0(path, ' is neither a datapath nor a mcmcpath.'))
  }
  if (eems.output) {
    outer <- scan(mcmcpath.outer, what = numeric(), quiet = TRUE)
  } else {
    outer <- scan(datapath.outer, what = numeric(), quiet = TRUE)
  }
  outer <- matrix(outer, ncol = 2, byrow = TRUE)
  ## "Close" the outline if the first row is not the same as the last row
  if (sum(head(outer, 1) != tail(outer, 1))) {
    outer <- rbind(outer, head(outer, 1))
  }
  if (!longlat) { outer <- outer[, c(2, 1)] }
  xlim <- range(outer[, 1])
  ylim <- range(outer[, 2])
  aspect <- abs((diff(ylim) / diff(xlim)) / cos(mean(ylim) * pi / 180))
  ## Choose the number of interpolation in each direction
  if (is.null(nxmrks) && is.null(nymrks)) {
    if (aspect > 1) {
      nxmrks <- 100
      nymrks <- round(nxmrks * aspect)
    } else {
      nymrks <- 100
      nxmrks <- round(nymrks / aspect)
    }
  }
  ## The interpolation points are equally spaced
  xmrks <- seq(from = xlim[1], to = xlim[2], length = nxmrks)
  ymrks <- seq(from = ylim[1], to = ylim[2], length = nymrks)
  marks <- cbind(rep(xmrks, times = nymrks), rep(ymrks, each = nxmrks))
  return(list(nxmrks = nxmrks, xmrks = xmrks, xlim = xlim, xspan = diff(xlim),
              nymrks = nymrks, ymrks = ymrks, ylim = ylim, yspan = diff(ylim),
              marks = marks, nmrks = c(nxmrks, nymrks), aspect = aspect,
              outer = outer))
}

read_voronoi <- function(mcmcpath,longlat,is.mrates) {
  if (is.mrates) {
    rates <- scan(paste(mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
    tiles <- scan(paste(mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xseed <- scan(paste(mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    yseed <- scan(paste(mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
  } else {
    rates <- scan(paste(mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
    tiles <- scan(paste(mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xseed <- scan(paste(mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    yseed <- scan(paste(mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    rates <- 1/(2*rates) ## N = 1/2q
  }
  if (!longlat) {
    tempi <- xseed
    xseed <- yseed
    yseed <- tempi
  }
  return(list(rates=rates,tiles=tiles,xseed=xseed,yseed=yseed))
}

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

average_over_mcmcpaths <- function(mcmcpath, dimns, is.mrates, longlat=F){
  niter <- 0
  mean.rate <- matrix(0,dimns$nxmrks,dimns$nymrks)
  for (path in mcmcpath) {
    rslt <- standardize_rates(path,dimns,longlat,is.mrates)
    niter <- niter + rslt$niter
    mean.rate <- mean.rate + rslt$vals
  }
  mean.rate <- mean.rate/niter
  return(list(mean=mean.rate))
}

standardize_rates <- function(mcmcpath,dimns,longlat,is.mrates) {
  voronoi <- read_voronoi(mcmcpath,longlat,is.mrates)
  rates <- voronoi$rates
  tiles <- voronoi$tiles
  xseed <- voronoi$xseed
  yseed <- voronoi$yseed
  vals <- matrix(0,dimns$nxmrks,dimns$nymrks)
  niter <- length(tiles)
  count <- 0
  means  <- matrix(0, dimns$nxmrks, dimns$nymrks)
  for (i in 1:niter) {
    now.tiles <- tiles[i]
    now.rates <- rates[(count+1):(count+now.tiles)]
    now.xseed <- xseed[(count+1):(count+now.tiles)]
    now.yseed <- yseed[(count+1):(count+now.tiles)]
    now.seeds <- cbind(now.xseed,now.yseed)
    rslts <- compute_contour_vals(dimns$marks, now.rates, now.seeds)
    vals <- vals + matrix(rslts,  nrow = dimns$nxmrks, ncol = dimns$nymrks)
    count <- count + now.tiles
  }
  return(list(vals=vals,niter=niter))
}

plot_contour <- function(mcmcpath,dimns,
                                         is.mrates, ...) {
  if (is.mrates) {
    message('Plotting effective migration rates m : posterior mean and variance')
    files <- c('/mcmcmtiles.txt','/mcmcmrates.txt', 
               '/mcmcxcoord.txt','/mcmcycoord.txt')
  } else {
    message('Plotting effective diversity rates q : posterior mean and variance')
    files <- c('/mcmcqtiles.txt','/mcmcqrates.txt', 
               '/mcmcwcoord.txt','/mcmczcoord.txt')
  }
  
  mcmcpath <- check_files_at_path(files, mcmcpath)
  n_runs <- length(mcmcpath)
  ## return an error insteadif (n_runs==0) { return(0) }
  averages <- average_over_mcmcpaths(mcmcpath, dimns, is.mrates, T)
  avg <- averages[[1]]
  return(add.contour(mcmcpath, dimns, avg,
                                       is.mrates, ...))
}

add.contour <- function(mcmcpath, dimns, avg, is.mrates, ...){
  y <- rep(dimns$ymrks, each=length(dimns$xmrks))       
  x <- rep(dimns$xmrks, length(dimns$ymrks))
  
  df <- data.frame(x=x, y=y, avg=c(avg))
  
  mu = mean(log10(avg))
  a = 10^(mu - 1)
  b = 10^(mu + 1)
  limits = c(min(c(avg, a)),
             max(c(avg, b)))
  
  if(is.mrates){
    eems.colors <- scale_fill_gradientn(colours=default_eems_colors(),
                                        name="m", limits=limits, trans = "log10")
  } else {
    eems.colors <- scale_fill_gradientn(colours=default_eems_colors(),
                                        name="N", limits=limits, trans = "log10")
  }
  contour_plot <- ggplot() + 
    geom_tile(data=df, aes(x=y, y=x, fill=avg)) + eems.colors 
  
  return(contour_plot)
}
