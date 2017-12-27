#' @title 
#' @param 
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


#' @title 
#' @param 
#' @export
read_voronoi <- function(mcmcpath,params) {
  if (params$is.mrates) {
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
  if (!params$longlat) {
    tempi <- xseed
    xseed <- yseed
    yseed <- tempi
  }
  return(list(rates=rates,tiles=tiles,xseed=xseed,yseed=yseed))
}

#' @title 
#' @param 
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
  edges <- read.table(paste0(mcmcpath, '/edges.txt'), colClasses = numeric())
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

#' @title 
#' @param 
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


#' @title 
#' @param 
#' @export
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
              marks = marks, nmrks = c(nxmrks, nymrks),
              outer = outer))
}




