#' @import mapproj
#' @import maps
#' @import ggplot2
#' @import dplyr


library(ggplot2)
library(maps)
library(mapproj)
library(dplyr)


add_graph <- function(g, color="black"){
  xstart <- g$demes[g$edges[,1],1]
  xend <- g$demes[g$edges[,2],1]
  ystart <- g$demes[g$edges[,1],2]
  yend <- g$demes[g$edges[,2],2]
  grid <- data.frame(xstart, xend, ystart, yend)
  geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend), data=grid, color=color, alpha=0.3)
}

average_over_mcmcpaths <- function(mcmcpath, dimns, params, longlat=F){
  niter <- 0
  mean.rate <- matrix(0,dimns$nxmrks,dimns$nymrks)
  for (path in mcmcpath) {
    rslt <- standardize_rates(path,dimns,longlat,params)
    niter <- niter + rslt$niter
    mean.rate <- mean.rate + rslt$vals
  }
  mean.rate <- mean.rate/niter
  return(list(mean=mean.rate))
}

standardize_rates <- function(mcmcpath,dimns,longlat,params) {
  voronoi <- read_voronoi(mcmcpath,longlat,params)
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

#' @title adds contours
#' @useDynLib plotmaps
#' @param mcmcpath 
#' @param dimns
#' @param g a ggplot2 object
#' @param params list of parameters
#'
#' @return ggplot2 object
#' @export
## We need to include the useDynLib somewhere in our documentation
## see: https://stackoverflow.com/questions/36605955/c-function-not-available
add_contours <- function(mcmcpath,dimns, g, params) {
  if (params$is.mrates) {
    message('Plotting migration rates m : posterior mean')
    files <- c('/mcmcmtiles.txt','/mcmcmrates.txt', 
               '/mcmcxcoord.txt','/mcmcycoord.txt')
  } else {
    message('Plotting population sizes N : posterior mean')
    files <- c('/mcmcqtiles.txt','/mcmcqrates.txt', 
               '/mcmcwcoord.txt','/mcmczcoord.txt')
  }
  
  mcmcpath <- check_files_at_path(files, mcmcpath)
  n_runs <- length(mcmcpath)
  ## return an error insteadif (n_runs==0) { return(0) }
  averages <- average_over_mcmcpaths(mcmcpath, dimns, params, T)
  avg <- averages[[1]]
  return(add_contour(mcmcpath, dimns, avg,
                                       g, params))
}

add_pts <- function(g, color="#efefefdd", const_size=T){
  tbl <- table(g$ipmap)
  ind <- as.numeric(names(tbl))
  sizes <- as.vector(tbl)
  df <- data.frame(x=g$demes[ind,1], y=g$demes[ind,2], sizes=sizes)
  if(const_size) {
    pts <- geom_point(aes(x=x, y=y), data=df, color=color, size=1.5)
  } else {
    pts <- geom_point(aes(x=x, y=y, size=sizes), data=df, color=color)
  }
}

add_contour <- function(mcmcpath, dimns, avg, g, params){
  x <- dimns$marks[,1]
  y <- dimns$marks[,2]
  
  df <- data.frame(x=x, y=y, avg=c(avg))
  
  mu = mean(log10(avg))
  a = 10^(mu - 1)
  b = 10^(mu + 1)
  limits = c(min(c(avg, a)),
             max(c(avg, b)))
  
  if (params$is.mrates) {
    eems.colors <- scale_fill_gradientn(colours=default_eems_colors(),
                                        name="m", limits=limits, trans = "log10")
  } else {
    eems.colors <- scale_fill_gradientn(colours=default_eems_colors(),
                                        name="N", limits=limits, trans = "log10")
  }
  
  g <- g + geom_tile(data=df, aes(x=x, y=y, fill=avg), alpha = 0.4) + 
    eems.colors + coord_fixed()
  
  graph <- read_graph(mcmcpath, longlat)
  
  if (params$add.graph){
    g <- g + add_graph(graph)
  }
  
  if (params$add.pts){
    g <- g + add_pts(graph, col = "black")
  }
  return(g)
}

get_boundary_map <- function(bbox){
  xlim=bbox[c('left', 'right')]
  ylim=bbox[c('bottom', 'top')]
  m1 <- map(interior=F, plot=F, 
            xlim=xlim, 
            ylim=ylim,
             resolution=0)
  m1$group <- cumsum(is.na(m1$x))
  m <- data.frame(long=m1$x, lat=m1$y,
                  group=m1$group)
  m <- m %>% filter(!is.na(long)) %>% filter(long > xlim[1]) %>% 
    filter(long < xlim[2]) %>% filter(lat > ylim[1]) %>% filter(lat < ylim[2])
  #m$long[m$long< -30] <- m$long[m$long< -30] +360   
  #m$lat[m$lat< -38] <- -38
  return(m)
}

#' make_map
#' @param 
#' @export
#' 
make_map <- function(dimns, params){
  boundary <- dimns$outer
  bbox <- c(left=min(boundary[,1]), right=max(boundary[,1]),
            bottom=min(boundary[,2]), top=max(boundary[,2]))
  bbox['top'] <- pmin(bbox['top'], 83)
  
  g <- ggplot()

  g=g+theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=12))         
  g=g+theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=12))
  g=g + theme_classic()
  
  if (params$add.countries){
    m_boundary <- get_boundary_map(bbox)
    g = g + coord_map("mercator", parameters=NULL,  xlim=bbox[c('left', 'right')],
                      ylim=bbox[c('bottom', 'top')]) + 
      xlim(bbox[c('left', 'right')])+ 
      ylim(bbox[c('bottom', 'top')])

    g = g + geom_path(data=m_boundary, aes(x=long, y=lat, group=group),  color='black', size=0.5, 
                      alpha = 1)
  }
  
  return(g)
}

