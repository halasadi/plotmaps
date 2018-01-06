#' @import mapproj
#' @import maps
#' @import ggplot2
#' @import dplyr

## constructs the dimensions of the main plot as a matrix representing a grid of pixels
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

## reads each voronoi tile
read_voronoi <- function(params) {
  if (params$is.mrates) {
    rates <- scan(paste(params$mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
    tiles <- scan(paste(params$mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xseed <- scan(paste(params$mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    yseed <- scan(paste(params$mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
  } else {
    rates <- scan(paste(params$mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
    tiles <- scan(paste(params$mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xseed <- scan(paste(params$mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    yseed <- scan(paste(params$mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    rates <- 1/(2*rates) ## N = 1/2q
  }
  if (!params$longlat) {
    tempi <- xseed
    xseed <- yseed
    yseed <- tempi
  }
  return(list(rates=rates,tiles=tiles,xseed=xseed,yseed=yseed))
}


add_graph <- function(g, color="black"){
  xstart <- g$demes[g$edges[,1],1]
  xend <- g$demes[g$edges[,2],1]
  ystart <- g$demes[g$edges[,1],2]
  yend <- g$demes[g$edges[,2],2]
  grid <- data.frame(xstart, xend, ystart, yend)
  geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend), data=grid, color=color, alpha=0.3)
}

compute_summary_statistic <- function(params, dimns){
  rslts <- compute_rates_each_pixel(params,dimns)
  
  l <- compute_scaling(params$mcmcpath)
  
  if (params$is.mrates & params$add.countries){
    rslts <- sqrt(rslts * l$m.scalingfactor)
  } else if (!params$is.mrates & params$add.countries){
    rslts <- rslts * l$N.scalingfactor
  }
  
  mean.rate = NA
  med.rate = NA
  upper.ci = NA
  
  # mean 
  if (params$plot.mean){
    mean.rate <- apply(rslts, c(2,3), mean)
  }

  # median
  if (params$plot.median | params$plot.sign){
    med.rate  <- apply(rslts, c(2,3), median)
  }

  if (params$plot.sign) {
    # measure of credibility of a barrier
    base.rate <- mean(med.rate)
    upper.ci  <- apply(rslts, c(2,3), function(x){ mean(x > base.rate)})
  }
  
  return(list(avg=mean.rate, med = med.rate, upper.ci=upper.ci))
}

compute_rates_each_pixel <- function(params,dimns) {
  voronoi <- read_voronoi(params)
  rates <- voronoi$rates
  tiles <- voronoi$tiles
  xseed <- voronoi$xseed
  yseed <- voronoi$yseed
  vals <- matrix(0,dimns$nxmrks,dimns$nymrks)
  niter <- length(tiles)
  count <- 0
  means  <- matrix(0, dimns$nxmrks, dimns$nymrks)
  pixelated.tiles <- array(dim = c(niter, dimns$nxmrks, dimns$nymrks))
  for (i in 1:niter) {
    now.tiles <- tiles[i]
    now.rates <- rates[(count+1):(count+now.tiles)]
    now.xseed <- xseed[(count+1):(count+now.tiles)]
    now.yseed <- yseed[(count+1):(count+now.tiles)]
    now.seeds <- cbind(now.xseed,now.yseed)
    rslts <- compute_contour_vals(dimns$marks, now.rates, now.seeds)
    pixelated.tiles[i,,] <- matrix(rslts,  nrow = dimns$nxmrks, ncol = dimns$nymrks)
    count <- count + now.tiles
  }
  return(pixelated.tiles)
}

#' @title adds contours
#' @useDynLib plotmaps
#' @param mcmcpath path to the mcmc output files
#' @param dimns grid layout (where each entry is a pixel) to build the picture
#' @param g a ggplot2 object
#'
#' @return ggplot2 object
#' @export
## We need to include the useDynLib somewhere in our documentation
## see: https://stackoverflow.com/questions/36605955/c-function-not-available
add_contours <- function(params,dimns, g) {

  n_runs <- length(params$mcmcpath)
  if (n_runs==0) { stop(paste0(params$mcmcpath, ' is not a valid mcmcpath')) }
  if (n_runs>1)  {  stop(paste0(params$mcmcpath, ' the MAPS plotting package plots only one run')) }
  
  
  summary_stats <- compute_summary_statistic(params, dimns)
  
  if (sum(c(params$plot.median, params$plot.mean, params$plot.sign)) != 1){
    stop("you must plot only one of the following: mean, median, or probability of positive sign")
  }
  
  return(add_contour(params, dimns, summary_stats, g))
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

add_contour <- function(params, dimns, summary_stats, g){
  x <- dimns$marks[,1]
  y <- dimns$marks[,2]
  
  if (params$plot.median){
    summary_stat <- summary_stats$med
  } else if (params$plot.mean) {
    summary_stat <- summary_stats$avg
  } else if (params$plot.sign) {
    summary_stat <- summary_stats$upper.ci
  }
  
  df <- data.frame(x=x, y=y, ss=c(summary_stat))
  
  mu <- mean(log10(summary_stat))
  a <- 10^(mu - 1)
  b <- 10^(mu + 1)
  colours <- default_eems_colors()
  
  if (params$plot.sign){
    
    legend.title <- "p(x > mean)"
    limits <- c(0, 1)
    trans <- "identity"
    
  } else {

    limits <- c(min(c(df$ss, a)),
               max(c(df$ss, b)))
    trans <- "log10"
    
    if (params$is.mrates) {
      
      if (params$add.countries){
        legend.title <- expression(paste(frac(km, sqrt(gen))))
      } else{
        legend.title <- "m"
      }
      
    } else {
      
      if (params$add.countries){
        legend.title <- expression(paste(frac(N, km^2)))
      } else{
        legend.title <- "N"
        trans <- "identity"
        
      }
    }
    
    if (params$plot.difference){
      legend.title <- paste("+log", legend.title)
      colours <- inferno_colors()
      df$ss <- log10(df$ss)
      mu <- mean(df$ss)
      trans <- "identity"
      
      if (!params$is.mrates){
        # need to offset by log10(2) because MAPS MCMC parameterized on log10 scale q coalescent
        # and here we plot, N = log(1/2q) => log10(N) = -(log10(2) + log10(eps))
        # where eps \approx 0
        df$ss  <- df$ss + log10(2)
        mu <- mu + log10(2)
      }
      
      a <- mu - 1
      b <- mu + 1
      limits <- c(min(c(df$ss, a)),
                  max(c(df$ss, b)))
    }
    
  }
  eems.colors <- scale_fill_gradientn(colours=colours,
                                      name=legend.title, limits = limits, 
                                      trans = trans)
  
  g <- g + geom_tile(data=df, aes(x=x, y=y, fill=ss), alpha = 1) + 
    eems.colors + theme(legend.key.width=unit(0.75, 'cm')) + 
    theme(legend.key.height=unit(2.25, 'cm')) + 
    theme(legend.text=element_text(size=15)) + coord_fixed() +
    theme(legend.title =element_text(size=15)) 
  
  graph <- read_graph(params$mcmcpath, params$longlat)
  
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

#' make_base
#' @param dimns the grid layout (nrow*ncol=npixels) to build the picture
#' @return ggplot2 pobject
#' @export
#' 
make_base <- function(dimns, params, g){
  
  g<- ggplot() + theme_classic() + theme(axis.line = element_blank(), axis.ticks=element_blank(),
                                axis.text.y=element_blank(), axis.text.x=element_blank(),
                                axis.title.x=element_blank(), axis.title.y=element_blank())
  
  if (params$add.countries){
    boundary <- dimns$outer
    bbox <- c(left=min(boundary[,1]), right=max(boundary[,1]),
              bottom=min(boundary[,2]), top=max(boundary[,2]))
    bbox['top'] <- pmin(bbox['top'], 83)
    m_boundary <- get_boundary_map(bbox)
    g = g + coord_map("mercator", parameters=NULL,  xlim=bbox[c('left', 'right')],
                      ylim=bbox[c('bottom', 'top')]) + 
      xlim(bbox[c('left', 'right')])+ 
      ylim(bbox[c('bottom', 'top')])
  }
  
  return(g)
}


#' add_map
#' @param dimns the grid layout (nrow*ncol=npixels) to build the picture
#' @return ggplot2 pobject
#' @export
#' 
add_map <- function(dimns, params, g){
  boundary <- dimns$outer
  bbox <- c(left=min(boundary[,1]), right=max(boundary[,1]),
            bottom=min(boundary[,2]), top=max(boundary[,2]))
  bbox['top'] <- pmin(bbox['top'], 83)
  m_boundary <- get_boundary_map(bbox)
  g = g + geom_path(data=m_boundary, aes(x=long, y=lat, group=group),  color='black', size=0.5, 
                    alpha = 1)
  return(g)
}

#' plot contour plot
#' @param params a list containing the plot options
#' @export
#' 
plot_contour <- function(params){
  dimns <- read_dimns(params$mcmcpath, params$longlat)
  
  
  g <- make_base(dimns, params) 
  g <- add_contours(params, dimns, g) 
  if (params$add.countries){
    g <- add_map(dimns, params, g)
  }

  
  filename <- params$outpath
  if (params$is.mrates){
    filename <- paste0(filename, "/mrates-")
  } else{
    filename <- paste0(filename, "/Nsizes-")
  }
  
  if (params$plot.mean){
    filename <- paste0(filename, "mean.pdf")
  } else if (params$plot.median) {
    filename <- paste0(filename, "median.pdf")
  } else{
    filename <- paste0(filename, "signplot.pdf")
  }

  ggsave(filename, width = params$width, height = params$height)
}


#' plot_maps maining plotting function
#' @param add.pts mark demes that are sampled (boolean)
#' @param add.graph overlay the underlying graph (boolean) 
#' @param add.countries add country borders (boolean)
#' @param plot.mean plots the median if TRUE and mean if FALSE (boolean)
#' @param longlat order of longitude/lattitude in the .coords and .outer file (boolean)
#' @param mcmcpath path to mcmc output files
#' @param outpath files will be written to this path
#' @param width width of the main MAPS plot
#' @param height height of the main MAPS plot
#' @param plot.difference set TRUE if user ran MAPS with the olderpath parameter set, this way
#'                        MAPS only estimates the difference (boolean)
#' @export
#' 
plot_maps <- function(add.pts = TRUE, add.graph = TRUE, add.countries = FALSE,
                     plot.median = TRUE, longlat, mcmcpath, outpath,
                     width = 10, height = 6, plot.difference=FALSE){
  
  dir.create(file.path(outpath), showWarnings = FALSE)
  
  files <- c('/mcmcmtiles.txt','/mcmcmrates.txt', 
             '/mcmcxcoord.txt','/mcmcycoord.txt',
             '/mcmcqtiles.txt','/mcmcqrates.txt',
             '/mcmcpilogl.txt', '/rdistoDemes.txt',
             '/rdistJtDobsJ.txt', '/rdistJtDhatJ.txt')
  mcmcpath <- check_files_at_path(files, mcmcpath)
  
  if (plot.median){
    plot.mean= FALSE
  } else {
    plot.mean = TRUE
  }

  params <- list(mcmcpath = mcmcpath, outpath = outpath, longlat = longlat,
                 is.mrates = TRUE, plot.mean = plot.mean, plot.median = plot.median,
                 plot.sign = FALSE, width = width, height = height, 
                 add.countries = add.countries, add.graph = add.graph, add.pts = add.pts, 
                 plot.difference = plot.difference)
  
  
  message('plotting migration surface')
  plot_contour(params)
  

  message('plotting migration surface sign plot')
  params$is.mrates = TRUE
  params$plot.sign = TRUE
  params$plot.mean = FALSE
  params$plot.median = FALSE
  plot_contour(params)
  
  message('plotting population-size surface')
  params$is.mrates = FALSE
  params$plot.sign = FALSE
  params$plot.mean = plot.mean
  params$plot.median = plot.median
  plot_contour(params)
  
  message('plotting population-size surface sign plot')
  params$is.mrates = FALSE
  params$plot.sign = TRUE
  params$plot.mean = FALSE
  params$plot.median = FALSE
  plot_contour(params)
  
  message('plotting diagonostics of model fit and MCMC convergence')
  plot_trace(mcmcpath, outpath) 
  plot_fit_data(mcmcpath, outpath, params$longlat)
}
