#' @import mapproj
#' @import maps
#' @import ggplot2
#' @import dplyr
#' @import sp
#' @import abind
#' @import tidyr

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
  
  pip <- point.in.polygon(point.x=marks[,1], point.y=marks[,2], pol.x = outer[,1], pol.y = outer[,2])
  filter <- pip == 1
  
  return(list(nxmrks = nxmrks, xmrks = xmrks, xlim = xlim, xspan = diff(xlim),
              nymrks = nymrks, ymrks = ymrks, ylim = ylim, yspan = diff(ylim),
              marks = marks, nmrks = c(nxmrks, nymrks),
              outer = outer, filter = filter))
}

## reads each voronoi tile
read_voronoi <- function(params, path) {
  if (params$is.mrates) {
    rates <- scan(paste(path,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
    tiles <- scan(paste(path,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xseed <- scan(paste(path,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    yseed <- scan(paste(path,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
    
    # old verion of runs
    #if (params$plot.difference){
    #  rates <- log10(rates)
    #}
    
  } else {
    rates <- scan(paste(path,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
    tiles <- scan(paste(path,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xseed <- scan(paste(path,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    yseed <- scan(paste(path,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    
    #if (params$plot.difference){
    #  rates <- -log10(rates)
    #} else{
    #  rates <- 1 / (2 * rates)
    #}
    if (params$plot.difference){
      rates <- -rates
    } 

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
  #geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend), data=grid, color=color, alpha=0.3)
  geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend), data=grid, color=color, alpha=1)
}

compute_summary_statistic <- function(params, dimns){
  
  rslts <- compute_rates_each_pixel(params,dimns)

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
  
  #if (params$is.mrates){
  #  old_rates = mean.rate
  #  save(old_rates, file = "1_5_old.rda")
  #  load(file = "1_5_old.rda")
  #  mean.rate = mean.rate - old_rates
  #}
  
  if (!params$plot.difference){
    mean.rate <- 10^mean.rate
    med.rate  <- 10^med.rate
    
    if (!params$is.mrates){
      mean.rate <- 1 / (2 * mean.rate)
      med.rate <- 1 / (2 * med.rate)
    }
    
  }
  
  l <- compute_scaling(params$mcmcpath[1])
  if (params$is.mrates & params$add.countries & !params$plot.difference){
    mean.rate <- sqrt(mean.rate * l$m.scalingfactor)
    med.rate <- sqrt(med.rate * l$m.scalingfactor)
  } else if (!params$is.mrates & params$add.countries & !params$plot.difference){
    mean.rate <- mean.rate * l$N.scalingfactor
    med.rate <-  med.rate * l$N.scalingfactor
  }
  
  return(list(avg=mean.rate, med = med.rate, upper.ci=upper.ci))
}

compute_rates_each_pixel <- function(params, dimns) {
  for (p in 1:length(params$mcmcpath)){
    path = params$mcmcpath[p]
    voronoi <- read_voronoi(params, path)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    vals  <- matrix(0, dimns$nxmrks, dimns$nymrks)
    niter <- length(tiles)
    count <- 0
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

    if (p == 1){
      ret <- pixelated.tiles
    } else{
      ret <- abind(ret, pixelated.tiles, along = 1)
    }
    
  } # end loop over paths

  return(ret)
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

  summary_stats <- compute_summary_statistic(params, dimns)

  return(add_contour(params, dimns, summary_stats, g))
}

add_pts <- function(g, const_size=F){
  tbl <- table(g$ipmap)
  ind <- as.numeric(names(tbl))
  sizes <- as.vector(tbl)
  
  #coords <- read.table('/Users/halasadi/eems2/data/Ralph+Coop/runs/overall/r1/popres_1_Inf/popres_1_Inf.coord')
  #coords <- unique(coords)
  #df <- data.frame(x = coords[,2], y = coords[,1])
  #pts <- geom_point(aes(x=x, y=y), data=df, color = "red", size = 3)
  
  df <- data.frame(x=g$demes[ind,1], y=g$demes[ind,2], sizes=sizes)
  if(const_size) {
    pts <- geom_point(aes(x=x, y=y), data=df, color="red", fill = "black", size=1.5)
  } else {
    pts <- geom_point(aes(x=x, y=y, size=sizes), data=df, colour = "red",
                      fill = "black", pch = 21, show.legend=FALSE)
  }
}

#' plot_voronoi_samples
#' @param is.mrates plot mvoronoi (TRUE) or qvoronoi (FALSE)
#' @param add.countries add country borders (boolean)
#' @param longlat order of longitude/lattitude in the .coords and .outer file (boolean)
#' @param mcmcpath path to mcmc output files
#' @param outpath files will be written to this path
#' @param width width of the main MAPS plot
#' @param height height of the main MAPS plot
#' @param nsamples number of samples to plot from posterior
#' @export
#' 
plot_voronoi_samples <- function(longlat, mcmcpath, outpath, is.mrates=TRUE, 
                                 nsamples = 10, add.countries = TRUE, width = 5, height = 5){
  
  dir.create(file.path(outpath), showWarnings = FALSE)
  files <- c('/mcmcmtiles.txt','/mcmcmrates.txt', 
             '/mcmcxcoord.txt','/mcmcycoord.txt',
             '/mcmcqtiles.txt','/mcmcqrates.txt',
             '/mcmcpilogl.txt', '/rdistoDemes.txt',
             '/rdistJtDobsJ.txt', '/rdistJtDhatJ.txt')
  mcmcpath <- check_files_at_path(files, mcmcpath)
  
  
  params <- list(mcmcpath = mcmcpath, outpath = outpath, longlat = longlat,
                 is.mrates = is.mrates, add.countries = add.countries,
                 plot.difference=FALSE)
  
  dimns <- read_dimns(params$mcmcpath[1], params$longlat)
  x <- dimns$marks[,1]
  y <- dimns$marks[,2]

  rslts <- compute_rates_each_pixel(params,dimns)
  niter <- dim(rslts)[1]
  samples <- sample(1:niter, nsamples)
  graph <- read_graph(params$mcmcpath[1], params$longlat)
  eems.colors <- scale_fill_gradientn(colours=default_eems_colors(), 
                                      trans = "log10")
  
  for (i in 1:length(samples)){
    df <- data.frame(x=x, y=y, rate=c(rslts[samples[i],,]))
    g<- ggplot() + theme_classic() + 
      theme(axis.line = element_blank(), axis.ticks=element_blank(),
            axis.text.y=element_blank(), axis.text.x=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank()) +
      geom_raster(data=df, aes(x=x, y=y, fill=rate), alpha = 1) + eems.colors +
      add_graph(graph) #+ add_pts(graph)
    if (params$add.countries){
      g <- add_map(dimns, params, g)
    }
    
    prefix <- 'm'
    if (!is.mrates){
      prefix <- 'q'
    }
    
    g
    ggsave(paste0(outpath, "/", prefix, "draw-", i, ".pdf"), width = width, height = height)
  }
}

get_title <- function(params){
  
  m          <- expression(m)
  N          <- expression(N)
  diff.m     <- bquote(paste(log10, bgroup("(", frac(m^"'", m),  ")")))
  diff.N     <- bquote(paste(log10, bgroup("(", frac(N^"'", N),  ")")))
  sigma      <- expression(sigma)
  D          <- expression('D'[e])
  diff.sigma <- bquote(paste(log10, bgroup("(", frac(sigma^"'", sigma),  ")")))
  diff.D     <- bquote(paste(log10, bgroup("(", frac('D'[e]^"'", 'D'[e]),  ")")))
  
  
  if (!params$plot.difference){
    
    if (params$add.countries){
      
      if (params$is.mrates){
        return(sigma)
      } else{
        return(D)
      }
    } else{
      if (params$is.mrates){
        return(m)
      } else{
        return(N)
      }
    
    }
  } else{
    
    if (params$add.countries){
      
      if (params$is.mrates){
        return(diff.sigma)
      } else{
        return(diff.D)
      }
    } else{
      if (params$is.mrates){
        return(diff.m)
      } else{
        return(diff.N)
      }
    }
    
  }

}


get_limits <- function(df, params){
  if (!params$plot.difference){
    mu <- mean(log10(df$ss))
    a <- 10^(mu - 2)
    b <- 10^(mu + 2)
    limits <- c(min(c(df$ss, a)),
                max(c(df$ss, b)))
  }
  
  if (params$plot.difference){
    mu <- 0
    a <- mu - 2
    b <- mu + 2
    limits <- c(min(c(df$ss, a)),
               max(c(df$ss, b)))
  }
  
  return(limits)
}
 
get_trans <- function(params){
  
  if (!is.na(params$trans)){
    trans = params$trans
    return(trans)
  }
  
  if (params$plot.difference){
    return("identity")
  }
  return("log10")  
}

get_color_gradient <- function(params, legend.title, limits, trans){

  if (params$set.range){
    if (params$is.mrates){
      my.limits <- params$m.limits
    } else{
      my.limits <- params$N.limits
    }
    color.gradient <- scale_fill_gradientn(colours=default_eems_colors(),
                                        name=legend.title, 
                                       trans = trans, na.value = "lightgray",
                                       limits = my.limits) 
  } else{
    color.gradient <- scale_fill_gradientn(colours=default_eems_colors(),
                                        name=legend.title, limits = limits, 
                                        trans = trans, na.value = "lightgray")
  }
  
  return(color.gradient)
}

add_contour <- function(params, dimns, summary_stats, g){
  x <- dimns$marks[,1]
  y <- dimns$marks[,2]
  
  if (params$plot.mean){
    summary_stat <- summary_stats$avg
  } else{
    summary_stat <- summary_stats$med
  }
  
  
  
  df <- data.frame(x=x, y=y, ss=c(summary_stat), filter = dimns$filter)
  print(summary(df$ss))
  df <- df[df$filter,]
  
  
  legend.title   <- get_title(params)
  trans          <- get_trans(params)
  limits         <- get_limits(df, params)
  color.gradient <- get_color_gradient(params, legend.title, limits, trans)
  
  if (params$plot.sign){
    ci  <- data.frame(x=x, y=y, ss=c(summary_stats$upper.ci), filter = dimns$filter)
    ci  <- ci[ci$filter,]
    df$ss[ci$ss > params$alpha  & ci$ss < (1-params$alpha)] <- NA
  }
  
  g <- g + geom_raster(data=df, aes(x=x, y=y, fill=ss), alpha = 1) + 
    color.gradient + theme(legend.key.width=unit(0.75, 'cm')) + 
    theme(legend.key.height=unit(2.25, 'cm')) + 
    theme(legend.text=element_text(size=15)) + coord_fixed() +
    theme(legend.title =element_text(size=15)) 
  
  graph <- read_graph(params$mcmcpath[1], params$longlat)
  
  if (params$add.graph){
    g <- g + add_graph(graph)
  }
  
  if (params$add.pts){
    g <- g + add_pts(graph)
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
  g = g + geom_path(data=m_boundary, aes(x=long, y=lat, group=group),  color='black', size=1, 
                    alpha = 1)
  return(g)
}

#' plot contour plot
#' @param params a list containing the plot options
#' @export
#' 
plot_contour <- function(params){
  dimns <- read_dimns(params$mcmcpath[1], params$longlat)
  
  
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
    if (params$plot.sign){
      filename <- paste0(filename, "mean-sign.pdf")
    } else{
      filename <- paste0(filename, "mean.pdf")
    }
  } else {
    if (params$plot.sign){
      filename <- paste0(filename, "median-sign.pdf")
    } else{
      filename <- paste0(filename, "median.pdf")
    }
  }
  
  ggsave(filename, width = params$width, height = params$height)
}


#' plot_maps maining plotting function
#' @param add.pts mark demes that are sampled (boolean)
#' @param add.graph overlay the underlying graph (boolean) 
#' @param add.countries add country borders (boolean)
#' @param plot.mean plots the mean if TRUE and mean if FALSE (boolean)
#' @param longlat order of longitude/lattitude in the .coords and .outer file (boolean)
#' @param mcmcpath path to mcmc output files
#' @param outpath files will be written to this path
#' @param width width of the main MAPS plot
#' @param height height of the main MAPS plot
#' @param plot.difference set TRUE if user ran MAPS with the olderpath parameter set, this way
#'                        MAPS only estimates the difference (boolean)
#' @param set.range TRUE if range is set as custom
#' @param m.limits limits of the migration surface plot (used if set.range==TRUE)
#' @param N.limits limits of the population-size surface plot (used if set.range==TRUE)
#' @param alpha the significance level, such that rates exceeding that level will be shown in the sign plot
#' @export
#' 
plot_maps <- function(add.pts = TRUE, add.graph = TRUE, add.countries = FALSE,
                     plot.mean = TRUE, longlat, mcmcpath, outpath,
                     width = 10, height = 6, plot.difference=FALSE, set.range = FALSE,
                     alpha = 0.1, m.limits = NA, N.limits = NA, trans = NA){
  
  dir.create(file.path(outpath), showWarnings = FALSE)
  
  files <- c('/mcmcmtiles.txt','/mcmcmrates.txt', 
             '/mcmcxcoord.txt','/mcmcycoord.txt',
             '/mcmcqtiles.txt','/mcmcqrates.txt',
             '/mcmcpilogl.txt', '/rdistoDemes.txt',
             '/rdistJtDobsJ.txt', '/rdistJtDhatJ.txt')
  mcmcpath <- check_files_at_path(files, mcmcpath)
  
  if (plot.mean){
    plot.median= FALSE
  } else {
    plot.median = TRUE
  }

  params <- list(mcmcpath = mcmcpath, outpath = outpath, longlat = longlat,
                 is.mrates = TRUE, plot.mean = plot.mean, plot.median = plot.median,
                 plot.sign = FALSE, width = width, height = height, 
                 add.countries = add.countries, add.graph = add.graph, add.pts = add.pts, 
                 plot.difference = plot.difference, set.range = set.range, alpha = alpha,
                 m.limits = m.limits, N.limits = N.limits, trans = trans)
  
  
  if (!is.na(params$trans)){
    if (!(params$trans %in% c("identity", "log10"))){
      stop("trans incorrectly specified, trans = identity or log10")
    }
  }
  
  if (params$set.range) {
    if (length(params$m.limits) < 2 | length(params$N.limits) < 2){
      stop("set.range = TRUE but m.limits or N.limits not specified")
    }
  }
  
  message('plotting migration surface')
  plot_contour(params)
  
  message('plotting population-size surface')
  params$is.mrates = FALSE
  params$plot.sign = FALSE
  plot_contour(params)
  
  if (!params$plot.difference){
    
    message('plotting migration surface sign plot')
    params$is.mrates = TRUE
    params$plot.sign = TRUE
    plot_contour(params)

    message('plotting population-size surface sign plot')
    params$is.mrates = FALSE
    params$plot.sign = TRUE
    plot_contour(params)
  }
  
  message('plotting diagonostics of model fit and MCMC convergence')
  plot_trace(mcmcpath, outpath) 
  #plot_parameter_distribution(mcmcpath, outpath)
  plot_fit_data(mcmcpath, outpath, params$longlat)
}
