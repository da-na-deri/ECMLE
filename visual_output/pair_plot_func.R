packages <- c("IDPmisc", "RColorBrewer", "colorRamps")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


# Custom plotting functions from your provided code
panel.hist2 <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  
  h <- hist(x, breaks = nB*4, plot = FALSE) #use 4 times the breaks!!!!!!!!
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; 
  if(max(y) > 0) y <- y/max(y)  # avoid division by zero
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

myImage = function (x, y = NULL, pixs = 1, zmax = NULL, ztransf = function(x) {x}, colramp = IDPcolorRamp, factors = c(FALSE, FALSE), matrix = FALSE, add = FALSE) 
{
  ## function hacked from IDPmisc, as I needed the values and not the histogram ##
  
  xy <- NaRV.omit(getXY(x, y))
  factors <- factors | sapply(xy, is.factor)
  xy <- sapply(xy, as.numeric)
  pixs <- (pixs/10)/2.54
  usr <- par("usr")
  if (factors[1]) {
    bx <- seq(min(xy[, 1] - 0.25), max(xy[, 1] + 0.25), length.out = 2 * diff(range(xy[, 1])) + 2) #bx and by are bins!
  }
  else {
    bx <- seq(usr[1], usr[2], length.out = round(par("pin")/pixs)[1] + 1)
  }
  if (factors[2]) {
    by <- seq(min(xy[, 2] - 0.25), max(xy[, 2] + 0.25), length.out = 2 * diff(range(xy[, 2])) + 2)
  }
  else {
    by <- seq(usr[3], usr[4], length.out = round(par("pin")/pixs)[2] + 1)
  }
  zz <- ztransf(table(cut(xy[, 1], breaks = bx), cut(xy[, 2], breaks = by)))
  zzmax <- ceiling(max(zz))
  if (is.null(zmax)) zmax <- zzmax
  else zmax <- ceiling(zmax)
  if (zmax < 1 || is.null(zmax)) {stop("zmax must be >= 1 and\n          plot(x,y,...) must have been called before calling this function!\n")}
  if (zzmax > zmax) stop("zmax too small! Densiest aereas are out of range!", call. = FALSE)
  zmax <- max(zmax, 2)
  lbx <- length(bx)
  lby <- length(by)
  xx <- 0.5 * (bx[-1] + bx[-lbx])
  yy <- 0.5 * (by[-1] + by[-lby])
  
  # Add the image plotting part like ipairs does
  if (!add) {
    image(x = xx, y = yy, zz, col = colramp(zmax), breaks = seq(0.5, zmax + 1, 1), xaxs = "r", yaxs = "r", add = TRUE)
    box()
  } else {
    image(x = xx, y = yy, zz, col = colramp(zmax), breaks = seq(0.5, zmax + 1, 1), xaxs = "r", yaxs = "r", add = TRUE)
  }
  
  if (matrix) invisible(list(xx,yy,zz,bx,by)) ######outputs list with bins as well
  else invisible(zzmax)
}

conf_contour2 = function(hist2d_obj, alpha, normalise=F, add=F, cl="black", lab=T, Xlim, Ylim) { #flattest, not edge
  ## based on drawcontour.R, written by J.D. Forester, 17 March 2008 ##
  znorm <- hist2d_obj[[3]] / ifelse(normalise,sum(hist2d_obj[[3]]),1)
  zdens <- sort(znorm, decreasing=T)
  Czdens <- cumsum(zdens)
  Czdens <- (Czdens/Czdens[length(zdens)]) #normalise levels
  
  for(cont.level in 1:length(alpha)){
    tmpp = max(which(Czdens<=alpha[cont.level]))
    if( is.infinite(tmpp) ) {
      cat("skipped level",alpha[cont.level],"\n", sep=" ")
      next
    }
    crit.val <- zdens[tmpp]
    contour(x=hist2d_obj[[1]], y=hist2d_obj[[2]], z=znorm, levels = crit.val, add=ifelse(add,T,ifelse(cont.level==1,F,T) ), method="flattest", labels = alpha[cont.level], drawlabels = lab, col = cl[cont.level], xlim=Xlim, ylim=Ylim) #edge
  }
}

pair_plot <- function(YourData, pixs=1, lab=NULL){ #like ipairs but no upper panel
  if ( is.null(lab) ){
    return(pairs(YourData, lower.panel=function(...) {par(new=TRUE); myImage(matrix(c(...),dim(YourData)[1],2), pixs=pixs, add = T)}, diag.panel=panel.hist2, upper.panel=NULL))
  } else {
    return(pairs(YourData, lower.panel=function(...) {par(new=TRUE); myImage(matrix(c(...),dim(YourData)[1],2), pixs=pixs, add = T)}, diag.panel=panel.hist2, upper.panel=NULL, labels = lab))
  }
}