
# Plot an arbitrary number of 2D contours. The first function is used as the background also.
ggplotMultipte2DContour <- function(..., xlim = c(0, 2*pi), ylim = c(0, 5), res = 100) {

  # Obtain a list of functions to plot.
  funlist <- list(...)

  nfun <- length(funlist)

  # Save the names of the arguments to use in the legend.
  thisCall <- as.character(match.call())
  funNames <- thisCall[2:(1 + nfun)]

  names(funlist) <- funNames

  # Grid of evaluated values.
  xyGrid <- expand.grid(x = seq(xlim[1], xlim[2], length.out = res),
                        y = seq(ylim[1], ylim[2], length.out = res))

  # Height (z) of the function that are provided.
  zMatrix <- sapply(funlist, function(x) apply(xyGrid, 1, x))

  # Combined matrix.
  xyzMatrix <- cbind(xyGrid, zMatrix)

  # Melt to provide ggplot a long format.
  moltenGrid <- melt(xyzMatrix, id.vars = c("x", "y"))

  # Return the plot
  p <- ggplot(moltenGrid, mapping = aes(x = x, y = y)) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    geom_raster(data = xyzMatrix, mapping = aes(x = x, y = y, fill = zMatrix[, 1])) +
    geom_contour(aes(z = value, group = variable, col = variable)) + theme_bw()

  p
}


plotIS <- function(posteriorKernel, importancePdf, xlim, ylim, res = 100) {

  # Grid of evaluated values.
  testGrid <- expand.grid(mu = seq(xlim[1], xlim[2], length.out = res),
                          kp = seq(ylim[1], ylim[2], length.out = res))

  # Height of the contours for the target and the IS density.
  posteriorProb  <- apply(testGrid, 1, posteriorKernel)
  importanceProb <- apply(testGrid, 1, importancePdf)


  # Dataset for the plot.
  gridDat <- cbind(testGrid, posteriorProb, importanceProb)

  ggplot(gridDat, mapping = aes(x = mu, y = kp)) +
    lims(x = xlim, y = ylim) +
    geom_raster(aes(fill = posteriorProb)) +
    # geom_point(data = data.frame(postSample), aes(x = mu, y = kp), inherit.aes = FALSE) +
    geom_contour(aes(z = importanceProb)) +
    geom_contour(aes(z = posteriorProb), col = "tomato")
}
