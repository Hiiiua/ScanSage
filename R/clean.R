#' Get binarized cimg
#'
#' @description
#' Convert the x3p and its mask to a cimg
#'
#' @param x3p_path A path to x3p
#' @param mask_path A path to the mask corresponding to the x3p file
#' @param downsample_m The downsample size
#' @param color color of the cimg, default to "W", which is non-na as white
#'
#' @importFrom terra as.raster
#' @importFrom x3ptools x3p_read x3p_add_mask x3p_sample
#' @importFrom png readPNG
#' @importFrom imager imfill
#'
#' @export
boundary_cimg <- function (x3p_path, mask_path, downsample_m, color="W"){
  land <- x3p_read(x3p_path)
  mask <- png::readPNG(mask_path)
  overlay <- x3p_add_mask(land, mask = terra::as.raster(mask))
  down0 <- overlay %>% x3p_sample(m=20)

  if (color=="W"){
    boundaryW = imfill(x=down0$header.info$sizeX, y=down0$header.info$sizeY, z=1, val=1)
    boundaryW[is.na(down0$surface.matrix)]=0
    return(boundaryW)
  }
  else{
    boundaryB = imfill(x=down0$header.info$sizeX, y=down0$header.info$sizeY, z=1, val=0)
    boundaryB[is.na(down0$surface.matrix)]=1
    return(boundaryB)
  }
}

#' Find highlighted regions
#'
#' @description
#' A helper function to determine the number of pixels in highlight regions
#'
#' @param bound_cimg A cimg from boundary_cimg output
#' @param maxY maximum???
#' @param px number of pixels
#'
#' @importFrom imager highlight fill
#' @export
outRegion <- function(bound_cimg, maxY=1.5, px = NA){
  numElement <- function(element){
    max(lengths(element))
  }

  if (is.na(px) == TRUE){ # automatically find the lines
    px = 1
    find = 0
    while (find == 0 & px <=5) {
      highlight = bound_cimg %>% fill(px) %>% highlight
      numHighlights = sapply(highlight, numElement)
      # The top 2 number of pixels of highlight regions
      numMax = sort(numHighlights, decreasing = TRUE)[1:3]

      if (sum(numMax > dim(bound_cimg)[2])==2){# top 2 has all y axis covered
        if (sum(numMax[1:2] <= maxY*dim(bound_cimg)[2]) <= 2){
          ms = order(numHighlights, decreasing = TRUE)[1:2]
          find = 1
          h1 = highlight[[ms[1]]]
          h2 = highlight[[ms[2]]]
        }
        else{
          find = 0}
      }
      else{
        find = 0
      }
      px = px+1
    }
  }

  else{
    find=1
    highlight = bound_cimg %>% fill(px) %>% highlight
    numHighlights = sapply(highlight, numElement)
    # The top 2 number of pixels of highlight regions
    numMax = sort(numHighlights, decreasing = TRUE)[1:3]
    ms = order(numHighlights, decreasing = TRUE)[1:2]
    h1 = highlight[[ms[1]]]
    h2 = highlight[[ms[2]]]
  }

  if(find==1){
    plot(bound_cimg)
    with(h1, lines(x,y, col='red', lwd=2))
    with(h2, lines(x,y, col='blue', lwd=2))
    return(list(px=px, h1=data.frame(h1$x, h1$y), h2=data.frame(h2$x, h2$y)))}
  else{
    return("Not available for px<=5")
  }
}


#' Create Alpha hull to identify grooves
#'
#'
#' @param boundaryW boundary_cimg in white mode
#' @param outRegionBoundary one set of the boundary points from outRegion function
#' @param alpha alpha used for generating alpha hull
#' @param pointsCol color to display boundary points
#' @param hullCol color to display the hull
#'
#' @importFrom alphahull ahull
#' @importFrom graphics points
#' @importFrom stats runif
#'
#' @export

boundaryAlphaHull <- function(boundaryW, outRegionBoundary, alpha=10,
                              pointsCol="red", hullCol="blue"){
  blue = outRegionBoundary

  # find the max and min of x values in the points set, so we know which boundary
  # and how to include corner points
  x = outRegionBoundary[,1]
  xmax = max(x)
  if (xmax < 0.5*dim(boundaryW)[1]){ # left boundary
    inside = data.frame(x=sample.int(round(min(blue[,1])), 30, replace = TRUE),
                        y = sample.int(round(max(blue[,2])), 30, replace = TRUE))
    colnames(inside) = colnames(blue)
    blue = rbind(blue, inside)
    blue = rbind(c(1, 0), blue, c(1, dim(boundaryW)[2])) # add corner points, 1 for visual purpose
  }
  else{ # right boundary, normally plot as red
    inside = data.frame(x=sample(seq(round(max(blue[,1])), dim(boundaryW)[1]), 30, replace = TRUE),
                        y = sample.int(round(max(blue[,2])), 30, replace = TRUE))
    colnames(inside) = colnames(blue)
    blue = rbind(blue, inside)
    blue = rbind(c(xmax, 0), blue, c(xmax, dim(boundaryW)[2]))
  }

  # jitter for alpha hull function
  jitterBlue = blue
  jitterBlue[,1] = blue[,1]+runif(dim(blue)[1], 0, 0.05)
  jitterBlue[,2] = blue[,2]+runif(dim(blue)[1], 0, 0.05)

  bluealphahull <- ahull(jitterBlue, alpha = alpha)
  ahullPoints = bluealphahull$arcs[,"end1"]
  print(paste("There are", length(ahullPoints), "alpha hull points, reduced ", dim(outRegionBoundary)[1],"Reduced ratio is ",
              length(ahullPoints)/dim(outRegionBoundary)[1])) # print how many points are taken
  plot(boundaryW)
  points(jitterBlue, col='green', pch=19, cex=0.5)
  plot(bluealphahull, add=TRUE, col='red', pch="1", cex=0.1)

  res = list(bluealphahull, bluealphahull$xahull[ahullPoints,])
  names(res) = c("ahull.obj", "points")
  return(res)
}
