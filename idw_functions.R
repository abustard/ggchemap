#' Grid irregular data and perform Inverse Distance Weighting Interpolation
#'
#' @param data A data.frame with columns x, y, and z with x and y being in projected map units, and z being the element to be interpolated
#' @param spacing The desired grid spacing for the resulting interpolated grid.
#' @param crop_result Logical. Should the resulting interpolation be cropped? If TRUE, the resulting grid will be cropped by the convex hull of the data points with a buffer of one grid spacing. If FALSE, returns the entire grid
#' @param idwp Power argument to be used by the inverse distance weighting algorithm. Defaults to 2
#'
idw_irregulargg <- function(data, spacing, crop_result = TRUE, idwp = 2){

  x.range <- range(data$x)
  y.range <- range(data$y)

  length_out <- max(c(x.range, y.range)/spacing)

  x<-seq(x.range[1], x.range[2], length.out= (x.range[2] - x.range[1])/spacing)
  y<-seq(y.range[1], y.range[2], length.out= (y.range[2] - y.range[1])/spacing)

  grd <- expand.grid(x,y)

  data_sp <- data
  coordinates(data_sp) = ~ x + y
  coordinates(grd) <- ~ Var1 + Var2
  gridded(grd) <- TRUE

  proj4string(data_sp) <- CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  proj4string(grd) <- CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  dat_idw <- gstat::idw(fill~ 1, locations = data_sp, newdata = grd, idp = idwp)

  if (crop_result == TRUE) {
    ch_data <- chull(data[, c("x", "y")])
    coords <- data[c(ch_data, ch_data[1]), ]
    crop_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords[,c("x", "y")])), ID=1)))
    crop_poly <- rgeos::gBuffer(crop_poly, width = spacing)
    proj4string(crop_poly) <- CRS("+proj=utm +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    dat_idw <- raster::crop(dat_idw, crop_poly)
  }

  dat_idw <- as.data.frame(dat_idw)

  data.frame(
    x = dat_idw$Var1,
    y = dat_idw$Var2,
    fill = dat_idw$var1.pred
  )
}

#'
#'
#'
StatIDW <- ggproto("StatIDW", Stat,
                   required_aes = c("x", "y"),
                   default_aes = aes(fill = ..fill..),

                   compute_group = function(data, scales,
                                            spacing = NULL,
                                            crop_result = NULL,
                                            idwp = NULL) {

                     print(spacing)
                     if (is.null(spacing)) {
                       spacing <- min(diff(range(data$x)), diff(range(data$y)))/250
                       warning(paste0("Grid size defaulting to", " ", spacing,
                                      "(minimum range/250), use 'spacing' to set value manually"), call. = FALSE)
                       print(spacing)
                     }

                     if (is.null(crop_result)) {
                       crop_result <- TRUE
                     }

                     if (is.null(idwp)) {
                       idwp <- 2
                       warning("Using idp value of 2", call. = FALSE)
                     }

                     idw_irregulargg(data, spacing, crop_result, idwp)

                   }
)

#'
#'
#'
#'
#'
stat_idw <- function(mapping = NULL, data = NULL, geom = "tile",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, ...) {
  layer(
    stat = StatIDW, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

