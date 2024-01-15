#' SWF+CLC2018
#'
#' A \code{SpatRaster} (EPSG: 4326).
#' 
#' This dataset provides a ?
#'
#' @format A \code{SpatRaster} containing the following elements:
#' \describe{
#'   \item{CLC}{SWF+CLC}
#' }
#' @usage load_swf_clc()
#' @source \url{https://land.copernicus.eu/global/products/clc}
#' @references \url{https://land.copernicus.eu/global/products/clc}
#' @name swf_clc
#' @docType data
#' @keywords datasets
#' @examples
#' clc_swf <- terra::unwrap(readRDS(system.file("extdata", "clc_swf.rds", package = "SWFmolder")))
NULL