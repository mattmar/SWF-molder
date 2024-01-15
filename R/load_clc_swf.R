#' Load SWF+CLC data
#'
#' This function loads and returns the Copernicus SWF+CLC stored within the package.
#'
#' @return A `SpatRaster` object representing the SWF+CLC data.
#' @export

load_clc_swf <- function() {
  # Set the path to the RDS file
  rds_file <- system.file("extdata", "clc_swf.rds", package = "SWFmolder")
  
  # Check if the file exists
  if (file.exists(rds_file)) {
    # Load the SpatRaster object from the RDS file
    clc_swf <- terra::unwrap(readRDS(rds_file))
  } else {
    stop("Data file 'clc_swf.rds' not found in the package.")
  }
  
  return(clc_swf)
}
