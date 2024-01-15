#' SWF Molder Function
#'
#' This function modifies a habitat matrix to increase the cover of a specified category
#' (swfCat) by clumping it within a defined kernel size, considering neighbor preferences
#' and optionally reducing the selection by a given factor.
#'
#' @param Hmatrix A matrix representing the initial habitat state.
#' @param swfCover The desired cover proportion for the swfCat category.
#' @param swfCat The category value in the matrix to be increased.
#' @param agriCat The category value in the matrix considered as non-habitat.
#' @param Q The number of cells to be potentially modified in each iteration.
#' @param fDiagonal Function to determine the prioritization of diagonal neighbors.
#' @param reduceTo A factor by which to reduce the selection of cells during processing.
#' @param iterations The number of iterations to run the modification process.
#' @param kernelCl The number of columns in the kernel (clumping window).
#' @param kernelRw The number of rows in the kernel (clumping window).
#' @param NNeighbors The number of neighbor cells to consider for potential habitat clumping.
#' @param maxDistance The maximum distance to look for neighbor cells.
#' @param queensCase Logical; if TRUE, considers all 8 directions for neighbors; if FALSE, only orthogonal.
#' @param np Number of cores for parallel processing.
#' @param deBug Logical; if TRUE, prints debugging information during processing.
#'
#' @return A list of matrices representing the state of the habitat matrix after each iteration.
#' @export
#' @examples
#' initial_matrix <- matrix(c(1,2,1,3,1,2,2,3), ncol=2)
#' result <- swf.molder(Hmatrix = initial_matrix, swfCover=0.75, swfCat=2, 
#'                      agriCat=1, Q=1, fDiagonal=TRUE, reduceTo=0, iterations=20, 
#'                      kernelCl=2, kernelRw=2, NNeighbors=0, maxDistance=1, 
#'                      queensCase=FALSE, np=1, deBug=FALSE)

swf.molder <- function(Hmatrix, swfCover=0.10, swfCat, agriCat, Q, fDiagonal, reduceTo=0, iterations = 20, kernelCl=20, kernelRw=20, NNeighbors=0, maxDistance = 1, queensCase=FALSE, np=1, deBug=FALSE) {
matrices.list <- list()
ncolumns <- ncol(Hmatrix)
nrows <- nrow(Hmatrix)
nkernels = (kernelCl*kernelRw) / (kernelCl+kernelRw)
counter <- 0
iteration = 0

SWFarea = length(which(Hmatrix%in%c(swfCat))) / length(which(Hmatrix%in%c(agriCat,swfCat)))

while(iteration < iterations && SWFarea<swfCover) {
iteration=iteration+1

temp.cl <- lapply(seq(1, ncolumns, kernelCl), function(cl) {
cl_end = min(cl + kernelCl - 1, ncolumns)

temp.rw <- lapply(seq(1, nrows, kernelRw), function(rw) {
rw_end = min(rw + kernelRw - 1, nrows)

this.window = Hmatrix[rw:rw_end, cl:cl_end]

# Finds coords with agriCat
agri.pos = as.data.frame(which(this.window == agriCat, arr.ind = TRUE))

# Finds coords with swfCat
swf.ed <- as.data.frame(which(this.window == swfCat, arr.ind = TRUE))

# Adds swf only if there is at least one swfCat pixel in the kernel
if ( (nrow(swf.ed) != 0 ) ) {
if (deBug) cat(paste("Col:",cl, "; Row:",rw,"; SWF pixel: ", nrow(swf.ed),"\n"))
# Chose a gravity cell
gravity.pos <- findSwfCells(this.window, swfCat, agriCat, NNeighbors, TRUE, maxDistance, np)

if(nrow(gravity.pos)>0) {

if( reduceTo!=0 ) {
sampledGravitypos <- floor(sample(1:nrow(gravity.pos)*reduceTo))
if(length(sampledGravitypos)>0) gravity.pos <- gravity.pos[sampledGravitypos,]
}

# Finds the closest Q number of neighbourns
agri.fate <- do.call(rbind.data.frame, ClosestDiagonalAgriCell(this.window, agri.pos, gravity.pos, swfCat, agriCat, fDiagonal, Q, np=np))

# Adds it to the new position
this.window[as.matrix(agri.fate)] <- swfCat
} else{counter=counter+1; iteration=iteration+1; return(this.window) }
} else{counter=counter+1; iteration=iteration+1; return(this.window) }
if(counter>nkernels) stop(paste("Nothing more to move, stopping at iteration: ", iter))
return(this.window)
})

do.call(rbind, temp.rw)
} )

outmatrix <- do.call(cbind, temp.cl)
matrices.list[[iteration]] <- outmatrix
Hmatrix <- outmatrix
SWFarea = length(which(Hmatrix%in%c(swfCat))) / length(which(Hmatrix%in%c(agriCat,swfCat)))
message(paste("Iteration: ", iteration, "; SWF cover:", round(SWFarea,2)))

}

return(matrices.list)
}