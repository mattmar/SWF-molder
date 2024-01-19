# Helper functions
#' @export
findSwfCells <- function(matrix, swfCat, agriCat, AgriNeighbors, queensCase = TRUE, maxDistance = 1, fDiagonal = FALSE, np = 1) {
	rows <- nrow(matrix)
	cols <- ncol(matrix)

    # Generate directions for neighbor checking
    if (queensCase) {
    	directions <- expand.grid(row = -maxDistance:maxDistance, col = -maxDistance:maxDistance)
    	directions <- subset(directions, !(row == 0 & col == 0))
    	} else {
    		directions <- do.call(rbind, lapply(1:maxDistance, function(dist) {
    			rbind(c(-dist, 0), c(dist, 0), c(0, -dist), c(0, dist))
    			}))
    	}

    # Function to count agriCat neighbors
    countAgriNeighbors <- function(r, c, matrix, agriCat, directions) {
    	sum(sapply(1:nrow(directions), function(dir) {
    		i <- r + directions[dir, "row"]
    		j <- c + directions[dir, "col"]
    		i >= 1 && i <= rows && j >= 1 && j <= cols && matrix[i, j] == agriCat
    		}))
    }

    # Define a function to check each cell
    checkCell <- function(r, c, matrix, swfCat, agriCat, AgriNeighbors, directions) {
    	if (matrix[r, c] == swfCat) {
    		agriNeighbors <- countAgriNeighbors(r, c, matrix, agriCat, directions)
    		if (agriNeighbors >= AgriNeighbors) {
    			return(c(r, c))
    		}
    	}
    	NULL
    }

    # Create a matrix of row and column indices
    indices <- as.matrix(expand.grid(row = 1:rows, col = 1:cols))

    # Use parallel processing
    resultCoords <- do.call(rbind, parallel::mclapply(1:nrow(indices), function(i) {
    	checkCell(indices[i, "row"], indices[i, "col"], matrix, swfCat, agriCat, AgriNeighbors, directions)
    	}, mc.cores = np))

    # Convert to logical matrix
    resultMatrix <- matrix(FALSE, nrow = rows, ncol = cols)
    if ( !is.null(nrow(resultCoords)) && nrow(resultCoords) > 0) {
    	resultMatrix[cbind(resultCoords[, 1], resultCoords[, 2])] <- TRUE
    }

    which(resultMatrix, arr.ind = TRUE)
}
