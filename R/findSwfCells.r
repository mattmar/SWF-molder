# Helper functions
#' @export
findSwfCells <- function(matrix, swfCat, agriCat, AgriNeighbors, queensCase = TRUE, maxDistance = 1, fDiagonal=FALSE, np = 1) {
	rows <- nrow(matrix)
	cols <- ncol(matrix)

  # Initialize a matrix to store the results
  resultMatrix <- matrix(FALSE, nrow = rows, ncol = cols)

  # Function to count agriCat neighbors
  countAgriNeighbors <- function(r, c, matrix, agriCat, queensCase, maxDistance) {
  	neighborCount <- 0

  	if (queensCase) {
  		directions <- expand.grid(row = -maxDistance:maxDistance, col = -maxDistance:maxDistance)
  		directions <- subset(directions, !(row == 0 & col == 0))
  		} else {
  			directions <- list()
  			for (dist in 1:maxDistance) {
  				directions <- c(directions, list(c(-dist, 0)), list(c(dist, 0)), 
  					list(c(0, -dist)), list(c(0, dist)))
  			}
  			directions <- do.call(rbind, directions)
  			colnames(directions) <- c("row", "col")
  		}

  		for (dir in 1:nrow(directions)) {
  			i <- r + directions[dir, "row"]
  			j <- c + directions[dir, "col"]
  			if (i >= 1 && i <= rows && j >= 1 && j <= cols && matrix[i, j] == agriCat) {
  				neighborCount <- neighborCount + 1
  			}
  		}
  		return(neighborCount)
  	}

# Create a grid of row and column indices
indices <- expand.grid(row = 1:rows, col = 1:cols)

# Define a function to check each cell
checkCell <- function(index, matrix, swfCat, agriCat, AgriNeighbors) {
	r <- index$row
	c <- index$col
	if (matrix[r, c] == swfCat) {
		agriNeighbors <- countAgriNeighbors(r, c, matrix, agriCat, queensCase, maxDistance)
      if (agriNeighbors >= AgriNeighbors) {  # Change here if you want "at least" AgriNeighbors
      return(list(row = r, col = c))
    }
  }
  return(NULL)
}

# Use mclapply for parallel processing
resultList <- parallel::mclapply(seq_len(nrow(indices)), function(i) checkCell(indices[i, ], matrix, swfCat, agriCat, AgriNeighbors), mc.cores = np)

# Filter out NULL results and combine
resultCoords <- if (length(Filter(Negate(is.null), resultList)) > 0) {
	as.matrix(do.call(rbind, Filter(Negate(is.null), resultList)))
	} else {
    matrix(numeric(0), ncol = 2)  # An empty matrix with 2 columns
  }

# Initialize resultMatrix
resultMatrix <- matrix(FALSE, nrow = rows, ncol = cols)

# Populate resultMatrix
if (nrow(resultCoords) > 0) {
	resultMatrix[cbind(as.numeric(resultCoords[,1]), as.numeric(resultCoords[,2]))] <- TRUE
}

return(which(resultMatrix,arr.ind=T))
}
