#' Main Function for Clumping Habitat Categories in a Binary Matrix
#'
#' This function modifies a discrete matrix to increase aggregation of SWF habitat categories, 
#' starting from an initial allocation of habitat (specified by swfCat) and non-habitat 
#' (specified by agriCat). The function follows a growing approach, where habitat seeds 
#' are clumped in two dimensions, while keeping the total habitat size constant.
#'
#' @param Hmatrix A matrix representing the initial habitat state.
#' @param swfCat Integer representing habitat category in the matrix.
#' @param agriCat Integer representing non-habitat category in the matrix.
#' @param foreCat Integer representing fores (habitat) category in the matrix.
#' @param Q Integer representing the number of cells to be moved in each kernel per iteration.
#' @param iterations Number of iterations for the aggregation process.
#' @param kernelCl Vertical size of the kernel for processing.
#' @param kernelRw Horizontal size of the kernel for processing.
#' @param NNeighbors Threshold of neighbours in a kernel below which habitat pixels are moved in a kernel.
#' @param maxDistance Searching radius for habitat pixel aggregation.
#' @param Density Strategy for selecting the cell as a gravity center. "H" to choose area with high habitat cover, "L" from low, "M" for median density.
#' @param queensCase If TRUE, all 8 directions are considered for neighbors; if FALSE, only orthogonal neighbors are considered.
#' @param np Number of cores for parallel processing.
#' @param deBug If TRUE, debugging information is provided.

swf.clumper.regular <- function(Hmatrix = null.mt, swfCat, agriCat, foreCat=3, Q, sigma, iterations = 20, kernelCl=20, kernelRw=20, NNeighbors=0, maxDistance = 5, Density="H", queensCase=FALSE, np=1, deBug=FALSE) {
	matrices.list <- list()
	ncolumns <- ncol(Hmatrix)
	nrows <- nrow(Hmatrix)
	nkernels = (kernelCl*kernelRw) / (kernelCl+kernelRw)
	counter <- 0
	numCores <- min(np, length(seq(1, ncolumns, kernelCl)))

	for (iteration in 1:iterations) {

		# This defines a matrix made only of swfCat and agriCat so to locate edges of interest
		Bmatrix <- Hmatrix
		Bmatrix[Bmatrix%in%c(foreCat)]<-swfCat
		Bmatrix[!Bmatrix%in%c(agriCat, swfCat)] <- swfCat  
		Hedges <- as.matrix(imager::cannyEdges(imager::as.cimg(Bmatrix), sigma=sigma))
		# image(as.matrix(imager::cannyEdges(as.cimg(as.cimg(Bmatrix)),sigma=2)))

		temp.cl <- mclapply(seq(1, ncolumns, kernelCl), function(cl) {
			cl_end = min(cl + kernelCl - 1, ncolumns) # if windows exceed columns then columns

			temp.rw <- lapply(seq(1, nrows, kernelRw), function(rw) {
				rw_end = min(rw + kernelRw - 1, nrows)
				
				this.window = Hmatrix[rw:rw_end, cl:cl_end]
				this.edges = Hedges[rw:rw_end, cl:cl_end]

				# Finds coords with agriCat
				agri.pos = as.data.frame(which(this.window == agriCat, arr.ind = TRUE))

				# Finds coords with swfCat and on edges: This approach first uses canny edges detector, then selects only habitat cells
				swf.ed <- as.data.frame(which(this.window == swfCat & this.edges, arr.ind = TRUE))
				# Finds coords with swfCat and not on edges: This approach first uses canny edges detector, then selects only habitat cells
				swf.NOed <- as.data.frame(which(this.window == swfCat & !this.edges, arr.ind = TRUE))

				# if( nrow(swf.ed) > 0) {
				# swf.ed = as.data.frame(which(this.window == swfCat, arr.ind = TRUE))
				# Reshaffle only if there is at least one 1 and
				# Reshaffle only if there is at least one non-neighbour 1 (maxDistance is the search radius around each 1)
				if ( (nrow(swf.ed) != 0 && (nrow(swf.ed) > 1 && nrow(agri.pos) > 1)) && areThereLonely(this.window, NNeighbors, queensCase, maxDistance) ) {
					if (deBug) cat(paste("Col: ",cl, " Row: ",rw," Nrow swf.ed: ", nrow(swf.ed),"\n"))
					# Chose a gravity cell which is in a high/low density part of the matrix, this is done using findCell which counts the neighborn in a radius for each cell and chose the one with most/least neighbourns
					gravity.pos <- findCell(swf.ed, ifelse(kernelCl>kernelRw,kernelCl,kernelRw)/2, Density)
					# Orders the 1's coordinates based on distance from gravity center 
					swf.NOed$dist = as.numeric(proxy::dist(gravity.pos, swf.NOed))
					swf.NOed = swf.NOed[order(swf.NOed$dist), 1:2]
					# Orders the 0's coordinates based on distance from gravity center 
					# agri.pos$dist = as.numeric(proxy::dist(gravity.pos, agri.pos))

					# Draw coins to chose a 0 cell, coins are biased with p equal to an exponential distance decay
					# roc <- rowSums(rmultinom(Q, 1:length(agri.pos$dist), 1/agri.pos$dist^3))
					# agri.fate = agri.pos[which(roc>=1), 1:2]
					# Otherwise draw a coin to chose between the closest two 0's
					# agri.fate <- head(agri.pos[order(agri.pos$dist), 1:2], Q)
					# Find the closest Q number of neighbourns, but give priority to diagonal neigbourns (to avoid straight line)
					agri.fate <- ClosestDiagonalAgriCell(this.window, agri.pos, gravity.pos, swfCat, agriCat, Q)
					# Decide which SWF is going to be moved, now the nth farthest from gravity center and not an edge
					particle.pos = tail(swf.NOed, nrow(agri.fate))

					# Remove the 1 that moves
					this.window[as.matrix(particle.pos)] <- agriCat
					# Add it to the new position
					this.window[as.matrix(agri.fate)] <- swfCat
					# } else{counter=counter+1; return(this.window) }
				} else{counter=counter+1; return(this.window) }
				if(counter>nkernels) stop(paste("Nothing more to move, stopping at iteration: ", iter))
				return(this.window)
				})

			do.call(rbind, temp.rw)
			}, mc.cores=numCores )

		outmatrix <- do.call(cbind, temp.cl)
		matrices.list[[iteration]] <- outmatrix
		Hmatrix <- outmatrix
		message(paste("Iteration: ", iteration))
	}

	return(matrices.list)
}

# Helper functions
#' Check if a Specific Number of Neighbors Exist for Habitat Cells
#'
#' This function checks if each habitat cell (specified by swfCat) in a matrix has a specific 
#' number of neighboring habitat cells within a given distance.
#'
#' @param matrix A binary matrix to check.
#' @param NNeighbors Exact number of neighbors to check for each habitat cell.
#' @param queensCase If TRUE, checks all 8 directions for neighbors; if FALSE, checks only 4.
#' @param maxDistance Maximum distance to look for neighbors.
#' @return TRUE if any habitat cell has the specified number of neighbors within the given distance, FALSE otherwise.

areThereLonely <- function(matrix, NNeighbors, queensCase = TRUE, maxDistance = 5) {
	rows <- nrow(matrix)
	cols <- ncol(matrix)

  # Function to count neighbors
  countNeighbors <- function(r, c) {
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
  			if (i >= 1 && i <= rows && j >= 1 && j <= cols && matrix[i, j] == 1) {
  				neighborCount <- neighborCount + 1
  			}
  		}
  		return(neighborCount)
  	}

  # Iterate through the matrix
  for (r in 1:rows) {
  	for (c in 1:cols) {
  		if (matrix[r, c] == 1) {
  			neighbors <- countNeighbors(r, c)
  			if (neighbors == NNeighbors) {
  				return(TRUE)
  			}
  		}
  	}
  }
  return(FALSE)
}

#' Calculate Euclidean Distance Between Two Points
#'
#' This function calculates the Euclidean distance between two points, typically in a matrix.
#'
#' @param x1 X-coordinate of the first point.
#' @param y1 Y-coordinate of the first point.
#' @param x2 X-coordinate of the second point.
#' @param y2 Y-coordinate of the second point.
#' @return Euclidean distance between the two points.

distance <- function(x1, y1, x2, y2) {
	sqrt((x2 - x1)^2 + (y2 - y1)^2)
}

#' Find a Cell Based on Local Density
#'
#' This function finds a cell in a matrix based on the local density of a specified category 
#' (habitat or non-habitat). It can be set to find cells in high, low, or median density areas.
#'
#' @param cells A matrix of cells to search from.
#' @param radius Radius within which to calculate local density.
#' @param Density Density selection mode: "H" for high, "L" for low, "M" for median density.
#' @return Coordinates of the chosen cell based on the specified density.

findCell <- function(cells, radius=5, Density="M", deBug=FALSE) {
    localDensities <- sapply(1:nrow(cells), function(i) countNeighbors(cells[i,], cells, radius))

    chosen <- switch(Density,
        "H" = which(localDensities == max(localDensities)),
        "L" = which(localDensities == min(localDensities)),
        "M" = which.min(abs(localDensities - median(localDensities)))  # Closest to median
    )

    # Handle case where 'chosen' is empty or multiple indices
    if (length(chosen) == 0) {
        if (deBug) cat("No cell found for Density:", Density, "\n")
        return(NULL)
    } else if (length(chosen) > 1) {
        chosen <- sample(chosen, 1)
    }

    Cell <- cells[chosen, ]
    if (deBug) cat("Chosen cell for Density:", Density, " - ", Cell, "\n")
    return(Cell)
}

#' Count Neighbors of a Cell Within a Radius
#'
#' Counts the number of neighbors a given cell has within a specified radius, 
#' useful for assessing habitat connectivity.
#'
#' @param cell The cell for which to count neighbors.
#' @param cells The matrix of cells to check against.
#' @param radius The radius within which to count neighbors.
#' @return The count of neighbors within the specified radius.

countNeighbors <- function(cell, cells, radius) {
	count <- 0
	for (i in 1:nrow(cells)) {
		if (distance(cell$row, cell$col, cells$row[i], cells$col[i]) <= radius && !(cell$row == cells$row[i] && cell$col == cells$col[i])) {
			count <- count + 1
		}
	}
	return(count)
}
#' Resample Elements from a Vector
#'
#' Resamples elements from a given vector, typically used for selecting random elements.
#'
#' @param x Vector to resample from.
#' @return A resampled element from the vector.
resample <- function(x, ...) x[sample.int(length(x), ...)]

# Find the Q closest neighbourns
ClosestDiagonalAgriCell <- function(matrix, agriCatList, targetPos, swfCat, agriCat, numCells) {
    nrows <- nrow(matrix)
    ncols <- ncol(matrix)
    
    r <- as.numeric(targetPos[1])
    c <- as.numeric(targetPos[2])

    # Initialize variables
    selectedNeighbors <- data.frame(row = integer(), col = integer())
    radius <- 1

    while(nrow(selectedNeighbors) < numCells && radius <= max(nrows, ncols)) {
        # Define neighbor positions within the current radius
        rowRange <- max(1, r-radius):min(nrows, r+radius)
        colRange <- max(1, c-radius):min(ncols, c+radius)
        
        # Find agriCat cells within the current radius
        inRadius <- agriCatList$row %in% rowRange & agriCatList$col %in% colRange
        validNeighbors <- agriCatList[inRadius, ]

        # Combine with previously selected neighbors and remove duplicates
        selectedNeighbors <- unique(rbind(selectedNeighbors, validNeighbors))

        radius <- radius + 1
    }

    # Order selected neighbors by distance
    orderedNeighbors <- orderCellsByDistance(selectedNeighbors, r, c)
    orderedNeighbors <- orderedNeighbors[1:min(nrow(orderedNeighbors), numCells), ]

    return(orderedNeighbors)
}

orderCellsByDistance <- function(cellList, targetRow, targetCol) {
    if (nrow(cellList) == 0) {
        return(cellList)
    }
    distances <- sqrt((cellList$row - targetRow)^2 + (cellList$col - targetCol)^2)
    orderedIndices <- order(as.matrix(distances))
    return(cellList[orderedIndices, , drop = FALSE])
}

# safeCannyEdges <- function(window) {
#   tryCatch({
#     # Check for homogeneity to avoid unnecessary edge detection
#     if (all(window == window[1,1])) {
#       return(matrix(FALSE, nrow = nrow(window), ncol = ncol(window)))
#     } else {
#       edges <- imager::cannyEdges(as.cimg(window))
#       return(as.matrix(edges) > 0)
#     }
#   }, error = function(e) {
#     # In case of error, return a matrix of FALSE (no edges detected)
#     return(matrix(FALSE, nrow = nrow(window), ncol = ncol(window)))
#   })
# }