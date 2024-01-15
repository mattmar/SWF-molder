#' @export
orderCellsByDistance <- function(cellList, targetRow, targetCol, fDiagonal) {
    if (nrow(cellList) == 0) {
        return(cellList)
    }

# Calculate Manhattan distance for each cell from the target
cellList$manhattanDist <- (abs(cellList$row - targetRow) + abs(cellList$col - targetCol))
# Report NN in the same order to the same score
cellList$manhattanDist[which(cellList$manhattanDist%%2==0)] <- cellList$manhattanDist[which(cellList$manhattanDist%%2==0)]-1

# Identify diagonal cells by checking if the absolute differences are the same
cellList$isDiagonal <- abs(cellList$row - targetRow) == abs(cellList$col - targetCol)

# Create a ranking score that gives a slight advantage to diagonal neighbors
# Diagonal cells will have their Manhattan distance decreased by 0.5
cellList$rankScore <- cellList$manhattanDist - 0.5 * cellList$isDiagonal

# Order by Manhattan distance first, then prioritize diagonals within each level
if(fDiagonal) {
    cellList <- cellList[order(cellList$manhattanDist, -cellList$isDiagonal),1:2]
    } else {
        cellList <- cellList[order(cellList$manhattanDist, cellList$isDiagonal),1:2]
    }

return(cellList)
}

