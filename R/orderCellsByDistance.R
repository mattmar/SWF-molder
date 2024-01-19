orderCellsByDistance <- function(cellList, targetRow, targetCol, ExpPriority, ExpDirection) {
    if (nrow(cellList) == 0) {
        return(cellList)
    }

    # Calculate Manhattan distance for each cell from the target
    cellList$manhattanDist <- abs(cellList$row - targetRow) + abs(cellList$col - targetCol)
    # Adjust Manhattan distance for even values
    even_indices <- cellList$manhattanDist %% 2 == 0
    cellList$manhattanDist[even_indices] <- cellList$manhattanDist[even_indices] - 1

    # Identify diagonal cells
    cellList$isDiagonal <- abs(cellList$row - targetRow) == abs(cellList$col - targetCol)

    # Assign a random tiebreaker score
    cellList$tiebreaker <- runif(nrow(cellList))
    if (ExpPriority == "mixed") {
        if (ExpDirection == "orthogonal") {
            cellList$tiebreaker[cellList$isDiagonal] <- 0
        } else if (ExpDirection == "diagonal") {
            cellList$tiebreaker[!cellList$isDiagonal] <- 0
        } else if (ExpDirection == "mixed") {
            cellList$tiebreaker[!cellList$isDiagonal] <- sample(c(0,1),1)
        }
    } else if (ExpPriority == "vertical") {
        cellList$tiebreaker <- 0
    }

    # Decide the ordering direction
    if (ExpDirection == "mixed") {
        ExpDirection <- sample(c("orthogonal", "diagonal"), 1)
    }

    # Order by Manhattan distance, diagonal, and tiebreaker
    if (ExpDirection == "diagonal") {
        cellList <- cellList[order(cellList$manhattanDist, -cellList$isDiagonal, cellList$tiebreaker), 1:2]
    } else if (ExpDirection == "orthogonal") {
        cellList <- cellList[order(cellList$manhattanDist, cellList$isDiagonal, cellList$tiebreaker), 1:2]
    }

    return(cellList)
}
