# # Find the Q closest neighbourns
# test <- matrix(c(2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1), nrow=4)
# gravity.pos <- findSwfCells(test,swfCat=2,agriCat=1,AgriNeighbors=4,TRUE,maxDistance=1)
# agri.pos <- as.data.frame(which(test==1,arr.ind=T))
# do.call(rbind.data.frame, ClosestDiagonalAgriCell(test, agri.pos, gravity.pos, 1, 2, 2, np=1))

# ClosestDiagonalAgriCell(test, agri.pos, gravity.pos, 1, 2, 2, np=1)
#' @export

ClosestDiagonalAgriCell <- function(matrix, agriCatList, targetPosMatrix, swfCat, agriCat, fDiagonal, numCells, np) {
    nrows <- nrow(matrix)
    ncols <- ncol(matrix)

    parallel::mclapply( 1:nrow(targetPosMatrix), function(index) {
        targetPos <- targetPosMatrix[index, ]
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

# Order selected neighbors by (distance)
orderedNeighbors <- orderCellsByDistance(selectedNeighbors, r, c, fDiagonal)
allNeighbors <- orderedNeighbors[1:min(nrow(orderedNeighbors), numCells), ]
return(allNeighbors)
},mc.cores=np)

}