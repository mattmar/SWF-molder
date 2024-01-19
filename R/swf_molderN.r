#' @title SWF Molder Network Function
#' @description Attemps tomodifies a habitat graph to increase the cover of a specified category
#' @param g Graph object representing the habitat.
#' @param swfCat Category value in the graph to be increased.
#' @param agriCat Category value in the graph to allocate.
#' @param iterations Number of iterations to run the modification process.
#' @param Q Number of nodes to be potentially modified in each iteration.
#' @param NNeighbors Number of agriCat neighbor nodes to consider for allocation.
#' @param swfCover Desired cover proportion for the swfCat category.
#' @param max_radius Maximum search radius for neighbours.
#' @param np Number of cores for parallel processing.
#' @export

swf_molderN <- function(g, swfCat, agriCat, iterations, Q, NNeighbors = 1, swfCover=0.3, max_radius=1, np=1) {
	iteration=0
	SWFarea = length(which(V(g)$value == swfCat))/length(which(V(g)$value %in% c(swfCat,agriCat)))
	graph.list <- list()
	while( iteration < iterations && SWFarea < swfCover ) {
# Identify nodes for potential modification
# pippo<-1
nNN <- integer(0)
radius=1
#Used for in the past, while seems to be more appropriat for growing searches
while( length(nNN) == 0 && radius <= max_radius ) {
	AgriV <- which(V(g)$value==swfCat)
	nNN <- unlist(parallel::mclapply(V(g)[AgriV], check_neighbors, g, agriCat, NNeighbors, radius, mc.cores = np))
	if(length(nNN) == 0) {
		radius <- radius + 1
	}
}
if( length(nNN)>0 ) {
nodes_to_allocate <- unlist(lapply(V(g)[nNN], function(nei) {
	nei_nodes <- neighbors(g, nei)
agri_nodes <- nei_nodes[V(g)[nei_nodes]$value == agriCat]
if (length(agri_nodes) >= Q) {
	return(sample(agri_nodes, Q))
	} else {
		return(agri_nodes)
	}
	}))
}

nodes_to_allocate <- nodes_to_allocate[!V(g)[nodes_to_allocate]$value %in% c(swfCat, 3)]

V(g)[nodes_to_allocate]$value <- swfCat
iteration = iteration+1
graph.list[[iteration]] <- g
SWFarea = length(which(V(g)$value == swfCat))/length(which(V(g)$value %in% c(swfCat,agriCat)))
message(paste("Iteration: ", iteration, "; SWF cover:", round(SWFarea,2)))
}
if(length(nNN)>0 && length(graph.list)>2) {
	if( identical(g, graph.list[[length(graph.list)-1]]) ) {
		break
	}
}
return(graph.list)
}

#' @title Check Neighbours Function
#' @description Checksfor the number of agriCat neighbors for each node within a specified radius.
#' @param v Node for check neighbors.
#' @param graph Graph object representing the habitat.
#' @param agriCat Category value in the graph to allocate.
#' @param NNeighbors Number of neighbors to check for each AgriCat node.
#' @param radius Radius within which to check for neighbors.
#' @return Node if it meets the NNeighbors condition.
#' @export
check_neighbors <- function(v, graph, agriCat, NNeighbors, radius) {
	nbrs <- neighborhood(graph, order = radius, nodes = v, mode = "all")[[1]]
	agri_neighbors <- sum(V(graph)[nbrs]$value == agriCat)
	if(agri_neighbors >= NNeighbors) {
		return(v)
	}
}

#' @title Graph from Matrix Function
#' @description Attempts to create a graph object from a matrix, with nodes representing cells and edges representing adjacency. I've taken this partially from StackExchange but it seems very unefficient 
#' @param mat Matrix representing the initial habitat state.
#' @param np Number of cores for parallel processing.
#' @return Graph object corresponding to the input matrix.
#' @export
GfM <- function(mat, np = 1) {
	nrows <- nrow(mat)
	ncols <- ncol(mat)
	nodes <- expand.grid(row = 1:nrows, col = 1:ncols)
	node_values <- as.vector(mat)

	get_edges <- function(node) {
		r <- node$row
		c <- node$col
		node_id <- (r - 1) * ncols + c

# Define the potential neighbor positions 
# (including diagonals)? I think it does but check
neighbors <- c(-1, 0, 1)
edges <- c()

for (dr in neighbors) {
	for (dc in neighbors) {
# Skip the current cell to avoid self-loops, hard ti solve for now works
if (dr == 0 && dc == 0) next
nr <- r + dr
nc <- c + dc
# Check if neighbor is within bounds
if (nr >= 1 && nr <= nrows && nc >= 1 && nc <= ncols) {
	neighbor_id <- (nr - 1) * ncols + nc
	edges <- c(edges, c(node_id, neighbor_id))
}
}
}
return(edges)
}

# Use parallel processing to get edges for all nodes (unefficient sometumes)
edge_list <- parallel::mclapply(seq_len(nrow(nodes)), function(i) get_edges(nodes[i,]), mc.cores = np)

# Combine all edges into a single vector
all_edges <- unlist(edge_list, recursive = FALSE)

# Create the graph from the edge list
g <- graph(edges = all_edges, directed = FALSE)
V(g)$value <- node_values

# Remove self-loops and simplify the graph (again?)
g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)

return(g)
}


#' @title Matrix from Graph Function
#' @description Creates a matrix from a graph object, with cells representing nodes.
#' @param g Graph object to convert into a matrix.
#' @param nrows Number of rows for the resulting matrix.
#' @param ncols Number of columns for the resulting matrix.
#' @export
MfG <- function(g, nrows, ncols) {
# Create an empty matrix of defined dimensions
mat <- matrix(0, nrow = nrows, ncol = ncols)
# Iterate through the nodes and set matrix values based on graph edges
# Could be done with lapply to improve speed?
for (node_id in 1:vcount(g)) {
	r <- ((node_id - 1) %/% ncols) + 1
	c <- ((node_id - 1) %% ncols) + 1
	mat[r,c] = V(g)[node_id]$value
}
return(mat)
}