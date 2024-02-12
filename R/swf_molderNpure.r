#' @title SWF Molder Network Function
#' @description Attemps tomodifies a habitat graph to increase the cover of a specified category
#' @param Hmatrix Matrix object representing the habitat.
#' @param swfCat Category value in the graph to be increased.
#' @param agriCat Category value in the graph to allocate.
#' @param iterations Number of iterations to run the modification process.
#' @param Q Number of nodes to be potentially modified in each iteration.
#' @param NNeighbors Number of agriCat neighbor nodes to consider for allocation.
#' @param swfCover Desired cover proportion for the swfCat category.
#' @param max_radius Maximum search radius for neighbours.
#' @param np Number of cores for parallel processing.
#' @export

swf_molderN <- function(Hmatrix, swfCat, agriCat, iterations, Q, NNeighbors = 1, swfCover=0.3, max_radius=1, np=1) {
	
	iteration=1
	nr <- nrow(Hmatrix)
	nc <- ncol(Hmatrix)
	graph.list <- list()

	g <- make_lattice(c(nr, nc), nei=2)
	V(g)$value <- ''
	V(g)$id <- 1:(nr*nc)

	swfCoverObs = length(which(Hmatrix%in%swfCat)) / length(which(Hmatrix%in%c(agriCat,swfCat)))
	all_nei <- neighborhood(g, max_radius)

	while( iteration < iterations && swfCoverObs < swfCover ) {

		agriPx <- which(Hmatrix==agriCat)
		swfPx <- which(Hmatrix==swfCat)
		
		neighbors_list <- mclapply(all_nei[swfPx],intersect,agriPx,mc.cores=np)

		# If there are nodes with the required number of neighbors, allocate them
		if( length(neighbors_list) > 0 ) {
			nodes_to_allocate <- unlist(lapply(neighbors_list, function(agri_nodes) {
				if(length(agri_nodes) >= Q) {
					return(sample(agri_nodes, Q))
					} else {
						return(agri_nodes)
					}
					}))

			Hmatrix[nodes_to_allocate] <- swfCat
			} else break

			graph.list[[iteration]] <- Hmatrix

			swfCoverObs = length(which(Hmatrix%in%swfCat)) / length(which(Hmatrix%in%c(agriCat,swfCat)))
			if( swfCoverObs >= swfCover ) break

			iteration = iteration+1
			message(paste("Iteration: ", iteration, "; SWF cover:", round(swfCoverObs,2)))

			if( iteration>2 ) {
				if( identical(Hmatrix, graph.list[[iteration-2]]) ) {
					break
				}
			}
		}
		return(graph.list)
	}