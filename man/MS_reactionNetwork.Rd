\name{MS_reactionNetwork}
\alias{MS_reactionNetwork}
\title{Build reaction-compound network}
\description{
This function generates a directed reaction-compound network. The network is formalized
as a three-column matrix, where each row represents an edge connecting two nodes
(from source to target).
}
\usage{
MS_reactionNetwork(metabo_paths)
}
\arguments{
  \item{metabo_paths}{character vector containing the KEGG IDs of the metabolic
  pathways of interest. See functions "MS_keggFinder( )" and "MS_getPathIds( )".
  }
  }
\value{
A three-column matrix where each row represents an edge between two nodes.
}
\note{
Reaction directionality reported in KEGG has been cross-validated with published
literature (Duarte et al., 2007).
}
\examples{\donttest{
reaction_network <- MS_reactionNetwork(metabo_paths = c("rno00010", "rno00562"))
}
}
