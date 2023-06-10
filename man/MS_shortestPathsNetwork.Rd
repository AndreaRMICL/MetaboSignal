\name{MS_shortestPathsNetwork}
\alias{MS_shortestPathsNetwork}
\title{Build shortest-path subnetwork}
\description{
This function allows calculating the shortest paths from a set of source nodes to
a set of target nodes, and representing them as a network. By default, the function
exports a network file and two attribute files ("NodesType.txt", "TargetNodes.txt"),
which can be imported into Cytoscape to visualize the network. The first attribute
file allows customizing the nodes of the network based on the molecular entity
they represent: signaling-gene, metabolic-gene, reaction or compound. The second
attribute file allows highlighting the source and target nodes.
}
\usage{
MS_shortestPathsNetwork(network_table, organism_code, source_nodes, target_nodes,
                        mode = "out", type = "first", distance_th = Inf, names = TRUE,
                        export_cytoscape = TRUE, file_name = "MS")
}
\arguments{
  \item{network_table}{three-column matrix where each row represents an edge
  between two nodes. See function "MS_keggNetwork( )".
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_keggFinder( )".
  }
  \item{source_nodes}{character vector containing the node IDs (typically genes)
  from which the shortest paths will be calculated. When using gene IDs make sure
  that they are consistent with the format of the network (i.e. organism-specific gene
  IDs or orthology IDs). Remember that Entrez IDs and gene symbols can be transformed
  into KEGG IDs with the function "MS_convertGene( )".
  }
  \item{target_nodes}{character vector containing the nodes IDs (typically compounds)
  to which the shortest paths will be calculated. Compound KEGG IDs can be obtained using
  the function "MS_keggFinder( )".
  }
  \item{mode}{character constant indicating whether a directed (mode = "out") or
  semi-directed (mode = "SP") network will be considered. "out" indicates that
  all the edges of the network will be considered as directed. "SP" indicates that
  all network will be considered as directed except the edges linked to target_node,
  which will be considered as undirected. The difference between the "out" and the
  "SP" options, is that the latter aids reaching target metabolites that are
  substrates of irreversible reactions.
  }
  \item{type}{character constant indicating whether all shortest paths or a
  single shortest path will be considered when there are several shortest paths
  between a source node and a target node. If type = "all", all shortest
  paths will be considered. If type = "first" a single path will be considered.
  If type = "bw" the path with the highest betweenness score will be considered.
  The betweenness score is calculated as the average betweenness of the gene
  nodes of the path. Note that using type = "bw" increases the time required to
  compute this function.
  }
  \item{distance_th}{establishes a shortest path length threshold. Only shortest
  paths with length below this threshold will be included in the network.
  }
  \item{names}{logical scalar indicating whether metabolite or gene KEGG IDs
  will be transformed into common metabolite names or gene symbols. Reaction IDs
  remain unchanged.
  }
  \item{export_cytoscape}{logical scalar indicating whether network and attribute
  Cytoscape files will be generated and exported.}
  \item{file_name}{character vector that allows customizing the name of the
  exported files.}
}
\value{
A matrix where each row represents an edge between two nodes. By
default, the function also generates a network file ("MS_Network.txt") and
two attribute files ("MS_NodesType.txt", "MS_TargetNodes.txt"), which can be
imported into Cytoscape to visualize the network.
}
\note{
The network-table generated with this function can be also visualized in R using
the igraph package. The network-table can be transformed into an igraph object
using the function "graph.data.frame( )" from igraph.
}
\references{
Csardi, G. & Nepusz, T. (2006). The igraph software package for complex network
research. InterJournal, Complex Systems, 1695.

Shannon, P., et al. (2003). Cytoscape: a software environment for integrated
models of biomolecular interaction networks. Genome Research, 13, 2498-2504.
}
\examples{\donttest{
data(MetaboSignal_table)

# Shortest paths from G6PC ("K01084") to pyruvate ("cpd:C00022") and
# to a-D-Glucose ("cpd:C00267")

subnet_first <- MS_shortestPathsNetwork(MetaboSignal_table, organism_code = "rno",
                                        source_nodes = "K01084",
                                        target_nodes = c("cpd:C00022", "cpd:C00267"),
                                        mode = "SP", type = "first")

subnet_all <- MS_shortestPathsNetwork(MetaboSignal_table, organism_code = "rno",
                                      source_nodes = "K01084",
                                      target_nodes = c("cpd:C00022", "cpd:C00267"),
                                      mode = "SP", type = "all")
}
}
