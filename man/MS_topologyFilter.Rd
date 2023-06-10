\name{MS_topologyFilter}
\alias{MS_topologyFilter}
\title{Filter network based on distances or betweenness}
\description{
This function allows reducing the dimensionality of a network, by removing nodes
that do not meet the established distance and/or node betweenness criteria.
}
\usage{
MS_topologyFilter(network_table, mode = "all", type, target_node, distance_th, bw_th)
}
\arguments{
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. See function "MS_keggNetwork( )".
  }
  \item{mode}{character constant indicating whether a directed ("out") or
  undirected ("all") network will be considered.}
  \item{type}{character constant used to establish the criteria for filtering
  the network. "bw" indicates that edges (i.e. rows of the network_table)
  containing at least one node with betweenness below bw_th will be neglected.
  "distance" indicates edges containing at least one node with shortest path
  length to the target_node above distance_th will be neglected. "all" indicates
  that edges containing at least one node with either betweenness below bw_th or
  distance above distance_th, will be neglected.}
  \item{target_node}{character vector containing the ID of the node to which the
  distances will be calculated.}
  \item{distance_th}{numeric value corresponding to the distance threshold.
  Nodes with shortest path length to the target_node above this threshold will
  be removed from the network-table.}
  \item{bw_th}{numeric value corresponding to the normalized-betweenness
  threshold. Nodes with betweenness below this threshold will be removed from
  the network-table. See also "MS_nodeBW( )".}
}
\value{
A three-column matrix where each row represents an edge between two nodes.
}
\references{
Csardi, G. & Nepusz, T. (2006). The igraph software package for complex network
research. InterJournal, Complex Systems, 1695.
}
\examples{
data(MetaboSignal_table)
# Remove edges containing nodes with distance to D-glucose ("cpd:C00031") > 2

network_filtered1 <- MS_topologyFilter(MetaboSignal_table, type = "distance",
                                       target_node = "cpd:C00031",
                                       distance_th = 2)

# Remove edges containing nodes with distance to D-glucose ("cpd:C00031") > 2 or
# normalized-betweenness < 0.00005

network_filtered2 <- MS_topologyFilter(MetaboSignal_table, type = "all",
                                       target_node = "cpd:C00031",
                                       distance_th = 2, bw_th = 0.00005)

# Note below that network_filtered1 has one edge more than network_filtered2. This is
# because "cpd:C00031" has betweenness = 0, and therefore it is removed in network_filtered2:

setdiff(as.vector(network_filtered1[, 1:2]),as.vector(network_filtered2[, 1:2]))
}

