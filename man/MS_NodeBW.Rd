\name{MS_NodeBW}
\alias{MS_NodeBW}
\title{Get distribution of node betweeness}
\description{
This function calculates the betweenness of each node of the network.
}
\usage{
MS_NodeBW(network_table, mode = "all", normalized = TRUE)
}
\arguments{
  \item{network_table}{two-column matrix where each row represents and edge
  between two nodes. See function "MetaboSignal_matrix ( )".
  }
  \item{mode}{character constant indicating whether a directed ("out") or
  undirected ("all") network will be considered. By default, mode = "all".
  }
  \item{normalized}{logical scalar indicating whether to normalize the betweeness
  scores. If TRUE, normalized betweenness scores will be returned. If FALSE, raw
  betweenness scores will be returned. By default, normalized = TRUE.
  }
}
\value{
A numeric vector containing the betweenness of each node of the network. The
function also produces and histogram showing the distribution of node betweenness.
}
\references{
Csardi, G. & Nepusz, T. (2006). The igraph software package for complex network
research. InterJournal, Complex Systems, 1695.
}
\examples{
data(MetaboSignal_table)
MS_NodeBW(MetaboSignal_table)
}

