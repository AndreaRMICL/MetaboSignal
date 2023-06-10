\name{MS_nodeBW}
\alias{MS_nodeBW}
\title{Get distribution of node betweeness}
\description{
This function calculates the betweenness of each node of the network.
}
\usage{
MS_nodeBW(network_table, mode = "all", normalized = TRUE)
}
\arguments{
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. See function "MS_keggNetwork( )".
  }
  \item{mode}{character constant indicating whether a directed ("out") or
  undirected ("all") network will be considered.
  }
  \item{normalized}{logical scalar indicating whether to normalize the betweeness
  scores. If TRUE, normalized betweenness scores will be returned. If FALSE, raw
  betweenness scores will be returned.
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
MS_nodeBW(MetaboSignal_table)
}

