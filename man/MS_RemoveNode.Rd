\name{MS_RemoveNode}
\alias{MS_RemoveNode}
\title{Remove undesired nodes from the network}
\description{
This function allows removing undesired nodes of the network-table.
}
\usage{
MS_RemoveNode(nodes, network_table)
}
\arguments{
  \item{nodes}{character vector containing the node IDs to be removed.
  }
  \item{network_table}{two-column matrix where each row represents and edge
  between two nodes. See function "MetaboSignal_matrix ( )".
  }
}
\value{
A two-column matrix corresponding to the input network-table without the undesired nodes.
}
\examples{
data(MetaboSignal_table)

# Remove glucose nodes

GlucoseRemoved <- MS_RemoveNode(nodes = c("cpd:C00267", "cpd:C00221", "cpd:C00031"),
                                MetaboSignal_table)
}
