\name{MS_removeNode}
\alias{MS_removeNode}
\title{Remove undesired nodes from the network}
\description{
This function allows removing undesired nodes from the network-table.
}
\usage{
MS_removeNode(nodes, network_table)
}
\arguments{
  \item{nodes}{character vector containing the node IDs to be removed.
  }
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. See function "MS_keegNetwork( )".
  }
}
\value{
A three-column matrix corresponding to the input network-table without the undesired nodes.
}
\examples{
data(MetaboSignal_table)

# Remove glucose nodes

glucoseRemoved <- MS_removeNode(nodes = c("cpd:C00267", "cpd:C00221", "cpd:C00031"),
                                MetaboSignal_table)
}
