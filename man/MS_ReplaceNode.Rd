\name{MS_replaceNode}
\alias{MS_replaceNode}
\title{Replace nodes of the network}
\description{
This function allows replacing node IDs of a network-table. It can be used to
cluster the IDs of chemical isomers (e.g. alpha-D-glucose ("cpd:C00267"),
D-glucose ("cpd:C00031"), and beta-D-glucose  ("cpd:C00021")) into a single ID.
}
\usage{
MS_replaceNode(node1, node2, network_table)
}
\arguments{
  \item{node1}{character vector containing the node IDs to be replaced.
  }
  \item{node2}{character vector containing the ID that will be used as a
  replacement.
  }
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. See function "MS_keggNetwork( )".
  }
}
\value{
A three-column matrix corresponding to the input network-table with replaced nodes.
}
\examples{
data(MetaboSignal_table)

# Cluster D-glucose isomers ("cpd:C00267","cpd:C00221","cpd:C00031")

glucoseClustered <- MS_replaceNode(node1 = c("cpd:C00267", "cpd:C00221"),
                                   node2 = "cpd:C00031", MetaboSignal_table)
}

