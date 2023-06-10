\name{MS2_mergeNetworks}
\alias{MS2_mergeNetworks}
\title{Merge networks}
\description{
This function allows merging two network-tables of interest.
}
\usage{
MS2_mergeNetworks(network_table1, network_table2)
}
\arguments{
  \item{network_table1}{three-column matrix where each row represents an edge
  between two nodes. See functions "MS_keggNetwork()" and "MS2_ppiNetwork()".
  }
  \item{network_table2}{three-column matrix where each row represents an edge
  between two nodes. See functions "MS_keggNetwork()" and "MS2_ppiNetwork()".
  }
}
\value{
A three-column matrix where each row represents an edge between two nodes.
}
\examples{
data(keggNet_example)
data(ppiNet_example)

# Fast example using subsets
global_network1 <- MS2_mergeNetworks(keggNet_example[1:10, ],
                                     ppiNet_example[1:10, ])

\donttest{ # Example using full datasets
global_network2 <- MS2_mergeNetworks(keggNet_example, ppiNet_example)
}
}


