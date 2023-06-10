\name{MS_exportCytoscape}
\alias{MS_exportCytoscape}
\title{Export network in cytoscape format}
\description{
This function generates a network file and two attribute files ("NodesType.txt",
"TargetNodes.txt"), which can be imported into Cytoscape to visualize the network.
The first attribute file allows customizing the nodes of the network based on
the molecular entity they represent: compound, reaction, metabolic-gene or
signaling-gene. The second attribute file allows highlighting a set of nodes of
interest.
}
\usage{
MS_exportCytoscape(network_table, organism_code, names = TRUE,
                   targets = NULL, file_name = "MS")
}
\arguments{
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. Nodes must be KEGG IDs, not common names. See function
  "MS_keggNetwork()". For human networks, Entrez gene IDs are also allowed.
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See function
  "MS_keggFinder( )".
  }
  \item{names}{logical scalar indicating whether metabolite or gene KEGG IDs
  will be transformed into common metabolite names or gene symbols. Reaction IDs
  remain unchanged.
  }
  \item{targets}{optional character vector containing the IDs of the target
  nodes to be discriminated from the other nodes of the network.
  }
  \item{file_name}{character vector that allows customizing the name of the
  exported files.
  }
}
\value{
A data frame where each row represents an edge between two nodes (from source to
target). The function also generates and exports a network file ("MS_Network.txt")
and two attribute files ("MS_NodesType.txt", "MS_TargetNodes.txt"), which can be
imported into Cytoscape to visualize the network.
}
\references{
Shannon P et al. (2003). Cytoscape: a software environment for integrated models of
biomolecular interaction networks. Genome Research, 13, 2498-2504.
}
\examples{\donttest{
data(MetaboSignal_table)
MS_exportCytoscape(MetaboSignal_table, organism_code = "rno", names = FALSE)
}
}
