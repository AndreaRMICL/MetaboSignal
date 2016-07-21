\name{MS_ToCytoscape}
\alias{MS_ToCytoscape}
\title{Export network in cytoscape format}
\description{
The function exports a network file ("CytoscapeNetwork.txt") and two attribute
files ("CytoscapeAttributesType.txt", "CytoscapeAttributesTarget.txt"), which
can be imported into cytoscape to visualize the network. The first attribute file
allows customizing the nodes of the network based on the molecular entity they
represent: metabolites, metabolic-genes, or signaling-genes. The second attribute
file allows discriminating a set of nodes of interest ("target") from any other
node ("untarget") of the network.
}
\usage{
MS_ToCytoscape (network_table, organism_code, names = TRUE, target_nodes = NULL,
                file_name = "Cytoscape")
}
\arguments{
  \item{network_table}{two-column matrix where each row represents and edge
  between two nodes. See function "MetaboSignal_matrix ( )".
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_FindKEGG( )".
  }
  \item{names}{logical scalar indicating whether the metabolite or gene KEGG IDs
  will be transformed into common metabolite names or gene symbols. Reaction IDs
  remain unchanged. By default, names = TRUE.
  }
  \item{target_nodes}{character vector containing the IDs of the target nodes to
  be discriminated from the other nodes of the network.This argument is optional.
  }
  \item{file_name}{character vector that allows customizing the name of the
  exported files. By default, the file_name = "Cytoscape".
  }
}
\value{
A data frame where each row represents an edge between two nodes. The function also
generates and exports a network file ("CytoscapeNetwork.txt") and two attribute
files ("CytoscapeAttributesType.txt", "CytoscapeAttributesTarget.txt"), which
can be imported into cytoscape to visualize the network.
}
\references{
Shannon, P., Markiel, A., Ozier, O., Baliga, N.S., Wang, J.T., Ramage, D., Amin, N.
& Ideker, B.S.T. (2003). Cytoscape: a software environment for integrated models of
biomolecular interaction networks. Genome Research, 13, 2498-2504.
}
\examples{\donttest{
Glucolysis <- MetaboSignal_matrix(metabo_paths = "mmu00010", organism_name = "mouse")
MS_ToCytoscape(Glucolysis, organism_code = "mmu")
}
}
