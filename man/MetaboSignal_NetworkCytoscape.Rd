\name{MetaboSignal_NetworkCytoscape}
\alias{MetaboSignal_NetworkCytoscape}
\title{Build shortest-path subnetwork}
\description{
This function allows calculating the shortest paths from a set of genes to a set
of metabolites, and representing them as a network-table (i.e. two-column matrix).
By default, the function exports a network file ("CytoscapeNetwork.txt") and two
attribute files ("CytoscapeAttributesType.txt", "CytoscapeAttributesTarget.txt"),
which can be imported into cytoscape to visualize the network. The first attribute
file allows customizing the nodes of the network based on the molecular entity
they represent: metabolic-genes, signaling-genes, or metabolites. The second
attribute file allows discriminating the source_genes and the target_metabolites
("target") from any other node ("untarget") of the network.

The network-table generated with this function can be further customized based
on different criteria. For instance, undesired nodes can be removed or replaced
using the functions "MS_RemoveNode( )" or "MS_ReplaceNode( )" respectively. The
final version of the network-table can be used to generate new cytoscape files
using the function "MS_ToCytoscape( )".
}
\usage{
MetaboSignal_NetworkCytoscape(network_table, organism_code, organism_name, source_genes,
                              target_metabolites, mode = "SP", type = "first",
                              distance_th = Inf, collapse_genes = FALSE, names = TRUE,
                              export_cytoscape = TRUE, file_name = "Cytoscape")
}
\arguments{
  \item{network_table}{two-column matrix where each row represents an edge
  between two nodes. See function "MetaboSignal_matrix ( )".
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_FindKEGG( )".
  }
  \item{organism_name}{character vector containing the common name of the organism of
  interest (e.g. "rat", "mouse", "human", "zebrafish") or taxonomy id. For more
  details, check: http://docs.mygene.info/en/latest/doc/data.html#species. This
  argument is only required when source_genes are gene symbols.
  }
  \item{source_genes}{character vector containing the genes from which the
  shortest paths will be calculated. All input genes need to have the same ID
  format. Possible ID formats are: entrez IDs, official gene symbols, or gene
  nodes of the network (i.e. KEGG orthology IDs or KEGG gene IDs). The latter
  option allows reducing the time required to compute this function. Entrez IDs
  or gene symbols can be transformed into KEGG IDs using the function
  "MS_GetKEGG_GeneID( )".
  }
  \item{target_metabolites}{character vector containing the KEGG IDs of the
  metabolites to which the shortest paths will be calculated. Compound KEGG IDs
  can be obtained using the function "MS_FindKEGG( )".
  }
  \item{mode}{character constant indicating whether a directed or an undirected
  network will be considered. "all" indicates that all the edges of the network
  will be considered as undirected. "out" indicates that all the edges of the
  network will be considered as directed. "SP" indicates that all network will
  be considered as directed except the edges linked to target metabolite, which
  will be considered as undirected. The difference between the "out" and the
  "SP" options, is that the latter aids reaching target metabolites that are
  substrates of irreversible reactions. By default, mode = "SP".
  }
  \item{type}{character constant indicating whether all shortest paths or a
  single shortest path will be considered when there are several shortest paths
  between a source_gene and a target_metabolite. If type = "all", all shortest
  paths will be considered. If type = "first" a single path will be considered.
  If type = "bw" the path with the highest betweenness score will be considered.
  The betweenness score is calculated as the average betweenness of the gene
  nodes of the path. Note that using type = "bw" increases the time required to
  compute this function. By default, type = "first".
  }
  \item{distance_th}{establishes a shortest path length threshold. Only shortest
  paths with length below this threshold will be included in the network. By
  default, distance_th = Inf.
  }
  \item{collapse_genes}{logical scalar indicating whether KEGG gene IDs will be
  transformed into orthology IDs. Since several gene isoforms are associated
  with the same orthology ID, this options leads to a dramatic decrease in the
  dimensionality of the network. This argument is ignored if the gene nodes of
  the network_table already represent orthology IDs. By default, collapse_genes
  = FALSE.
  }
  \item{names}{logical scalar indicating whether the metabolite or gene KEGG IDs
  will be transformed into common metabolite names or gene symbols. Reaction IDs
  remain unchanged. By default, names = TRUE.
  }
  \item{export_cytoscape}{logical scalar indicating whether network and attribute
  cytoscape files will be generated and exported. By default, export_cytoscape =
  TRUE.
  }
  \item{file_name}{character vector that allows customizing the name of the
  exported files. By default, file_name = "Cytoscape".
  }
}
\value{
A two-column matrix where each row represents an edge between two nodes. By
default, the function also generates a network file ("CytoscapeNetwork.txt") and
two attribute files ("CytoscapeAttributesType.txt", "CytoscapeAttributesTarget.txt"),
which can be imported into cytoscape to visualize the network.
}
\note{
The network-table generated with this function can be also visualized in R using
the igraph package. The network-table can be transformed into an igraph object
using the function "graph.data.frame( )" from igraph.
}
\references{
Csardi, G. & Nepusz, T. (2006). The igraph software package for complex network
research. InterJournal, Complex Systems, 1695.

Shannon, P., Markiel, A., Ozier, O., Baliga, N.S., Wang, J.T., Ramage, D., Amin,
N. & Ideker, B.S.T. (2003). Cytoscape: a software environment for integrated
models of biomolecular interaction networks. Genome Research, 13, 2498-2504.
}
\examples{
data(MetaboSignal_table)

# Shortest-path subnetwork from Foxo1 (84482), Ldha (24533) to alpha D-glucose
#("cpd:C00267") and lactate ("cpd:C00186"). Different source_gene formats are valid:

# 1) Source_genes as network IDs (in this case orthology IDs): fastest option.
# To get gene KEGG IDs use "MS_GetKEGG_GeneID( )", as shown below:
\donttest{
MS_GetKEGG_GeneID(c("foxo1", "ldha"), organism_code = "rno", organism_name = "rat")
}
subnet_KEGG <- MetaboSignal_NetworkCytoscape(MetaboSignal_table, organism_code="rno",
                                             source_genes = c("K07201", "K00016"),
                                             target_metabolites = c("cpd:C00267",
                                                                    "cpd:C00186"),
                                             names = FALSE)
\donttest{
# 2) Source_genes as entrez IDs

subnet_Entrez <- MetaboSignal_NetworkCytoscape(MetaboSignal_table, organism_code="rno",
                                               source_genes = c("84482", "24533"),
                                               target_metabolites = c("cpd:C00267",
                                                                    "cpd:C00186"),
                                               names = FALSE)

# 3) Source_genes as symbols

subnet_Symbol <- MetaboSignal_NetworkCytoscape(MetaboSignal_table,
                                               organism_code="rno", organism_name ="rat",
                                               source_genes = c("foxo1", "ldha"),
                                               target_metabolites = c("cpd:C00267",
                                                                      "cpd:C00186"),
                                               names = FALSE)
}
}


