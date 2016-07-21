\name{MetaboSignal_distances}
\alias{MetaboSignal_distances}
\title{Calculate gene-metabolite distance matrix}
\description{
This function generates a distance matrix containing the length of all shortest
paths from a set of genes (or reactions) to a set of metabolites. The shortest
path length between two nodes is defined as the minimum number of edges between
these two nodes.
}
\usage{
MetaboSignal_distances(network_table, organism_code, organism_name, mode = "SP",
                       source_genes = "all", target_metabolites = "all", names = FALSE)
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
  \item{mode}{character constant indicating whether a directed or an undirected
  network will be considered. "all" indicates that all the edges of the network
  will be considered as undirected. "out" indicates that all the edges of the
  network will be considered as directed. "SP" indicates that all network will
  be considered as directed except the edges linked to target metabolite, which
  will be considered as undirected. The difference between the "out" and the "SP"
  options, is that the latter aids reaching target metabolites that are substrates
  of irreversible reactions. By default, mode = "SP".
  }
  \item{source_genes}{character vector containing the genes from which the
  shortest paths will be calculated. All input genes need to have the same ID
  format. Possible ID formats are: entrez IDs, official gene symbols, or gene
  nodes of the network (i.e. KEGG orthology IDs or KEGG gene IDs). The latter
  option allows reducing the time required to compute this function. Entrez IDs
  or gene symbols can be transformed into KEGG IDs using the function
  "MS_GetKEGG_GeneID( )". By default, genes = "all" indicating that all genes or
  reactions of the network will be used.
  }
  \item{target_metabolites}{character vector containing the KEGG IDs of the
  metabolites to which the shortest paths will be calculated. Compound KEGG IDs
  can be obtained using the function "MS_FindKEGG( )". By default,
  metabolites = "all", indicating that all metabolites of the network will be
  used.
  }
  \item{names}{logical scalar indicating whether the metabolite IDs or gene KEGG
  IDs will be transformed into common metabolite names or gene symbols. Reaction
  IDs remain unchanged. By default, names = FALSE.
  }
 }
\value{
A matrix containing the shortest path length from the genes or reactions (in the
rows) to the metabolites (in the columns). For unreacheable metabolites Inf is
included.
}
\references{
Csardi, G. & Nepusz, T. (2006). The igraph software package for complex network
research. InterJournal, Complex Systems, 1695.
}
\examples{\donttest{
data(MetaboSignal_table)

# Distances from Ship2 (65038) and Ppp2r5b (309179) to D-glucose ("cpd:C00031")

distances_targets <- MetaboSignal_distances(MetaboSignal_table, organism_code = "rno",
                                            source_genes = c("65038", "309179"),
                                            target_metabolites = "cpd:C00031",
                                            names = TRUE)
}
# Distances from all genes to all metabolites of the network

distances_all <- MetaboSignal_distances(MetaboSignal_table, organism_code = "rno")
}

