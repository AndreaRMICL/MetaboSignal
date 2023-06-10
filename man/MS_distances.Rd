\name{MS_distances}
\alias{MS_distances}
\title{Calculate gene-metabolite distance matrix}
\description{
This function generates a distance matrix containing the length of all shortest
paths from a set of genes (or reactions) to a set of metabolites. The shortest
path length between two nodes is defined as the minimum number of edges between
these two nodes.
}
\usage{
MS_distances(network_table, organism_code, mode = "SP", source_genes = "all",
             target_metabolites = "all", names = FALSE)
}
\arguments{
  \item{network_table}{three-column matrix where each row represents an edge
  between two nodes. See function "MS_keggNetwork( )".}
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_keggFinder( )".}
  \item{mode}{character constant indicating whether a directed or an undirected
  network will be considered. "all" indicates that all the edges of the network
  will be considered as undirected. "out" indicates that all the edges of the
  network will be considered as directed. "SP" indicates that all network will
  be considered as directed except the edges linked to target metabolite, which
  will be considered as undirected. The difference between the "out" and the "SP"
  options, is that the latter aids reaching target metabolites that are substrates
  of irreversible reactions.}
  \item{source_genes}{character vector containing the genes from which the
  shortest paths will be calculated. Remember that Entrez IDs or gene symbols
  can be transformed into KEGG IDs using the function "MS_convertGene( )".
  By default, source_genes = "all", indicating that all the genes of the network
  will be used.}
  \item{target_metabolites}{character vector containing the KEGG IDs of the
  metabolites to which the shortest paths will be calculated. Compound KEGG IDs
  can be obtained using the function "MS_keggFinder( )". By default,
  target_metabolites = "all", indicating that all the metabolites of the network
  will be used.}
  \item{names}{logical scalar indicating whether metabolite or gene KEGG IDs will be
  transformed into common metabolite names or gene symbols. Reaction
  IDs remain unchanged.}
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

MS_convertGene(genes = c("65038","309179"), "rno", "rat", output = "matrix")

distances_targets <- MS_distances(MetaboSignal_table, organism_code = "rno",
                                  source_genes = c("K15909", "K11584"),
                                  target_metabolites = "cpd:C00031",
                                  names = TRUE)
}
# Distances from all genes to all metabolites of the network

distances_all <- MS_distances(MetaboSignal_table, organism_code = "rno")
}

