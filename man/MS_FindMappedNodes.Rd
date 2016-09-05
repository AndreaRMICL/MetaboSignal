\name{MS_FindMappedNodes}
\alias{MS_FindMappedNodes}
\title{Map gene IDs or metabolite IDs onto the network}
\description{
This function can be used to find out if a set of genes or metabolites of interest
can be mapped onto the MetaboSignal_network.
}
\usage{
MS_FindMappedNodes(nodes, network_table, organism_code, organism_name,
                   orthology = TRUE)
}
\arguments{
  \item{nodes}{character vector containing the IDs of the genes or the metabolites to
  be mapped onto the network. All IDs need to correspond to the same molecular entity
  (i.e. gene or metabolite). For metabolites, KEGG IDs are required (see function
  "MS_FindKEGG( )"). For genes, entrez IDs or official symbols can be used, but note
  that all genes need to be in the same ID format (i.e. entrez or symbols). It is
  preferable to use entrez IDs rather than gene symbols, since some gene symbols
  are not unique.
  }
  \item{network_table}{two-column matrix where each row represents and edge
  between two nodes. See function "MetaboSignal_matrix ( )".
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_FindKEGG( )". This argument is not required to map metabolites.
  }
  \item{organism_name}{character vector containing the common name of the organism of
  interest (e.g. "rat", "mouse", "human", "zebrafish") or taxonomy id. For more
  details, check: http://docs.mygene.info/en/latest/doc/data.html#species. This
  argument is not required to map metabolites.
  }
  \item{orthology}{logical scalar indicating whether the gene nodes of the
  network_table represent KEGG orthology IDs (FALSE) or specific KEGG gene IDs
  (TRUE). By default, orthology = TRUE.
  }
}
\value{
A list reporting which genes or metabolites can or cannot be mapped onto the network.
}
\references{
Carlson, M. org.Hs.eg.db: Genome wide annotation for Human.R package version 3.2.3.

Mark, A., Thompson, R. and Wu, C. (2014) mygene: Access MyGene.Info_ services.R package version 1.6.0.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{
data(MetaboSignal_table)

# Map D-glucose ("cpd:C00031") and taurine ("cpd:C00245") onto the network

MS_FindMappedNodes(nodes = c("cpd:C00031","cpd:C00245"), MetaboSignal_table,
                   orthology = TRUE)
\donttest{
# Map entrez IDs onto the network

MS_FindMappedNodes(nodes = c("303565", "24267", "11114"), MetaboSignal_table,
                   organism_code = "rno", organism_name = "rat", orthology = TRUE)

# Map gene symbols onto the network

MS_FindMappedNodes(nodes = c("G6pc3","Comt"), MetaboSignal_table, organism_code = "rno",
                   organism_name = "rat",orthology = TRUE)

}
}
