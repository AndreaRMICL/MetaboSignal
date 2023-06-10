\name{MS_findMappedNodes}
\alias{MS_findMappedNodes}
\title{Map gene IDs or metabolite IDs onto the network}
\description{
This function can be used to find out if a set of genes or metabolites of interest
can be mapped onto the network.
}
\usage{
MS_findMappedNodes(nodes, network_table)
}
\arguments{
  \item{nodes}{character vector containing the IDs of the genes or the metabolites
  to be mapped onto the network. Remember that Entrez IDs or gene symbols can be
  transformed into KEGG IDs using the function "MS_convertGene( )".}
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. See function "MS_keggNetwork( )".
  }
}
\value{
A list reporting which genes or metabolites can or cannot be mapped onto the network.
}
\references{
Carlson, M. org.Hs.eg.db: Genome wide annotation for Human.R package version >=
3.2.3.

Mark, A., et al.(2014) mygene: Access MyGene.Info_ services. R package version >=
1.6.0.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{
data(MetaboSignal_table)

# Map D-glucose ("cpd:C00031"), taurine ("cpd:C00245"), and aldh ("K00128") onto
# onto the network

MS_findMappedNodes(nodes = c("cpd:C00031","cpd:C00245", "K00128"), MetaboSignal_table)
}
