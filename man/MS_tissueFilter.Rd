\name{MS_tissueFilter}
\alias{MS_tissueFilter}
\title{Filter network based on tissue expression data}
\description{
This function allows filtering a network based on tissue expression data from the
Human Protein Atlas, by removing signaling genes that are not detected in the
target tissue(s) (reliability = "approved" or "supported"). This function can be
only used to filter human networks.
}
\usage{
MS_tissueFilter(network_table, tissue, input_format = "kegg", expand_genes = FALSE)
}
\arguments{
  \item{network_table}{three-column matrix where each row represents an edge
  between two nodes. The gene nodes of this network must be human specific gene
  IDS (not orthologies). For this, use the function "MS_keggNetwork( )" with
  expand_genes = TRUE.}
  \item{tissue}{character vector indicating the tissue(s) of interest. Signaling
  genes (i.e. non-enzymatic genes) not detected in the target tissue(s)
  (reliability = "approved" or "supported") will be removed from the network.
  Check all possible tissues in the "hpaNormalTissue" dataset.}
  \item{input_format}{character vector indicating the gene format in the input
  network_table ("entrez" or "kegg").}
  \item{expand_genes}{logical scalar indicating whether the gene nodes in the
  filtered network will represent orthology IDs (expand_genes = FALSE) or
  organism-specific gene IDs (expand_genes = TRUE).}
}
\value{
A three-column matrix where each row represents an edge between two nodes.
}

\references{
Gatto, L. hpar: Human Protein Atlas in R.R package version 1.12.0.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{\donttest{
# Build network

net <- MS_keggNetwork(metabo_paths = "hsa00010", signaling_paths = "hsa04014",
                      expand_genes = TRUE)

# Filter network by liver and cluster genes by orthology

net_filtered <- MS_tissueFilter(net, tissue = "liver")
}
}
