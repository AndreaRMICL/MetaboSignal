\name{MS_keggNetwork}
\alias{MS_keggNetwork}
\title{Build MetaboSignal network-table}
\description{
This function generates a directed network-table (i.e. three-column matrix), where
each row represents an edge connecting two nodes (from source to target). Nodes
represent different molecular entities: metabolic-genes (i.e. genes encoding
enzymes that catalyze metabolic reactions), signaling-genes (e.g. kinases),
reactions and compounds (metabolites, drugs or glycans). The third column of the
matrix indicates the interaction type. Compound-gene (or gene-compound) interactions
are designated as: "k_compound:reversible" or "kegg_compound:irreversible",
depending on the direction of the interaction. Other types of interactions
correspond to gene-gene interactions. When KEGG reports various types of interaction
for the same gene pair, the "interaction_type" is collapsed using "/".

The network-table generated with this function can be customized based on several
criteria. For instance, undesired nodes can be removed or replaced using the
functions "MS_removeNode( )" or "MS_replaceNode( )" respectively. Also, the
network can be filtered according to different topological parameters (e.g. node
betweenness) using the function "MS_topologyFilter( )".
}
\usage{
MS_keggNetwork(metabo_paths, signaling_paths, expand_genes = FALSE,
               convert_entrez = FALSE)
}
\arguments{
  \item{metabo_paths}{character vector containing the KEGG IDs of the metabolic
  pathways of interest (organism-specific). Pathway IDs take the form: "organism
  code + 5-digit number". For example, the ID of the rat "glycolysis/gluconeogenesis"
  pathway is "rno00010". See functions "MS_keggFinder( )" and "MS_getPathIds( )".
  }
  \item{signaling_paths}{character vector containing the KEGG IDs for the
  signaling pathways of interest (organism-specific). For example, the ID
  for the pathway "insulin signaling pathway" in the rat is "rno04910".
  See functions "MS_keggFinder( )" and "MS_getPathIds( )".
  }
  \item{expand_genes}{logical scalar indicating whether the gene nodes will
  represent orthology IDs (FALSE) or organism-specific gene IDs (TRUE).
  }
  \item{convert_entrez}{logical scalar indicating whether the KEGG gene IDs will
  be transformed into Entrez IDs. This argument will be ignored if expand_genes =
  FALSE, or if the input paths are not human-specific.}
}
\value{
A three-column matrix where each row represents an edge between two nodes.
}
\note{
Reaction directionality reported in KEGG has been cross-validated with published
literature (Duarte et al., 2007).
}
\references{
Davidovic, L., et al. (2011). A metabolomic and systems biology perspective on the
brain of the fragile X syndrome mouse model. Genome Research, 21, 2190-2202.

Duarte, N.C., et al. (2007). Global reconstruction of the human metabolic network
based on genomic and bibliomic data. Proceedings of the National Academy of Sciences,
104, 1777-1782.

Posma, J.M., et al.(2014). MetaboNetworks, an interactive Matlab-based toolbox for
creating, customizing and exploringsub-networks from KEGG. Bioinformatics, 30, 893-895.

Zhang, J.D. & Wiemann, S. (2009). KEGGgraph: a graph approach to KEGG PATHWAY in
R and Bioconductor. Bioinformatics, 25, 1470-1471.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{\donttest{
# MetaboSignal network-table with organism-specific gene nodes

MS_netIsoforms <- MS_keggNetwork(metabo_paths = c("rno00010", "rno00562"),
                                 signaling_paths = c("rno04910", "rno04151"),
                                 expand_genes = TRUE)

# MetaboSignal network-table with orthology gene nodes

MS_netK <- MS_keggNetwork(metabo_paths = c("rno00010", "rno00562"),
                         signaling_paths = c("rno04910", "rno04151"))

# MetaboSignal network-table with human Entrez gene IDs

MS_netEntrez <- MS_keggNetwork(metabo_paths = c("hsa00010", "hsa00562"),
                               signaling_paths = c("hsa04910", "hsa04151"),
                               expand_genes = TRUE, convert_entrez = TRUE)
}
}
