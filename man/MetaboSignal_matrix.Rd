\name{MetaboSignal_matrix}
\alias{MetaboSignal_matrix}
\title{Build MetaboSignal network-table}
\description{
This function generates a directed network-table (i.e. two-column matrix), where
each row represents an edge connecting two nodes (from node in column 1 to node
in column 2). Nodes represent four different molecular entities: metabolic-genes
(i.e. genes encoding enzymes that catalyze metabolic reactions), signaling-genes (e.g.
kinases), reactions or metabolites. It is possible to build a tissue-specific
network-table that excludes signaling genes that are not expressed in a given
tissue. Tissue expression data is obtained using the hpar package, which is
based on the The Human Protein Atlas database. The genes "non detected" in the
target tissue (reliability = supportive) are neglected.

The network-table generated with this function can be customized based on several
criteria. For instance, undesired nodes can be removed or replaced using the
functions "MS_RemoveNode( )" or "MS_ReplaceNode( )" respectively. Also, the
network can be filtered according to different topological parameters (e.g. node
betweenness) using the function "MS_FilterNetwork( )".
}
\usage{
MetaboSignal_matrix(metabo_paths, signaling_paths, organism_name, tissue = "all",
                    expand_genes = FALSE)
}
\arguments{
  \item{metabo_paths}{character vector containing the KEGG IDs of the metabolic
  pathways of interest (organism-specific). For example, the KEGG ID for the
  pathway "glycolysis/gluconeogenesis" in the rat is"rno00010". See the function
  "MS_FindKEGG( )".
  }
  \item{signaling_paths}{character vector containing the KEGG IDs for the
  signaling pathways of interest (organism-specific). For example, the KEGG ID
  for the pathway "insulin signaling pathway" in the rat is "rno04910".
  }
  \item{organism_name}{character vector containing the common name of the organism of
  interest (e.g. "rat", "mouse", "human", "zebrafish") or taxonomy id. For more
  details, check: http://docs.mygene.info/en/latest/doc/data.html#species. This
  argument is only required when filtering genes by tissue expression.
  }
  \item{tissue}{character vector containing the name(s) of the target tissue(s).
  By default, tissue = "all" indicating that signaling gene nodes will not be
  filtered by tissue expression. Otherwise, possible tissues are those included in
  the dataset hpaNormalTissue (see levels(hpaNormalTissue[,2])).
  }
  \item{expand_genes}{logical scalar indicating whether the gene nodes will
  represent orthology IDs (FALSE) or organim-specific gene IDs (TRUE). By default,
  expand_genes = FALSE.
  }
}
\value{
A two-column matrix where each row represents an edge between two nodes.
}
\note{
Reaction directionality reported in KEGG has been cross-validated with published
literature (Duarte et al., 2007).
}
\references{
Carlson, M. org.Hs.eg.db: Genome wide annotation for Human.R package version 3.2.3.

Davidovic, L., Navratil, V., Bonaccorso, C.M., Catania, M.V., Bardoni, B. & Dumas,
M.E. (2011). A metabolomic and systems biology perspective on the brain of the
fragile X syndrome mouse model. Genome Research, 21, 2190-2202.

Duarte, N.C., Becker, S.A., Jamshidi, N., Thiele, I., Mo, M.L., Vo, T.D., Srivas,
R. & Palsson, O.B. (2007). Global reconstruction of the human metabolic network
based on genomic and bibliomic data. Proceedings of the National Academy of Sciences,
104, 1777-1782.

Gatto, L. hpar: Human Protein Atlas in R.R package version 1.12.0.

Mark, A., Thompson, R. and Wu, C. (2014) mygene: Access MyGene.Info_ services.R
package version 1.6.0.

Posma, J.M., Robinette, S.L., Holmes, E. & Nicholson, J.K. (2014). MetaboNetworks,
an interactive Matlab-based toolbox for creating, customizing and exploring
sub-networks from KEGG. Bioinformatics, 30, 893-895.

Zhang, J.D. & Wiemann, S. (2009). KEGGgraph: a graph approach to KEGG PATHWAY in
R and Bioconductor. Bioinformatics, 25, 1470-1471.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{\donttest{
# MetaboSignal network-table with organism-specific gene nodes

MetaboSignal_tableIsoforms <- MetaboSignal_matrix(metabo_paths = c("rno00010", "rno00562"),
                                                  signaling_paths = c("rno04910",
                                                                    "rno04151"),
                                                  expand_genes = TRUE)

# MetaboSignal network-table with orthology gene nodes

MetaboSignal_table <- MetaboSignal_matrix(metabo_paths = c("rno00010", "rno00562"),
                                          signaling_paths = c("rno04910", "rno04151"))

# MetaboSignal network-table orthology gene nodes filtered by liver

MetaboSignal_tableLiver <- MetaboSignal_matrix(metabo_paths = "hsa00010",
                                               signaling_paths = "hsa04151",
                                               organism_name = "human",
                                               tissue = "liver")
}
}



