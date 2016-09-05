\name{MS_GetKEGG_GeneID}
\alias{MS_GetKEGG_GeneID}
\title{Transform entrez IDs or gene symbols into KEGG IDs}
\description{
This function allows transforming entrez gene IDs or official gene symbols into
KEGG IDs (orthology IDs or organism-specific gene IDs). The transformed KEGG IDs
can be stored and used as source_genes in the functions "MetaboSignal_distances( )"
or "MetaboSignal_NetworkCytoscape( )". This strategy allows reducing the time
required to compute these functions.
}
\usage{
MS_GetKEGG_GeneID(genes, organism_code, organism_name, output = "vector",
                  orthology = TRUE)
}
\arguments{
  \item{genes}{character vector containing the entrez IDs or official symbols of
  the genes of interest. All genes need to be in the same ID format (i.e. entrez
  or symbols). It is preferable to use entrez IDs rather than gene symbols, since
  some gene symbols are not unique.
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_FindKEGG( )".
  }
  \item{organism_name}{character vector containing the common name of the organism of
  interest (e.g. "rat", "mouse", "human", "zebrafish") or taxonomy id. For more
  details, check: http://docs.mygene.info/en/latest/doc/data.html#species. This
  argument is only required when gene symbols are used.
  }
  \item{output}{character constant indicating whether the function will return a
  vector containing mapped and transformed KEGG IDs (output = "vector"), or a
  matrix containing both mapped entrez IDs or gene symbols and their
  corresponding KEGG IDs (output = "matrix"). By default, output = "vector".
  }
  \item{orthology}{logical scalar indicating whether the gene IDs will be
  transformed into orthology IDs or into organism-specific gene IDs. By
  default, orthology = TRUE.
  }
}
\value{
A character vector containing mapped and transformed KEGG IDs or a matrix
containing both mapped entrez IDs or gene symbols and their corresponding KEGG
IDs.
}
\references{
Carlson M.(2016). org.Hs.eg.db: Genome wide annotation for Human, Bioconductor,
v3.2.3.

Mark A, et al.(2016). mygene: Access MyGene.Info_ services. Bioconductor, v1.8.0.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{
# Transform gene symbol Hoga1 (293949) into rat-specific KEGG ID

MS_GetKEGG_GeneID(genes = "Hoga1", organism_code = "rno", organism_name = "rat",
                  orthology = FALSE)
\donttest{
MS_GetKEGG_GeneID(genes = "Hoga1", "rno", "rat", output = "matrix", orthology = FALSE)

# Transform entrez ID 293949 into orthology KEGG ID

MS_GetKEGG_GeneID(genes = "293949", organism_code = "rno", output = "matrix")
}
}
