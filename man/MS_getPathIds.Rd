\name{MS_getPathIds}
\alias{MS_getPathIds}
\title{Get pathway identifiers of a given organism}
\description{
This function retrieves the identifiers (IDs) of all metabolic and signaling KEGG
pathways of a given organism. These pathway IDs can be used to build a MetaboSignal
network with the function "MS_keggNetwork( )".
}
\usage{
MS_getPathIds(organism_code)
}
\arguments{
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_keggFinder( )".}
}
\value{
This function returns a matrix, where each row contains the ID, description, category,
and type (i.e. "metabolic" or "signaling") of each pathway. This matrix is also
exported in a file named "organism-code_pathways.txt".
}

\references{
Tenenbaum, D. KEGGREST: Client-side REST access to KEGG. R package version
>= 1.17.0.
}
\examples{\donttest{
rat_paths <- MS_getPathIds(organism_code = "rno")
human_paths <- MS_getPathIds(organism_code = "hsa")
}
}
