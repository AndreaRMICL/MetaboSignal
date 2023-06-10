\name{MS_changeNames}
\alias{MS_changeNames}
\title{Transform KEGG IDs into common names}
\description{
This function allows transforming KEGG IDs of genes or compounds into their
corresponding common names (for compounds) or symbols (for genes).
}
\usage{
MS_changeNames(nodes, organism_code)
}
\arguments{
  \item{nodes}{character vector or matrix containing the KEEG IDs of either
  metabolites, genes (organism-specific or orthology), or reactions. It also converts
  human Entrez gene IDs into symbols.}
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_keggFinder( )". This argument is ignored when nodes are metabolites.}
}
\value{
A character string or a matrix containing the common metabolite names or gene
symbols corresponding to the input KEGG IDs. Reaction IDs remain unchanged.
}
\references{
http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{
\donttest{
MS_changeNames(c("rno:84482", "K01084", "cpd:C00267"), "rno")
MS_changeNames("K01082", organism_code = "rno")
}
}

