\name{MS_ChangeNames}
\alias{MS_ChangeNames}
\title{Transform KEGG IDs into common names}
\description{
This function allows transforming KEGG IDs of genes or metabolites into their
corresponding common names (for metabolites) or symbols (for genes).
}
\usage{
MS_ChangeNames(nodes, organism_code)
}
\arguments{
  \item{nodes}{character vector or matrix containing the KEEG IDs of either
  metabolites, genes (organism-specific or orthology), or reactions.
  }
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". See the function
  "MS_FindKEGG( )". This argument is ignored when nodes are metabolites.
  }
}
\value{
A character string or a matrix containing the common metabolite names or gene
symbols corresponding to the input KEGG IDs. Reaction IDs remain unchanged.
}
\references{
http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{
MS_ChangeNames("K01082", organism_code = "rno")
\donttest{
MS_ChangeNames(c("rno:84482", "K01084", "cpd:C00267"), "rno")
}
}

