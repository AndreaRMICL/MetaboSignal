\name{MS_keggFinder}
\alias{MS_keggFinder}
\title{Get KEGG IDs for compounds, organisms or pathways}
\description{
This function returns a list of entries corresponding to one of the following
KEGG databases: "compound", "organism", "pathway". It can also find entries with
matching query keywords in a given database.
}
\usage{
MS_keggFinder(KEGG_database, match = NULL, organism_code)
}
\arguments{
  \item{KEGG_database}{character vector containing the name of the KEGG database
  of interest: "compound", "organism", "pathway".}
  \item{match}{character vector containing one or more elements (i.e. key words)
  to be matched as compound names.}
  \item{organism_code}{character vector containing the KEGG code for the organism
  of interest. For example the KEGG code for the rat is "rno". This argument is only
  required for KEGG_database = "pathway".}
}
\value{
By default, a matrix where each row contains the KEGG entries of the database of
interest. When using the option "match" a list is returned, each list element
containing information of matched entries.
}
\examples{
MS_keggFinder(KEGG_database = "compound", match = "acetoacetic acid")
\donttest{
MS_keggFinder(KEGG_database = "organism", match = c("rattus","human"))

MS_keggFinder(KEGG_database = "pathway", match = c("glycol", "insulin signal", "akt"),
            organism_code = "rno")

}
}
