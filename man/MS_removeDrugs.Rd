\name{MS_removeDrugs}
\alias{MS_removeDrugs}
\title{Remove edges containing drug nodes}
\description{
This function allows removing edges containing drug ("dr:") nodes.
}
\usage{
MS_removeDrugs(network_table)
}
\arguments{
  \item{network_table}{three-column matrix where each row represents and edge
  between two nodes. See function "MS_keggNetwork( )".
  }
}
\value{
A three-column matrix corresponding to the input network-table without the drug
nodes.
}
\examples{
data(MetaboSignal_table)

# Remove drug nodes if present

drugsRemoved <- MS_removeDrugs(MetaboSignal_table)
}
