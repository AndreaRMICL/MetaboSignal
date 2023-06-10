\name{MS2_ppiNetwork}
\alias{MS2_ppiNetwork}
\title{Build signaling-transduction network}
\description{
This function generates a directed regulatory network by merging interactions
reported in two literature-curated resources: OmniPath and TRRUST. The network
is formalized as a three-column matrix, where each row represents an edge
connecting two nodes (from source to target). The third column indicates
the type of interaction, as well as the source of the interaction (OmniPath = "o_",
TRRUST = "t_").  Nodes represent gene Entrez IDs.
}
\usage{
MS2_ppiNetwork(datasets = "all")
}
\arguments{
  \item{datasets}{character vector indicating the datasets that will be used to
  build the network ("all", "omnipath","trrust"). It is also possible to select
  databases included within OmniPath (e.g. datasets = c("biogrid", "string"))
}
}
\value{
A three-column matrix where each row represents an edge between two nodes.
}
\note{
The dataset "regulatory_interactions" contains details regarding primary database
reference(s) as well as literature reference(s) of each of the regulatory interactions.
The users are fully responsible for respecting the terms of the these databases and for
citing them when required. The users can edit/update this dataset if needed.
}
\references{
Ceol, A., et al. (2007). DOMINO: a database of domain-peptide interactions. Nucleic
Acid Research, 35, D557-60.

Cui, Q., et al. (2007). A map of human cancer signaling. Molecular Systems Biology,
3:152.

Diella, F., et al. (2004). Phospho.ELM: a database of experimentally verified
phosphorylation sites in eukaryotic proteins. BMC Bioinformatics, 22, 5:79.

Dinkel, H., et al. (2012). ELM--the database of eukaryotic linear motifs. Nucleic
Acid Research, 40, D242-51.

Han, H., et al. (2015). TRRUST: a reference database of human transcriptional
regulatory interactions. Scientific Reports, 15, 11432.

Hornbeck, P.V., et al. (2012). PhosphoSitePlus: a comprehensive resource for investigating
the structure and function of experimentally determined post-translational modifications
in man and mouse. Nucleic Acid Research, 40, D261-70.

Korcsmaros, T., et al. (2010). Uniformly curated signaling pathways reveal
tissue-specific cross-talks and support drug target discovery. Bioinformatics,
26, 2042:2050.

Lynn, D.J., et al. (2008). InnateDB: facilitating systems-level analyses of the
mammalian innate immune response. Molecular Systems Biology, 4, 218.

Orchard, S., et al. (2014). The MIntAct project--IntAct as a common curation platform
for 11 molecular interaction databases. Nucleic Acid Research, 242, D358-63.

Pagel, P., et al. (2005). The MIPS mammalian protein-protein interaction database.
Bioinformatics, 21, 832-834.

Papp, D., et al. (2012). The NRF2-related interactome and regulome contain
multifunctional proteins and fine-tuned autoregulatory loops. FEBS Letters, 586,
1795-802.

Pawson, A.J., et al. (2014). The IUPHAR/BPS Guide to PHARMACOLOGY: an expert-driven
knowledgebase of drug targets and their ligands. Nucleic Acids Research, 42, D1098-106.

Peri, S., et al. (2003). Development of human protein reference database as an initial
platform for approaching systems biology in humans. Genome Research, 13, 2363-2371.

Turei, D., et al. (2015). Autophagy Regulatory Network - a systems-level bioinformatics
resource for studying the mechanism and regulation of autophagy. Autophagy, 11, 155-165.

Turei, D., et al. (2016). OmniPath: guidelines and gateway for literature-curated
signaling pathway resources. Nature methods, 13, 966-967.

Sarkar, D., et al. (2015). LMPID: a manually curated database of linear motifs
mediating protein-protein interactions. Database(Oxford), pii: bav014.

Shin, Y.C., et al. (2011). TRIP Database: a manually curated database of protein-protein
interactions for mammalian TRP channels. Nucleic Acids Research, 39, D356-61.

Snel, B., et al. (2000). STRING: a web-server to retrieve and display the repeatedly
occurring neighbourhood of a gene. Nucleic Acids Research, 28, 3442-3444.

Xenarios, I., et al. (2002). DIP, the Database of Interacting Proteins: a research
tool for studying cellular networks of protein interactions. Nucleic Acids Research,
30, 303-305.
}
\examples{
# Build regulatory network using the OmniPath dataset only
omnipath_net <- MS2_ppiNetwork(datasets = "omnipath")

# Build regulatory network using the TRRUST dataset only
trrust_net <- MS2_ppiNetwork(datasets = "trrust")

# Build regulatory network using interactions from STRING and BioGRID
biogridstring_net <- MS2_ppiNetwork(datasets = c("biogrid", "string"))
}
