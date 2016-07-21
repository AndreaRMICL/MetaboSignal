######################### TEST METABOSIGNAL ###################################

library(RUnit)

#### Test MetaboSignal_matrix ####

### Examples of situations were MetaboSignal_matrix gives an error

# Attempt to build a MetaboSignal network with pathIDs from different organisms
MS_DifferentOrganism = try (MetaboSignal_matrix(metabo_paths = "rno00010",
                                                signaling_paths = "hsa04922"), silent = TRUE)

# Attempt to build a network with non-existent path IDs
MS_IncorrectPath = try (MetaboSignal_matrix(metabo_paths = "X"), silent = TRUE)

## Tests

test.MetaboSignal_matrix <- function() {

  # Non-existent path ID an error
  checkException(MS_DifferentOrganism, silent = TRUE)

  # Path IDs from different organisms produces an error
  checkException(MS_IncorrectPath, silent = TRUE)
}

test.MetaboSignal_matrix()

################################################################################
################################################################################

#### Test MetaboSignal_distances ####

## Examples of distance_matrix generated with MetaboSignal_distances

# Distance matrix generated using different modes ("SP","all","out")
distances_SP = MetaboSignal_distances(MetaboSignal_table, organism_code = "rno",
                                      mode = "SP" )

distances_all = MetaboSignal_distances(MetaboSignal_table, organism_code = "rno",
                                       mode = "all" )

distances_out = MetaboSignal_distances(MetaboSignal_table, organism_code = "rno",
                                       mode = "out")

# Distance matrix from all genes to one metabolite ("cpd:C00267")
distances_M = MetaboSignal_distances(MetaboSignal_table, organism_code ="rno",
                                     target_metabolites = "cpd:C00267")

# Distance matrix from all genes to a set of metabolites, including both mapped
#("cpd:C00267") and non-mapped("X") metabolites
distances_IncorrectMetabo = MetaboSignal_distances(MetaboSignal_table, organism_code = "rno",
                                                   target_metabolites = c("cpd:C00267", "X"))

# Attempt to build a distance matrix when none of the target_metabolites can be
#onto the network
distances_IncorrectMetaboAll = try(MetaboSignal_distances(MetaboSignal_table,organism_code ="rno",
                                                          target_metabolites = "X"),
                                   silent = TRUE)

## Tests

test.MetaboSignal_distances <- function() {
  # MetaboSignal_distances from all genes to all metabolites generates a matrix
  checkTrue(is.matrix(distances_SP) & is.matrix (distances_all) & is.matrix
            (distances_out))

  # MetaboSignal_distances from all genes to one metabolite generates a matrix
  checkTrue(is.matrix(distances_M))

  # MetaboSignal_distances ignores source_nodes or target_metabolites not-mapped
  #onto the network
  checkTrue(is.matrix(distances_IncorrectMetabo))
  checkTrue(identical(distances_M, distances_IncorrectMetabo))

  # MetaboSignal_distances produces an error when none of the source_genes or
  #target_metabolites can be mapped onto the network
  checkException(distances_IncorrectMetaboAll, silent=TRUE)

}

test.MetaboSignal_distances()

################################################################################
################################################################################

#### Test MetaboSignal_NetworkCytoscape ####

## Examples of subnetworks generated with MetaboSignal_NetworkCytoscape

subnet_1st = MetaboSignal_NetworkCytoscape (MetaboSignal_table, organism_code = "rno",
                                            source_genes = c("K04456","K08074"),
                                            target_metabolites = c("cpd:C00267"),
                                            type = "first", names = FALSE,
                                            export_cytoscape = FALSE)

subnet_bw = MetaboSignal_NetworkCytoscape (MetaboSignal_table, organism_code = "rno",
                                           source_genes = c("K04456","K08074"),
                                           target_metabolites = c("cpd:C00267"),
                                           type = "bw", names = FALSE,
                                           export_cytoscape = FALSE)

subnet_all = MetaboSignal_NetworkCytoscape (MetaboSignal_table, organism_code = "rno",
                                            source_genes = c("K04456","K08074"),
                                            target_metabolites = c("cpd:C00267"),
                                            type = "all", names = FALSE,
                                            export_cytoscape = FALSE)


## Tests
test.MetaboSignal_NetworkCytoscape <- function() {
  # MetaboSignal_NetworkCytoscape generates a matrix
  checkTrue(is.matrix(subnet_1st) & is.matrix(subnet_bw) & is.matrix(subnet_all))
}
test.MetaboSignal_NetworkCytoscape()

################################################################################
################################################################################

#### Test MS_FilterNetwork####

test.MS_FilterNetwork <- function() {
  # MS_FilerNetwork returns a matrix
  checkTrue(is.matrix(MS_FilterNetwork(MetaboSignal_table, type = "all",
                                       target_node = "cpd:C00267", distance_th = 4,
                                       bw_th = 0.00005)))
}
test.MS_FilterNetwork()


################################################################################
################################################################################

#### Test MS_FindMappedNodes####

test.MS_FindMappedNodes <- function() {
  # MS_FindMappedNodes returns a list
  checkTrue(is.list(MS_FindMappedNodes(nodes = c("cpd:C00267","cpd:C00021"),
                                       MetaboSignal_table)))
}
test.MS_FindMappedNodes()

################################################################################
################################################################################

#### Test MS_GetShortestpaths####

test.MS_GetShortestpaths <- function() {
  # MS_GetShortestpaths returns a vector or a matrix
  checkTrue(is.vector(MS_GetShortestpaths(MetaboSignal_table, source_node = "K13952",
                                          target_node = "cpd:C00469")))
  checkTrue(is.matrix(MS_GetShortestpaths(MetaboSignal_table, source_node = "K13952",
                                          target_node = "cpd:C05125", mode ="all",
                                          type = "all")))
}
test.MS_GetShortestpaths()

################################################################################
################################################################################

#### Test MS_NodeBW ####

test.MS_NodeBW <- function() {
  # MS_NodeBW returns a vector
  checkTrue(is.vector(MS_NodeBW(MetaboSignal_table)))
}
test.MS_NodeBW()

################################################################################
################################################################################

#### Test MS_RemoveNode ####

test.MS_RemoveNode <- function() {
  # MS_RemoveNode returns a matrix
  checkTrue(is.matrix(MS_RemoveNode(nodes = "K04465", MetaboSignal_table)))

  # MS_RemoveNode removes the edges containing the undesired node
  checkTrue(nrow(MS_RemoveNode(nodes = "K04465", MetaboSignal_table)) <
              nrow (MetaboSignal_table))
}
test.MS_RemoveNode()

################################################################################
################################################################################

#### Test MS_ReplaceNode ####

test.MS_ReplaceNode <- function() {
  # MS_RemoveNode returns a matrix
  checkTrue(is.matrix(MS_ReplaceNode(node1 = "cpd:C00267", node2 = "cpd:C00031",
                                     MetaboSignal_table)))
}
test.MS_ReplaceNode()


