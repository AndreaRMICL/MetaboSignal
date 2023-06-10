######################### TEST METABOSIGNAL ###################################

library(RUnit)

#### Test MS_keggNetwork ####

### Examples of situations were MS_keggNetwork gives an error

# Attempt to build a MetaboSignal network with pathIDs from different organisms
MS_DifferentOrganism = try (MS_keggNetwork(metabo_paths = "rno00010",
                                           signaling_paths = "hsa04922"),
                            silent = TRUE)

# Attempt to build a network with non-existent path IDs
MS_IncorrectPath = try (MS_keggNetwork(metabo_paths = "X"), silent = TRUE)

## Tests

test.MS_keggNetwork <- function() {

  # Non-existent path ID an error
  checkException(MS_DifferentOrganism, silent = TRUE)

  # Path IDs from different organisms produces an error
  checkException(MS_IncorrectPath, silent = TRUE)
}

test.MS_keggNetwork()

################################################################################
################################################################################

#### Test MS_distances ####

## Examples of distance_matrix generated with MS_distances

# Distance matrix generated using different modes ("SP","all","out")
distances_SP = MS_distances(MetaboSignal_table, organism_code = "rno",
                            mode = "SP" )

distances_all = MS_distances(MetaboSignal_table, organism_code = "rno",
                             mode = "all" )

distances_out = MS_distances(MetaboSignal_table, organism_code = "rno",
                             mode = "out")

# Distance matrix from all genes to one metabolite ("cpd:C00267")
distances_M = MS_distances(MetaboSignal_table, organism_code ="rno",
                                     target_metabolites = "cpd:C00267")

# Distance matrix from all genes to a set of metabolites, including both mapped
#("cpd:C00267") and non-mapped("X") metabolites
distances_IncorrectMetabo = MS_distances(MetaboSignal_table, organism_code = "rno",
                                         target_metabolites = c("cpd:C00267", "X"))

# Attempt to build a distance matrix when none of the target_metabolites can be
#onto the network
distances_IncorrectMetaboAll = try(MS_distances(MetaboSignal_table,organism_code ="rno",
                                                          target_metabolites = "X"),
                                   silent = TRUE)

## Tests

test.MS_distances <- function() {
  # MS_distances from all genes to all metabolites generates a matrix
  checkTrue(is.matrix(distances_SP) & is.matrix (distances_all) & is.matrix
            (distances_out))

  # MS_distances from all genes to one metabolite generates a matrix
  checkTrue(is.matrix(distances_M))

  # MS_distances ignores source_nodes or target_metabolites not-mapped
  #onto the network
  checkTrue(is.matrix(distances_IncorrectMetabo))
  checkTrue(identical(distances_M, distances_IncorrectMetabo))

  # MS_distances produces an error when none of the source_genes or
  #target_metabolites can be mapped onto the network
  checkException(distances_IncorrectMetaboAll, silent=TRUE)

}

test.MS_distances()

################################################################################
################################################################################

#### Test MS_shortestPathsNetwork ####

## Examples of subnetworks generated with MS_shortestPathsNetwork

subnet_1st = MS_shortestPathsNetwork (MetaboSignal_table, organism_code = "rno",
                                      source_nodes = c("K04456","K08074"),
                                      target_nodes = c("cpd:C00267"),
                                      type = "first", names = FALSE,
                                      export_cytoscape = FALSE)

subnet_bw = MS_shortestPathsNetwork (MetaboSignal_table, organism_code = "rno",
                                     source_nodes = c("K04456","K08074"),
                                     target_nodes = c("cpd:C00267"),
                                     type = "bw", names = FALSE,
                                     export_cytoscape = FALSE)

subnet_all = MS_shortestPathsNetwork (MetaboSignal_table, organism_code = "rno",
                                      source_nodes = c("K04456","K08074"),
                                      target_nodes = c("cpd:C00267"),
                                      type = "all", names = FALSE,
                                      export_cytoscape = FALSE)


## Tests
test.MS_shortestPathsNetwork <- function() {
  # MS_shortestPathsNetwork generates a matrix
  checkTrue(is.matrix(subnet_1st) & is.matrix(subnet_bw) & is.matrix(subnet_all))
}
test.MS_shortestPathsNetwork()

################################################################################
################################################################################

#### Test MS_topologyFilter ####

test.MS_topologyFilter <- function() {
  # MS_topologyFilter returns a matrix
  checkTrue(is.matrix(MS_topologyFilter(MetaboSignal_table, type = "all",
                                       target_node = "cpd:C00267", distance_th = 4,
                                       bw_th = 0.00005)))
}
test.MS_topologyFilter()


################################################################################
################################################################################

#### Test MS_findMappedNodes####

test.MS_findMappedNodes <- function() {
  # MS_findMappedNodes returns a list
  checkTrue(is.list(MS_findMappedNodes(nodes = c("cpd:C00267","cpd:C00021"),
                                       MetaboSignal_table)))
}
test.MS_findMappedNodes()

################################################################################
################################################################################

#### Test MS_shortestPaths####

test.MS_shortestPaths <- function() {
  # MS_shortestPaths returns a vector or a matrix
  checkTrue(is.vector(MS_shortestPaths(MetaboSignal_table, source_node = "K13952",
                                       target_node = "cpd:C00469")))
  checkTrue(is.matrix(MS_shortestPaths(MetaboSignal_table, source_node = "K13952",
                                       target_node = "cpd:C05125", mode ="SP",
                                       type = "all")))
}
test.MS_shortestPaths()

################################################################################
################################################################################

#### Test MS_NodeBW ####

test.MS_nodeBW <- function() {
  # MS_nodeBW returns a vector
  checkTrue(is.vector(MS_nodeBW(MetaboSignal_table)))
}
test.MS_nodeBW()

################################################################################
################################################################################

#### Test MS_removeNode ####

test.MS_removeNode <- function() {
  # MS_removeNode returns a matrix
  checkTrue(is.matrix(MS_removeNode(nodes = "K04465", MetaboSignal_table)))

  # MS_removeNode removes the edges containing the undesired node
  checkTrue(nrow(MS_removeNode(nodes = "K04465", MetaboSignal_table)) <
              nrow (MetaboSignal_table))
}
test.MS_removeNode()

################################################################################
################################################################################

#### Test MS_replaceNode ####

test.MS_replaceNode <- function() {
  # MS_RemoveNode returns a matrix
  checkTrue(is.matrix(MS_replaceNode(node1 = "cpd:C00267", node2 = "cpd:C00031",
                                     MetaboSignal_table)))
}
test.MS_replaceNode()


