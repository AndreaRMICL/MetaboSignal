#################### MetaboSignal_distances ###############
MS_distances = function(network_table, organism_code, mode = "SP",
                        source_genes = "all", target_metabolites = "all",
                        names = FALSE) {

    ## Check mode and type
    check_mode_type(mode = mode)

    ## Network_table must be a 3-column
    network_table = unique(network_table)
    check_matrix_v2(network_table, n = 2)

    ## Define all nodes of the network
    all_nodes = unique(as.vector(network_table[, 1:2]))

    ## Define gene and metabolite nodes of the network
    cpd_ind = grep("cpd:|dr:|gl:", all_nodes)

    if (length(cpd_ind) == length(all_nodes)) { ## There are not genes
        to_print=paste("Could not build distance matrix. It seems that the",
                       "network_table does not have genes", sep = " ")
        stop(to_print)
    }
    if (length(cpd_ind) == 0) {
      to_print=paste("Could not build distance matrix. It seems that the",
                     "network_table does not have compounds", sep = " ")
      stop(to_print)
    }
    genes = all_nodes[-cpd_ind]
    metabolites = all_nodes[cpd_ind]

    ## Build distance_matrix: from all nodes to all nodes.
    if (mode == "all") {
        MetaboSignal_network_i = graph.data.frame(network_table[, 1:2], directed = TRUE)
        distance_matrix = distances(MetaboSignal_network_i, mode = mode)
    } else {
        MetaboSignal_network_i = graph.data.frame(network_table[, 1:2], directed = TRUE)
        distance_matrix = distances(MetaboSignal_network_i, mode = "out")
    }

    ## Build distanceGM_matrix: from genes to metabolites##
    distanceGM_matrix = distance_matrix[genes, metabolites]
    distanceGM_matrix = matrix(distanceGM_matrix, ncol = length(metabolites),
                               nrow = length(genes))
    colnames(distanceGM_matrix) = metabolites
    rownames(distanceGM_matrix) = genes
    distanceGM_final = distanceGM_matrix

    ## Correct distances for SP mode
    if (mode == "SP") {
        DGM = distanceGM_matrix  # distance matrix from genes to metabolites
        D = distance_matrix  # distance matrix from all nodes to all nodes.
        DG = distance_matrix[genes, ]
        if (is.vector(DG)) { # DG needs to be in matrix
            Names = names(DG)
            DG = matrix(DG, ncol = ncol(distance_matrix))
            colnames(DG) = Names
        }
        DGM_corrected = DGM

        for (i in 1:ncol(DGM)) {
            metabolites = colnames(DGM)  ##all metabolites of the network_table
            metabolite = metabolites[i]
            all_nodes = colnames(D)  ## all nodes of the network_table

            all_substrates = network_table[, 1]
            index_substrates = which(all_substrates == metabolite)
            all_reactions = network_table[index_substrates, 2]

            if (length(index_substrates) >= 1) {# the metabolite is a substrate
                indexS = unlist(lapply(all_reactions, find_node_index, all_nodes))

                submatrix = DG[, c(indexS)]
                submatrix = submatrix + 1
                submatrix = cbind(submatrix, DGM[, i])

                min = apply(submatrix, 1, min)
                DGM_corrected[, i] = min
            }
        }
        distanceGM_final = DGM_corrected
    }

    ## Select submatrix#
    if (source_genes[1] == "all") {
        mapped_sources = genes
        unmapped_sources = NULL
    } else {
        source_genes = unique(gsub(" ", "", source_genes))  #remove potential white spaces
        mapped_sources = intersect(source_genes, genes)
        if (length(mapped_sources) == 0) {
            stop("None of the source_genes was mapped onto the network")
        }
        unmapped_sources = setdiff(source_genes, genes)
    }
    if (target_metabolites[1] == "all") {
        mapped_targets = metabolites
        unmapped_targets = NULL
    } else {
        target_metabolites = gsub(" ", "", target_metabolites)
        target_metabolites = unique(target_metabolites)
        mapped_targets = intersect(target_metabolites, metabolites)
        if (length(mapped_targets) == 0) {
            stop("None of the target_metabolites was mapped onto the network")
        }
        unmapped_targets = setdiff(target_metabolites, metabolites)
    }

    submatrixGM = distanceGM_final[as.character(mapped_sources), mapped_targets]
    submatrixGM = matrix(submatrixGM, ncol = length(mapped_targets),
                         nrow = length(mapped_sources))
    rownames(submatrixGM) = mapped_sources
    colnames(submatrixGM) = mapped_targets

    if (names == TRUE) {
      rownames(submatrixGM) = MS_changeNames(rownames(submatrixGM),
                                             organism_code)
      colnames(submatrixGM) = MS_changeNames(colnames(submatrixGM),
                                             organism_code)
    }

    all_unmapped = c(unmapped_sources, unmapped_targets)

    if (length(all_unmapped) > 0) {
        message("Note: some source_nodes or target_nodes were not mapped onto the network")
    }
    return(submatrixGM)
}
