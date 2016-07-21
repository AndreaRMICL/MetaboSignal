#################### EXTERNAL FUNCTION ####################

#################### MetaboSignal_distances ###############

MetaboSignal_distances = function(network_table, organism_code,
                                  organism_name, mode = "SP", source_genes = "all",
                                  target_metabolites = "all", names = FALSE) {

    ## Check mode and type
    check_mode_type(mode = mode)

    ## Force network_table to be a unique 2-column matrix
    network_table = check_matrix(network_table)

    ## Define all nodes of the network
    all_nodes = unique(as.vector(network_table))

    ## Define gene nodes of the network
    genes = c()
    for (node in all_nodes) {
        if (grepl("rn", node) == "TRUE" | grepl("K", node) == "TRUE" |
            grepl(organism_code, node) == "TRUE") {
            genes = c(genes, node)
        }
    }
    if (length(genes) == 0) {
        # This could because the network_table has specific organism
        #IDs (e.g. 'rno:') that are not consistent with the
        #organism_code(e.g.'mmu')
        to_print=paste("Could not build distance matrix. It might be that",
                       "the network_table does not have gene nodes or that the",
                       "organism_code is incorrect")
        stop(to_print)
    }

    ## Define metabolite nodes of the network
    metabolites = c()
    for (node in all_nodes) {
        if (grepl("cpd", node) == "TRUE") {
            metabolites = c(metabolites, node)
        }
    }
    if (length(metabolites) == 0) {
        stop("The input network does not have metabolite nodes")
    }

    ## Build distance_matrix: from all nodes to all nodes.
    if (mode == "all") {
        MetaboSignal_network_i = graph.data.frame(network_table, directed = TRUE)
        distance_matrix = distances(MetaboSignal_network_i, mode = mode)
    } else {
        MetaboSignal_network_i = graph.data.frame(network_table, directed = TRUE)
        distance_matrix = distances(MetaboSignal_network_i, mode = "out")
    }

    ## Build distanceGM_matrix: from genes to metabolites##
    metabolite_index = c()
    for (metabolite in metabolites) {
        index = which(colnames(distance_matrix) == metabolite)
        metabolite_index = c(metabolite_index, index)
    }
    gene_index = c()
    for (gene in genes) {
        index = which(rownames(distance_matrix) == gene)
        gene_index = c(gene_index, index)
    }
    distanceGM_matrix = distance_matrix[c(gene_index), c(metabolite_index)]
    distanceGM_matrix = matrix(distanceGM_matrix,
                               ncol = length(metabolite_index),
                               nrow = length(gene_index))
    colnames(distanceGM_matrix) = metabolites
    rownames(distanceGM_matrix) = genes
    distanceGM_final = distanceGM_matrix

    ## Correct distances for SP mode
    if (mode == "SP") {
        DGM = distanceGM_matrix  # distance matrix from genes to metabolites
        D = distance_matrix  # distance matrix from all nodes to all nodes.
        DG = distance_matrix[c(gene_index), ]
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
                indexS = c()
                for (s in 1:length(all_reactions)) {
                  indexS = c(indexS, which(all_nodes == all_reactions[s]))
                }
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
    all_genes = rownames(distanceGM_final)
    all_metabolites = colnames(distanceGM_final)
    source_genes = gsub(" ", "", source_genes)  #remove potential white spaces
    source_genes = unique(source_genes)
    target_metabolites = gsub(" ", "", target_metabolites)
    target_metabolites = unique(target_metabolites)
    answer_genes = source_genes[1] == "all"
    answer_metabolites = target_metabolites[1] == "all"

    if (answer_genes == TRUE) {
        source_genes = all_genes
        source_genes_backup = source_genes
    } else {
        # Find out if the genes represent KEGG orthology IDs
        ortho_answer = unique(grepl(organism_code, all_genes))
        ortho_answer = TRUE %in% ortho_answer

        if (ortho_answer == TRUE) {# Network genes represent KEGG specific IDs
            # Check if source genes are already organism specific KEGG IDs
            if (grepl(organism_code, source_genes[1]) == FALSE) {
                source_genes_backup = source_genes
                source_genes = try(MS_GetKEGG_GeneID(genes = source_genes,
                                                     organism_name = organism_name,
                                                     orthology = FALSE),
                                   silent = TRUE)
            } else {
                source_genes_backup = source_genes
            }
        } else {# Network genes represent KEGG orthology IDs
            # Check if source genes are already orthology IDs
            if ((substr(source_genes[1], 1, 1) == "K" & nchar(source_genes[1])
                 == 6) == FALSE) {
                source_genes_backup = source_genes
                source_genes = try(MS_GetKEGG_GeneID(genes = source_genes,
                                                     organism_name = organism_name,
                                                     orthology = TRUE),
                                   silent = TRUE)
            } else {
                source_genes_backup = source_genes
            }
        }
    }
    if (answer_metabolites == TRUE) {
        target_metabolites = all_metabolites
    } else {
        target_metabolites_backup = target_metabolites
    }

    if(grepl("Error",source_genes)[1] == TRUE){
      stop (" None of the genes is reported in KEGG. Could be that the organism_name is incorrect")
    }

    index_genes = c()
    for (g in 1:length(source_genes)) {
        if (source_genes[g] %in% all_genes == TRUE) {
            index = which(all_genes == source_genes[g])
            index_genes = c(index_genes, index)
        }
    }
    index_metabolites = c()
    for (m in 1:length(target_metabolites)) {
        target_metabolites_backup = target_metabolites

        if (target_metabolites[m] %in% all_metabolites == TRUE) {
            index = which(all_metabolites == target_metabolites[m])
            index_metabolites = c(index_metabolites, index)
        }
    }
    if (length(index_metabolites) >= 1 & length(index_genes) >= 1) {
        if (answer_genes == FALSE) {
          message ("Building distance matrix", "\n")
        }
        gene_names = rownames(distanceGM_final)
        metabolite_names = colnames(distanceGM_final)
        distanceGM_final = distanceGM_final[index_genes, index_metabolites]

        if (is.vector(distanceGM_final) == TRUE) {
            distanceGM_final = matrix(distanceGM_final,
                                      nrow = length(index_genes),
                                      ncol = length(index_metabolites))
            colnames(distanceGM_final) = metabolite_names[index_metabolites]
            rownames(distanceGM_final) = gene_names[index_genes]
        }

        if (names == TRUE) {
            rownames(distanceGM_final) =
                MS_ChangeNames(rownames(distanceGM_final),organism_code)
            colnames(distanceGM_final) =
                MS_ChangeNames(colnames(distanceGM_final),organism_code)
        }

        ## Final report
        check_unmapped(target_metabolites_backup, index_metabolites,
                       source_genes_backup, index_genes)

        return(distanceGM_final)

    } else {
        to_print=paste("Impossible to build a distance matrix with these genes",
                       "and/or metabolites. None of the",
                       "source_genes and/or none of the target_metabolites was",
                       "mapped onto the network")
        stop(to_print)
    }
}
