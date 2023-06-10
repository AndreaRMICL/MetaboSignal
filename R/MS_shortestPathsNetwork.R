############### gene_SP_v2 #############
# This function generates a shortest-path (SP) subnetwork from a gene to a metabolite
gene_SP_v2 = function(gene, metabolite, network_i, network_table,
                      distance_table, distance_th, mode, type, BW_matrix,
                      networkBW_i) {

    index_gene = which(rownames(distance_table) == gene)
    index_metabolite = which(colnames(distance_table) == metabolite)
    distanceGM = distance_table[index_gene, index_metabolite]

    if (distanceGM < distance_th) {
        # If distance is not Inf,

        if (type == "first") {
            if (mode == "all") {
                ASP = get.shortest.paths(network_i, gene, metabolite,
                                         mode = "all")
            } else {
                ASP = get.shortest.paths(network_i, gene, metabolite,
                                         mode = "out")
            }
            ASP = ASP[[1]]
        } else {
            if (mode == "all") {
                ASP = get.all.shortest.paths(network_i, gene,
                                             metabolite, mode = "all")

            } else {
                ASP = get.all.shortest.paths(network_i, gene,
                                             metabolite, mode = "out")
            }
            ASP = ASP$res
        }

        all_paths = lapply(ASP, ASP_paths)
        all_paths = unique(do.call(rbind, all_paths))
        rownames(all_paths) = NULL

        ## Select bw-ranked shortest path of multiple shortest paths
        if (type == "bw" & nrow(all_paths) > 1) {
            all_paths = BW_ranked_SP(all_paths, BW_matrix, networkBW_i, mode)
        }

        ## Transform the shortest path matrix into a network table
        all_pathsList = split(all_paths, row(all_paths))
        all_pathsnetworkList = lapply(all_pathsList, path_as_network)
        subnetwork = unique(do.call(rbind, all_pathsnetworkList))
        rownames(subnetwork) = NULL  # subnetwork is a network-table

        ## Correct directionality in case of SP mode
        if (mode == "SP") {
            for (i in 1:nrow(subnetwork)) {
                all_edges = paste(network_table[, 1], network_table[, 2], sep = "_")
                my_edge = paste(subnetwork[i, 1], subnetwork[i, 2], sep = "_")

                if(length(intersect(my_edge, all_edges)) == 0) {
                    subnetwork[i, ] = rev(subnetwork[i, ])
                }
            }
        }

        ## Add edges for reversible interactions
        reverseSubnetwork = cbind(subnetwork[, 2], subnetwork[, 1])
        if (mode == "all") {
            # Build an undirected network
            pathsGM = rbind(subnetwork, reverseSubnetwork)
            pathsGM = unique(pathsGM)
            rownames(pathsGM) = NULL
        } else {
            rev_edges = paste(reverseSubnetwork[, 1], reverseSubnetwork[, 2], sep = "_")
            all_edges = paste(network_table[, 1], network_table[, 2], sep = "_")
            edges_toadd = intersect(rev_edges, all_edges)

            pathsGM = subnetwork

            if(length(edges_toadd) > 0) {
                edges_toaddM = do.call(rbind, strsplit(edges_toadd, "_"))
                pathsGM = rbind(pathsGM, edges_toaddM)

            }
            pathsGM = unique(pathsGM)
            rownames(pathsGM) = NULL
        }
    } else {
        pathsGM = NULL
    }
    return(pathsGM)
}

############### metabolite_SP #############
metabolite_SP = function(metabolite, source_genes, network_table, distance_table,
                         distance_th, mode, type, BW_matrix, networkBW_i) {
    if (mode == "SP") {
        network_i = igraph_edited(network_table, metabolite)
    } else {
        network_i = graph.data.frame(network_table, directed = TRUE)
    }
    metabo_pathsGM = lapply(source_genes, gene_SP_v2, metabolite, network_i,
                            network_table, distance_table, distance_th, mode, type,
                            BW_matrix, networkBW_i)
    metabo_pathsGM = unique(do.call(rbind, metabo_pathsGM))
    return(metabo_pathsGM)
}

####################### map_interaction_type #######################
map_interaction_type = function (st, network) {
    collapsed_edges = paste(network[, 1], network[, 2], sep = "_")
    ind = which(collapsed_edges == st)
    new_edge = c(unlist(strsplit(st, "_")), network[ind, 3])
    return(new_edge)
}

############### MS_shortestPathsNetwork #############
MS_shortestPathsNetwork = function(network_table, organism_code, source_nodes,
                                   target_nodes, mode = "out", type = "first",
                                   distance_th = Inf, names = TRUE,
                                   export_cytoscape = TRUE, file_name = "MS") {

    ## Check mode
    if((mode %in% c("SP", "out")) == FALSE) {
        stop ("mode should be SP or out")
    }

    ## Check type
    if((type %in% c("first", "all", "bw")) == FALSE) {
        stop ("type should be first, all or bw")
    }

    ## Check that network is correct
    network = unique(network_table)
    check_matrix_v2(network, n = 3)

    ## Make source_nodes and target_nodes unique
    source_nodes = unique(source_nodes)
    target_nodes = unique(target_nodes)

    ## Iniciate node bw matrix in case it's needed
    BW_matrix = matrix(c("Node", 1), ncol = 2)
    colnames(BW_matrix) = c("node", "bw")
    networkBW_i = graph.data.frame(network, directed = TRUE)

    message("Calculating distances", "\n")

    if (mode == "SP") {
        ## Mode SP assummes that all target nodes must be metabolites
        ans_cpd_target = grep("cpd:", target_nodes)
        if (length(ans_cpd_target) < length(target_nodes)) {
            stop("All target nodes must be compounds to use the SP mode")
        }
        ans_cpd_source = grep("cpd", source_nodes)

        if(length(ans_cpd_source) > 0) {
            stop("All source nodes must be genes to use the SP mode")
        }
        distance_table = MS_distances(network, organism_code = organism_code,
                                      mode = mode, names = FALSE)
    } else {
        MetaboSignal_network_i = graph.data.frame(network[ , 1:2], directed = TRUE)
        distance_table = distances(MetaboSignal_network_i, mode = mode)
    }

    ## Check if the source_genes and target_metabolites can be
    ## mapped onto the network#
    all_nodes = unique(as.vector(network[, 1:2]))

    mapped_sources = intersect(source_nodes, all_nodes)
    non_mapped_sources = setdiff(source_nodes, all_nodes)
    mapped_targets= intersect(target_nodes, all_nodes)
    non_mapped_targets = setdiff(target_nodes, all_nodes)

    if(length(mapped_sources) == 0) {
        stop ("None of the source nodes was mapped onto the network")
    }
    if(length(mapped_targets) == 0) {
        stop ("None of the target nodes was mapped onto the network")
    }

    ## Calculate shortest path network
    message("Building shortest path network", "\n")

    all_pathsGM = lapply(mapped_targets, metabolite_SP, mapped_sources,
                         network_table = network[, 1:2], distance_table,
                         distance_th, mode, type,
                         BW_matrix, networkBW_i)
    all_pathsGM = unique(do.call(rbind, all_pathsGM))

    if (length(all_pathsGM) == 0) {
        to_print = paste("Could not connect any of the source_nodes with any",
                         "of the target_nodes")
        stop(to_print)
    }

    st_edges = paste(all_pathsGM[, 1], all_pathsGM[, 2], sep = "_")
    all_pathsGM = do.call(rbind, lapply(st_edges, map_interaction_type, network))
    colnames(all_pathsGM) = c("source", "target", "type")

    if (export_cytoscape == TRUE) {
        message("Creating cytoscape files")
        message()
        all_targets = unique(c(mapped_sources, mapped_targets))
        net_cyto = MS_exportCytoscape(all_pathsGM, organism_code = organism_code,
                                      names = names, file_name = file_name,
                                      targets = all_targets)
    }

    ## Final report
    non_mapped_all = c(non_mapped_sources, non_mapped_targets)

    if(length(non_mapped_all) > 0) {
      report_message = paste("Note that some source_nodes or some target_nodes",
                             "were not mapped onto the network")
      #message(report_message)
      #message()
    } else {
        report_message = paste("Note that all source_nodes and target_nodes",
                               "were successfully mapped onto the network")
    }

    if (names == FALSE) {
        message(report_message)
        message()
        return(all_pathsGM)
    }

    all_pathsGM_names = network_names(all_pathsGM, organism_code)
    message(report_message)
    message()
    return(all_pathsGM_names)
}
