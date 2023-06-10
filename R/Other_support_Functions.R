#################### MS_changeNames ########################
MS_changeNames = function(nodes, organism_code) {

    compoundM = c()
    ## Check if there are compound nodes##
    if (length(grep("cpd:", nodes) > 0)) {
        # Get compounds table
        file = "https://rest.kegg.jp/list/compound"
        response = getURL(file)
        metaboliteM = convertTable(response)
        compoundM = rbind(compoundM, metaboliteM)
    }
    ## Check if there are drug nodes##
    if (length(grep("dr:", nodes) > 0)) {
        # Get drugs table
        file = "https://rest.kegg.jp/list/drug"
        response = getURL(file)
        drugM = suppressWarnings(convertTable(response))
        compoundM = rbind(compoundM, drugM)
    }
    ## Check if there are glycan compounds
    if (length(grep("gl:", nodes) > 0)) {
        # Get glycan table
        file = "https://rest.kegg.jp/list/glycan"
        response = getURL(file)
        glycanM = convertTable(response)
        compoundM = rbind(compoundM, glycanM)
    }

    for (i in seq_along(nodes)) {
        node = nodes[i]  # Backup for last else of the loop
        if (grepl("cpd:|gl:|dr:", node) == TRUE) {
            index = which(compoundM[, 1] == nodes[i])
            if (length(index) > 0) {
                name = compoundM[index, 2]
                name = unlist(strsplit(name, "[;]"))
                #name = unlist(strsplit(name, "[,]"))
                nodes[i] = gsub(" ", "-", name[1])
            }
        } else if ((grepl("K", node) == TRUE | grepl(organism_code,
                                                     node) == TRUE)) {
            nodes[i] = From_geneID_to_symbol(node)

            if (nchar(nodes[i]) > 12 & grepl("", node) == FALSE) {
                nodes[i] = paste("LOC", substr(node, 5, nchar(node)), sep = "")
            }

            if (nchar(nodes[i]) > 12 & grepl("", node) == TRUE) {
                nodes[i] = node  #example dme:Dmel_CG10000
            }
        } else if (organism_code == "hsa" & !is.na (as.numeric(node))) {
            node_ent = mapIds(EnsDb.Hsapiens.v75, keys = node,
                              keytype = "ENTREZID", column = "SYMBOL")
            if (length(node_ent) == 0) {
                nodes[i] = node
            } else if (is.na(node_ent)) {
                nodes[i] = node
            } else {
                nodes[i] = node_ent
            }
        } else {
            nodes[i] = node
        }
    }
    return(nodes)
}

#################### MS_shortestPaths ######################
MS_shortestPaths = function(network_table, source_node, target_node,
                            mode = "out", type = "first") {

    ## Check that source_node and target_node are different
    if (source_node == target_node) {
        stop ("source_node and target_not must be different")
    }

    ## Network_table must be a matrix with at least two columns
    check_matrix_v2(network_table, n = 2)
    network_table = network_table[, 1:2]

    ## Check mode and type
    check_mode_type(mode = mode, type = type)

    ## Iniciate node bw matrix in case it's needed#
    BW_matrix = matrix(c("Node", 1), ncol = 2)
    colnames(BW_matrix) = c("node", "bw")
    networkBW_i = graph.data.frame(network_table, directed = TRUE)
    # SP mode does not make sense for bw.

    all_nodes = unique(as.vector(network_table))
    answer1 = source_node %in% all_nodes
    answer2 = target_node %in% all_nodes
    global_answer = unique(c(answer1, answer2))

    if ("FALSE" %in% global_answer) {
        stop("Source_node and/or target_node is not present in the network")
    } else {
        ## Get shortest paths
        if (mode == "SP") {
            network_i = igraph_edited(network_table, target_node)
        } else {
            network_i = graph.data.frame(network_table, directed = TRUE)
        }
        if (type == "first") {# Gets one shortest path
            if (mode == "all") {
                ASP = get.shortest.paths(network_i, source_node, target_node,
                                         mode = "all")
            } else {
                ASP = get.shortest.paths(network_i, source_node, target_node,
                                         mode = "out")
            }
            ASP = ASP[[1]]

        } else {# Gets all shortest paths (needed for both mode="all" or "bw"
            if (mode == "all") {
                ASP = get.all.shortest.paths(network_i, source_node,
                                             target_node, mode = "all")
            } else {
                ASP = get.all.shortest.paths(network_i, source_node,
                                             target_node, mode = "out")
            }
            ASP = ASP$res
        }

        all_paths = lapply(ASP, ASP_paths)
        all_paths = unique(do.call (rbind, all_paths))
        rownames(all_paths) = NULL

        if (type == "all") {# Selects all shortest paths
            output_path = all_paths

        } else {# Selects only one shortest path (bw-ranked or first)
            answer_vector = is.vector(all_paths)
            #Check if there is only one path or several paths
            if (answer_vector == TRUE) {
                path = as.character(all_paths)
            } else {
                if (type == "first") {
                    path = as.character(all_paths[1, ])
                } else { # type='bw'
                    path = as.character(BW_ranked_SP(all_paths, BW_matrix,
                                                     networkBW_i, mode))
                }
            }
            output_path = path
        }

        ## Path as network
        if(length(output_path) == 0) {
            return(NULL)
        }

        return(output_path)
    }
}

#################### MS_topologyFilter ####################
MS_topologyFilter = function(network_table, mode = "all", type, target_node = NULL,
                            distance_th, bw_th) {

    # Check mode and type
    check_mode_type(mode2 = mode, type2 = type)

    # Check if type is consistent
    if (type != "bw" & length(target_node) == 0) {
      stop("target_node is missing")
    }

    if (type != "bw" & length(target_node) > 0) {
      if (length(target_node) > 1) { # target_node must have length 1
        target_node = target_node[1]
        warning ("target_node has length > 1 and only the first element was used")
      }
      if (target_node %in% as.vector(network_table) == FALSE) {
      stop("target_node is not in the network")
      }
    }

    # Check network_table
    check_matrix_v2(network_table, n = 2)
    network_tableBU = network_table # in case it comes from MS2
    network_table = network_table[, 1:2] # in case it comes from MS2

    network_i = graph.data.frame(network_table, directed = TRUE)
    all_nodes = unique(as.vector(network_table))

    if (type != "bw") {
      distance_matrix = distances(network_i, mode = mode)
      index_target = which(colnames(distance_matrix) == target_node)
      Distances = distance_matrix[, index_target]
      Distances = as.numeric(Distances)
    }

    if (type != "distance") {
        if (mode == "all") {
            bw = as.numeric(betweenness(network_i, all_nodes, directed = FALSE,
                                        weights = NULL, nobigint = TRUE,
                                        normalized = TRUE))
        } else {
            bw = betweenness(network_i, all_nodes, directed = TRUE,
                             weights = NULL, nobigint = TRUE, normalized = TRUE)
        }
    }

    edges_response = vector(mode = "character", length = nrow(network_table))
    for (i in 1:nrow(network_table)) {
        node1 = as.character(network_table[i, 1])
        node2 = as.character(network_table[i, 2])

        if(type == "distance") {
            distance1 = Distances[which(rownames(distance_matrix) == node1)]
            distance2 = Distances[which(rownames(distance_matrix) == node2)]
        } else if (type == "bw") {
            bw1 = bw[which(all_nodes == node1)]
            bw2 = bw[which(all_nodes == node2)]
        } else { # type = "all"
            distance1 = Distances[which(rownames(distance_matrix) == node1)]
            distance2 = Distances[which(rownames(distance_matrix) == node2)]
            bw1 = bw[which(all_nodes == node1)]
            bw2 = bw[which(all_nodes == node2)]
        }

        if (type == "all") {
            if (distance1 > distance_th | distance2 > distance_th |
                bw1 < bw_th | bw2 < bw_th) {
                edges_response[i] = "unwanted"
            } else {
                edges_response[i] = "wanted"
            }
        } else if (type == "distance") {
            if (distance1 > distance_th | distance2 > distance_th) {
                edges_response[i] = "unwanted"
            } else {
                edges_response[i] = "wanted"
            }
        } else {
            if (bw1 < bw_th | bw2 < bw_th) {
                edges_response[i] = "unwanted"
            } else {
                edges_response[i] = "wanted"
            }
        }
    }

    edges_unwanted = which(edges_response == "unwanted")

    if (length(edges_unwanted) >= 1) {
        network_table2 = network_tableBU[-c(edges_unwanted), ]
        network_table2 = matrix(network_table2, ncol = ncol(network_tableBU))
        colnames(network_table2) = colnames(network_tableBU)
    } else {
      network_table2 = network_table
      warning ("Filtering was ignored: check filtering parameters")
    }

    ## Report features
    network_features(network_table2[, 1:2])

    return(network_table2)

}

#################### MS_nodeBW ####################
MS_nodeBW = function(network_table, mode = "all", normalized = TRUE) {

    ## Check mode and type
    check_mode_type(mode2 = mode)

    ## Check network_table
    check_matrix_v2(network_table, n = 2)
    network_table = network_table[, 1:2] # in case it comes from MS2.

    network_i = graph.data.frame(network_table, directed = TRUE)

    nodes = unique(as.vector(network_table))
    if (mode == "all") {
        BW = betweenness(network_i, nodes, directed = FALSE,
            weights = NULL, nobigint = TRUE, normalized = normalized)
    } else {
        BW = betweenness(network_i, nodes, directed = TRUE, weights = NULL,
             nobigint = TRUE, normalized = normalized)
    }
    bw = as.numeric(BW)
    hist(bw, col = "gray", main = paste("Node betweeness distribution"))
    return(BW)
}
