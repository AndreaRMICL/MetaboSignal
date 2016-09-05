#################### MS_ChangeNames ########################

MS_ChangeNames = function(nodes, organism_code) {

    ## Check if there are compound nodes##
    if (length(grep("cpd:", nodes) > 0)) {
        # Get compounds table
        file = "http://rest.kegg.jp/list/compound"
        response = getURL(file)
        compoundM = convertTable(response)
    }

    for (i in seq_along(nodes)) {
        node = nodes[i]  # Backup for last else of the loop
        if (grepl("cpd:", node) == TRUE) {
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
        } else {
            nodes[i] = node
        }
    }
    return(nodes)
}

#################### MS_GetShortestpaths ######################

MS_GetShortestpaths = function(network_table, source_node, target_node,
                               mode = "SP", type = "first") {

    ## Force network_table to be a unique 2-column matrix
    network_table = check_matrix(network_table)

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
                  path = as.character(BW_ranked_SP(all_paths, BW_matrix, networkBW_i, mode))
                }
            }
            output_path = path
        }
        return(output_path)
    }
}

#################### MS_FilterNetwork ####################

MS_FilterNetwork = function(network_table, mode = "all", type, target_node = NULL,
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

    # Force network_table to be a unique 2-column matrix
    network_table = check_matrix(network_table)

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
        network_table2 = network_table[-c(edges_unwanted), ]
        network_table2 = matrix(network_table2, ncol = 2)
        colnames(network_table2) = c("node1", "node2")
    } else {
      network_table2 = network_table
      warning ("Filtering was ignored: check filtering parameters")
    }
    colnames(network_table2) = c("node1", "node2")

    ## Report features
    network_features(network_table2)

    return(network_table2)

}

#################### MS_NodeBW ####################

MS_NodeBW = function(network_table, mode = "all", normalized = TRUE) {

    ## Check mode and type
    check_mode_type(mode2 = mode)

    ## Force network_table to be a unique 2-column matrix
    network_table = check_matrix(network_table)

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


#################### MS_ToCytoscape ####################

MS_ToCytoscape = function(network_table, organism_code, names = TRUE,
                          target_nodes = NULL, file_name = "Cytoscape") {

    ## Force network_table to be a unique 2-column matrix
    network_table = check_matrix(network_table)

    ## Build cytoscape network with interaction based on edge directionality
    network_tableCytoscape = cbind(network_table[, 1], rep("interaction",
        nrow(network_table)), network_table[, 2])
    #The column interaction will be replace by reversible or irreversible.
    rows = split(network_tableCytoscape, row(network_tableCytoscape))

    index_to_include = rep(1, times = length(rows))
    for (i in 1:nrow(network_tableCytoscape)) {
        a = c(network_tableCytoscape[i, 1], network_tableCytoscape[i, 2],
              network_tableCytoscape[i, 3])
        a_2 = rev(a)

        if (Match(rows, a_2) == TRUE) {
            network_tableCytoscape[i, 2] = "reversible"
            index = as.numeric(which(network_tableCytoscape[, 1] == a_2[1]
                & network_tableCytoscape[, 2] == a_2[2]
                & network_tableCytoscape[, 3] == a_2[3], arr.ind = FALSE))
            index_to_include[index] = 0
        } else (network_tableCytoscape[i, 2] = "irreversible")
    }
    network_table_interactions = network_tableCytoscape
    index_to_remove = which (index_to_include == 0)

    if (length(index_to_remove) >= 1) {
        network_table_interactions = network_table_interactions[-index_to_remove, ]
        network_table_interactions = matrix(network_table_interactions, ncol = 3)
        # In case the network originally had only two edges that were reversible.
        # we force vector to matrix
    }
    colnames(network_table_interactions) = c("source_node", "interaction",
                                             "target_node")
    rownames(network_table_interactions) = NULL

    if (names == TRUE) {
        all_nodes_network = unique(as.vector(network_table_interactions))
        all_nodes_names = MS_ChangeNames(all_nodes_network, organism_code)
        for (i in seq_along(all_nodes_names)) {
            network_table_interactions[network_table_interactions ==
                all_nodes_network[i]] = all_nodes_names[i]
        }
    }
    cytoscape = as.data.frame(network_table_interactions, rownames = NULL)
    file_nameN = paste(file_name, "Network.txt", sep = "")

    write.table(cytoscape, file_nameN, row.names = FALSE, sep = "\t",
                quote = FALSE, col.names = FALSE)

    ## Create network attributes
    nodes = unique(as.vector(network_table))
    node_type_all = as.character(sapply (nodes, get_molecule_type,
                                         organism_code = organism_code))

    if (names == TRUE) {
        for (i in seq_along(all_nodes_names)) {
            nodes[nodes == all_nodes_network[i]] = all_nodes_names[i]
        }
    }
    node_typeM = cbind(nodes, node_type_all)
    node_typeM = as.data.frame(node_typeM)
    file_nameAtype = paste(file_name, "AttributesType.txt", sep = "")

    write.table(node_typeM, file_nameAtype, row.names = FALSE, sep = "\t",
                quote = FALSE, col.names = FALSE)

    ## Create node attribute to HL genes of interest
    if (length(target_nodes) > 0) {
        nodes = unique(as.vector(network_table))
        node_HL_all = as.character(sapply (nodes, get_target_type,
                                           target_nodes = target_nodes))

        if (names == TRUE) {
            for (i in seq_along(all_nodes_names)) {
                nodes[nodes == all_nodes_network[i]] = all_nodes_names[i]
            }
        }
        node_HLM = cbind(nodes, node_HL_all)
        node_HLM = as.data.frame(node_HLM)
        file_nameAtarget = paste(file_name, "AttributesTarget.txt", sep = "")

        write.table(node_HLM, file_nameAtarget, row.names = FALSE, sep = "\t",
                    quote = FALSE, col.names = FALSE)
    }
    return(cytoscape)
}

