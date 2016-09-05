############### gene_SP #############
# This function generates a shortest-path (SP) subnetwork from a gene to a metabolite
gene_SP = function(gene, metabolite, network_i, network_table,
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
                supernetworkDirection = rbind(subnetwork, network_table)
                intersection_subnetwork = supernetworkDirection[duplicated(supernetworkDirection),]
                intersection_subnetwork = matrix(intersection_subnetwork, ncol = 2)
                rows = split(intersection_subnetwork, row(intersection_subnetwork))
                edge = subnetwork[i, ]

                if (Match(rows, edge) == FALSE) {
                  ## This will correct for the directionality adjustment.  For
                  ## example if the target node it's the sustrate of the
                  ## irreversible reaction, it should be ploted as substrate not
                  ## as product.
                  subnetwork[i, ] = rev(edge)
                }
            }
        }

        ## Add edges for reversible interactions
        reverseSubnetwork = cbind(subnetwork[, 2], subnetwork[, 1])
        supernetwork = rbind(reverseSubnetwork, network_table)
        if (mode == "all") {
            # Build an undirected network
            pathsGM = rbind(subnetwork, reverseSubnetwork)
            pathsGM = unique(pathsGM)
            rownames(pathsGM) = NULL
        } else {
            intersection_reverseSubnetwork = supernetwork[duplicated(supernetwork),]
            intersection_reverseSubnetwork = matrix(intersection_reverseSubnetwork,
                ncol = 2)

            if (nrow(intersection_reverseSubnetwork) > 0) {
                subnetwork = rbind(subnetwork, intersection_reverseSubnetwork)
            }
            pathsGM = unique(subnetwork)
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
    metabo_pathsGM = lapply(source_genes, gene_SP, metabolite, network_i,
                            network_table, distance_table, distance_th, mode, type,
                            BW_matrix, networkBW_i)
    metabo_pathsGM = unique(do.call(rbind, metabo_pathsGM))
    return(metabo_pathsGM)
}

############### MetaboSignal_NetworkCytoscape #############

MetaboSignal_NetworkCytoscape = function(network_table, organism_code,
    organism_name, source_genes, target_metabolites, mode = "SP",
    type = "first", distance_th = Inf, collapse_genes = FALSE,
    names = TRUE, export_cytoscape = TRUE, file_name = "Cytoscape") {

    ## Check mode and type
    check_mode_type(mode = mode, type = type)

    ## Force network_table to be a unique 2-column matrix
    network_table = check_matrix(network_table)

    ## Iniciate node bw matrix in case it's needed
    BW_matrix = matrix(c("Node", 1), ncol = 2)
    colnames(BW_matrix) = c("node", "bw")
    networkBW_i = graph.data.frame(network_table, directed = TRUE)

    ## Collapse the network by orthology
    if (collapse_genes == TRUE) {
        distance_table = MetaboSignal_distances(network_table,
            organism_code, organism_name, mode = mode, names = FALSE)
        all_genes = rownames(distance_table)
        ortho_answer = unique(grepl(organism_code, all_genes))
        ortho_answer = TRUE %in% ortho_answer

        if (ortho_answer == FALSE) {
            to_print = paste("Collapse was ignored because the gene nodes were",
                "already collapsed")
            warning(to_print, "\n")

        } else {
            network_table = orthology_clustering(network_table, organism_code, all_genes)
            networkBW_i = graph.data.frame(network_table, directed = TRUE)
        }
    }

    ## Get distance table#
    distance_table = MetaboSignal_distances(network_table, organism_code,
        organism_name, mode = mode, names = FALSE)

    ## Check if the input genes are present in the network_table
    all_genes = rownames(distance_table)
    ortho_answer = unique(grepl(organism_code, all_genes))
    ortho_answer = TRUE %in% ortho_answer  # are genes orthology IDs ?
    source_genes = gsub(" ", "", source_genes)  # remove potential white spaces.
    source_genes = unique(source_genes)
    source_genes_backup = source_genes
    target_metabolites = gsub(" ", "", target_metabolites)
    target_metabolites = unique(target_metabolites)
    target_metabolites_backup = target_metabolites

    if (ortho_answer == TRUE) {
        # gene nodes are organism specific IDs.  Check if the source
        # genes represent entrez IDs or symbols
        if (grepl(organism_code, source_genes[1]) == FALSE) {
            source_genes = try(MS_GetKEGG_GeneID(genes = source_genes,
                organism_code = organism_code, organism_name = organism_name,
                orthology = FALSE), silent = TRUE)
        }
    } else {
        # Gene nodes are orthology IDs.
        if ((substr(source_genes[1], 1, 1) == "K" & nchar(source_genes[1]) ==
            6) == FALSE) {
            # The input genes are entrezIDs or symbols
            source_genes = try(MS_GetKEGG_GeneID(genes = source_genes,
                organism_code = organism_code, organism_name = organism_name,
                orthology = TRUE), silent = TRUE)
        }
    }

    ## Check if it was possible to map source_genes onto KEGG IDs
    if (grepl("Error", source_genes)[1]) {
        stop("None of the genes is reported in KEGG. Could be that organism_name or organism_code are incorrect")
    }

    ## Check if the source_genes and target_metabolites can be
    ## mapped onto the network#
    all_nodes = unique(as.vector(network_table))
    all_metabolites = all_nodes[grep("cpd:", all_nodes)]

    mapped_genes = intersect(source_genes, all_nodes)
    mapped_metabolites = intersect(target_metabolites, all_metabolites)

    if (length(mapped_genes) > 0) {
        source_genes = mapped_genes
    } else {
        stop_source_genes = paste("None of the source_genes is present in",
            "the network. Could be that the organism_code is incorrect")
        stop(stop_source_genes)
    }

    if (length(mapped_metabolites) > 0) {
        target_metabolites = mapped_metabolites
    } else {
        stop("None of the target_metabolites is present in the network")
    }

    ## Calculate shortest path network
    message("Building shortest path network", "\n")

    all_pathsGM = lapply(target_metabolites, metabolite_SP, source_genes,
                         network_table, distance_table, distance_th, mode, type,
                         BW_matrix, networkBW_i)
    all_pathsGM = unique(do.call(rbind, all_pathsGM))
    colnames(all_pathsGM) = c("node1", "node2")

    if (length(all_pathsGM) >= 1) {
        all_pathsGM_names = all_pathsGM
        if (names == TRUE) {
            all_nodes_network = unique(as.vector(all_pathsGM))
            all_nodes_names = MS_ChangeNames(all_nodes_network,
                organism_code)

            for (i in seq_along(all_nodes_names)) {
                all_pathsGM_names[all_pathsGM_names == all_nodes_network[i]] = all_nodes_names[i]
            }

        } else (all_pathsGM_names = all_pathsGM)

        ## Export network in cytoscape format
        namescyto = names

        if (export_cytoscape == TRUE) {
            message("Creating cytoscape files", "\n")
            all_target_nodes = c(source_genes, target_metabolites)

            cytoscape = MS_ToCytoscape(all_pathsGM, organism_code,
                names = namescyto, target_nodes = all_target_nodes,
                file_name = file_name)
        }

        ## Final report
        check_unmapped(target_metabolites_backup, target_metabolites,
            source_genes_backup, source_genes)

        return(all_pathsGM_names)

    } else {
        to_print = paste("Could not connect any of the source_genes with any",
            "of the target_metabolites")
        stop(to_print)
    }
}

