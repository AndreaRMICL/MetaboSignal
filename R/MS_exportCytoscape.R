############### cyto_directionality #############
cyto_directionality = function(edge) {
    uncollapsed_edge = unlist(strsplit(edge, ">"))
    new_edge = c(sort(uncollapsed_edge[1:2]), uncollapsed_edge[3])
    return(new_edge)
}

############### get_source ############# Updated
get_source = function(edge) {
    ans_k = grepl("k_", edge)
    ans_o = grepl("o_", edge)
    ans_t = grepl("t_", edge)
    #ans_m = grepl("m_", edge)

   # if (ans_m) { #microRNA
        #return("mirtarbase")
    #}
    if (ans_k + ans_o + ans_t == 3) {
        return("kegg; omnipath;trrust")
    }
    if (ans_k + ans_o + ans_t == 0) {
        return("other")
    }
    if (ans_k == TRUE & ans_o == TRUE & ans_t == FALSE) {
        return("kegg;omnipath")
    }
    if (ans_k == TRUE & ans_o == FALSE & ans_t == TRUE) {
        return("kegg;trrust")
    }
    if (ans_k == FALSE & ans_o == TRUE & ans_t == TRUE) {
        return("omnipath;trrust")
    }
    if (ans_k == TRUE & ans_o == FALSE & ans_t == FALSE) {
        return("kegg")
    }
    if (ans_k == FALSE & ans_o == TRUE & ans_t == FALSE) {
        return("omnipath")
    }
    if (ans_k == FALSE & ans_o == FALSE & ans_t == TRUE) {
        return("trrust")
    }
}

#################### get_molecule_type ####################### Updated
get_molecule_type = function(node, organism_code) {
    nodeBU = node # new
    if( organism_code == "hsa" & !is.na(suppressWarnings(as.numeric(node)))) { # if it is entrez
        node = paste("hsa:", node, sep = "")
        #node = conv_entrez_kegg(node, source = "entrez", organism_code = "hsa")
    }
    if (grepl("cpd:|dr:|gl", node) == TRUE) {
        node_type = "compound"
    } else if (grepl("rn", node) == TRUE) {
        node_type = "reaction"
    } else if (grepl(organism_code, node) == TRUE | grepl("K", node) == TRUE) {
        file = paste("https://rest.kegg.jp/get/", node, sep = "")
        lines = try(readLines(file), silent = TRUE)
        if (grepl("Error", lines[1]) == FALSE) {
            enzyme_lines = grep("EC:", lines[1:5])
            metabo_lines = grep("Metabolism", lines)
            if (length(enzyme_lines) >= 1 & length(metabo_lines) > 0) {
                node_type = "metabolic-gene"
            } else (node_type = "signaling-gene")
        } else {
            if (nodeBU %in% regulatory_interactions[, c(1, 3)]) {
                node_type = "signaling-gene"
            } else {
                node_type = "other"
            }
        }
    } else if (grepl("hsa-miR", node)) {
        node_type = "microRNA"
    } else {
        node_type = "other"
    }

    return(node_type)
}

############### MS_exportCytoscape #############
MS_exportCytoscape = function(network_table, organism_code, names = TRUE,
    targets = NULL, file_name = "MS") {

    ## Check that network is correct
    network = unique(network_table)
    check_matrix_v2(network, n = 3)

    ## Check targets
    targets = intersect(targets, unique(as.vector(network)))

    if (length(targets) == 0) {
        targets = NULL
    }

    ## Arrange directionality
    rev_ind = grep("compound:reversible", network[, 3])

    if (length(rev_ind) > 0) {
        rev_edges = paste(network[rev_ind, 1], network[rev_ind, 2],
                          network[rev_ind, 3], sep = ">")
        new_edges = do.call(rbind, lapply(rev_edges, cyto_directionality))
        network = rbind(network[-rev_ind, ], unique(new_edges))
        rownames(network) = NULL
    }

    if (names == TRUE) {
        cytoscape_net = network_names(network, organism_code)
    } else {
        cytoscape_net = network
    }

    ## Get source
    all_net_edges = paste(cytoscape_net[, 1], cytoscape_net[, 2],
                          cytoscape_net[, 3], sep = "_")
    sources = as.character(sapply(all_net_edges, get_source))
    cytoscape_netDF = unique(as.data.frame(cbind(cytoscape_net, database = sources),
        rownames = NULL)) # Make it unique in case of duplicated common names

    file_nameN = paste(file_name, "Network.txt", sep = "_")
    write.table(cytoscape_netDF, file_nameN, row.names = FALSE,
        sep = "\t", quote = FALSE, col.names = TRUE)

    ## Get node type
    all_nodes = unique(as.vector(network[, 1:2]))
    types = sapply(all_nodes, get_molecule_type, organism_code = organism_code)
    node_type = cbind(all_nodes, types)
    colnames(node_type) = c("node", "type")
    rownames(node_type) = NULL

    if (names == TRUE) { # This has been updated
        on_nodes = unique(cbind(as.vector(network[, 1:2]),
                                as.vector(cytoscape_net[, 1:2])))
        # in case the several KEGG IDs are associated to the same symbol
        for (i in 1:nrow(node_type)) {
            ind = which(on_nodes[, 1] == node_type[i, 1])[1]
            node_type[i, 1] = on_nodes[ind, 2]
        }
        colnames(node_type) = c("node", "type")
        rownames(node_type) = NULL
    }

    node_typeDF = unique(as.data.frame(node_type, rownames = NULL))
    file_nameNT = paste(file_name, "NodesType.txt", sep = "_")

    write.table(node_typeDF, file_nameNT, row.names = FALSE,
        sep = "\t", quote = FALSE, col.names = TRUE)

    ## Write targets
    if (!is.null(targets)) { # This has been updated
        if (names == TRUE) {
            for (i in 1:length(targets)) {
                ind = which(on_nodes[, 1] == targets[i])[1]
                targets[i] = on_nodes[ind, 2]
            }
        }
        targetDF = unique(as.data.frame(cbind(targets, "target"), rownames = NULL))
        colnames(targetDF) = c("node", "target")
        file_nameT = paste(file_name, "TargetNodes.txt", sep = "_")
        write.table(targetDF, file_nameT, row.names = FALSE,
            sep = "\t", quote = FALSE, col.names = TRUE)

    }
    return(cytoscape_netDF)
}
