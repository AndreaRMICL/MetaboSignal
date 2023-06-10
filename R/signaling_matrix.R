#################### get_interactiontype ####################
get_interactiontype = function(path) {

    file = paste("https://rest.kegg.jp/get/", path, "/kgml", sep = "")
    pathway = try(getURL(file), silent = TRUE)
    interactions = parseKGML2DataFrame(pathway, reactions = FALSE)

    if(nrow(interactions) == 0) {
        return(NULL)
    } else {
        return(interactions)
    }
}

#################### collapse_interactions ####################
collapse_interactions = function(interaction, all_interactions) {
    ind = which(all_interactions[, 1] == interaction[1] &
                all_interactions[, 2] == interaction[2])
    subtypes = sort(unique(all_interactions[ind, "subtype"])) ## Updated on 5/09/19 - Thank you: Shilpa Harshan!
    subtype = paste(subtypes, collapse = "/")
    subtype = gsub(" ", "-", subtype)
    new_line = c(interaction[1], interaction[2], subtype)
    return(new_line)
}

#################### build_signal_edges ####################
build_signal_edges = function(node, edges) {
    linked_nodes = edges[[node]]

    if (length(linked_nodes) >= 1) {
      # Node with edges
        new_edge = cbind(rep(node, length(linked_nodes)), edges[[node]])
    } else {
      new_edge = NULL
    }
    return(new_edge)
}

#################### signaling_matrix ####################
signaling_matrix = function(global_network_all) {

    signaling_table = NULL

    ## Create network_edges matrix
    nodes = nodes(global_network_all)
    edges = KEGGgraph::edges(global_network_all)

    network_edges = lapply(nodes, build_signal_edges, edges)
    network_edges = unique(do.call(rbind, network_edges))

    if (length(network_edges) == 0) {
        to_print = ("Impossible to build a signaling network")
        warning(to_print, "\n")
    } else {
        network_edges = matrix(network_edges, ncol = 2)

        ## Remove the edges that contain paths
        all_lines_paths = c(grep("path", network_edges[, 1]),
            grep("path", network_edges[, 2]))
        # All edges contain paths.
        if (length(all_lines_paths) == nrow(network_edges)) {
            to_print = ("Impossible to build a signaling network")
            warning(to_print, "\n")
        } else {
            if (length(all_lines_paths) >= 1) {
                network_edges = network_edges[-(all_lines_paths), ]
                network_edges = matrix(network_edges, ncol = 2)
            }
            ## Remove non-biological nodes from KEGG#
            nonbio_nodes = paste("fibrates", "non-steroids", "agents",
                                 "Antiinflammatory", "BR:br08303",
                                 "Thiazolidinediones", sep = "|")
            nonbio_edges = c(grep(nonbio_nodes, network_edges[,1]),
                             grep(nonbio_nodes, network_edges[,2]))
            nonbio_edges = unique(nonbio_edges)

            if (length(nonbio_edges) >= 1) {
                network_edges = network_edges[-c(nonbio_edges), ]
            }
            signaling_table = unique(network_edges)
            signaling_table = matrix(signaling_table, ncol = 2)
            colnames(signaling_table) = c("node1", "node2")
        }
    }
    return(signaling_table)
}

#################### MS_interactionType ####################
MS_interactionType = function(signaling_paths, all_paths) {

    ## Get correct paths
    paths_org = intersect(all_paths, signaling_paths)

    paths_orgDF = do.call(rbind, lapply(paths_org, get_interactiontype))

    if(is.null(paths_orgDF)) {
        return (NULL)
    }

    interactionM = na.omit(unique(as.matrix(paths_orgDF, ncol = 3)))
    # This has been updated-July 2017
    rownames(interactionM) = NULL

    ## Collapse interactions ##
    IM_lines = split(interactionM, row(interactionM))
    interactionM_collapsed = do.call(rbind, lapply(IM_lines, collapse_interactions,
                                                   interactionM))
    interactionM_collapsed = as.matrix(unique(interactionM_collapsed), ncol = 3)
    interactionM_collapsed = gsub("compound", "indirect-compound", interactionM_collapsed)
    colnames(interactionM_collapsed) = c("source", "target", "interaction_type")
    rownames(interactionM_collapsed) = NULL

    return(interactionM_collapsed)
}

