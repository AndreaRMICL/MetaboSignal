#################### merge_interactions ####################
merge_interactions = function(st, network, target = 3) {
    #print(st)
    collapsed_edges = paste(network[, 1], network[, target], sep = "_")
    ind = which(collapsed_edges == st)
    type = paste(sort(network[ind, "interaction_type"]), collapse = ";")
    new_edge = c(unlist(strsplit(st, "_")), type)
    return(new_edge)
}

#################### MS2_ppiNetwork ####################
MS2_ppiNetwork = function(datasets = "all") {

    ## Check datasets
    datasets = unique(tolower(datasets))

    if(datasets[1] != "all") {
        all_databases = unique(unlist(strsplit(regulatory_interactions[, "database_reference"], ";")))
        ans = setdiff(datasets, tolower(all_databases))
        if(length(ans) > 0) {
            stop ("Invalid datasets")
        } else {
            datasets = paste(datasets, collapse = "|")
        }
    } else {
        datasets = paste(c("omnipath","trrust"), collapse = "|")
    }

    ## Select subset
    ind = grep(datasets, tolower(regulatory_interactions[, "database_reference"]))
    regulatory_interactions = regulatory_interactions[ind, ]
    collapsed_edges = paste(regulatory_interactions[, "source_entrez"],
                            regulatory_interactions[, "target_entrez"], sep = "_")
    rownames(regulatory_interactions) = collapsed_edges
    if (length(collapsed_edges) != length(unique(collapsed_edges))) {
        ## Needs to collapse interactions
        duplicated_edges = collapsed_edges[duplicated(collapsed_edges)]
        unique_edges = regulatory_interactions[setdiff(collapsed_edges,
            duplicated_edges), c("source_entrez", "target_entrez", "interaction_type")]
        new_edges = do.call(rbind, lapply(duplicated_edges, merge_interactions,
                                          regulatory_interactions))
        regulatory_interactions_updated = unique(rbind(unique_edges, new_edges))
    } else {
        regulatory_interactions_updated =
          regulatory_interactions[, c("source_entrez", "target_entrez", "interaction_type")]
    }
        rownames(regulatory_interactions_updated) = NULL
        colnames(regulatory_interactions_updated) = c("source", "target", "interaction_type")

        network_features(regulatory_interactions_updated[, 1:2])

        return(regulatory_interactions_updated)

}

#################### MS2_mergeNetworks ####################
MS2_mergeNetworks = function(network_table1, network_table2) {

    ## Check networks
    check_matrix_v2(network_table1)
    check_matrix_v2(network_table2)

    network1 = network_table1[, 1:3]
    network2 = network_table2[, 1:3]

    ## Merge networks
    merged_network = unique(rbind(network1, network2))

    ## Merge interactions if necessary
    collapsed_edges = paste(merged_network[, 1], merged_network[, 2],
                            sep = "_")
    rownames(merged_network) = collapsed_edges
    duplicated_edges = collapsed_edges[duplicated(collapsed_edges)]

    if(length(duplicated_edges) == 0) {
        network_features(merged_network[, 1:2])
        rownames(merged_network) = NULL
        return(merged_network)
    }

    unique_edges = merged_network[setdiff(collapsed_edges, duplicated_edges), 1:3]

    new_edges = do.call(rbind, lapply(duplicated_edges, merge_interactions,
                                       merged_network, target = 2))

     merged_network_updated = rbind(unique_edges, new_edges)
     rownames(merged_network_updated) = NULL
     colnames(merged_network_updated) = colnames(merged_network)[1:3]

     network_features(merged_network_updated[, 1:2])

     return(merged_network_updated)

}
