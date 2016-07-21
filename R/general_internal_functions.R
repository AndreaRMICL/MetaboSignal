##This script contains internal functions that are used in several of the other
#scripts

#################### convertTable ####################
convertTable = function(res) {
    if (nchar(res) == 0) {
        print("no result")
        result = NULL
    } else {
        rows = strsplit(res, "\n")
        rows.len = length(rows[[1]])
        result = matrix(unlist(lapply(rows, strsplit, "\t")), nrow = rows.len,
                        byrow = TRUE)
    }
    return(result)
}

####################### igraph_edited ######################
# This function makes all the edges linked to the metabolite as
#reversible if the metabolite is a substrate. It is used to calculate
#shortest paths with SP mode.
igraph_edited = function(MetaboSignal_table, metabolite) {

    # Chek if the metabolite is a acts as a substrate (first column)
    index = which(MetaboSignal_table == metabolite, arr.ind = TRUE)[, 1]

    if (length(index) == 0) {
        MetaboSignal_table = MetaboSignal_table
    } else {
        linked_nodes = unique(c(MetaboSignal_table[index, 1], MetaboSignal_table[index, 2]))
        index_metabolite = grep(metabolite, linked_nodes)
        linked_nodes = linked_nodes[-index_metabolite]
        new_edges = rbind(cbind(rep(metabolite, length(linked_nodes)), linked_nodes),
                          cbind(linked_nodes, rep(metabolite, length(linked_nodes))))

        MetaboSignal_table = unique(rbind(MetaboSignal_table, new_edges))

    }
    MetaboSignal_networkEdited_i = graph.data.frame(MetaboSignal_table,
                                                    directed = TRUE)

    return(MetaboSignal_networkEdited_i)
}

######################### match ############################
Match = function(x, want) {
    out = sapply(x, function(x, want) isTRUE(all.equal(x, want)), want)
    any(out)
}

####################### check_matrix #######################
# The goal of this function is to make sure that the network_table is a 2 column
#matrix
check_matrix = function(network_table) {
    network_table = unique(network_table)
    network_table = try(matrix(network_table, ncol = 2), silent = TRUE)
    rows = try(nrow(network_table), silent = TRUE)  # example network_table = list()
    if (rows == 0) {
        stop("network_table needs to be a 2-column matrix", call. = FALSE)
    }
    if (grepl("Error", network_table)[1]) {
        stop("network_table needs to be a 2-column matrix", call. = FALSE)
    }
    return(network_table)
}

####################### match_KEGG #######################
match_KEGG = function(match, target_column, target_matrix) {
    matchM = list()
    match = tolower(match)
    match = unique(match)
    length(matchM) = length(match)
    names(matchM) = match

    for (key_word in match) {
        key_wordS = unlist(strsplit(key_word, " "))

        if (length(key_wordS) > 1) {
            main_index = grep(key_wordS[1], tolower(target_column))
            if (length(main_index) == 0) {
                to_print = paste("Could not match", key_word, sep = " ")
                message(to_print)
            } else {
                secondary_index = c()
                index_to_intersect = main_index
                for (p in 2:length(key_wordS)) {
                  indexS = grep(key_wordS[p], tolower(target_column))
                  if (length(intersect(main_index, indexS)) == 0) {
                    indexS = NULL
                  }
                  if (length(indexS) > 0) {
                    secondary_index = c(secondary_index, indexS)
                    secondary_index = unique(secondary_index)
                    index_to_intersect = intersect(index_to_intersect, indexS)
                  }
                }
                if (length(secondary_index) > 0) {
                  if (length(index_to_intersect) > 0) {
                    index = index_to_intersect
                  } else {
                    intersection = intersect(main_index, secondary_index)
                    if (length(intersection) > 0) {
                      index = intersection
                    } else {
                      index = main_index
                    }
                  }
                } else {
                  index = main_index
                }
                key_wordLines = target_matrix[index, ]
                matchM[[key_word]] = key_wordLines
            }
        } else {
            key_word = key_wordS
            index = grep(key_word, tolower(target_column))
            if (length(index) == 0) {
                to_print = paste("Could not match", key_word, sep = " ")
                warning(to_print)
            } else {
                key_wordLines = target_matrix[index, ]
                matchM[[key_word]] = key_wordLines
            }
        }
    }
    return(matchM)
}

####################### check_mode_type #######################
check_mode_type = function(mode = "all", type = "all", mode2 = "all", type2 = "all") {
    if ( mode[1] %in% c("all", "out", "SP") == FALSE) {
      stop ("Invalid mode. Possible values for mode are: all, SP, out", call. = FALSE)
    }
    if ( mode2[1] %in% c("all", "out") == FALSE) {
      stop ("Invalid mode. Possible values for mode are: all, out", call. = FALSE)
    }
    if ( type[1] %in% c("bw", "first", "all") == FALSE){
      stop ("Invalid type. Possible values for type are: bw, first, all", call. = FALSE)
    }
    if ( type2[1] %in% c("bw", "distance", "all") == FALSE){
      stop ("Invalid type. Possible values for type are: bw, distance, all", call. = FALSE)
    }
}

####################### network_features #######################
network_features = function(network_table) {
    all_nodes = unique(as.vector(network_table))
    nodes_number = length(all_nodes)
    message("Network features:")
    to_print = paste("Number of nodes:", nodes_number, sep = "")
    message(to_print)
    to_print = paste("Number of edges:", nrow(network_table), sep = "")
    message(to_print,"\n")
}

####################### check_unmapped #######################
check_unmapped = function (metabolites_backup, metabolites, genes_backup, genes) {
    if (length(metabolites_backup) > length(metabolites)) {
      difference = length(metabolites_backup) - length(metabolites)
      to_print = paste("n =", difference, "target_metabolites not mapped onto",
                       "the network. To check them, use the function",
                       "MS_FindMappedNodes()")
      warning(to_print, call. = FALSE)
    }
    if (length(genes_backup) > length(genes)) {
      difference = length(genes_backup) - length(genes)
      to_print = paste("n =", difference, "source_genes were not reported in",
                       "KEGG or were not mapped into the network.",
                       "To check them, use the function",
                       "MS_FindMappedNodes()")
      warning(to_print, call. = FALSE)
    }
}
