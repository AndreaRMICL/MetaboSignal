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
#This function makes all the edges linked to the metabolite as reversible if the
#metabolite is a substrate. It is used to calculate shortest paths with SP mode.
igraph_edited = function(MetaboSignal_table, metabolite) {

    # Chek if the metabolite is a acts as a substrate (first column)
    index = which(MetaboSignal_table == metabolite, arr.ind = TRUE)[, 1]

    if (length(index) == 0) {
        MetaboSignal_table = MetaboSignal_table
    } else {
        linked_nodes = unique(c(MetaboSignal_table[index, 1],
                                MetaboSignal_table[index, 2]))
        index_metabolite = grep(metabolite, linked_nodes)
        linked_nodes = linked_nodes[-index_metabolite]
        new_edges = rbind(cbind(rep(metabolite, length(linked_nodes)),
            linked_nodes), cbind(linked_nodes, rep(metabolite,
            length(linked_nodes))))

        MetaboSignal_table = unique(rbind(MetaboSignal_table, new_edges))

    }
    MetaboSignal_networkEdited_i = graph.data.frame(MetaboSignal_table, directed = TRUE)

    return(MetaboSignal_networkEdited_i)
}

######################### Match ############################
Match = function(x, want) {
    out = sapply(x, function(x, want) isTRUE(all.equal(x, want)), want)
    any(out)
}

####################### check_matrix_v2 #######################
#The goal of this function is to make sure that the network_table is a 2 column matrix
check_matrix_v2 = function(network_table, n = 3) {
    if (!is.matrix(network_table)) {
        stop ("network_table must be a 3-column matrix")
    }
    if (ncol(network_table) < n) {
        stop ("network_table must be a 3-column matrix")
    }
}

####################### match_KEGG #######################
match_KEGG = function(match, target_column, target_matrix) {
    match = tolower(match)
    match = unique(match)
    matchM = vector (mode = "list", length = length(match))
    names (matchM) = match

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
                to_print = paste("Could not match", key_word,
                  sep = " ")
                warning(to_print)
            } else {
                key_wordLines = target_matrix[index, ]
                matchM[[key_word]] = key_wordLines
            }
        }
    }
    return(matchM)
}

####################### ASP_paths #######################
ASP_paths = function (ASP) {
  shortpath = rownames(as.matrix(unlist(ASP)))
  return(shortpath)
}

####################### get_bw_score #######################
get_bw_score = function (node, BW_matrix) {
    index = which(BW_matrix[, 1] == node)
    if (length(index) > 0) {
        node_bw = as.numeric(BW_matrix[index, 2])
    } else {
        node_bw = -1 # This has been changed (July 2017)
    }
    return(node_bw)
}

####################### get_global_BW_score #######################
get_global_BW_score = function (row, BW_matrix) {
    path_individual = as.character(row)
    BW = sapply(path_individual, get_bw_score, BW_matrix)
    score_BW = as.numeric(BW)
    score_BW = score_BW[score_BW >= 0] ## This has been changed (July 2017)
    #score_BW = sum(BW)/(length(BW)-sum(BW==0))
    score_BW = sum(score_BW)/length(score_BW)
    return(score_BW)
}

####################### BW_ranked_SP #######################
BW_ranked_SP = function (all_paths, BW_matrix, networkBW_i, mode) {

    all_nodes = unique(as.vector(all_paths))
    index_cpd = grep("cpd:|dr:|gl:", all_nodes) # This has ben changed:July 2017

    if (length(index_cpd) > 0) {
    # This should be always true because if there are not
    # compounds in the network the function MetaboSignal_distances wont work.
        gene_nodes = all_nodes[-index_cpd]
    } else{
        gene_nodes = all_nodes # This was changed: July 2017
    }

    gene_nodes = unique(gene_nodes)

    for (gene in gene_nodes) {

        if (gene %in% BW_matrix[, 1] == FALSE) {
            if (mode == "all") {
                bw = betweenness(networkBW_i, gene, directed = FALSE,
                                 weights = NULL, nobigint = TRUE, normalized = TRUE)
            } else {
                bw = betweenness(networkBW_i, gene, directed = TRUE,
                                 weights = NULL, nobigint = TRUE, normalized = TRUE)
            }
            bw_line = c(gene, as.numeric(bw))
            BW_matrix = rbind(BW_matrix, bw_line)
            BW_matrix = unique(BW_matrix)
        }
  }
    ## Get global BW score for each path
    Global_BW_score = sapply (split(all_paths, row(all_paths)), get_global_BW_score,
                              BW_matrix)

    maxBW = which(Global_BW_score == max(Global_BW_score))[1]
    path = all_paths[maxBW, ]
    path = as.character(path)
    all_paths = matrix(path, ncol = length(path))

    return(all_paths)
}

####################### check_mode_type #######################
check_mode_type = function(mode = "all", type = "all", mode2 = "all", type2 = "all") {
    if (mode[1] %in% c("all", "out", "SP") == FALSE) {
        stop("Invalid mode. Possible values for mode are: all, SP, out",
            call. = FALSE)
    }
    if (mode2[1] %in% c("all", "out") == FALSE) {
        stop("Invalid mode. Possible values for mode are: all, out",
            call. = FALSE)
    }
    if (type[1] %in% c("bw", "first", "all") == FALSE) {
        stop("Invalid type. Possible values for type are: bw, first, all",
             call. = FALSE)
    }
    if (type2[1] %in% c("bw", "distance", "all") == FALSE) {
        stop("Invalid type. Possible values for type are: bw, distance, all",
            call. = FALSE)
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
    message(to_print, "\n")
}

#################### find_node_index ####################
find_node_index = function(node, target) {
    index = which(target == node)
    return(index)
}

#################### conv_entrez_kegg ####################
conv_entrez_kegg = function(genes, source = "entrez", organism_code) {
    if (source == "kegg") {
        res = keggConv("ncbi-geneid", genes)
    } else {
        ncbi_genes = paste("ncbi-geneid:", genes, sep = "")
        res = keggConv(organism_code, ncbi_genes)
    }
    if(length(res) == 0) {
        res = genes
    }
    return(res)
}

#################### conv_entrez_symbol ####################
conv_entrez_symbol = function(gene, organism_name, source = "entrez") {
    response = query(q = gene, size = 1, species = organism_name)
    if (NROW(response$hits) != 0L) {
        ID = as.matrix(unlist(response))
        rows = rownames(ID)
        ID = as.character(ID)
        if (source == "entrez") {
            index = which(rows == "hits.symbol")
        } else {
            index = which(rows == "hits.entrezgene")
        }
        if (length(index) > 0) {
            res = ID[index]
        } else {
            res = "NF"  # not found
        }
    } else {
        res = "NF"
    }
    return(res)
}

#################### link_kogene ####################
link_kogene = function(gene, koTable) {
    index_gene = which(koTable[, 1] == gene)
    if (length(index_gene) > 0) {
        ko_line = koTable[index_gene[1], 2]
    } else {
        ko_line = "NF"
    }
    return(ko_line)
}

##################### path_as_network ######################
path_as_network = function(path) {
    all_edges = c()
    for (i in 1:(length(path) - 1)) {
        edge = c(path[i], path[i + 1])
        all_edges = rbind(all_edges, edge)
    }
    return(all_edges)
}

#################### From_geneID_to_symbol ################
From_geneID_to_symbol = function(ID) {
    file = file = paste("https://rest.kegg.jp/get/", ID, sep = "")
    find = try(readLines(file), silent = TRUE)
    if (grepl("Error", find[1]) == FALSE) {
        find = find[2]
        all_names = substr(find, 13, nchar(find))
        name = unlist(strsplit(all_names, "[,]"))
        name = unlist(strsplit(name, "[;]"))
        if (grepl("E", name[1]) == TRUE) {
            # Avoids taking weird names
            name = sort(name, decreasing = FALSE)
        }
        name = name[1]  # Selects only the first gene name
        #name = toupper(name)
        name = gsub(" ", "", name)
    } else {
        name = ID
    }
    return(name)
}

#################### MS_FindCompound ####################
MS_FindCompound = function(match = NULL) {
    file = "https://rest.kegg.jp/list/compound"
    response = getURL(file)
    compoundM = convertTable(response)
    colnames(compoundM) = c("KEGG compound", "common names")
    rownames(compoundM) = NULL
    if (length(match) >= 1) {
        target_matrix = compoundM
        target_column = compoundM[, 2]
        matchM = match_KEGG(match, target_column, target_matrix)
        return(matchM)
    } else (return(compoundM))
}

#################### MS_FindOrganism ####################
MS_FindOrganism = function(match = NULL) {
    file = "https://rest.kegg.jp/list/organism"
    response = getURL(file)
    organismM = convertTable(response)
    colnames(organismM) = c("T", "organism_code", "organism_name",
        "description")
    rownames(organismM) = NULL
    if (length(match) >= 1) {
        target_matrix = organismM
        target_column = organismM[, 3]
        matchM = match_KEGG(match, target_column, target_matrix)
        return(matchM)
    } else (return(organismM))
}

#################### MS_FindPathway ####################
MS_FindPathway = function(match = NULL, organism_code = NULL) {
    file = paste("https://rest.kegg.jp/list/pathway/", organism_code, sep = "")
    response = try(getURL(file), silent = TRUE)
    if (nchar(response) == 0) {
        stop("A valid organism_code is required for KEGG_entry = pathway")
    }
    pathM = convertTable(response)
    colnames(pathM) = c("path_ID", "path_Description")
    rownames(pathM) = NULL
    if (length(match) >= 1) {
        target_matrix = pathM
        target_column = pathM[, 2]
        matchM = match_KEGG(match, target_column, target_matrix)
        return(matchM)
    } else (return(pathM))
}

############### network_names #############
network_names = function(all_pathsGM, organism_code) {
    all_pathsGM_names = all_pathsGM
    all_nodes = unique(as.vector(all_pathsGM[, 1:2]))
    all_names = MS_changeNames(all_nodes, organism_code)

    for (i in seq_along(all_names)) {
        all_pathsGM_names[all_pathsGM_names == all_nodes[i]] = all_names[i]
    }
    return(all_pathsGM_names)
}
