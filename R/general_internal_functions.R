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

####################### check_matrix #######################
#The goal of this function is to make sure that the network_table is a 2 column matrix
check_matrix = function(network_table) {
    network_table = unique(network_table)
    network_table = try(matrix(network_table, ncol = 2), silent = TRUE)
    rows = try(nrow(network_table), silent = TRUE)  # example network_table = list()
    if (rows == 0) {
        stop("network_table needs to be a 2-column matrix", call. = FALSE)
    }
    if (grepl("Error", network_table[1])) {
        stop("network_table needs to be a 2-column matrix", call. = FALSE)
    }
    return(network_table)
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
        node_bw = 0
    }
    return(node_bw)
}

####################### get_global_BW_score #######################
get_global_BW_score = function (row, BW_matrix) {
    path_individual = as.character(row)
    BW = sapply(path_individual, get_bw_score, BW_matrix)
    score_BW = as.numeric(BW)
    score_BW = sum(BW)/(length(BW)-sum(BW==0))

    return(score_BW)
}

####################### BW_ranked_SP #######################
BW_ranked_SP = function (all_paths, BW_matrix, networkBW_i, mode) {

    all_nodes = unique(as.vector(all_paths))
    index_cpd = grep("cpd:", all_nodes)

    if (length(index_cpd) > 0) {
    # This should be always true because if there are not
    # compounds in the network the function MetaboSignal_distances wont work.
        gene_nodes = all_nodes[-index_cpd]
    } else (gene_nodes = all_nodes)

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

####################### check_unmapped #######################
check_unmapped = function(metabolites_backup, metabolites, genes_backup, genes) {

    if (length(metabolites_backup) > length(metabolites)) {
        difference = length(metabolites_backup) - length(metabolites)
        to_print = paste("n =", difference, "target_metabolites not mapped onto",
            "the network. To check them, use the function", "MS_FindMappedNodes()")
        warning(to_print, call. = FALSE)
    }
    if (length(genes_backup) > length(genes)) {
        difference = length(genes_backup) - length(genes)
        to_print = paste("n =", difference, "source_genes not reported in",
            "KEGG or not mapped onto the network.", "To check them, use the function",
            "MS_FindMappedNodes()")
        warning(to_print, call. = FALSE)
    }
}

#################### find_unwanted_edge ####################
find_unwanted_edge = function(edge, organism_code, expand_genes) {
    if (expand_genes == TRUE) {
        pattern = paste(organism_code, "cpd:|rn:", sep = "|")
    } else {
        # if(grepl(organism_code, edge[1]) | grepl(organism_code, edge[2])) {
        #     return ("unwanted")
        # }
        pattern = "K|cpd:|rn:"
    }
    if (grepl(pattern, edge[1]) == FALSE | grepl(pattern, edge[2]) ==
        FALSE | edge[1] == edge[2]) {
        return("unwanted")
    } else {
        return("wanted")
    }
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

#################### find_gene_metabo #######################
find_gene_metabo = function(node, organism_code) {
    pattern_gene = paste(organism_code, "K|rn", sep = "|")
    if (grepl(pattern_gene, node)) {
        return("gene")
    } else if (grepl("cpd:", node)) {
        return("metabolite")
    } else {
        return("other")
    }
}

#################### get_molecule_type #######################
get_molecule_type = function(node, organism_code) {
    if (grepl("cpd:", node) == TRUE) {
        node_type = "metabolite"
    } else if (grepl(organism_code, node) == TRUE | grepl("K", node) == TRUE) {
        file = paste("http://rest.kegg.jp/get/", node, sep = "")
        lines = try(readLines(file), silent = TRUE)
        if (grepl("Error", lines[1]) == FALSE) {
            enzyme_lines = grep("EC:", lines[1:5])
            metabo_lines = grep("Metabolism", lines)
            if (length(enzyme_lines) >= 1 & length(metabo_lines) > 0) {
                node_type = "metabolic_gene"
            } else (node_type = "signaling_gene")
        } else (node_type = "other")
    } else {
        node_type = "other"
    }
    return(node_type)
}

#################### get_target_type #######################
get_target_type = function(node, target_nodes) {
    if (node %in% target_nodes == TRUE) {
        type = "target"
    } else {
        type = "non-target"
    }
    return(type)
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

##################### orthology_clustering ######################
#This function clusters organism-specific KEGG gene IDs into orthology IDs
orthology_clustering = function(network_table, organism_code, all_genes) {

    to_print = ("Clustering gene IDs into orthology IDs")
    message(to_print, "\n")
    index_genes = grep(organism_code, all_genes)
    gene_nodes = unique(all_genes[index_genes])

    file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code, sep = "")
    response_ko = getURL(file_ko)
    koTable = convertTable(response_ko)
    koTable[, 2] = substr(koTable[, 2], 4, 9)

    for (i in 1:nrow(network_table)) {
        for (z in 1:ncol(network_table)) {
            index = which(koTable[, 1] == network_table[i, z])
            if (length(index) > 0) {
                network_table[i, z] = koTable[index, 2]
            }
        }
    }

    # Remove edges without ko orthology
    unwanted_idx = c(grep(organism_code, network_table[, 1]),
        grep(organism_code, network_table[, 2]))

    if (length(unwanted_idx) > 0) {
        network_table = network_table[-unwanted_idx, ]
    }

    network_table = unique(network_table)
    return(network_table)
}

#################### From_geneID_to_symbol ################
From_geneID_to_symbol = function(ID) {
    file = file = paste("http://rest.kegg.jp/get/", ID, sep = "")
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
        name = toupper(name)
        name = gsub(" ", "", name)
    } else {
        name = ID
    }
    return(name)
}

#################### MS_FindCompound ####################
MS_FindCompound = function(match = NULL) {
    file = "http://rest.kegg.jp/list/compound"
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
    file = "http://rest.kegg.jp/list/organism"
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
    file = paste("http://rest.kegg.jp/list/pathway/", organism_code, sep = "")
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

#################### MS_FindMappedGenes ####################
MS_FindMappedGenes = function(genes, organism_code, organism_name,
                              network_table, orthology = TRUE) {

    List_mappedID = vector (mode = "list", length = 4)
    names(List_mappedID) = c("Input gene IDs mapped onto the network",
        "All input gene IDs not mapped onto the network", "Input gene IDs not found in KEGG",
        "KEGG gene IDs not mapped onto the network")

    genes = gsub(" ", "", genes)  #remove potential white spaces
    genes = tolower(unique(genes))

    transformed_genesM = try(MS_GetKEGG_GeneID(genes, organism_code,
        organism_name, output = "matrix", orthology = orthology), silent = TRUE)

    if (grepl("Error", transformed_genesM[1]) == TRUE) {
        List_mappedID[[2]] = genes
        List_mappedID[[3]] = genes
        return(List_mappedID)
    }

    if (is.matrix(transformed_genesM) == TRUE) {
        if (ncol(transformed_genesM) == 3) { # Input genes are symbols.
            transformed_IDs = as.character(transformed_genesM[, 3])
        } else (transformed_IDs = as.character(transformed_genesM[, 2]))
        # Input genes are entrez IDs.
    } else {
        List_mappedID[[2]] = genes
        List_mappedID[[3]] = genes
        return(List_mappedID)
    }

    ## Find input genes that are not mapped in KEGG#
    IDs_not_inKEGG = setdiff(genes, transformed_IDs)

    ## Find ko or organism specific genes (transformed genes) that
    ## are mapped in KEGG but not in the network
    transformed_genes = as.character(transformed_genesM[, 1])
    all_nodes = unique(as.vector(network_table))
    genes_not_mapped = setdiff(transformed_genes, all_nodes)

    # Find gene input IDs (symbols or entrez) that of the
    # transformed genes that were not mapped onto the network.
    if (length (genes_not_mapped) > 0) {
        not_mapped_idx = unlist(lapply(genes_not_mapped, find_node_index,
                                       transformed_genes))
        IDs_not_mapped = transformed_IDs[not_mapped_idx]
    } else {
        IDs_not_mapped = NULL
    }

    IDs_not_found_all = c(IDs_not_inKEGG, IDs_not_mapped)
    IDs_found = setdiff(genes, IDs_not_found_all)

    if (length(IDs_not_found_all) > 0) {
        List_mappedID[[1]] = IDs_found
        List_mappedID[[2]] = IDs_not_found_all
        List_mappedID[[3]] = IDs_not_inKEGG
        List_mappedID[[4]] = IDs_not_mapped
    } else {
        to_print = ("All genes were mapped onto the network")
        message(to_print)
        List_mappedID[[1]] = IDs_found
    }
    return(List_mappedID)
}

################## MS_FindMappedMetabolites ###################
MS_FindMappedMetabolites = function(metabolites, network_table) {
    List_metabo = vector(mode = "list", length = 2)
    names(List_metabo) = c("metabolites mapped onto the network",
        "metabolites not mapped onto the network")
    all_nodes = unique(as.vector(network_table))
    metabolites_not_found = setdiff(metabolites, all_nodes)
    metabolites_found = setdiff(metabolites, metabolites_not_found)

    if (length(metabolites_not_found) > 0) {
        List_metabo[[1]] = metabolites_found
        List_metabo[[2]] = metabolites_not_found
    } else {
        to_print = ("All metabolites were mapped onto the network")
        message(to_print)
        List_metabo[[1]] = metabolites_found
    }
    return(List_metabo)
}
