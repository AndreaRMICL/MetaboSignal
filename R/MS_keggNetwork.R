#################### get_metabonet ###################
get_metabonet = function(path, all_paths, organism_code) {
    message(path)
    if (path %in% all_paths) {
        if (substr(path, (nchar(organism_code) + 1), nchar(path)) == "01100") {
          # remove metabolic pathways map
            parsed_path = NULL
        } else {
            # Check that the input path exists
            file = paste("https://rest.kegg.jp/get/", path, "/kgml", sep = "")
            pathway = try(getURL(file), silent = TRUE)
            reactions = try(getReactions(parseKGML(pathway)), silent = TRUE)

            if (grepl("Error", reactions[1]) == TRUE) {
                to_print = paste(path, "-path ID without XML:path removed", sep = "")
                message(to_print)
                parsed_path = NULL
            } else {
                parsed_path = capture.output(reactions, file = NULL)
            }
        }
    } else {
        to_print = paste(path, "-incorrect path ID:path removed", sep = "")
        message(to_print)
        parsed_path = NULL
    }
    return(parsed_path)
}

#################### get_signalnet ###################
get_signalnet = function(path, all_paths) {
    message(path)
    if (path %in% all_paths) {
        file = paste("https://rest.kegg.jp/get/", path, "/kgml", sep = "")
        pathway = try(getURL(file), silent = TRUE)
        path_parsed = try(parseKGML(pathway), silent = TRUE)
        path_network = try(KEGGpathway2Graph(path_parsed, genesOnly = FALSE,
            expandGenes = TRUE), silent = TRUE)
        urlcheck = try(edges(path_network), silent = TRUE)
        if (grepl("Error", urlcheck)[1] == TRUE) {
            to_print = paste(path, "-path ID without XML:path removed", sep = "")
            message(to_print)
            parsed_path = "bad_path"
        } else {
            parsed_path = path_network
        }
    } else {
        to_print = paste(path, "-incorrect path ID:path removed", sep = "")
        message(to_print)
        parsed_path = "bad_path"
    }
    return(parsed_path)
}

#################### get_reaction_type ####################
get_reaction_type = function(edge_compound, edge_compound_rev) {
    reversible_edges = intersect(edge_compound, edge_compound_rev)

    if (length(reversible_edges) > 0) {
        irreversible_edges = setdiff(edge_compound, reversible_edges)
        rev_edges = do.call(rbind, strsplit(reversible_edges, "_"))
        edges_direc = cbind(rev_edges, "compound:reversible")

        if (length(irreversible_edges) > 0) {
            irev_edges = do.call(rbind, strsplit(irreversible_edges, "_"))
            irev_edges = cbind(irev_edges, "compound:irreversible")
            edges_direc = rbind(edges_direc, irev_edges)
        }
    } else { # no reversible edges
        irev_edges = do.call(rbind, strsplit(edge_compound, "_"))
        edges_direc = cbind(irev_edges, "compound:irreversible")
    }
    colnames(edges_direc) = c("source", "target", "interaction_type")
    rownames(edges_direc) = NULL

    return(edges_direc)
}

#################### get_signaling_type ####################
get_signaling_type = function(st, interaction_type) {
    collapsed_edges = paste(interaction_type[, 1], interaction_type[, 2],
                            sep = "_")
    ind = which(collapsed_edges == st)

    if(length(ind) == 0) {
        #print(st)
        type = "unknown"
    } else {
        type = gsub(";", "/", interaction_type[ind, 3])
    }
    return(type)
}

#################### map_interactions ####################
map_interactions = function(MetaboSignal_table, interaction_type) {
    edge_net = paste(MetaboSignal_table[, 1], MetaboSignal_table[, 2], sep = "_")
    edge_net_rev = paste(MetaboSignal_table[, 2], MetaboSignal_table[, 1], sep = "_")
    cpd_ind = grep("cpd:|gl:|dr:", edge_net) # compounds drugs or glycans

    if (length(cpd_ind) > 0 & (length(cpd_ind) < length(edge_net))) {
      # There are compounds and genes
        edge_gene = edge_net[-cpd_ind]
        signaling_interactions = sapply(edge_gene, get_signaling_type,
                                        interaction_type)
        edge_compound = edge_net[grep("cpd:|gl:|dr:", edge_net)] # compounds drugs or glycans
        edge_compound_rev = edge_net_rev[grep("cpd:|gl:|dr:", edge_net)]
        metabo_net = get_reaction_type (edge_compound, edge_compound_rev)
        gene_net = cbind(MetaboSignal_table[-cpd_ind, c(1:2)], signaling_interactions)
        gene_net = gsub("_", "-", gene_net)
        #metabo_net = cbind(MetaboSignal_table[cpd_ind, c(1:2)], metabolic_interactions)
        merged_net = rbind(metabo_net, gene_net)

    } else if (length(cpd_ind) == length(edge_net)) { # no signaling interactions
        edge_compound = edge_net[grep("cpd:|gl:|dr:", edge_net)] # compounds drugs or glycans
        edge_compound_rev = edge_net_rev[grep("cpd:|gl:|dr:", edge_net)]
        metabo_net = get_reaction_type (edge_compound, edge_compound_rev)
        merged_net = metabo_net
    } else { # there are not compound interactions
        edge_gene = edge_net
        signaling_interactions = sapply(edge_gene, get_signaling_type,
                                        interaction_type)
        gene_net = cbind(MetaboSignal_table[ , c(1:2)], signaling_interactions)
        merged_net = gsub("_", "-", gene_net)
    }
    colnames(merged_net) = c("source", "target", "interaction_type")
    rownames(merged_net) = NULL
    merged_net[, 3] = paste("k", merged_net[, 3], sep = "_")
    return(merged_net)
}

#################### find_autoedges ####################
find_autoedges = function(edge) {
    if(edge[1] == edge[2]) {
        return("unwanted")
    } else {
        return("wanted")
    }
}

#################### MS_keggNetwork ###################
MS_keggNetwork = function(metabo_paths = NULL, signaling_paths = NULL,
                          expand_genes = FALSE, convert_entrez = FALSE) {

    ########## 0)Preparatory steps###########

    ## Check if metabo_paths or signaling_paths exists
    if (is.vector(metabo_paths) == FALSE & is.vector(signaling_paths) == FALSE) {
        stop("At least one metabo_path or one signaling_path is required")
    }

    ## Make metabo_paths and signaling_paths unique
    metabo_paths = unique(metabo_paths)
    signaling_paths = unique(signaling_paths)

    ## Check that all the paths belong to the same organism
    input_paths = c(metabo_paths, signaling_paths)

    ## Account for organism code with 4 digits
    org_ans = suppressWarnings(is.na(as.numeric(unlist(strsplit(input_paths[1], ""))[4])))

    if(org_ans == TRUE) {
        organism_code = substr(input_paths, 1, 4)
    } else {
        organism_code = substr(input_paths, 1, 3)
    }

    organism_code = unique(organism_code)
    if (length(organism_code) > 1) {
        stop("All paths have to belong to the same organism: check path IDs")
    }

    ## Get all paths of the organism of interest
    lines = try (keggList("pathway", organism_code), silent = TRUE)

    if (grepl("Error", lines)[1]) { # example: metabo_paths = 'X'
        stop("Incorrect path IDs")
    }
    all_paths = names(lines)

    ## Check if input paths are included
    if(!is.null(metabo_paths)) {
        if (length(intersect(metabo_paths, all_paths)) == 0) {
            stop ("All metabo_paths seem to be incorrect")
        }
    }

    if(!is.null(signaling_paths)) {
        if (length(intersect(signaling_paths, all_paths)) == 0) {
            stop ("All signaling_paths seem to be incorrect")
        }
    }

    ## Initiate paths_included
    paths_included = rep(1, times = length(input_paths))
    names(paths_included) = input_paths

    ## Iniciate metabolic and signaling networks
    metabolic_table_RG = NULL
    signaling_table = NULL

    ########## 1)Build metabolic table###########

    if (length(metabo_paths) >= 1) {

        to_print = ("Building metabolic table")
        message(to_print)
        to_print = ("Reading paths:")
        message(to_print)

        ### Get KGML files and transform them into reaction files####
        list_parsed_paths = lapply(metabo_paths, get_metabonet, all_paths, organism_code)
        names(list_parsed_paths) = metabo_paths
        path_names = metabo_paths

        if (length(list_parsed_paths) == 0) {
            to_print = ("Impossible to build a metabolic network")
            warning(to_print, "\n")

        } else {
            ## Remove empty paths: for example oxidative phosphorylation#
            length_paths = sapply(list_parsed_paths, length)
            empty_paths = which(length_paths <= 1)

            if (length(empty_paths) >= 1) {
                list_parsed_paths = list_parsed_paths[-c(empty_paths)]
                paths_included[path_names[empty_paths]] = 0
                path_names = path_names[-c(empty_paths)]
            }
            if (length(list_parsed_paths) == 0) {
                to_print = ("Impossible to build a metabolic network")
                warning(to_print, "\n")
            } else {

                ### Create metabolic_table_RG ####
                metabolic_table_RG = metabolic_matrix(path_names,
                  list_parsed_paths, organism_code, expand_genes = TRUE)

                ## Collapse metabolic_table_RG
                idx = grepl("_", metabolic_table_RG)
                nodes = metabolic_table_RG[idx]
                metabolic_table_RG[idx] = substr(nodes, 11, nchar(nodes))
            }
        }
    }
    ############ End Metabo #####################

    ########## 2)Build signaling table###########

    if (length(signaling_paths) >= 1) {
        message()
        to_print = ("Building signaling table")
        message(to_print)
        to_print = ("Reading paths:")
        message(to_print)

        network_list = lapply(signaling_paths, get_signalnet, all_paths)

        if (length(network_list) == 0) {
            to_print = ("Impossible to build a signaling network")
            warning(to_print, "\n")

        } else {
            ## Remove empty paths
            length_paths = sapply(network_list, is.vector)
            empty_paths = which(length_paths == 1)

            if (length(empty_paths) >= 1) {
                network_list = network_list[-c(empty_paths)]
                paths_included[signaling_paths[empty_paths]] = 0
            }
            if (length(network_list) == 0) {
                to_print = ("Impossible to build a signaling network")
                warning(to_print, "\n")

            } else {
                ## Print ("hola")
                global_network_all = mergeGraphs(network_list)

                ### Create signaling_table ####
                signaling_table = signaling_matrix(global_network_all)

                ### Get interactions
                interaction_type = MS_interactionType(signaling_paths,
                                                      all_paths)
            }
        }
    }
    ################## End signal###################

    ########## 3)Build MetaboSignal table ###########

    #MetaboSignal_table = NULL  #initiate the variable
    answer_metabo = is.matrix(metabolic_table_RG)
    answer_signal = is.matrix(signaling_table)

    if (answer_metabo == FALSE & answer_signal == FALSE) {
        stop ("Impossible to build a MetaboSignal network with these paths")
    }
    if (answer_metabo == TRUE & answer_signal == TRUE) {
        message()
        to_print = ("Building MetaboSignal table")
        message(to_print)
        MetaboSignal_table = unique(rbind(metabolic_table_RG, signaling_table))

    } else if (answer_metabo == TRUE & answer_signal == FALSE) {
        MetaboSignal_table = unique(metabolic_table_RG)

    } else if (answer_metabo == FALSE & answer_signal == TRUE) {
        MetaboSignal_table = unique(signaling_table)
    }

    ## Remove auto-edges
    edges_response = sapply(split(MetaboSignal_table, row(MetaboSignal_table)),
                            find_autoedges)
    edges_unwanted_global = grep("unwanted", edges_response)

    if (length(edges_unwanted_global) >= 1) {
        MetaboSignal_table = MetaboSignal_table[-c(edges_unwanted_global), ]
        MetaboSignal_table = matrix(MetaboSignal_table, ncol = 2)
    }

    ## Retrieve interaction types
    message()
    to_print = ("Retrieving interactions type")
    message(to_print)
    MetaboSignal_interactions = map_interactions(MetaboSignal_table, interaction_type)

    ## Transform kegg IDs into orthology IDs if required
    if (expand_genes == FALSE) {
        message()
        to_print = ("Transforming gene IDs into orthology IDs")
        message(to_print)

        file_ko = paste("https://rest.kegg.jp/link/ko/", organism_code, sep = "")
        response_ko = getURL(file_ko)
        koTable = convertTable(response_ko)
        koTable[, 2] = substr(koTable[, 2], 4, 9)

        for (r in 1:nrow(MetaboSignal_interactions)) {
            for (c in 1:ncol(MetaboSignal_interactions[, 1:2])) {
                index = which(koTable[, 1] == MetaboSignal_interactions[r, c])
                if (length(index) > 0) {
                    MetaboSignal_interactions[r, c] = koTable[index, 2]
                }
            }
        }
    }

    ## Transform KEGG genes into entrez IDs
    if (convert_entrez & expand_genes & organism_code == "hsa") {
        message()
        message("Transforming kegg ids into entrez ids")
        message()

        net_nodes = unique(as.vector(MetaboSignal_interactions[, 1:2]))
        hsa_nodes = net_nodes[grep("hsa", net_nodes)]
        genes_list = split(hsa_nodes, ceiling(seq_along(hsa_nodes)/100))
        entrez_res = unlist(lapply(genes_list, conv_entrez_kegg, source = "kegg"))
        entrez_ids = as.character(substr(entrez_res, 13, nchar(entrez_res)))
        start = unlist(gregexpr(pattern = "hsa", names(entrez_res)))
        kegg_ids = substr(names(entrez_res), start, nchar(names(entrez_res)))

        for (i in seq_along(entrez_res)) {
            MetaboSignal_interactions[MetaboSignal_interactions == kegg_ids[i]] = entrez_ids[i]
        }
    }

    ## Report network features#
    MetaboSignal_interactions = unique(MetaboSignal_interactions)
    message()
    network_features(MetaboSignal_interactions[, 1:2])

    if (sum(paths_included == 0)) {
        paths_removed = names(paths_included[paths_included == 0])
        path_line = paste(paths_removed, collapse = ",")
        message()
        to_print = paste("Some path IDs were not used:", path_line, sep = "")
        warning(to_print, "\n")
    }
    return(MetaboSignal_interactions)
}
