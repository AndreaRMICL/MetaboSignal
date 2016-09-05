#################### get_metabonet ###################
get_metabonet = function(path, all_paths) {
    message(path)
    if (path %in% all_paths) {
        if (substr(path, 4, nchar(path)) == "01100") { # remove metabolic pathways map
            parsed_path = NULL
        } else {
            # Check that the input path exists
            file = paste("http://rest.kegg.jp/get/", path, "/kgml", sep = "")
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
        file = paste("http://rest.kegg.jp/get/", path, "/kgml", sep = "")
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

#################### MetaboSignal_matrix ###################

MetaboSignal_matrix = function(metabo_paths = NULL, signaling_paths = NULL,
                               organism_name = NULL, tissue = "all", expand_genes = FALSE) {

    ########## 0)Preparatory steps###########

    ## Check if metabo_paths or signaling_paths exists
    if (is.vector(metabo_paths) == FALSE & is.vector(signaling_paths) == FALSE) {
        stop("At least one metabo_path or one signaling_path is required")
    }

    ## Check if organism_name exits
    if (length(organism_name) == 0 & tissue[1] != "all") {
        stop("An organism_name is required for tissue filtering")
    }

    ## Make metabo_paths and signaling_paths unique
    metabo_paths = unique(metabo_paths)
    signaling_paths = unique(signaling_paths)

    ## Check that all the paths belong to the same organism
    input_paths = c(metabo_paths, signaling_paths)
    organism_code = substr(input_paths, 1, 3)
    organism_code = unique(organism_code)
    if (length(organism_code) > 1) {
        stop("All paths have to belong to the same organism: check path IDs")
    }

    ## Get all paths of the organism of interest
    lines = try (keggList("pathway", organism_code), silent = TRUE)

    if (grepl("Error", lines)[1]) { # example: metabo_paths = 'X'
        stop("Incorrect path IDs")
    }
    all_paths = substr(names(lines), 6, 13)

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
        list_parsed_paths = lapply(metabo_paths, get_metabonet, all_paths)
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
                  list_parsed_paths, organism_code, expand_genes)
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
                global_network_all = mergeGraphs(network_list)

                ### Create signaling_table ####
                signaling_table = signaling_matrix(global_network_all,
                  tissue, organism_code, organism_name, expand_genes)
            }
        }
    }

    ################## End signal###################

    ########## 3)Build MetaboSignal table ###########

    MetaboSignal_table = NULL  #initiate the variable
    answer_metabo = is.matrix(metabolic_table_RG)
    answer_signal = is.matrix(signaling_table)

    if (answer_metabo == TRUE & answer_signal == TRUE) {
        message()
        to_print = ("Building MetaboSignal table")
        message(to_print)
        MetaboSignal_table = rbind(metabolic_table_RG, signaling_table)
        MetaboSignal_table = unique(MetaboSignal_table)

    } else if (answer_metabo == TRUE & answer_signal == FALSE) {
        MetaboSignal_table = unique(metabolic_table_RG)

    } else if (answer_metabo == FALSE & answer_signal == TRUE) {
        MetaboSignal_table = unique(signaling_table)
    }

    # Collapse network
    if (is.matrix(MetaboSignal_table) == FALSE) {
        stop("Impossible to build a MetaboSignal_table with these paths",
            "\n")
    } else {
        idx = grepl("_", MetaboSignal_table)
        nodes = MetaboSignal_table[idx]
        MetaboSignal_table[idx] = substr(nodes, 11, nchar(nodes))

        ## Final changes: if an edge contains node that donnt contain
        ## k, cpd, rn or organism_code, or both nodes of the same egde
        ## are duplicated,it's removed.

        edges_response = sapply(split(MetaboSignal_table, row(MetaboSignal_table)),
            find_unwanted_edge, organism_code = organism_code,
            expand_genes = expand_genes)
        edges_unwanted_global = grep("unwanted", edges_response)

        if (length(edges_unwanted_global) >= 1) {
            MetaboSignal_table = MetaboSignal_table[-c(edges_unwanted_global), ]
            MetaboSignal_table = matrix(MetaboSignal_table, ncol = 2)
        }

        if (is.matrix(MetaboSignal_table) == FALSE) {
            stop("Impossible to build a MetaboSignal_table with these paths",
                "\n")
        } else {
            if (nrow(MetaboSignal_table) == 0) {
                stop("Impossible to build a MetaboSignal_table with these paths")
            } else {
                MetaboSignal_table = unique(MetaboSignal_table)
                MetaboSignal_table = matrix(MetaboSignal_table, ncol = 2)
                colnames(MetaboSignal_table) = c("node1", "node2")
                rownames(MetaboSignal_table) = NULL
            }
        }

        ## Report network features#
        message()
        network_features(MetaboSignal_table)

        if (sum(paths_included == 0)) {
            paths_removed = names(paths_included[paths_included == 0])
            path_line = paste(paths_removed, collapse = ",")
            message()
            to_print = paste("Some path IDs were not used:",
                path_line, sep = "")
            warning(to_print, "\n")
        }
    }
    return(MetaboSignal_table)
}

