#################### MetaboSignal_matrix ###################

MetaboSignal_matrix = function(metabo_paths = NULL, signaling_paths = NULL,
                               organism_name = NULL, tissue = "all",
                               expand_genes = FALSE) {

    ########## 0)Preparatory steps###########

    ## Check if metabo_paths or signaling_paths exists
    if (is.vector(metabo_paths) == FALSE & is.vector(signaling_paths) == FALSE) {
        stop("At least one metabo_path or one signaling_path is required")
    }

    ## Check that all the paths belong to the same organism
    input_paths = c(metabo_paths, signaling_paths)
    organism_code = substr(input_paths, 1, 3)
    organism_code = unique(organism_code)
    if (length(organism_code) > 1) {
        stop("All paths have to belong to the same organism: check path IDs")
    }

    ## Check if organism_name exits
    if (length(organism_name) == 0 & tissue[1] != "all"){
      stop ("An organism_name is required for tissue filtering")
    }

    # Iniciate paths_removed
    paths_removed = c()

    ## Get all paths of the organism of interest
    file = paste("http://rest.kegg.jp/list/pathway/", organism_code, sep = "")
    options(warn = -1)
    lines = try(readLines(file), silent = TRUE)
    options(warn = 0)

    if (grepl ("Error", lines)[1]) { # example: metabo_paths = "X"
      stop("Incorrect path IDs")
    }
    all_paths = substr(lines, 6, 13)

    ## Iniciate metabolic and signaling networks
    metabolic_table_RG = c()
    signaling_table = c()

    ########## 1)Build metabolic table###########

    if (length(metabo_paths) >= 1) {

        to_print = ("Building metabolic table")
        message(to_print)
        to_print = ("Reading paths:")
        message(to_print)

        ### Get KGML files and transform them into reaction files####
        path_names = c()
        paths_removed = c()
        list_parsed_paths = list()

        for (path in metabo_paths) {
            message(path)
            if (path %in% all_paths) {
                if (substr(path, 4, nchar(path)) == "01100") { #remove metabolic pathways map
                  paths_removed=c(paths_removed, path)
                } else {
                  # Check that the input path exists
                  file = paste("http://rest.kegg.jp/get/", path, "/kgml", sep = "")
                  pathway = try(getURL(file), silent = TRUE)
                  reactions = try(getReactions(parseKGML(pathway)), silent = TRUE)
                  Org_Name = try(parseKGML(pathway)@pathwayInfo@org, silent = TRUE)
                  pathway_name = try(parseKGML(pathway)@pathwayInfo@title, silent = TRUE)
                  pathway_name_cleaned = try(gsub(" ", "_", gsub("[^[:alnum:] ]","",
                                                                 pathway_name)),
                                             silent = TRUE)
                  name = paste(Org_Name, pathway_name_cleaned, sep = "_")
                  parsed_path = capture.output(reactions, file = NULL)
                  path_names = c(path_names, name)

                  if (grepl("Error", parsed_path)[1] == TRUE) {
                    to_print = paste(path, "-path ID without XML:path removed",
                                     sep = "")
                    message(to_print)
                    paths_removed = c(paths_removed, path)
                    list_parsed_paths[[name]] = c()
                  } else {
                    if(length(parsed_path) <= 1) {
                      paths_removed = c(paths_removed,path)
                    }
                    list_parsed_paths[[name]] = parsed_path
                  }
                }
            } else {
                to_print = paste(path, "-incorrect path ID:path removed",
                                 sep = "")
                message(to_print)
                paths_removed = c(paths_removed, path)
            }
        }
        if (length(list_parsed_paths) == 0) {
            to_print = ("Impossible to build a metabolic network")
            warning(to_print, "\n")

        } else {
            ## Remove empty paths: for example oxidative phosphorylation#
            length_paths = sapply(list_parsed_paths,length)
            empty_paths = which(length_paths <= 1)

            if (length(empty_paths) >= 1) {
                list_parsed_paths = list_parsed_paths[-c(empty_paths)]
                path_names = path_names[-c(empty_paths)]
            }
            if (length(list_parsed_paths) == 0) {
                to_print = ("Impossible to build a metabolic network")
                warning(to_print, "\n")
            } else {

                ### Create metabolic_table_RG ####
                metabolic_table_RG = metabolic_matrix(path_names,
                                                      list_parsed_paths,
                                                      organism_code,
                                                      expand_genes)
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

        ## Remove bad paths
        response_paths = sapply (signaling_paths, find_bad_path, all_paths = all_paths)
        bad_path = grep("bad", response_paths)

        if (length(bad_path) >= 1) {
            paths_removed = c(paths_removed, signaling_paths [ c(bad_path)])
            signaling_paths = signaling_paths[-c(bad_path)]
        }
        if (length(signaling_paths) == 0) {
            to_print = ("Impossible to build a signaling network")
            warning(to_print, "\n")
        } else {
            network_list = list()
            for (path in signaling_paths) {
                file = paste("http://rest.kegg.jp/get/", path, "/kgml", sep = "")
                pathway = try(getURL(file), silent = TRUE)
                path_parsed = try(parseKGML(pathway), silent = TRUE)
                path_network = try(KEGGpathway2Graph(path_parsed,
                                                     genesOnly = FALSE,
                                                     expandGenes = TRUE),
                                   silent = TRUE)
                urlcheck = try(edges(path_network), silent = TRUE)

                if (grepl("Error", urlcheck)[1] == TRUE) {
                  paths_removed = c(paths_removed, path)
                  to_print = paste(path, "-path ID without XML:path removed",
                    sep = "")
                  message(to_print)
                } else {
                  message(path)
                  network_list[[path]] = path_network
                }
            }
            if (length(network_list) == 0) {
                to_print = ("Impossible to build a signaling network")
                warning(to_print, "\n")
            } else {
                global_network_all = mergeGraphs(network_list)
                global_network_all

                ### Create signaling_table ####
                signaling_table = signaling_matrix(global_network_all,
                  tissue, organism_code, organism_name, expand_genes)
            }
        }
    }

    ################## End signal###################

    ########## 3)Build MetaboSignal table ###########

    MetaboSignal_table = c()
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
        stop("Impossible to build a MetaboSignal_table with these paths","\n")
    } else {
        for (i in 1:nrow(MetaboSignal_table)) {
            for (z in 1:ncol(MetaboSignal_table)) {
                node = MetaboSignal_table[i, z]
                if (grepl("_", node) == TRUE) {
                  geneID = substr(node, 11, nchar(node))
                  MetaboSignal_table[i, z] = geneID
                }
            }
        }

        ## Final changes: if an edge contains node that donnt contain k, cpd, rn or
        ## organism_code, or both nodes of the same egde are duplicated,it's removed.

        edges_response = sapply(split(MetaboSignal_table, row(MetaboSignal_table)),
                                find_unwanted_edge, organism_code = organism_code, 
                                expand_genes = expand_genes)
        edges_unwanted_global = grep("unwanted", edges_response)

        if (length(edges_unwanted_global) >= 1) {
            MetaboSignal_table = MetaboSignal_table[-c(edges_unwanted_global), ]
            MetaboSignal_table = matrix(MetaboSignal_table, ncol = 2)
        }

        if (is.matrix(MetaboSignal_table) == FALSE) {
            stop("Impossible to build a MetaboSignal_table with these paths","\n")
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

        if (length(paths_removed) > 0) {
            path_line = c()
            for (path in paths_removed) {
                path_line = paste(path_line, path, sep = ",")
            }
            path_line = substr(path_line, 2, nchar(path_line))
            message()
            to_print = paste("Some path IDs were not used:", path_line,
                             sep = "")
            warning(to_print, "\n")
        }
    }
    return(MetaboSignal_table)
}

