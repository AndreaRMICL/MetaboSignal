#################### INTERNAL FUNCTIONS ####################

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
# This function clusters organism-specific KEGG gene IDs into orthology IDs
orthology_clustering=function(network_table, organism_code, all_genes) {

    to_print = ("Clustering gene IDs into orthology IDs")
    message(to_print,"\n" )
    index_genes = grep(organism_code, all_genes)
    gene_nodes = unique(all_genes[index_genes])

    file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code, sep = "")
    response_ko = getURL(file_ko)
    koTable = convertTable(response_ko)
    koTable[, 2] = substr(koTable[, 2], 4, 9)

    ko_M = c()
    for (gene in gene_nodes) {
        index_gene = which(koTable[, 1] == gene)
        if (length(index_gene) > 0) {
            ko_line = c(gene, koTable[index_gene[1], 2])
            ko_M = rbind(ko_M, ko_line)
        }
        rownames(ko_M) = NULL
        ko_M = unique(ko_M)
    }
    for (i in 1:nrow(network_table)) {
        for (z in 1:ncol(network_table)) {
            index = which(ko_M[, 1] == network_table[i, z])
            if (length(index) > 0) {
                network_table[i, z] = ko_M[index, 2]
            }
        }
    }
    network_table=unique(network_table)
    return(network_table)
}


#################### EXTERNAL FUNCTION ####################

############### MetaboSignal_NetworkCytoscape #############

MetaboSignal_NetworkCytoscape = function(network_table, organism_code,
                                         organism_name, source_genes,
                                         target_metabolites, mode = "SP",
                                         type = "first", distance_th = Inf,
                                         collapse_genes = FALSE, names = TRUE,
                                         export_cytoscape = TRUE,
                                         file_name = "Cytoscape") {

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
        distance_table = MetaboSignal_distances(network_table, organism_code,
                                                organism_name, mode = mode,
                                                names = FALSE)
        all_genes = rownames(distance_table)
        ortho_answer = unique(grepl(organism_code, all_genes))
        ortho_answer = TRUE %in% ortho_answer

        if (ortho_answer == FALSE) {
            to_print = paste("Collapse was ignored because the gene nodes were",
                             "already collapsed")
            warning(to_print, "\n")

        } else {
            network_table=orthology_clustering(network_table, organism_code,
                                               all_genes)
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

    if (ortho_answer == TRUE) {# gene nodes are organism specific IDs.
        # Check if the source genes represent entrez IDs or symbols
        if (grepl(organism_code, source_genes[1]) == FALSE) {
            source_genes = try(MS_GetKEGG_GeneID(genes = source_genes,
                organism_name = organism_name, orthology = FALSE), silent = TRUE)
        }
    } else {# Gene nodes are orthology IDs.
        if ((substr(source_genes[1], 1, 1) == "K" & nchar(source_genes[1]) ==
            6) == FALSE) { # The input genes are entrezIDs or symbols
            source_genes = try(MS_GetKEGG_GeneID(genes = source_genes,
                organism_name = organism_name, orthology = TRUE), silent = TRUE)
        }
    }

    ##  Check if it was possible to map source_genes onto KEGG IDs
    if (grepl("Error",source_genes)[1]) {
        stop ("None of the source_genes is reported in KEGG. Could be that the organism_name is incorrect")
    }

    ## Check if the source_genes and target_metabolites can be mapped onto the
    # network#
    all_nodes = unique(as.vector(network_table))
    index_unwanted = c()

    for (i in 1:length(source_genes)) {
        answer = source_genes[i] %in% all_nodes
        if (answer == FALSE) {
            index_unwanted = c(index_unwanted, i)
        }
    }
    if (length(index_unwanted) >= 1) {
        source_genes = source_genes[-c(index_unwanted)]

        if (length(source_genes) == 0) {
            stop_source_genes = paste ("None of the source_genes is present in",
              "the network. Could be that the organism_code is incorrect")
            stop(stop_source_genes)
        }
    }
    index_unwanted = c()
    for (i in 1:length(target_metabolites)) {
        answer = target_metabolites[i] %in% all_nodes
        if (answer == FALSE | grepl("cpd:",target_metabolites[i]) == FALSE) {
            index_unwanted = c(index_unwanted, i)
        }
    }
    if (length(index_unwanted) >= 1) {
        metabolites_not_found = target_metabolites[index_unwanted]
        target_metabolites = target_metabolites[-c(index_unwanted)]

        if (length(target_metabolites) == 0) {
            stop("None of the target_metabolites is present in the network")
        }
    }

    ## Calculate shortest path network
    message("Building shortest path network","\n")
    all_pathsGM = c()
    for (metabolite in target_metabolites) {

        if (mode == "SP") {
            network_i = igraph_edited(network_table, metabolite)
        } else {
            network_i = graph.data.frame(network_table, directed = TRUE)
        }

        for (gene in source_genes) {
            index_gene = which(rownames(distance_table) == gene)
            index_metabolite = which(colnames(distance_table) ==
                metabolite)
            distanceGM = distance_table[index_gene, index_metabolite]

            if (distanceGM < distance_th) {# If distance is not Inf,
                empty_matrix = c()

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
                    ASP = get.all.shortest.paths(network_i, gene, metabolite,
                                                 mode = "all")

                  } else {
                    ASP = get.all.shortest.paths(network_i, gene, metabolite,
                                                 mode = "out")
                  }
                  ASP = ASP$res
                }

                all_paths = c()
                for (i in 1:length(ASP)) {
                  shortpath = rownames(as.matrix(unlist(ASP[[i]])))
                  all_paths = rbind(all_paths, shortpath)
                  all_paths = unique(all_paths)
                  rownames(all_paths) = NULL
                }

                ## Select bw-ranked shortest path of multiple shortest paths
                if (type == "bw" & nrow(all_paths) > 1) {
                  all_nodes = unique(as.vector(all_paths))
                  index_cpd = grep("cpd:", all_nodes)

                  if (length(index_cpd) > 0) {
                    # This should be always true because if there are not
                    # compounds in the network the function
                    # MetaboSignal_distances wont work.
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

                  Global_BW_score = c()
                  for (i in 1:nrow(all_paths)) {
                    BW = c()
                    path_individual = as.character(all_paths[i,])

                    for (z in 1:length(path_individual)) {
                      index = which(BW_matrix[, 1] == path_individual[z])
                      if (length(index) > 0) {
                        BW = c(BW, as.numeric(BW_matrix[index, 2]))
                      }
                    }
                    BW = as.numeric(BW)
                    score_BW = as.numeric(BW)
                    score_BW = sum(BW)/length(BW)
                    Global_BW_score = c(Global_BW_score, score_BW)
                  }
                  maxBW = which(Global_BW_score == max(Global_BW_score))[1]
                  path = all_paths[maxBW, ]
                  path = as.character(path)
                  all_paths = matrix(path, ncol = length(path))
                }

                ## Transform the shortest path matrix into a network table
                all_pathsList = split(all_paths, row(all_paths))
                all_pathsnetworkList = lapply(all_pathsList, path_as_network)
                subnetwork = unique(do.call(rbind, all_pathsnetworkList))
                rownames(subnetwork) = NULL  # subnetwork is a network-table

                ## Correct directionality in case of SP mode
                if (mode == "SP") {
                  for (i in 1:nrow(subnetwork)) {
                    supernetworkDirection = rbind(subnetwork,
                      network_table)
                    intersection_subnetwork =
                        supernetworkDirection[duplicated(supernetworkDirection), ]
                    intersection_subnetwork = matrix(intersection_subnetwork, ncol = 2)
                    rows = split(intersection_subnetwork, row(intersection_subnetwork))
                    edge = subnetwork[i, ]

                    if (Match(rows, edge) == FALSE) {
                      ## This will correct for the directionality adjustment.
                      ## For example if the target node it's the sustrate of the
                      ## irreversible reaction, it should be ploted as substrate
                      ## not as product.
                      subnetwork[i, ] = rev(edge)
                    }
                  }
                }

                ## Add edges for reversible interactions
                reverseSubnetwork = cbind(subnetwork[, 2], subnetwork[, 1])
                supernetwork = rbind(reverseSubnetwork, network_table)
                if (mode == "all") {# Build an undirected network
                  subnetwork = rbind(subnetwork, reverseSubnetwork)
                  all_pathsGM = rbind(all_pathsGM, subnetwork)
                  all_pathsGM = unique(all_pathsGM)
                  rownames(all_pathsGM) = NULL
                } else {
                  intersection_reverseSubnetwork =
                      supernetwork[duplicated(supernetwork),]
                  intersection_reverseSubnetwork =
                      matrix(intersection_reverseSubnetwork, ncol = 2)

                  if (nrow(intersection_reverseSubnetwork) > 0) {
                    subnetwork = rbind(subnetwork, intersection_reverseSubnetwork)
                  }
                  all_pathsGM = rbind(all_pathsGM, subnetwork)
                  all_pathsGM = unique(all_pathsGM)
                  colnames(all_pathsGM) = c("node1", "node2")
                  rownames(all_pathsGM) = NULL
                }
            }
        }
    }

    if (length(all_pathsGM) >= 1) {
        all_pathsGM_names = all_pathsGM
        if (names == TRUE) {
            all_nodes_network = unique(as.vector(all_pathsGM))
            all_nodes_names = MS_ChangeNames(all_nodes_network, organism_code)

            for (i in 1:length(all_nodes_names)) {
                all_pathsGM_names[all_pathsGM_names == all_nodes_network[i]] =
                    all_nodes_names[i]
            }

        } else (all_pathsGM_names = all_pathsGM)

        ## Export network in cytoscape format
        namescyto = names

        if (export_cytoscape == TRUE) {
            message("Creating cytoscape files","\n")
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

