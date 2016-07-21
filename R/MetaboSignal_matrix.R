globalVariables("directionality_reactions")

#################### INTERNAL FUNCTIONS ####################

#################### metabolic_matrix ####################
metabolic_matrix = function(path_names, list_parsed_paths, organism_code,
                            expand_genes) {

    metabolic_table_RG = c()
    metabolic_table = c()

    for (path in path_names) {
        lines = list_parsed_paths[[path]]

        ## Select reactions
        global_lines_reactions = intersect(grep("rn",lines), grep("Name",lines))
        Reactions = lines[c(global_lines_reactions)]
        Reactions = substr(Reactions, 11, nchar(Reactions))

        ## Select substrates
        global_lines_substrates = grep("Substrate Name", lines)
        Substrates = lines[c(global_lines_substrates)]
        Substrates = substr(Substrates, 21, nchar(Substrates))

        #Select products
        global_lines_products = grep("Product Name", lines)
        Products = lines[c(global_lines_products)]
        Products = substr(Products, 19, nchar(Products))

        #Select type
        global_lines_type = grep("Type", lines)
        Type = lines[c(global_lines_type)]
        Type = substr(Type, 11, nchar(Type))

        ## Build network edges
        all_edges = c()
        for (i in 1:length(Reactions)) {
            reactions = unlist(strsplit(Reactions[i], " "))
            substrates = unlist(strsplit(Substrates[i], ";"))
            substrates = unlist(strsplit(substrates, " "))
            products = unlist(strsplit(Products[i], ";"))
            products = unlist(strsplit(products, " "))
            type = Type[i]

            # Correct for direction based on Duarte et al., 2007
            for (reaction in reactions) {
                a = which(directionality_reactions[,1] == reaction)
                if (length(a) >= 1) {
                  type = directionality_reactions[,2][a]
                }
            }

            # Build edges
            for (reaction in reactions) {
                substrate_lines = cbind(substrates, rep(reaction,
                                                        length(substrates)))
                product_lines = cbind(rep(reaction, length(products)), products)
                edges_reaction = rbind(substrate_lines, product_lines)

                if (type == "reversible") {
                  reverse = cbind(edges_reaction[, 2], edges_reaction[, 1])
                  edges_reaction = rbind(edges_reaction, reverse)
                }
                all_edges = rbind(all_edges, edges_reaction)
            }
        }
        metabolic_table = rbind(metabolic_table, all_edges)
        metabolic_table = unique(metabolic_table)
        metabolic_table = matrix(metabolic_table, ncol = 2)
    }

    ## Link reactions to genes (organism specific or orthology IDs)
    metabolic_table_RG = c()

    file_enzyme = paste("http://rest.kegg.jp/link/enzyme/", organism_code, sep = "")
    response_enzyme = getURL(file_enzyme)
    enzymeTable = convertTable(response_enzyme)

    file_reaction = "http://rest.kegg.jp/link/reaction/enzyme"
    response_reaction = getURL(file_reaction)
    reactionTable = convertTable(response_reaction)

    file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code, sep="")
    response_ko = getURL(file_ko)
    koTable = convertTable(response_ko)
    koTable[, 2] = substr(koTable[, 2], 4, 9)

    if (is.matrix(enzymeTable) & is.matrix(reactionTable) & is.matrix(koTable)) {
        to_print = ("Linking reactions to genes")
        message(to_print)

        genes = enzymeTable[, 1]
        enzymes = enzymeTable[, 2]
        ko = vector()
        length(ko) = length(genes)

        for (i in 1:length(genes)) {
            index = which(koTable[, 1] == genes[i])
            if (length(index) > 0) {
                ko[i] = koTable[index[1], 2]
            }
        }
        reactionM = c()
        for (z in 1:length(enzymes)) {
            index = which(reactionTable[, 1] == enzymes[z])
            if (length(index) > 0) {
                reactions = reactionTable[index, 2]
                line = cbind(reactions, rep(enzymes[z], length(reactions)),
                             rep(genes[z], length(reactions)), rep(ko[z],
                                                        length(reactions)))
                reactionM = rbind(reactionM, line)
            }
        }
        reactionko = unique(na.omit(reactionM[, c(1, 4)]))
        reactiongene = unique(na.omit(reactionM[, c(1, 3)]))

        for (y in 1:nrow(metabolic_table)) {
            reaction_genes_all = c()
            edge = as.character(metabolic_table[y, ])
            index_rn = grep("rn", edge)
            indexM = c()
            if (index_rn == 1) {
                indexM = 2
            } else {
                indexM = 1
            }

            reaction = edge[index_rn]
            metabolite = edge[indexM]

            if (expand_genes == FALSE) {# Genes represent orthology IDs.
                index_reaction = which(reactionko[, 1] == reaction)

                if (length(index_reaction) > 0) {
                  reaction_genes = paste(reaction, reactionko[index_reaction, 2],
                                         sep = "_")
                } else {
                  reaction_genes = reaction
                }
                reaction_genes_all = c(reaction_genes_all, reaction_genes)
            } else {# Genes represent organism specific IDs.
                index_reaction = which(reactiongene[, 1] == reaction)

                if (length(index_reaction) > 0) {
                  reaction_genes = paste(reaction, reactiongene[index_reaction, 2],
                                         sep = "_")
                } else {
                  reaction_genes = reaction
                }
                reaction_genes_all = c(reaction_genes_all, reaction_genes)
            }

            if (index_rn == 1) {
                new_edges = cbind(reaction_genes_all, rep(metabolite,
                  length(reaction_genes_all)))
            } else {
                new_edges = cbind(rep(metabolite, length(reaction_genes_all)),
                  reaction_genes_all)
            }
            metabolic_table_RG = rbind(metabolic_table_RG, new_edges)
        }

        metabolic_table_RG = unique(metabolic_table_RG)
        metabolic_table_RG = matrix(metabolic_table_RG, ncol = 2)
        colnames(metabolic_table_RG) = c("node1", "node2")
        rownames(metabolic_table_RG) = NULL
    }
    return(metabolic_table_RG)
}


#################### signaling_matrix ####################
signaling_matrix = function(global_network_all, tissue, organism_code,
                            organism_name, expand_genes) {

    signaling_table = c()

    ## Create network_edges matrix
    nodes = nodes(global_network_all)
    edges = KEGGgraph::edges(global_network_all)
    network_edges = c()

    for (i in 1:length(nodes)) {
        linked_nodes = edges[[i]]
        if (length(linked_nodes) >= 1) {# Node with edges
            new_edge = cbind(rep(nodes[i], length(linked_nodes)), edges[[i]])
            network_edges = rbind(network_edges, new_edge)
        }
    }
    if (length(network_edges) == 0) {
        to_print = ("Impossible to build a signaling network")
        warning(to_print, "\n")
    } else {
      network_edges=matrix(network_edges, ncol = 2)
        ## Remove the edges that contain paths
        all_lines_paths = c(grep("path", network_edges[,1]),
                            grep("path", network_edges[,2]))

        if (length(all_lines_paths) == nrow(network_edges)) {
            # All edges contain paths.
            to_print = ("Impossible to build a signaling network")
            warning(to_print, "\n")
        } else {
            if (length(all_lines_paths) >= 1) {
                network_edges = network_edges[-(all_lines_paths), ]
                network_edges = matrix(network_edges, ncol = 2)
            }

            ## Remove non-biological nodes from KEGG#
            nonbio_nodes = c("fibrates)", "non-steroids)", "agents,",
                "Antiinflammatory", "BR:br08303(A10BG", "BR:br08303(C10AB",
                "BR:br08303(S01BC", "Thiazolidinediones)")
            nonbio_edges = c()
            for (a in 1:length(nonbio_nodes)) {
                answer = which(network_edges == nonbio_nodes[a], arr.ind = TRUE)[, 1]
                if (length(answer) >= 1) {
                  nonbio_edges = c(nonbio_edges, answer)
                }
            }
            nonbio_edges = unique(nonbio_edges)

            if (length(nonbio_edges) >= 1) {
                network_edges = network_edges[-c(nonbio_edges), ]
            }

            ## Define gene nodes of the network
            nodes_network = as.vector(network_edges)
            nodes_network = unique(nodes_network)
            nodes_network = sort(nodes_network, decreasing = FALSE)
            answer = grepl(organism_code, nodes_network)
            index_rno = which(answer == "TRUE")[1]

            nodes_network_rno = nodes_network[index_rno:length(nodes_network)]

            ## Filter the network by tissue expression

            if (tissue[1] != "all") {
                to_print = paste("Preparing signaling genes to be",
                  "filtered by tissue: transforming", "gene IDs into",
                  "human ensembl IDs", sep = " ")
                message(to_print)

                ## Transform KEGG gene IDs into entrez IDs
                if (organism_code == "hsa") {
                  hsa_id = as.character(substr(nodes_network_rno, 5,
                                               nchar(nodes_network_rno)))
                  genes_found = nodes_network_rno

                  ## Tranform entrez IDs into ENSEMBL IDs
                  X = org.Hs.egENSEMBL
                  mapped_genes = mappedkeys(X)
                  xx = as.list(X[mapped_genes])

                  gene_enSall = c()
                  for (s in 1:length(hsa_id)) {
                    index = which(names(xx) == hsa_id[s])
                    ensembl = unlist(xx[index])[1]
                    # If there are several ensembl IDs only one is selected
                    ensembl = as.character(ensembl)

                    if (length(ensembl) >= 1) {
                      gene_enS = c(genes_found[s], hsa_id[s], ensembl)
                      gene_enSall = rbind(gene_enSall, gene_enS)
                    }
                  }
                  if (is.matrix(gene_enSall) == FALSE) {# Filtering ignored
                    gene_enSall = cbind(nodes_network_rno, nodes_network_rno,
                                        nodes_network_rno)
                  }
                } else {# If organism is not human
                    ## Transform KEGG gene IDs into symbols
                  symbols = c()
                  genes_found = c()
                  for (x in 1:length(nodes_network_rno)) {
                    gene = nodes_network_rno[x]
                    gene = substr(gene, 5, nchar(gene))
                    ID = try(as.matrix(unlist(query(q = gene,
                      size = 1, species = organism_name)$hits)), silent = TRUE)
                    if (grepl("Error", ID)[1] == FALSE) {
                      rows = rownames(ID)
                      ID = as.character(ID)
                      index = which(rows == "symbol")
                      if (length(index) > 0) {
                        symbols = c(symbols, ID[index])
                        genes_found = c(genes_found, nodes_network_rno[x])
                      }
                    }
                  }
                  ## Transform gene symbols into ENSEMBL IDs
                  if (length(genes_found) > 0) {
                    symbols = toupper(symbols)
                    gene_enSall = try(suppressMessages(select(org.Hs.eg.db,
                      keys = symbols, columns = c("ENSEMBL"),
                      keytype = "SYMBOL")), silent = TRUE)

                    if (grepl("Error", gene_enSall)[1] == FALSE) {
                      vector = as.character(nrow(gene_enSall))
                      gene_enSall = cbind(vector, gene_enSall)
                      gene_enSall = as.matrix(gene_enSall)

                      for (i in 1:length(symbols)) {
                        index = which(gene_enSall[, 2] == symbols[i])
                        gene_enSall[index, 1] = genes_found[i]
                      }
                      ## Remove duplicated ensembl IDs#
                      index_wanted = c()
                      for (symbol in symbols) {
                        index = which(gene_enSall[, 2] == symbol)[1]
                        index_wanted = c(index_wanted, index)
                      }
                      gene_enSall = gene_enSall[index_wanted, ]
                      gene_enSall = na.omit(gene_enSall)
                      gene_enSall = matrix(gene_enSall, ncol = 3)
                    } else {# Filtering will be ignored
                      gene_enSall = cbind(nodes_network_rno, nodes_network_rno,
                                          nodes_network_rno)
                    }
                  } else {# Filtering will be ignored
                    gene_enSall = cbind(nodes_network_rno, nodes_network_rno,
                                        nodes_network_rno)
                  }
                }

                ## Find genes that are not expressed in the target tissue
                ## We only remove non-detected genes (reliability=suportive)
                ensembl = as.character(gene_enSall[, 3])
                if (grepl(organism_code, ensembl)[1] == TRUE) {
                  to_print = paste("Filtering by genes was ignored.",
                    "Check if the organism is", "consistent with path IDs")
                  message(to_print)

                } else {
                  index_wanted = c()
                  index_unwanted = c()
                  edges_not_in_tissue = c()

                  message("Filtering genes by tissue:")
                  pr_value = round(0.5 * length(ensembl))

                  for (t in 1:length(ensembl)) {
                    tissues_detected = c()
                    # Report progress
                    if (t == pr_value) {
                      message(" -Progress:50% completed")
                    }
                    if (t == length(ensembl)) {
                      message(" -Progress:100% completed")
                    }
                    tissues = as.matrix(getHpa(ensembl[t], hpadata = "NormalTissue"))

                    if ("Supportive" %in% as.character(tissues[, 6]) == FALSE) {
                      tissues_detected = tissue # Filtering ignored
                    } else {
                      indexSup = grep("Supportive", as.character(tissues[, 6]))
                      tissues = tissues[indexSup, ]
                      undetected = c(grep("Not detected", tissues[, 4]))

                      if (length(undetected) >= 1) {
                        # if the gene is undetected in some tissues
                        if (nrow(tissues) - length(undetected) > 1) {
                          tissues_detected = tissues[-c(undetected), 2]
                        } else if (nrow(tissues) - length(undetected) == 1) {
                          tissues_detected = tissues[-c(undetected), ][2]
                        } else if (nrow(tissues) - length(undetected) == 0) {
                          tissues_detected = "undetected"
                        }
                      } else { # If the gene is detected in all tissues
                        tissues_detected = tissues[, 2]
                      }
                    }
                    tissues_detected = unique(as.character(tissues_detected))

                    if (length(intersect(tissue, tissues_detected)) >= 1) {
                      index_wanted = c(index_wanted, t)

                    } else {
                      index_unwanted = c(index_unwanted, t)
                    }
                  }
                  if (length(index_unwanted) > 0) { # there are undetected genes
                    genes_not_in_tissue = as.character(gene_enSall[index_unwanted, 1])

                    ## Remove edges that contain undetected genes
                    edges_not_in_tissue = c()

                    for (gene in genes_not_in_tissue) {
                      index_unwanted_edge = which(network_edges == gene, arr.ind = TRUE)[,1]
                      edges_not_in_tissue = c(edges_not_in_tissue, index_unwanted_edge)
                    }

                    edges_not_in_tissue = unique(edges_not_in_tissue)

                    if (nrow(network_edges) > length(edges_not_in_tissue)) {
                      network_edges = network_edges[-c(edges_not_in_tissue), ]

                    } else { #  None of the edges can be kept after filtering
                      to_print = paste("Filtering was ignored because",
                        "none of the signaling genes is", "detected in the",
                        "target tissue")
                      warning(to_print, "\n")
                    }
                  }
                }
            }
            network_edges_new_2 = matrix(network_edges, ncol = 2)

            ## Redefine nodes of the network
            nodes_network = as.vector(network_edges_new_2)
            nodes_network = unique(nodes_network)
            nodes_network = sort(nodes_network, decreasing = FALSE)
            answer = grepl(organism_code, nodes_network)
            index_rno = which(answer == "TRUE")[1]

            nodes_network_rno = nodes_network[index_rno:length(nodes_network)]

            ## Add orthology IDs
            if (expand_genes == FALSE) {
                message()
                to_print = ("Transforming gene IDs into orthology IDs")
                message(to_print)

                file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code,
                                sep = "")
                response_ko = getURL(file_ko)
                koTable = convertTable(response_ko)
                koTable[, 2] = substr(koTable[, 2], 4, 9)

                ko_genesM = c()
                for (gene in nodes_network_rno) {
                  index_gene = which(koTable[, 1] == gene)
                  if (length(index_gene) > 0) {
                    ko_line = c(gene, koTable[index_gene[1], 2])
                    ko_genesM = rbind(ko_genesM, ko_line)
                  }
                  rownames(ko_genesM) = NULL
                  ko_genesM = unique(ko_genesM)
                }

                for (r in 1:nrow(network_edges_new_2)) {
                  for (c in 1:ncol(network_edges_new_2)) {
                    index = which(ko_genesM[, 1] == network_edges_new_2[r, c])
                    if (length(index) > 0) {
                      network_edges_new_2[r, c] = ko_genesM[index, 2]
                    }
                  }
                }
            }
            network_edges_new_2 = unique(network_edges_new_2)
            signaling_table = network_edges_new_2
            signaling_table = matrix(signaling_table, ncol = 2)
            colnames(signaling_table) = c("node1", "node2")
        }
    }

    return(signaling_table)
}

#################### find_duplicated_nodes ####################
#Remove edges with node1=node2
find_duplicated_nodes = function(network) {
  duplicated_nodes = c()
  for (r in 1:nrow(network)) {
    if (network[r, 1] == network[r, 2]) {
      duplicated_nodes = c(duplicated_nodes, r)
    }
  }
  return(duplicated_nodes)
}

###########################################################
#################### EXTERNAL FUNCTION ####################

#################### MetaboSignal_matrix ###################

MetaboSignal_matrix = function(metabo_paths = NULL, signaling_paths = NULL,
    organism_name = NULL, tissue = "all", expand_genes = FALSE) {

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
            empty_paths = c()
            for (p in 1:length(list_parsed_paths)) {
                if (length(list_parsed_paths[[p]]) <= 1) {
                  empty_paths = c(empty_paths, p)
                }
            }
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
        bad_path = c()
        for (p in 1:length(signaling_paths)) {
            answer = signaling_paths[p] %in% all_paths
            if (answer == FALSE) {
                bad_path = c(bad_path, p)
                paths_removed = c(paths_removed, signaling_paths[p])
                to_print = paste(signaling_paths[p], "-incorrect path ID:",
                                 "path removed", sep = "")
                message(to_print)
            }
        }
        if (length(bad_path) >= 1) {
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

        ## Final changes if a node doesnt contain k, cpd, rn or
        ## organism_code, it's removed.
        node1 = MetaboSignal_table[, 1]
        node2 = MetaboSignal_table[, 2]
        all_edges_unwanted_1 = c()
        all_edges_unwanted_2 = c()
        edges_unwanted_global = c()

        for (p in 1:length(node1)) {
            if (grepl("K", node1[p]) == FALSE & grepl("cpd",
                node1[p]) == FALSE & grepl(organism_code, node1[p]) ==
                FALSE & grepl("rn", node1[p]) == FALSE) {
                all_edges_unwanted_1 = c(all_edges_unwanted_1, p)
            }
        }
        for (q in 1:length(node2)) {
            if (grepl("K", node2[q]) == FALSE & grepl("cpd",
                node2[q]) == FALSE & grepl(organism_code, node2[q]) ==
                FALSE & grepl("rn", node2[q]) == FALSE) {
                all_edges_unwanted_2 = c(all_edges_unwanted_2, q)
            }
        }

        duplicated_nodes = find_duplicated_nodes(MetaboSignal_table)
        edges_unwanted_global = unique(c(all_edges_unwanted_1,
                                         all_edges_unwanted_2, duplicated_nodes))

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

