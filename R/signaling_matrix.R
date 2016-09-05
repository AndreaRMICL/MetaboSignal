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

#################### link_hsagene_ensembl ####################
link_hsagene_ensembl = function(hsa_id, gene_found, xx) {
    index = which(names(xx) == hsa_id)
    ensembl = unlist(xx[index])[1]
    # If there are several ensembl IDs only one is selected
    ensembl = as.character(ensembl)
    if (length(ensembl) >= 1) {
        gene_enS = c(gene_found, hsa_id, ensembl)
        return(gene_enS)
    }
}

#################### unique_symbol_ensembl ###################
unique_symbol_ensembl = function (symbol, gene_enSall){
  index = which(gene_enSall[, 2] == symbol)[1]
  return(index)
}


#################### filter_genes_tissue ####################
filter_genes_tissue = function(ensembl, tissue, m_value, f_value) {

    if (ensembl == m_value) {
        message(" -Progress:50% completed")
    }
    if (ensembl == f_value) {
        message(" -Progress:100% completed")
    }
    tissues = as.matrix(getHpa(ensembl, hpadata = "NormalTissue"))

    if ("Supportive" %in% as.character(tissues[, 6]) == FALSE) {
        tissues_detected = tissue  # Filtering ignored
    } else {
        indexSup = grep("Supportive", as.character(tissues[,6]))
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
        } else {
            # If the gene is detected in all tissues
            tissues_detected = tissues[, 2]
        }
    }
    tissues_detected = unique(as.character(tissues_detected))

    if (length(intersect(tissue, tissues_detected)) >= 1) {
        return("gene_wanted")

    } else {
        return("gene_unwanted")
    }
}

#################### signaling_matrix ####################
signaling_matrix = function(global_network_all, tissue, organism_code,
    organism_name, expand_genes) {

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
            nonbio_nodes = paste("fibrates", "non-steroids", "agents", "Antiinflammatory",
                                 "BR:br08303","Thiazolidinediones", sep = "|")
            nonbio_edges = c(grep(nonbio_nodes, network_edges[,1]),
                             grep(nonbio_nodes, network_edges[,2]))
            nonbio_edges = unique(nonbio_edges)

            if (length(nonbio_edges) >= 1) {
                network_edges = network_edges[-c(nonbio_edges), ]
            }

            ## Define gene nodes of the network
            nodes_network = unique(as.vector(network_edges))
            index_rno = grep(organism_code, nodes_network)
            nodes_network_rno = nodes_network[index_rno]

            ## Filter the network by tissue expression

            if (tissue[1] != "all") {
                to_print = paste("Preparing signaling genes to be",
                  "filtered by tissue: transforming", "gene IDs into",
                  "human ensembl IDs", sep = " ")
                message(to_print)

                ## Split genes in list of max 100 genes (for multiple query in KEGG API)
                genes_list = split(nodes_network_rno, ceiling(seq_along(nodes_network_rno)/100))

                ## Transform KEGG IDs into entrez IDs #
                entrez_res = unlist(lapply(genes_list, conv_entrez_kegg, source = "kegg"))
                entrez_ids = as.character(substr(entrez_res, 13, nchar(entrez_res)))
                start = unlist(gregexpr(pattern = organism_code, names(entrez_res)))
                kegg_ids = substr(names(entrez_res), start, nchar(names(entrez_res)))

                ## Transform KEGG gene IDs into entrez IDs
                if (organism_code == "hsa") {
                  hsa_id = entrez_ids
                  genes_found = kegg_ids

                  ## Tranform entrez IDs into ENSEMBL IDs
                  X = org.Hs.egENSEMBL
                  mapped_genes = mappedkeys(X)
                  xx = as.list(X[mapped_genes])
                  NL = list()
                  NL[["xx"]] = xx
                  hsa_gene_ensembl = mapply(link_hsagene_ensembl,
                    hsa_id, genes_found, MoreArgs = NL, SIMPLIFY = FALSE)
                  gene_enSall = unique(do.call(rbind, hsa_gene_ensembl))

                  if (is.matrix(gene_enSall) == FALSE) {
                    # Filtering ignored
                    gene_enSall = cbind(nodes_network_rno, nodes_network_rno,
                      nodes_network_rno)
                  }
                } else {
                  # If organism is not human Transform KEGG gene IDs into
                  # symbols
                  symbols = sapply(entrez_ids, conv_entrez_symbol,
                                   organism_name, source = "entrez")
                  index_NF = grep("NF", symbols)

                  if (length(index_NF) > 0) {
                    symbols = symbols[-index_NF]
                    genes_found = kegg_ids[-index_NF]
                  } else {
                    genes_found = kegg_ids
                  }

                  ## Transform gene symbols into ENSEMBL IDs
                  if (length(genes_found) > 0) {
                    symbols = toupper(symbols)
                    gene_enSall = try(suppressMessages(select(org.Hs.eg.db,
                      keys = symbols, columns = c("ENSEMBL"),
                      keytype = "SYMBOL")), silent = TRUE)
                    # if keys are not valid the function returns an error
                    # message.

                    if (grepl("Error", gene_enSall)[1] == FALSE) {
                      vector = as.character(nrow(gene_enSall))
                      gene_enSall = cbind(vector, gene_enSall)
                      gene_enSall = as.matrix(gene_enSall)

                      for (i in 1:length(symbols)) {
                        index = which(gene_enSall[, 2] == symbols[i])
                        gene_enSall[index, 1] = genes_found[i]
                      }
                      ## Remove duplicated ensembl IDs#
                      index_wanted = sapply(symbols, unique_symbol_ensembl,gene_enSall)

                      gene_enSall = gene_enSall[index_wanted, ]
                      gene_enSall = na.omit(gene_enSall)
                      gene_enSall = matrix(gene_enSall, ncol = 3)

                    } else {
                      # Filtering will be ignored
                      gene_enSall = cbind(nodes_network_rno, nodes_network_rno,
                                          nodes_network_rno)
                    }
                  } else {
                    # Filtering will be ignored
                    gene_enSall = cbind(nodes_network_rno, nodes_network_rno,
                      nodes_network_rno)
                  }
                }
                ## Find genes that are not expressed in the target tissue We
                ## only remove non-detected genes (reliability=suportive)
                ensembl = as.character(gene_enSall[, 3])

                if (grepl(organism_code, ensembl)[1] == TRUE) {
                  to_print = paste("Filtering by genes was ignored.",
                    "Check if the organism is", "consistent with path IDs")
                  message(to_print)

                } else {
                  m_value = ensembl[round(0.5 * length(ensembl))]
                  f_value = ensembl[length(ensembl)]
                  message("Filtering genes by tissue:")
                  response_filter = sapply(ensembl, filter_genes_tissue,
                    tissue, m_value = m_value, f_value = f_value)
                  index_unwanted = as.numeric(which("gene_unwanted" == response_filter))

                  if (length(index_unwanted) > 0) {
                    # there are undetected genes
                    genes_not_in_tissue = as.character(gene_enSall[index_unwanted, 1])

                    ## Remove edges that contain undetected genes
                    edges_not_in_tissueL = sapply(lapply(split(network_edges, row(network_edges)),
                                                         intersect, y = genes_not_in_tissue), length)
                    edges_not_in_tissueL = as.vector (edges_not_in_tissueL)
                    edges_not_in_tissue = unique(which(edges_not_in_tissueL >= 1))

                    if (nrow(network_edges) > length(edges_not_in_tissue)) {
                      network_edges = network_edges[-c(edges_not_in_tissue), ]

                    } else {
                      # None of the edges can be kept after filtering
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
            nodes_network = unique(as.vector(network_edges_new_2))
            index_rno = grep(organism_code, nodes_network)
            nodes_network_rno = nodes_network[index_rno]

            ## Add orthology IDs
            if (expand_genes == FALSE) {
                message()
                to_print = ("Transforming gene IDs into orthology IDs")
                message(to_print)

                file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code, sep = "")
                response_ko = getURL(file_ko)
                koTable = convertTable(response_ko)
                koTable[, 2] = substr(koTable[, 2], 4, 9)

                for (r in 1:nrow(network_edges_new_2)) {
                  for (c in 1:ncol(network_edges_new_2)) {
                    index = which(koTable[, 1] == network_edges_new_2[r, c])
                    if (length(index) > 0) {
                      network_edges_new_2[r, c] = koTable[index, 2]
                    }
                  }
                }
            }
            signaling_table = unique(network_edges_new_2)
            signaling_table = matrix(signaling_table, ncol = 2)
            colnames(signaling_table) = c("node1", "node2")
        }
    }
    return(signaling_table)
}

