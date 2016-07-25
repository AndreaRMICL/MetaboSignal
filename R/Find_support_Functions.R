################ MS_FindKEGG ####################

MS_FindKEGG = function(KEGG_database, match = NULL, organism_code = NULL) {
    if (KEGG_database == "compound") { ### Find_Compound ####
      KEGG_value = MS_FindCompound (match = match)

    } else if (KEGG_database == "organism") { ### Find_Organism ####
      KEGG_value = MS_FindOrganism (match = match)

    } else if (KEGG_database == "pathway"){### Find_pathway ####
      KEGG_value = MS_FindPathway(match = match, organism_code = organism_code)

    } else {
      stop ("Possible KEGG_entry values are: compound, organism, or pathway")
    }
  message()
  return(KEGG_value)
  }

################ MS_GetKEGG_GeneID ################

MS_GetKEGG_GeneID = function(genes, organism_name = NULL, output = "vector", orthology = TRUE) {

    message("Transforming input gene IDs into gene nodes:")

    genes = gsub(" ", "", genes)  #remove potential white spaces
    symbols = c()

    ## Check if the genes are symbols or entrezIDs.
    options(warn = -1)
    entrez_answer = !is.na(as.numeric(genes[1]))
    options(warn = 0)

    if (entrez_answer == FALSE) {# The input genes are symbols
        if (length(organism_name) == 0) {
            stop ("An organism_name is required when genes represent symbols")
        }
        symbols = tolower(genes)
        entrezgenes = c()
        symbolsfound = c()
        entrezID = c()
        message(" -Transforming gene symbols into entrez IDs")

        for (symbol in symbols) {
            entrezID = try(as.matrix(unlist(query(q = symbol,
                size = 1, species = organism_name)$hits)), silent = TRUE)
            if (grepl("Error", entrezID)[1] == FALSE) {
                rows = rownames(entrezID)
                entrezID = as.character(entrezID)
                index = which(rows == "entrezgene")
                if (length(index) > 0) {
                  entrezgenes = c(entrezgenes, entrezID[index])
                  symbolsfound = c(symbolsfound, symbol)
                }
            }
        }
        if (length(entrezgenes) == 0) {
            stop("None of the gene symbols was found")
        }
    } else (entrezgenes = genes)

    if (length(entrezgenes) > 0) { # if the function continues, this is TRUE
        KEGG_GeneID_all = c()
        KEGG_GeneID_allM = c()
        message(" -Transforming entrez IDS into gene nodes")
        message()

        for (i in seq_along(entrezgenes)) {
            file = paste("http://rest.kegg.jp/conv/genes/ncbi-geneid:",
                entrezgenes[i], sep = "")
            KEGG_gene = readLines(file)

            if (nchar(KEGG_gene) >= 1) {
                # if the file does not really exist, KEGG_gene=''
                KEGG_GeneID = unlist(strsplit(KEGG_gene, "\t",
                  fixed = FALSE, perl = FALSE, useBytes = FALSE))[2]

                if (orthology == TRUE) {
                  filename = paste("http://rest.kegg.jp/get/", KEGG_GeneID, sep = "")
                  ko_file = try(readLines(filename), silent = TRUE)

                  if (grepl("Error", ko_file[1]) == FALSE) {
                    index_ko = grep("ORTHOLOGY", ko_file)
                    if (length(index_ko) > 0) {
                      ko = substr(ko_file[index_ko], 13, 18)
                      KEGG_GeneID_all = c(KEGG_GeneID_all, ko)
                      line = c(ko, entrezgenes[i])
                      if (length(symbols) > 0) {
                        line = c(line, symbolsfound[i])
                      }
                      KEGG_GeneID_allM = rbind(KEGG_GeneID_allM, line)
                      if (length(symbols) > 0) {
                        colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID",
                                                       "symbol")
                      } else {
                        colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID")
                      }
                      rownames(KEGG_GeneID_allM) = NULL
                    }
                  }
                } else {# We return the actual organism specific genes
                  KEGG_GeneID_all = c(KEGG_GeneID_all, KEGG_GeneID)
                  line = c(KEGG_GeneID, entrezgenes[i])
                  if (length(symbols) > 0) {
                    line = c(line, symbolsfound[i])
                  }
                  KEGG_GeneID_allM = rbind(KEGG_GeneID_allM, line)
                  if (length(symbols) > 0) {
                    colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID",
                                                   "symbol")
                  } else {
                    colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID")
                  }
                  rownames(KEGG_GeneID_allM) = NULL
                }
            }
        }
        if (length(KEGG_GeneID_all) > 0) {
            if (output == "vector") {
                return(KEGG_GeneID_all)
            } else {
                return(KEGG_GeneID_allM)
            }
        } else {
            stop ("None of the genes is reported in KEGG")
        }
    }
}

################ MS_FindMappedNodes ################

MS_FindMappedNodes = function(nodes, network_table, organism_name = NULL, orthology = TRUE){

  if(grepl("cpd:", nodes[1])) {
    message("Mapping metabolites onto the network", "\n")
    response = MS_FindMappedMetabolites(metabolites = nodes, network_table)

  } else {
    message("Mapping genes onto the network", "\n")
    if (length(organism_name) == 0){
      stop("An organism_name is required to map genes")
    }
    response = MS_FindMappedGenes(genes = nodes, organism_name = organism_name,
                                  network_table = network_table, orthology = orthology)
  }
  return(response)
}

################ MS_ReplaceNode ################

MS_ReplaceNode = function(node1, node2, network_table) {
    for (node in node1) {
        network_table[network_table == node] = node2
    }
    network_table = unique(network_table)
    return(network_table)
}

################ MS_RemoveNode ################

MS_RemoveNode = function(nodes, network_table) {
    for (node in nodes) {
        network_table[network_table == node] = NA
    }
    network_table = na.omit(network_table)
    network_table = unique(network_table)
    return(network_table)
}
