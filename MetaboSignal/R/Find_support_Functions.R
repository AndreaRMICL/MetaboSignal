#################### INTERNAL FUNCTIONS ####################

##### MS_FindCompound ####
MS_FindCompound = function(match = NULL) {
  file = "http://rest.kegg.jp/list/compound"
  response = getURL(file)
  compoundM = convertTable(response)
  colnames(compoundM) = c("KEGG compound", "common names")
  rownames(compoundM) = NULL
  if (length(match) >= 1) {
    target_matrix = compoundM
    target_column = compoundM[,2]
    matchM = match_KEGG (match, target_column , target_matrix)
    return(matchM)
  } else (return(compoundM))
}

##### MS_FindOrganism ####
MS_FindOrganism = function(match = NULL) {
  file = "http://rest.kegg.jp/list/organism"
  response = getURL(file)
  organismM = convertTable(response)
  colnames(organismM) = c("T","organism_code","organism_name","description")
  rownames(organismM) = NULL
  if (length(match) >= 1) {
    target_matrix = organismM
    target_column = organismM[,3]
    matchM = match_KEGG (match, target_column , target_matrix)
    return(matchM)
  } else (return(organismM))
}

##### MS_FindPathway ####
MS_FindPathway = function(match = NULL, organism_code = NULL) {
  file = paste("http://rest.kegg.jp/list/pathway/", organism_code, sep = "")
  response = try(getURL(file), silent = TRUE)
  if (nchar(response) == 0){
    stop ("A valid organism_code is required for KEGG_entry = pathway")
  }
  pathM = convertTable(response)
  colnames(pathM) = c("path_ID", "path_Description")
  rownames(pathM) = NULL
  if (length(match) >= 1) {
    target_matrix = pathM
    target_column = pathM[,2]
    matchM = match_KEGG (match, target_column , target_matrix)
    return(matchM)
  } else (return(pathM))
}

##### MS_FindMappedGenes ####
MS_FindMappedGenes = function(genes, organism_name, network_table, orthology = TRUE) {
  List_mappedID = list()
  length(List_mappedID) = 4
  names(List_mappedID) = c("Input gene IDs mapped onto the network",
                           "All input gene IDs not mapped onto the network",
                           "Input gene IDs not found in KEGG", "KEGG gene IDs not mapped onto the network")
  genes = gsub(" ", "", genes)  #remove potential white spaces
  genes = unique(genes)
  genes = tolower(genes)
  transformed_genesM = c()
  transformed_genesM = try(MS_GetKEGG_GeneID(genes, organism_name, output = "matrix",
                                             orthology = orthology), silent = TRUE)

  if (grepl("Error",transformed_genesM)[1]==TRUE) {
    List_mappedID[[2]] = genes
    List_mappedID[[3]] = genes
    return(List_mappedID)
  }
  if (is.matrix(transformed_genesM) == TRUE) {
    if (ncol(transformed_genesM) == 3) {# Input genes are symbols.
      transformed_IDs = as.character(transformed_genesM[,3])
    } else (transformed_IDs = as.character(transformed_genesM[, 2]))
    # Input genes are entrez IDs.
  } else {
    List_mappedID[[2]] = genes
    List_mappedID[[3]] = genes
    return(List_mappedID)
  }

  ## Find input genes that are not mapped in KEGG#
  IDs_not_inKEGG = c()
  IDs_not_inKEGG = setdiff(genes, transformed_IDs)

  ## Find ko or organism specific genes (transformed genes) that
  ## are mapped in KEGG but not in the network
  transformed_genes = as.character(transformed_genesM[, 1])
  all_nodes = unique(as.vector(network_table))
  genes_not_mapped = setdiff(transformed_genes, all_nodes)

  # Find gene input IDs (symbols or entrez) that of the
  # transformed genes that were not mapped onto the network.
  IDs_not_mapped = c()
  for (gene in genes_not_mapped) {
    index = which(transformed_genes == gene)
    IDs_not_mapped = c(IDs_not_mapped, (transformed_IDs[index]))
  }
  IDs_not_found_all = c()
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

##### MS_FindMappedMetabolites ####
MS_FindMappedMetabolites = function(metabolites, network_table) {
  List_metabo = list()
  length(List_metabo) = 2
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

#################### EXTERNAL FUNCTIONS ####################

##### MS_FindKEGG ####

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

##### MS_GetKEGG_GeneID ####

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

        for (i in 1:length(entrezgenes)) {
            file = paste("http://rest.kegg.jp/conv/genes/ncbi-geneid:",
                entrezgenes[i], sep = "")
            KEGG_gene = readLines(file)

            if (nchar(KEGG_gene) >= 1) {
                # if the file does not really exist, KEGG_gene=''
                KEGG_GeneID = unlist(strsplit(KEGG_gene, "\t",
                  fixed = FALSE, perl = FALSE, useBytes = FALSE))[2]

                if (orthology == TRUE) {
                  filename = paste("http://rest.kegg.jp/get/", KEGG_GeneID,
                                   sep = "")
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

##### MS_FindMappedNodes ####

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

##### MS_ReplaceNode ####

MS_ReplaceNode = function(node1, node2, network_table) {
    for (node in node1) {
        network_table[network_table == node] = node2
    }
    network_table = unique(network_table)
    return(network_table)
}

##### MS_RemoveNode ####

MS_RemoveNode = function(nodes, network_table) {
    for (node in nodes) {
        network_table[network_table == node] = NA
    }
    network_table = na.omit(network_table)
    network_table = unique(network_table)
    return(network_table)
}
