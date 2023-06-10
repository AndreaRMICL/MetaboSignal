################ MS_keggFinder ####################
MS_keggFinder = function(KEGG_database, match = NULL, organism_code = NULL) {
    if (KEGG_database == "compound") {
        ### Find_Compound ####
        KEGG_value = MS_FindCompound(match = match)

    } else if (KEGG_database == "organism") {
        ### Find_Organism ####
        KEGG_value = MS_FindOrganism(match = match)

    } else if (KEGG_database == "pathway") {
        ### Find_pathway ####
        KEGG_value = MS_FindPathway(match = match, organism_code = organism_code)

    } else {
        stop("Possible KEGG_entry values are: compound, organism, or pathway")
    }
    message()
    return(KEGG_value)
}

################ MS_convertGene ################
MS_convertGene = function(genes, organism_code, organism_name = NULL,
                          output = "vector", orthology = TRUE) {

    message("Transforming input gene IDs into gene nodes:")

    genes = gsub(" ", "", genes)  #remove potential white spaces
    genes = unique(genes)

    ## Check if the genes are symbols or entrezIDs.
    entrez_answer = !is.na(suppressWarnings((as.numeric(genes[1]))))

    if (entrez_answer == FALSE) {
        # The input genes are symbols
        if (length(organism_name) == 0) {
            stop("An organism_name is required when genes represent symbols")
        }
        message(" -Transforming gene symbols into entrez IDs")
        symbols = tolower(genes)

        entrezgenes = sapply(symbols, conv_entrez_symbol, organism_name, source = "symbol")
        entrezgenes = as.character(entrezgenes)
        index_NF = grep("NF", entrezgenes)

        if (length(index_NF) > 0) {
            entrezgenes = entrezgenes[-index_NF]
            symbolsfound = symbols[-index_NF]
        } else {
            symbolsfound = symbols
        }

        if (length(entrezgenes) == 0) {
            stop("None of the gene symbols was found")
        }
    } else {
        entrezgenes = genes
    }

    message(" -Transforming entrez IDS into gene nodes")
    message()

    ## Split genes in list of max 100 genes (for multiple query in KEGG API)
    genes_list = split(entrezgenes, ceiling(seq_along(entrezgenes)/100))

    ## Get KEGG IDs
    KEGG_res = unlist(lapply(genes_list, conv_entrez_kegg, source = "entrez",
                             organism_code = organism_code))
    kegg_ids = unique(as.character(KEGG_res))

    if (length(kegg_ids) == 0) {
        stop("None of the genes is reported in KEGG")
    }
    start = unlist(gregexpr(pattern = "ncbi", names(KEGG_res))) + 12
    entrez_mapped = substr(names(KEGG_res), start, nchar(names(KEGG_res)))
    entrez_mapped = unique(entrez_mapped)

    ## Get index symbols
    if (entrez_answer == FALSE) {
        index_mapped = numeric(length = length(entrez_mapped))
        for (i in seq_along(entrez_mapped)) {
            index_mapped[i] = which(entrezgenes == entrez_mapped[i])[1] #
            #in case there are several symbols associated to the same entrez
        }
        symbols_mapped = symbolsfound[index_mapped]
    }

    ## Transform the KEGG IDs into orthology IDs
    if (orthology == TRUE) {
        file_ko = paste("https://rest.kegg.jp/link/ko/", organism_code, sep = "")
        response_ko = getURL(file_ko)
        if(nchar(response_ko) == 0) {
          stop ("Could not get orthology IDs. Check that organism_code is correct")
        }
        koTable = convertTable(response_ko)
        koTable[, 2] = substr(koTable[, 2], 4, 9)

        ko_genes = as.character(sapply(kegg_ids, link_kogene, koTable))

        if (entrez_answer == TRUE) {
            KEGG_GeneID_allM = cbind(ko_genes, entrez_mapped)
            colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID")
        } else {
            KEGG_GeneID_allM = cbind(ko_genes, entrez_mapped, symbols_mapped)
            colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID", "symbol")
        }

        ## Remove NF
        index_NF = grep("NF", ko_genes)

        if (length(index_NF) > 0) {
            if (length(index_NF) == length(ko_genes)) {
                stop("None of the genes is reported in KEGG")
            } else {
                back_up = KEGG_GeneID_allM
                KEGG_GeneID_allM = KEGG_GeneID_allM[-index_NF,]
                KEGG_GeneID_allM = matrix(KEGG_GeneID_allM, ncol = ncol(back_up))
                colnames(KEGG_GeneID_allM) = colnames(back_up)
            }
        }
    } else { # organism specific genes
        if (entrez_answer == TRUE) {
            KEGG_GeneID_allM = cbind(kegg_ids, entrez_mapped)
            colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID")
        } else {
            KEGG_GeneID_allM = cbind(kegg_ids, entrez_mapped, symbols_mapped)
            colnames(KEGG_GeneID_allM) = c("KEGG_ID", "entrez_ID", "symbol")
        }
    }

    if (output == "vector") {
        return(as.character(KEGG_GeneID_allM[, 1]))
    } else {
        return(KEGG_GeneID_allM)
    }
}

################ MS_findMappedNodes ################
MS_findMappedNodes = function(nodes, network_table) {
    #check_matrix_v2(network_table)
    all_nodes = unique(as.vector(network_table))
    mapped_nodes = intersect(nodes, all_nodes)
    unmapped_nodes = setdiff(nodes, all_nodes)

    res_list = list(mapped_nodes = mapped_nodes, unmapped_nodes = unmapped_nodes)

    return(res_list)
}

################ MS_replaceNode ################
MS_replaceNode = function(node1, node2, network_table) {
    for (node in node1) {
        network_table[network_table == node] = node2
    }
    network_table = unique(network_table)
    return(network_table)
}

################ MS_removeNode ################
MS_removeNode = function(nodes, network_table) {
    for (node in nodes) {
        network_table[network_table == node] = NA
    }
    network_table = na.omit(network_table)
    network_table = unique(network_table)
    return(network_table)
}

################ MS_removeDrugs ################
MS_removeDrugs = function(network_table) {
    check_matrix_v2(network_table, n = 2)
    ind = unique(c(grep("dr:", network_table[, 1]), grep("dr:", network_table[, 2])))

    if (length(ind) > 0) {
          network_table = network_table[-ind, ]
          network_features(network_table[, 1:2])
          return(network_table)
    }
    message("Drug nodes not found")
    return(network_table)
}

################ MS_getPathIds ################
MS_getPathIds = function(organism_code) {
    path_ids = keggLink("pathway", organism_code)
    path_ids = unique(gsub("path:", "", path_ids))
    metabo_pathways = vector(mode = "list", length = length(path_ids))
    signaling_pathways = vector(mode = "list", length = length(path_ids))

    for(path in path_ids) {
        #print(path)
        res = keggGet(path)
        class = res[[1]]$CLASS

        if (!is.null(class)) {
            if (grepl("Metabolism", class) & class != "Metabolism; Overview") {
                metabo_pathways[[path]] = c(path, res[[1]]$NAME, class, "metabolic")
            }
            if (grepl("Metabolism", class) == FALSE &
                grepl("Genetic Information Processing;", class) == FALSE) {
                signaling_pathways[[path]] = c(path, res[[1]]$NAME, class, "signaling")
            }
        }
    }

    metabo_pathways = do.call(rbind, metabo_pathways)
    signaling_pathways = do.call(rbind, signaling_pathways)

    all_pathways = rbind(metabo_pathways, signaling_pathways)

    colnames(all_pathways) = c("Path_id", "Path_description", "Path_category", "Path_type")
    rownames(all_pathways) = NULL

    file_name = paste(organism_code, "pathways.txt", sep = "_")

    all_pathwaysDF = as.data.frame(all_pathways)

    write.table(all_pathwaysDF, file_name, row.names = FALSE,
                sep = "\t", quote = FALSE, col.names = TRUE)

    return(all_pathways)
}

