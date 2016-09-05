################ MS_FindKEGG ####################

MS_FindKEGG = function(KEGG_database, match = NULL, organism_code = NULL) {
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

################ MS_GetKEGG_GeneID ################

MS_GetKEGG_GeneID = function(genes, organism_code, organism_name = NULL,
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
        file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code, sep = "")
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

################ MS_FindMappedNodes ################

MS_FindMappedNodes = function(nodes, network_table, organism_code = NULL,
                              organism_name = NULL, orthology = TRUE) {

    nodes = unique(nodes)

    if (grepl("cpd:", nodes[1])) {
        message("Mapping metabolites onto the network", "\n")
        response = MS_FindMappedMetabolites(metabolites = nodes, network_table)

    } else {
        message("Mapping genes onto the network", "\n")
        if (length(organism_name) == 0 | length(organism_code) == 0) {
            stop("An organism_name and an organism_code are required to map genes")
        }
        response = MS_FindMappedGenes(genes = nodes, organism_code = organism_code,
                                      organism_name = organism_name,
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

