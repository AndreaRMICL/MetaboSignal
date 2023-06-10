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

#################### get_tissue_level ###################
get_tissue_level = function(ensembl_gene, my_hpar_data, alternative_levels) {

    ind = which(my_hpar_data[, "Gene"] == ensembl_gene)
    levels_gene = unique(as.character(my_hpar_data[ind, "Level"]))

    if (length (intersect(alternative_levels, levels_gene)) == 0) {
      return("undetected")
    }
    return("detected")
}

#################### MS_tissueFilter  ###################
MS_tissueFilter = function(network_table, tissue, input_format = "kegg",
                           expand_genes = FALSE) {

    ## Make network unique
    network = unique(network_table)

    ## Check input network
    check_matrix_v2(network, n = 2)

    ## Input format should be kegg or entrez
    if ((input_format %in% c("kegg", "entrez")) == FALSE) {
        stop("Input format must be kegg or entrez")
    }

    if (input_format == "kegg") {
        if(length(grep("hsa", as.vector(network_table))) == 0) {
            stop ("Gene nodes must represent human specific kegg ids")
        }
    }

    ## The filtering will only be apply to genes not linked to compounds
    cpd_ind = c(grep("cpd:|dr:|gl:", network[, 1]),
                grep("cpd:|dr:|gl:", network[, 2]))

    if(length(cpd_ind) > 0) {
        if (length(cpd_ind) == nrow(network)) {
            stop("network_table seems not to have signaling genes")
        }
        cpd_edges = network[cpd_ind, ]
        gene_edges = network[ -cpd_ind, ]
        gene_ids = unique(as.vector(gene_edges[, c(1, 2)]))
    } else {
        cpd_edges = NULL
        gene_edges = network
        gene_ids = unique(as.vector(gene_edges[, c(1, 2)]))
    }

    ## If genes are in kegg ids get entrez ids
    if (input_format == "kegg") {
        genes_list = split(gene_ids, ceiling(seq_along(gene_ids)/100))
        entrez_res = unlist(lapply(genes_list, conv_entrez_kegg, source = "kegg"))
        entrez_ids = as.character(substr(entrez_res, 13, nchar(entrez_res)))
        start = unlist(gregexpr(pattern = "hsa", names(entrez_res)))
        kegg_ids = substr(names(entrez_res), start, nchar(names(entrez_res)))
    } else {
        kegg_ids = gene_ids
        entrez_ids = gene_ids
    }

    ## Tranform entrez IDs into ENSEMBL IDs
    X = org.Hs.egENSEMBL
    mapped_genes = mappedkeys(X)
    xx = as.list(X[mapped_genes])
    NL = list()
    NL[["xx"]] = xx
    hsa_gene_ensembl = mapply(link_hsagene_ensembl,
                              entrez_ids, kegg_ids, MoreArgs = NL, SIMPLIFY = FALSE)
    gene_enSall = unique(do.call(rbind, hsa_gene_ensembl))

    ensembl = as.character(gene_enSall[, 3])

    ## Get HPA data, and get the data with good reliability for the tissues of interest
    hpar_data = hpaNormalTissue
    tissue = paste(tissue, collapse = "|")
    tissue_index = grep(tissue, hpar_data[, "Tissue"])
    level_index = grep("Not detected", hpar_data[, "Level"])
    reliability_index =  grep("Supportive|Approved|Supported",
                              hpar_data[, "Reliability"])
    common_index = intersect(tissue_index, reliability_index)

    if(length(common_index) == 0) {
        stop("Tissue-filtering failed: tissue seems to be invalid")
    }

    hpar_data_filtered = hpar_data[common_index, ]

    ## Map ensembl ids of interest in HPA data
    ensembl_in_hpar = intersect(ensembl, hpar_data_filtered[, "Gene"])

    if(length(ensembl_in_hpar) == 0) { ## none of the genes was mapped
        stop("Tissue-filtering failed: none of the genes seems to be in HPA")
    }

    find_index_hpar = function(name, all_names) {
        return(which(all_names == name))
    }

    ## Get HPA subset containing only the mapped ensembl ids
    my_hpar_ind = unique(unlist(lapply(ensembl_in_hpar, find_index_hpar,
                                       hpar_data_filtered[, "Gene"])))
    my_hpar_data = hpar_data_filtered[my_hpar_ind, ]

    ## Get unwanted genes: genes that are non detected in all input tissue(s)
    alternative_levels = setdiff(levels(hpar_data[, "Level"]), "Not detected")
    tissue_ans = sapply(ensembl_in_hpar, get_tissue_level,
                        my_hpar_data, alternative_levels)

    if (("undetected" %in% tissue_ans) == FALSE) {
      return("All genes satisfy filtering parameters")
    }
    unwanted_ensembl = ensembl_in_hpar[tissue_ans == "undetected"]
    rownames(gene_enSall) = gene_enSall[, 3]
    genes_not_in_tissue = as.character(gene_enSall[unwanted_ensembl, 1])

    edges_not_in_tissueL = sapply(lapply(split(gene_edges, row(gene_edges)),
                                         intersect, y = genes_not_in_tissue), length)
    edges_not_in_tissueL = as.vector (edges_not_in_tissueL)
    edges_not_in_tissue = unique(which(edges_not_in_tissueL >= 1))

    if (nrow(gene_edges) > length(edges_not_in_tissue)) {
      gene_edges_final = gene_edges[-c(edges_not_in_tissue), ]
      network_final = rbind(gene_edges_final, cpd_edges)
      colnames(network_final) = colnames(network)
      rownames(network_final) = NULL
    } else {
        return("None of the genes satisfies the filtering parameters")
    }

    ## Transform kegg IDs into orthology IDs if required
    if (expand_genes == FALSE) {
        message()
        to_print = ("Transforming gene IDs into orthology IDs")
        message(to_print)

        file_ko = paste("https://rest.kegg.jp/link/ko/", "hsa", sep = "")
        response_ko = getURL(file_ko)
        koTable = convertTable(response_ko)
        koTable[, 2] = substr(koTable[, 2], 4, 9)

        for (r in 1:nrow(network_final)) {
            for (c in 1:ncol(network_final[, 1:2])) {
                index = which(koTable[, 1] == network_final[r, c])
                if (length(index) > 0) {
                    network_final[r, c] = koTable[index, 2]
                }
            }
        }
    }

    ## Report network features#
    network_final = unique(network_final)
    message()
    network_features(network_final[, 1:2])

    return(network_final)
}
