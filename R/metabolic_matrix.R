globalVariables("directionality_reactions")

#################### link_reaction_gene ####################
link_reaction_gene = function(enzyme, gene, ko, reactionTable) {
    index = which(reactionTable[, 1] == enzyme)
    if (length(index) > 0) {
        reactions = reactionTable[index, 2]
        line = cbind(reactions, rep(enzyme, length(reactions)),
                     rep(gene, length(reactions)),
                     rep(ko, length(reactions)))
        line = matrix(line, ncol = 4, nrow = length(reactions))
    } else {
      line = NULL
    }
    return(line)
}

#################### new_edge ####################
new_edge = function (reaction, substrates, products, Type) {
    substrate_lines = cbind(substrates, rep(reaction, length(substrates)))
    product_lines = cbind(rep(reaction, length(products)), products)
    edges_reaction = rbind(substrate_lines, product_lines)

    if (Type == "reversible") {
        reverse = cbind(edges_reaction[, 2], edges_reaction[, 1])
        edges_reaction = rbind(edges_reaction, reverse)
    }
    return(edges_reaction)
}

#################### add_reaction_edge ####################
add_reaction_edge = function(Reaction, Substrate, Product, Type) {

    reactions = unlist(strsplit(Reaction, " "))
    substrates = unlist(strsplit(Substrate, ";"))
    substrates = unlist(strsplit(substrates, " "))
    products = unlist(strsplit(Product, ";"))
    products = unlist(strsplit(products, " "))

    # Correct for direction based on Duarte et al., 2007
    for (reaction in reactions) {
        a = which(directionality_reactions[, 1] == reaction)
        if (length(a) >= 1) {
            Type = directionality_reactions[, 2][a]
        }
    }

    # Build edges
    all_edges = lapply(reactions, new_edge, substrates, products, Type)
    all_edges = unique(do.call(rbind, all_edges))
    return(all_edges)
}

#################### network_from_path ####################
network_from_path = function (path, list_parsed_paths) {
    lines = list_parsed_paths[[path]]

    ## Select reactions
    global_lines_reactions = intersect(grep("rn", lines), grep("Name", lines))
    Reactions = lines[c(global_lines_reactions)]
    Reactions = substr(Reactions, 11, nchar(Reactions))

    ## Select substrates
    global_lines_substrates = grep("Substrate Name", lines)
    Substrates = lines[c(global_lines_substrates)]
    Substrates = substr(Substrates, 21, nchar(Substrates))

    # Select products
    global_lines_products = grep("Product Name", lines)
    Products = lines[c(global_lines_products)]
    Products = substr(Products, 19, nchar(Products))

    # Select type
    global_lines_type = grep("Type", lines)
    Type = lines[c(global_lines_type)]
    Type = substr(Type, 11, nchar(Type))

    ## Build network edges
    all_edges = mapply(add_reaction_edge, Reactions, Substrates, Products,
                       Type, SIMPLIFY = FALSE)
    all_edges = unique(do.call(rbind, all_edges))
    return(all_edges)
}

##################### bind_reaction_gene ####################################
bind_reaction_gene = function (reaction, matrix) {
    index_reaction = which(matrix[, 1] == reaction)
    if (length(index_reaction) > 0) {
      reaction_genes = paste(reaction, matrix[index_reaction, 2], sep = "_")
    } else {
      reaction_genes = reaction
    }
    return(reaction_genes)
}

##################### build_RG_edge ####################################
build_RG_edge = function(edge, reactionko, reactiongene, expand_genes) {
    index_rn = grep("rn", edge)
    indexM = numeric(length = 1)
    if (index_rn == 1) {
        indexM = 2
    } else {
        indexM = 1
    }

    reaction = edge[index_rn]
    metabolite = edge[indexM]

    if (expand_genes == FALSE) { # Genes represent orthology IDs.
        reaction_genes_all = bind_reaction_gene(reaction, matrix = reactionko)

    } else { # Genes represent organism specific IDs.
        reaction_genes_all = bind_reaction_gene(reaction, matrix = reactiongene)
    }

    if (index_rn == 1) {
        new_edges = cbind(reaction_genes_all, rep(metabolite,
            length(reaction_genes_all)))
    } else {
        new_edges = cbind(rep(metabolite, length(reaction_genes_all)), reaction_genes_all)
    }
    return(new_edges)
}

#################### metabolic_matrix ####################
metabolic_matrix = function(path_names, list_parsed_paths, organism_code,
                            expand_genes) {
    metabolic_table_RG = NULL # initiate the variable
    metabolic_table = lapply(path_names, network_from_path, list_parsed_paths )
    metabolic_table = unique(do.call(rbind, metabolic_table))
    metabolic_table = matrix(metabolic_table, ncol = 2)

    ## Link reactions to genes (organism specific or orthology IDs)

    file_enzyme = paste("https://rest.kegg.jp/link/enzyme/", organism_code, sep = "")
    response_enzyme = getURL(file_enzyme)
    enzymeTable = convertTable(response_enzyme)

    file_reaction = "https://rest.kegg.jp/link/reaction/enzyme"
    response_reaction = getURL(file_reaction)
    reactionTable = convertTable(response_reaction)

    file_ko = paste("https://rest.kegg.jp/link/ko/", organism_code, sep = "")
    response_ko = getURL(file_ko)
    koTable = convertTable(response_ko)
    koTable[, 2] = substr(koTable[, 2], 4, 9)

    if (is.matrix(enzymeTable) & is.matrix(reactionTable) & is.matrix(koTable)) {
        to_print = ("Linking reactions to genes")
        message(to_print)

        genes = enzymeTable[, 1]
        enzymes = enzymeTable[, 2]
        ko = vector(length = length(genes))

        for (i in seq_along(genes)) {
            index = which(koTable[, 1] == genes[i])
            if (length(index) > 0) {
                ko[i] = koTable[index[1], 2]
            } else {
                ko[i] = genes[i]
            }
        }
        ## Link reactions to genes
        RT = list()
        RT[["reactionTable"]] = reactionTable
        reactions_genes = mapply(link_reaction_gene, enzymes, genes, ko,
                                 MoreArgs = RT, SIMPLIFY = FALSE)
        reactionM = unique(do.call(rbind, reactions_genes))

        reactionko = unique(na.omit(reactionM[, c(1, 4)]))
        reactiongene = unique(na.omit(reactionM[, c(1, 3)]))

        all_edges_RG = lapply(split(metabolic_table, row(metabolic_table)),
            build_RG_edge, reactionko = reactionko, reactiongene = reactiongene,
            expand_genes = expand_genes)

        metabolic_table_RG = unique(do.call(rbind, all_edges_RG))
        metabolic_table_RG = matrix(metabolic_table_RG, ncol = 2)
        colnames(metabolic_table_RG) = c("node1", "node2")
        rownames(metabolic_table_RG) = NULL
    }
    return(metabolic_table_RG)
}

################## get_metabonetR ####################
get_metabonetR <- function(path) {
    message(path)
    if (substr(path, 4, nchar(path)) == "01100") { # remove metabolic pathways map
        parsed_path <- NULL
    } else {
        # Check that the input path exists
        file <- paste("https://rest.kegg.jp/get/", path, "/kgml", sep = "")
        pathway <- try(getURL(file), silent = TRUE)
        reactions <- try(getReactions(parseKGML(pathway)), silent = TRUE)

        if (grepl("Error", reactions[1]) == TRUE) {
            to_print <- paste(path, "-path ID without XML:path removed", sep = "")
            message(to_print)
            parsed_path <- NULL
        } else {
            parsed_path <- capture.output(reactions, file = NULL)
        }
    }
    return(parsed_path)
}


