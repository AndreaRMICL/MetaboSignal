#################### MS_reactionNetwork ###################
MS_reactionNetwork <- function(metabo_paths) {

    ### Get KGML files and transform them into reaction files####
    list_parsed_paths <- lapply(metabo_paths, get_metabonetR)
    names(list_parsed_paths) <- metabo_paths
    path_names <- metabo_paths

    if (length(list_parsed_paths) == 0) {
        to_print <- ("Impossible to build a metabolic network")
        warning(to_print, "\n")

    } else {
        ## Remove empty paths: for example oxidative phosphorylation#
        length_paths <- sapply(list_parsed_paths, length)
        empty_paths <- which(length_paths <= 1)

        if (length(empty_paths) >= 1) {
            list_parsed_paths = list_parsed_paths[-c(empty_paths)]
            path_names = path_names[-c(empty_paths)]
        }
        if (length(list_parsed_paths) == 0) {
            to_print <- ("Impossible to build a metabolic network")
            warning(to_print, "\n")
        } else {

            ### Create metabolic_table ####
            metabolic_table <- lapply(path_names, network_from_path, list_parsed_paths )
            metabolic_table <- unique(do.call(rbind, metabolic_table))
            metabolic_table <- matrix(metabolic_table, ncol = 2)

        }
    }
    ## Add directionality
    if(all(metabo_paths == "rn01100")) {
      metabolic_table_final <- cbind(metabolic_table, type = "k_compound:reversible")
    } else {
      rownames(metabolic_table) <- paste(metabolic_table[, 1], metabolic_table[, 2], sep = "_")
      check_rev <- c(rownames(metabolic_table), paste(metabolic_table[, 2], metabolic_table[, 1], sep = "_"))
      id_dup <- which(duplicated(check_rev))

      if (length(id_dup) > 0) { ## there are reversible reactions
        dup_names <- paste(check_rev[id_dup], collapse = "|")
        irrev_rn <- cbind(metabolic_table[-grep(dup_names, rownames(metabolic_table)), ],
                          type = "k_compound:irreversible")
        sort_dup <- lapply(check_rev[id_dup], function(x) sort(unlist(strsplit(x, "_"))))
        rev_rn <- cbind(unique(do.call(rbind, sort_dup)), type = "k_compound:reversible")

        metabolic_table_final <- rbind(irrev_rn, rev_rn)
      } else {
        metabolic_table_final <- cbind(metabolic_table, type = "k_compound:irreversible")
      }
    }
    colnames(metabolic_table_final) <- c("source", "target", "type")
    rownames(metabolic_table_final) <- NULL

    return(metabolic_table_final)
}
