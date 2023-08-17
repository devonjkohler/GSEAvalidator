#' Prepares MS label-free data for GSEA correlation analysis
#'
#' Takes the ProteinLevel output of the MSstats dataProcess function as input.
#'
#' @param data The ProteinLevelData object returned from the MSstats dataProcess
#' function.
#' @param gene_map table which maps uniprot protein names to gene names, if
#' required.
#' @param parse_genes if gene name has a "_" (e.g. SPTBN2_HUMAN), setting this
#' to TRUE will remove the string after the underscore.
#'
#' @import tidyverse
#' @import stringr
#' @import dplyr
#' @import tidyr
#'
#' @export
prep_msstats_data = function(data, gene_map=NULL, parse_gene=FALSE){

  if (parse_gene){
    data$Protein = unlist(
      lapply(data$Protein, function(x){str_split(x, "_")[[1]][[1]]})
      )
  }
  if (!is.null(gene_map)){
    data = merge(data, gene_map, all.x=TRUE, all.y=FALSE,
                 by.x="Protein", by.y="From")
    data$Protein = data$To
  }


  parse_data = data %>% select(Protein, originalRUN, LogIntensities)
  wide_data = parse_data %>% pivot_wider(names_from=Protein,
                                         values_from=LogIntensities)
  wide_data = wide_data[2:length(colnames(wide_data))]
  return(wide_data)
}

#' Generates list of different correlation matrices
#'
#' Takes the output of the prep_msstats_data function as input
#'
#' @param data output of the prep_msstats_data function.
#' @param methods list of correlation methods to use (in cor function)
#'
#' @export
gen_correlation_matrix = function(
    data,
    methods = c("pearson", "spearman")){

  correlations = list()
  for (m in methods){
    cor_mat = suppressWarnings(cor(data, method=m, use="pairwise.complete.obs"))
    correlations[m] = list(abs(cor_mat))
  }

  return(correlations)
}

#' Checks correlations between genes in GSEA pathways
#'
#' takes as input correlation list and gsea file path
#'
#' @param correlation_data output of gen_correlation_matrix function
#' @param gsea_path file path to GSEA JSON file
#' @param threshold correlation threshold to count as significant correlation
#'
#' @import rjson
#' @import combinat
#' @import data.table
#'
#' @export
test_GSEA_pathways = function(correlation_data, gsea_path, threshold=.5){

  gsea = fromJSON(file=gsea_path)
  pathways = names(gsea)

  measured_genes = colnames(correlation_data[[1]])

  result_list = list()

  for (path in seq_along(pathways)){
    current_path = gsea[[path]]
    genes_in_path = current_path$geneSymbols

    measured_genes_in_path = intersect(genes_in_path, measured_genes)

    if (length(measured_genes_in_path) > 2){
      combination_mat = combn(measured_genes_in_path, 2)
      total_combs = dim(combination_mat)[[2]]

      corr_types = names(correlation_data)

      for (i in seq_along(corr_types)){
        sig_tracker = 0
        for (node in seq_len(total_combs)){

          node_idx1 = match(combination_mat[1,node],
                            colnames(correlation_data[[i]]))
          node_idx2 = match(combination_mat[2,node],
                            colnames(correlation_data[[i]]))

          check = correlation_data[[i]][node_idx1, node_idx2]
          if (is.na(check)){
            check = 0
          }
          if (check > threshold){
            sig_tracker = sig_tracker + 1
          }
        }
        temp_data = data.table("pathway"=c(pathways[[path]]),
                               "correlation"=c(corr_types[[i]]),
                               "total_genes"=c(length(genes_in_path)),
                               "measured_genes"=c(length(measured_genes_in_path)),
                               "percent_measured"=c(
                                 length(measured_genes_in_path)/length(genes_in_path)),
                               "total_tests"=c(total_combs),
                               "sig_corrs"=c(sig_tracker),
                               "percent"=c(sig_tracker/total_combs))

        result_list <- append(result_list, list(temp_data))
      }
    } else {
      temp_data = data.table("pathway"=c(pathways[[path]]),
                             "correlation"=c("no measured nodes"),
                             "total_genes"=c(length(genes_in_path)),
                             "measured_genes"=c(0),
                             "percent_measured"=c(0),
                             "total_tests"=c(0),
                             "sig_corrs"=c(0),
                             "percent"=c(0))
    result_list <- append(result_list, list(temp_data))
    }
  }

  results = rbindlist(result_list)
  return(results)
}

