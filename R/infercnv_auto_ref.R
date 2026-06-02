# Helper function
.run_infercnv <- function(
    expr_mat,
    annot_df,
    ref_names,
    out_dir,
    gene_pos_file,
    cutoff,
    ncores,
    no_plot = FALSE
) {
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = expr_mat,
    gene_order_file   = gene_pos_file,
    annotations_file  = annot_df,
    delim             = "\t",
    ref_group_names   = ref_names,
    chr_exclude       = "chrM"
  )
  infercnv::run(
    obj,
    cutoff=cutoff,
    out_dir=out_dir,
    denoise=TRUE,
    HMM=FALSE,
    cluster_by_groups = FALSE,
    write_phylo = TRUE,
    analysis_mode = "cells",
    num_threads = ncores,
    no_plot = no_plot
  )
}

#' Run inferCNV with automatic same-population reference selection
#'
#' Runs an inferCNV workflow on a selected query population while automatically
#' selecting a cleaner reference subset from a user-defined reference annotation.
#' If fewer than 100 reference spots are available after filtering, all reference
#' spots are used directly. If at least 100 reference spots are available, the
#' function first runs reference-free inferCNV on the reference spots, performs
#' k-means clustering across increasing k values, and selects the lowest-CNV
#' cluster at a local CNV-variance minimum as the final reference.
#'
#' The reference annotation must be included in `query_anno`, because the method
#' assumes that reference candidates are selected from within the same queried
#' cell type or tissue compartment.
#'
#' @param obj A Seurat object, typically prepared with semla.
#' @param output_base Character. Base output directory where inferCNV results,
#' reference-selection files, and selected reference barcodes are written.
#' @param gene_pos_file Character. Path to the gene position file used by
#' `infercnv::CreateInfercnvObject()`.
#' @param hist_column Character. Metadata column containing the annotations used
#' for query and reference selection. Defaults to `"Histology"`.
#' @param query_anno Character vector. Annotation labels in `hist_column` to keep
#' for the final inferCNV query.
#' @param ref_anno Character vector. Annotation label(s) in `hist_column` used as
#' candidate reference spots. Must be included in `query_anno`.
#' @param min_feats Numeric. Minimum number of detected features required for a
#' spot to be retained. Default: `500`.
#' @param min_counts Numeric. Minimum number of counts required for a spot to be
#' retained. Default: `1000`.
#' @param cutoff Numeric. Cutoff passed to `infercnv::run()`. Default: `0.1`.
#' @param ncores Integer. Number of threads passed to inferCNV. Default: `10`.
#' @param prop_min Numeric. Minimum selected reference cluster size as a fraction
#' of all candidate reference spots. Default: `0.10` (suitable for 10x).
#' @param k_window Integer. Number of k values before and after a candidate k
#' used when detecting a local CNV-variance minimum. Default: `2`.
#' @param plot Logical. Whether inferCNV should generate its default plots.
#' Default: `TRUE`.
#'
#' @return
#' A named list is returned to R and can be saved as an object,
#' for example `res <- run_infercnv_auto_ref(...)`.
#' \describe{
#'   \item{`ref_out_dir`}{
#'   Path to the output directory containing the final inferCNV results.
#'   }
#'
#'   \item{`ref_cells`}{
#'   Character vector containing the spot barcodes selected as the final
#'   reference population.
#'   }
#'
#'   \item{`k_best`}{
#'   Optimal k selected during automatic reference selection. Returns `NA`
#'   when fewer than 100 candidate reference spots are available and all
#'   reference spots are used directly.
#'   }
#'
#'   \item{`n_reference_candidates`}{
#'   Number of candidate reference spots available after filtering.
#'   }
#' }
#'
#' @details
#' Candidate reference spots are first selected using `ref_anno` within the
#' filtered query population. If fewer than 100 candidate reference spots are
#' available, all reference spots are used directly. Otherwise, inferCNV is
#' first run without a reference on the candidate reference population.
#'
#' The resulting inferCNV expression matrix is used to calculate a per-spot
#' CNV score:
#'
#' \deqn{CNV\ score = mean((expr - 1)^2)}
#'
#' k-means clustering is then performed across increasing values of k. For each
#' k, the cluster with the lowest average CNV score is identified. The optimal
#' k is selected as a local minimum in CNV score, requiring the selected cluster
#' to contain at least `prop_min` of all candidate reference spots.
#'
#' The final selected reference spots are relabeled as `"pure_ref"` and used as
#' the inferCNV reference population for the full query dataset.
#'
#' Output directories:
#'
#' \describe{
#'   \item{`ref_no-ref`}{
#'   Reference-only inferCNV run used for automatic reference selection.
#'   }
#'
#'   \item{`ref_no-ref/k_selection`}{
#'   Intermediate k-means clustering results and reference-selection diagnostics.
#'   Contains:
#'   \itemize{
#'     \item `k_by_cnvscore_k_results.rds`
#'     \item `cnvVar_up_to_bestKplusX.csv`
#'     \item `cnvVar_by_k_windowX.png`
#'   }
#'   }
#'
#'   \item{`query_with_ref`}{
#'   Final inferCNV run when fewer than 100 candidate reference spots are
#'   available and all reference spots are used directly.
#'   }
#'
#'   \item{`query_with_ref_kX`}{
#'   Final inferCNV run using the automatically selected reference cluster,
#'   where `X` corresponds to the selected value of k.
#'   }
#' }
#'
#' The selected reference barcodes are written to:
#'
#' \itemize{
#'   \item `raw_ref_cells.txt` when all candidate reference spots are used.
#'   \item `ref_cells.txt` when automatic reference selection is performed.
#' }
#'
#' @examples
#' \dontrun{
#' res <- run_infercnv_auto_ref(
#'   obj = se_obj,
#'   output_base = "results/infercnv/pt13",
#'   gene_pos_file = "ref/gencode_v38_gene_pos.txt",
#'   hist_column = "Histology",
#'   query_anno = c("Benign", "Cancer", "GG3", "GG4", "GG5"),
#'   ref_anno = "Benign",
#'   min_feats = 500,
#'   min_counts = 1000,
#'   cutoff = 0.1,
#'   ncores = 10,
#'   plot = FALSE
#' )
#'
#' res$ref_out_dir
#' res$k_best
#' length(res$ref_cells)
#' head(res$ref_cells)
#' }
#'
#' @import readr
#' @import dplyr
#' @import ggplot2
#' @importFrom stats kmeans
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#'
#' @export
run_infercnv_auto_ref <- function(
    obj,
    output_base,
    gene_pos_file,
    hist_column = "Histology",
    query_anno,
    ref_anno,
    min_feats = 500,
    min_counts = 1000,
    cutoff = 0.1,
    ncores = 10,
    prop_min = 0.10,
    k_window = 2,
    plot = TRUE
) {

  if (missing(query_anno)) {
    stop("You must provide 'query_anno'")
  }

  if (missing(ref_anno)) {
    stop("You must provide 'ref_anno'")
  }

  if (!hist_column %in% colnames(obj@meta.data)) {
    stop(
      paste0(
        "Column '",
        hist_column,
        "' not found in object metadata"
      )
    )
  }

  if (!all(ref_anno %in% query_anno)) {
    stop(
      paste0(
        "'ref_anno' must also be included in 'query_anno'. ",
        "Automatic reference selection currently assumes the reference ",
        "population originates from the same queried cell type population."
      )
    )
  }

  if (!is.logical(plot) || length(plot) != 1 || is.na(plot)) {
    stop("'plot' must be either TRUE or FALSE")
  }

  SeuratObject::Idents(obj) <- hist_column
  query_obj <- semla::SubsetSTData(obj,
                            idents = query_anno,
                            expression = nFeature_Spatial >= min_feats &
                              nCount_Spatial >= min_counts
  )
  ref_obj <- semla::SubsetSTData(query_obj,
                          idents = ref_anno
  )
  mat_ref <- Seurat::GetAssayData(ref_obj,
                          assay = "Spatial",
                          layer = "counts")
  n_ref <- ncol(mat_ref)
  if (n_ref < 100) {
    message(
      "Fewer than 100 reference spots found. ",
      "Using all available reference spots directly."
    )

    ref_cand <- colnames(mat_ref)
    ref_out_dir <- file.path(output_base, "query_with_ref")
    dir.create(ref_out_dir, recursive = TRUE, showWarnings = FALSE)
    write_lines(ref_cand, file.path(ref_out_dir, "raw_ref_cells.txt"))

  } else {
    message(
      "Found ", n_ref, " reference spots. ",
      "Running automatic reference selection."
    )

    sel_dir   <- file.path(output_base, "ref_no-ref", "k_selection")
    dir.create(sel_dir, recursive=TRUE, showWarnings=FALSE)
    kres_file <- file.path(sel_dir, "k_by_cnvscore_k_results.rds")

    dir_ref <- file.path(output_base, "ref_no-ref")
    dir.create(
      dir_ref,
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!file.exists(file.path(dir_ref,"run.final.infercnv_obj"))) {
      message("   * Running reference free inferCNV on selected spots based selected annotation...")
      ann_ref <- data.frame(sample=ref_obj$section_name,
                            row.names=colnames(mat_ref)
      )
      .run_infercnv(mat_ref,
                    ann_ref,
                    NULL,
                    dir_ref,
                    gene_pos_file = gene_pos_file,
                    cutoff = cutoff,
                    ncores = ncores
      )
    }
    inf_ref <- readRDS(file.path(dir_ref,"run.final.infercnv_obj"))
    expr_data <- inf_ref@expr.data
    cnv_mat   <- t(expr_data)                 # cells × genes
    cnv_var   <- colMeans((expr_data - 1)^2)  # per-cell CNV score

    # k-means loop (reuse old results if available)
    kres   <- list()
    stat   <- tibble()  # will store k / size / cnv
    found  <- FALSE
    k_start <- 2

    # Load existing kres if present
    if (file.exists(kres_file)) {
      message("   * Loading existing k-results...")
      kres_loaded <- readRDS(kres_file)
      for (kn in names(kres_loaded)) {
        kres[[kn]] <- kres_loaded[[kn]]
        best_cl <- kres_loaded[[kn]]$best_cl
        sums    <- kres_loaded[[kn]]$sums
        stat <- bind_rows(stat,
                          tibble(k   = kres_loaded[[kn]]$k,
                                 size= sums$size[sums$cluster==best_cl],
                                 cnv = sums$avg [sums$cluster==best_cl]))
      }
      # check for valley in existing data
      if (nrow(stat) >= (1 + 2 * k_window)) {
        for (idx in (k_window+1):(nrow(stat)-k_window)) {
          cur   <- stat$cnv[idx]
          prevs <- stat$cnv[(idx-k_window):(idx-1)]
          nexts <- stat$cnv[(idx+1):(idx+k_window)]
          if (all(cur < prevs) && all(cur < nexts) &&
              stat$size[idx] >= prop_min * n_ref) {
            k_best <- stat$k[idx]
            found  <- TRUE
            message("   * Found best-k in existing data: k", k_best)
            break
          }
        }
      }
      k_start <- ifelse(found, NA, max(stat$k) + 1)
    }

    # if not found, extend k-means until valley found
    if (!found) {
      k <- k_start
      while (!found) {
        message("      Computing k = ", k)
        km <- kmeans(cnv_mat, centers=k)
        df <- tibble(
          barcode = rownames(cnv_mat),
          cluster = factor(km$cluster, levels=1:k),
          cnv     = cnv_var
        )
        sums    <- df %>% group_by(cluster) %>% summarize(
          avg  = mean(cnv),
          size = n(),
          .groups = "drop"
        )
        best_cl <- sums %>% slice_min(avg) %>% pull(cluster) %>% as.character()
        kres[[paste0("k",k)]] <- list(k=k, df=df, sums=sums, best_cl=best_cl)
        stat <- bind_rows(stat,
                          tibble(k=k,
                                 size = sums$size[sums$cluster==best_cl],
                                 cnv  = sums$avg [sums$cluster==best_cl]))
        # valley detection
        if (nrow(stat) >= (1+2*k_window)) {
          idx <- nrow(stat) - k_window
          if (idx > k_window) {
            cur   <- stat$cnv[idx]
            prevs <- stat$cnv[(idx-k_window):(idx-1)]
            nexts <- stat$cnv[(idx+1):(idx+k_window)]
            if (all(cur < prevs) && all(cur < nexts) &&
                stat$size[idx] >= prop_min * n_ref) {
              k_best <- stat$k[idx]
              found  <- TRUE
              message("   * Found best-k by extension: k", k_best)
            }
          }
        }
        k <- k + 1
        if (k > 40) {
          message("   ** Reached k=40 without finding valley")
          break
        }
      }

      # save updated kres
      saveRDS(kres, kres_file)
      if (!found) {
        stop("Reached k = 40 without finding a valid reference cluster.")
      }
    }

    # Plot normalized size vs CNV variance up to best_k + window
    stat_norm <- stat %>%
      mutate(
        cnv_norm  = cnv   / max(cnv),
        size_prop = size  / n_ref
      )
    stat_plot <- stat_norm %>% filter(k <= k_best + k_window)
    write_csv(stat_plot,
              file.path(sel_dir, paste0("cnvVar_up_to_bestKplus",k_window,".csv")))
    p <- ggplot(stat_plot, aes(x=k, y=cnv_norm, group=1)) +
      geom_line(color="royalblue") + geom_point(size=3) +
      geom_text(aes(label=k), vjust=-0.5) +
      xlab("k") + ylab("Normalized CNV variance") +
      ggtitle("Auto reference selection") +
      theme_minimal()
    ggsave(file.path(sel_dir, paste0("cnvVar_by_k_window",k_window,".png")),
           p, width=6, height=4)
    message("   * Size vs CNV-var plot saved")

    # Extract best-cluster barcodes as reference
    best_r     <- kres[[paste0("k",k_best)]]
    ref_cand      <- best_r$df$barcode[best_r$df$cluster==best_r$best_cl]
    ref_out_dir <- file.path(output_base,
                             paste0("query_with_ref_k",
                                    k_best))
    dir.create(ref_out_dir, recursive=TRUE, showWarnings=FALSE)
    write_lines(ref_cand, file.path(ref_out_dir,"ref_cells.txt"))
  }

  no_plot <- !plot

  # Run inferCNV with reference on provided query
  if (file.exists(file.path(ref_out_dir,'run.final.infercnv_obj'))) {
    message("   * inferCNV-with-ref already exists")
  } else {
    expr_all <- Seurat::GetAssayData(query_obj, assay="Spatial", layer="counts")

    infercnv_group_col <- paste0(hist_column, "_pure_ref")
    meta_all <- query_obj@meta.data %>%
      mutate(
        !!infercnv_group_col := ifelse(
          rownames(.) %in% ref_cand,
          "pure_ref",
          .data[[hist_column]]
        )
      )
    ann_all  <- data.frame(group=meta_all[[infercnv_group_col]],
                           row.names=rownames(meta_all))
    message("   * Running inferCNV on query with selected reference")
    .run_infercnv(expr_all,
                  ann_all,
                  "pure_ref",
                  ref_out_dir,
                  gene_pos_file = gene_pos_file,
                  cutoff = cutoff,
                  ncores = ncores,
                  no_plot = no_plot)
  }
  message("   * Finished: ", ref_out_dir)

  return(list(
    ref_out_dir = ref_out_dir,
    ref_cells = ref_cand,
    k_best = if (exists("k_best")) k_best else NA_integer_,
    n_reference_candidates = n_ref
  ))
}

