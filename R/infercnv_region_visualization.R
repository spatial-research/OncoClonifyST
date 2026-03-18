#' Plot inferCNV signal across a chromosome or chromosome arm
#'
#' This function visualizes inferCNV expression values for genes along a selected chromosome,
#' grouped by reference and clone labels, and optionally marks a centromere-based split between chromosome arms.
#'
#' @param infercnv_obj An inferCNV object containing expression data, gene ordering, and observation groups.
#' @param chromosome Character string naming the chromosome to plot. Example `"chr8"`.
#' @param split_by_arm Logical. If `TRUE`, draws a vertical line at the centromere position to separate chromosome arms.
#' @param centromere_pos Optional named vector or list of centromere positions, where names correspond to chromosomes.
#' Required if `split_by_arm = TRUE`.
#' @param clones_to_plot Optional character vector of clone names to include. If `NULL`, all available clones and references are shown.
#' @param plot_title Optional custom plot title.
#' @param x_lab Optional x-axis label.
#' @param y_lab Optional y-axis label.
#' @param line_width Numeric line width for the chromosome arm separator.
#' @param facet_ncol Number of columns to use when faceting by clone.
#' @param free_y Logical. If `TRUE`, each facet uses its own y-axis scale.
#' @param box_fill Fill color for the boxplots.
#'
#' @return A `ggplot2` object showing gene-level inferCNV signal across the selected chromosome for each reference or clone group.
#'
#' @details
#' - Genes are ordered by genomic start position using `infercnv_obj@gene_order`.
#' - Only genes present in both `infercnv_obj@gene_order` and `infercnv_obj@expr.data` are plotted.
#' - Reference spots are labeled as `"Reference"` and observed spots are grouped using inferCNV clone assignments.
#' - If `split_by_arm = TRUE`, a dashed line is added at the supplied centromere position.
#'
#' @examples
#' plot_infercnv_chr_boxplots(
#'   infercnv_obj = infercnv_final,
#'   chromosome = "chr8",
#'   split_by_arm = TRUE,
#'   centromere_pos = c(chr8 = 45600000),
#'   clones_to_plot = c("Clone.1", "Clone.2"),
#'   box_fill = "lightblue"
#' )
#'

plot_infercnv_chr_boxplots <- function(infercnv_obj,
                                       chromosome,
                                       split_by_arm = FALSE,
                                       centromere_pos = NULL,
                                       clones_to_plot=NULL,
                                       plot_title = NULL,
                                       x_lab = NULL,
                                       y_lab = NULL,
                                       line_width = 1.2,
                                       facet_ncol = 1,
                                       free_y = TRUE,
                                       box_fill = "lightblue") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  expr <- infercnv_obj@expr.data
  gene_order <- infercnv_obj@gene_order
  
  if (is.null(expr) || nrow(expr) == 0 || ncol(expr) == 0) {
    stop("infercnv_obj@expr.data is empty or NULL.")
  }
  
  if (is.null(gene_order) || nrow(gene_order) == 0) {
    stop("infercnv_obj@gene_order is empty or NULL.")
  }
  
  if (!all(c("chr", "start") %in% colnames(gene_order))) {
    stop("infercnv_obj@gene_order must contain at least the columns 'chr' and 'start'.")
  }
  
  if (is.null(rownames(gene_order))) {
    stop("infercnv_obj@gene_order must have gene names as rownames.")
  }
  
  # Build spot annotation table from infercnv object
  spot_annot_df <- data.frame(Spot = colnames(expr), stringsAsFactors = FALSE)
  spot_annot_df$Type <- "Reference"
  spot_annot_df$Clone <- "Reference"
  
  obs_groups <- infercnv_obj@observation_grouped_cell_indices
  
  if (!is.null(obs_groups) && length(obs_groups) > 0) {
    obs_idx <- unlist(obs_groups)
    spot_annot_df$Type[obs_idx] <- "Query"
    
    for (grp in names(obs_groups)) {
      spot_annot_df$Clone[obs_groups[[grp]]] <- grp
    }
  }
  
  # Optional clone filtering
  if (is.null(clones_to_plot)) {
    clone_levels <- c("Reference", names(obs_groups))
  } else {
    missing_clones <- setdiff(clones_to_plot, names(obs_groups))
    if (length(missing_clones) > 0) {
      warning(
        paste(
          "These clones were not found and will be ignored:",
          paste(missing_clones, collapse = ", ")
        )
      )
    }
    
    valid_clones <- clones_to_plot[clones_to_plot %in% names(obs_groups)]
    
    if (length(valid_clones) == 0) {
      stop("None of the requested clones in clones_to_plot were found in the infercnv object.")
    }
    
    clone_levels <- c("Reference", valid_clones)
    
    spot_annot_df <- spot_annot_df[spot_annot_df$Clone %in% clone_levels, , drop = FALSE]
  }
  
  # Subset gene order to chromosome of interest
  chr_info <- gene_order[gene_order$chr == chromosome, , drop = FALSE]
  
  if (nrow(chr_info) == 0) {
    stop(paste("No genes found in infercnv_obj@gene_order for", chromosome))
  }
  
  chr_info <- chr_info[order(chr_info$start), , drop = FALSE]
  chr_genes <- rownames(chr_info)
  
  # Keep only genes present in expression matrix
  genes_present <- chr_genes[chr_genes %in% rownames(expr)]
  
  if (length(genes_present) == 0) {
    stop(paste("No genes from", chromosome, "found in infercnv_obj@expr.data"))
  }
  
  expr_sub <- expr[genes_present, spot_annot_df$Spot, drop = FALSE]
  
  # Convert to long format using base R
  expr_long <- as.data.frame(as.table(expr_sub), stringsAsFactors = FALSE)
  colnames(expr_long) <- c("Gene", "Spot", "Expression")
  
  # Merge with annotations
  expr_long <- merge(expr_long, spot_annot_df, by = "Spot", all.x = TRUE)
  
  if (any(is.na(expr_long$Clone))) {
    warning("Some spots in the expression matrix were not matched to Clone labels.")
  }
  
  expr_long$Gene <- factor(expr_long$Gene, levels = genes_present)
  expr_long$Clone <- factor(expr_long$Clone, levels = clone_levels)
  
  if (is.null(plot_title)) {
    plot_title <- paste("Expression of", chromosome, "genes by Reference and Clone")
  }
  
  if (is.null(x_lab)) {
    x_lab <- paste0("Gene (", chromosome, ", ordered by position)")
  }
  
  facet_scales <- if (isTRUE(free_y)) "free_y" else "fixed"
  
  p <- ggplot2::ggplot(expr_long, ggplot2::aes(x = Gene, y = Expression)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = box_fill) +
    ggplot2::facet_wrap(~ Clone, ncol = facet_ncol, scales = facet_scales) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, size = 4, vjust = 0.5, hjust = 1),
      strip.text = ggplot2::element_text(size = 12, face = "bold")
    ) +
    ggplot2::labs(
      title = plot_title,
      x = x_lab,
      y = y_lab
    )
  
  if (split_by_arm) {
    if (is.null(centromere_pos)) {
      stop("split_by_arm = TRUE requires centromere_pos to be provided.")
    }
    
    if (!chromosome %in% names(centromere_pos)) {
      stop(paste("No centromere position found for", chromosome))
    }
    
    split_pos <- centromere_pos[[chromosome]]
    
    chr_info_present <- chr_info[match(genes_present, rownames(chr_info)), , drop = FALSE]
    
    left_idx <- which(chr_info_present$start < split_pos)
    right_idx <- which(chr_info_present$start > split_pos)
    
    if (length(left_idx) == 0 || length(right_idx) == 0) {
      warning(paste("Could not place arm split line for", chromosome,
                    "- centromere position falls outside plotted gene positions"))
    } else {
      i_left <- max(left_idx)
      i_right <- min(right_idx)
      arm_split_x <- mean(c(i_left, i_right))
      
      p <- p + ggplot2::geom_vline(
        xintercept = arm_split_x,
        linetype = "dashed",
        linewidth = line_width
      )
    }
  }
  
  return(p)
}
