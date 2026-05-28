remove_singletons <- function(
    object,
    label_col,
    sample_col = "sample_id",
    nNeighbors = 6,
    new_col = NULL,
    section_map = NULL
) {
  if (missing(label_col)) {
    
    stop("You must provide 'label_col'")
    
  }
  if (is.null(new_col)) {
    new_col <- paste0(label_col, "_nosingletons")
  }
  
  spatnet <- GetSpatialNetwork(object, nNeighbors = nNeighbors)
  
  meta_df <- object[[]] %>%
    rownames_to_column("barcode") %>%
    select(barcode, all_of(sample_col), all_of(label_col)) %>%
    rename(
      section = all_of(sample_col),
      label = all_of(label_col)
    )
  
  if (is.null(section_map)) {
    stop("Please provide a section_map, e.g. c('A1'='-1', 'A2'='-2', ...)")
  }
  
  clean_one_section <- function(net_df, meta_df_section) {
    edges <- net_df %>%
      select(from, to) %>%
      left_join(
        meta_df_section %>% select(barcode, from_label = label),
        by = c("from" = "barcode")
      ) %>%
      left_join(
        meta_df_section %>% select(barcode, to_label = label),
        by = c("to" = "barcode")
      ) %>%
      mutate(
        same_label_neighbor = !is.na(from_label) & !is.na(to_label) & from_label == to_label
      )
    
    support_tbl <- edges %>%
      group_by(from) %>%
      summarise(
        n_same_label_neighbors = sum(same_label_neighbor, na.rm = TRUE),
        .groups = "drop"
      )
    
    meta_df_section %>%
      left_join(support_tbl, by = c("barcode" = "from")) %>%
      mutate(
        n_same_label_neighbors = dplyr::coalesce(n_same_label_neighbors, 0L),
        cleaned_label = case_when(
          !is.na(label) & n_same_label_neighbors == 0 ~ NA,
          TRUE ~ label
        )
      ) %>%
      select(barcode, cleaned_label, n_same_label_neighbors)
  }
  
  cleaned_list <- lapply(names(spatnet), function(net_id) {
    sid <- names(section_map)[match(paste0("-", net_id), section_map)]
    
    if (is.na(sid)) {
      stop(paste("No sample_id mapping found for network section:", net_id))
    }
    
    net_df <- spatnet[[net_id]]
    meta_df_section <- meta_df %>% filter(section == sid)
    
    clean_one_section(net_df, meta_df_section)
  })
  
  cleaned_df <- bind_rows(cleaned_list)
  
  object[[new_col]] <- cleaned_df$cleaned_label[
    match(colnames(object), cleaned_df$barcode)
  ]
  
  object[[paste0(new_col, "_n_same_label_neighbors")]] <- cleaned_df$n_same_label_neighbors[
    match(colnames(object), cleaned_df$barcode)
  ]
  
  object
}


remove_small_spatial_islands <- function(
    object,
    label_col = "clone_na0",
    sample_col = "sample_id",
    nNeighbors = 6,
    min_island_size = 2,
    new_col = NULL,
    section_map = NULL
) {
  if (is.null(new_col)) {
    new_col <- paste0(label_col, "_no_small_islands")
  }
  
  if (is.null(section_map)) {
    stop("Please provide a section_map, e.g. c('A1'='-1', 'A2'='-2', ...)")
  }
  
  spatnet <- GetSpatialNetwork(object, nNeighbors = nNeighbors)
  
  meta_df <- object[[]] %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::select(
      barcode,
      dplyr::all_of(sample_col),
      dplyr::all_of(label_col)
    ) %>%
    dplyr::rename(
      section = dplyr::all_of(sample_col),
      label = dplyr::all_of(label_col)
    )
  
  clean_one_section <- function(net_df, meta_df_section, min_island_size) {
    same_label_edges <- net_df %>%
      dplyr::select(from, to) %>%
      dplyr::left_join(
        meta_df_section %>% dplyr::select(barcode, from_label = label),
        by = c("from" = "barcode")
      ) %>%
      dplyr::left_join(
        meta_df_section %>% dplyr::select(barcode, to_label = label),
        by = c("to" = "barcode")
      ) %>%
      dplyr::filter(
        !is.na(from_label),
        !is.na(to_label),
        from_label == to_label
      ) %>%
      dplyr::select(from, to)
    
    g <- igraph::graph_from_data_frame(
      same_label_edges,
      directed = FALSE,
      vertices = meta_df_section %>% dplyr::select(name = barcode)
    )
    
    comp <- igraph::components(g)
    
    comp_df <- tibble::tibble(
      barcode = names(comp$membership),
      component = as.integer(comp$membership)
    ) %>%
      dplyr::left_join(
        tibble::tibble(
          component = seq_along(comp$csize),
          island_size = as.integer(comp$csize)
        ),
        by = "component"
      )
    
    meta_df_section %>%
      dplyr::left_join(comp_df, by = "barcode") %>%
      dplyr::mutate(
        island_size = dplyr::case_when(
          is.na(label) ~ NA_integer_,
          TRUE ~ dplyr::coalesce(island_size, 1L)
        ),
        cleaned_label = dplyr::case_when(
          !is.na(label) & island_size < min_island_size ~ NA,
          TRUE ~ label
        )
      ) %>%
      dplyr::select(barcode, cleaned_label, island_size)
  }
  
  cleaned_list <- lapply(names(spatnet), function(net_id) {
    sid <- names(section_map)[match(paste0("-", net_id), section_map)]
    
    if (is.na(sid)) {
      stop(paste("No sample_id mapping found for network section:", net_id))
    }
    
    net_df <- spatnet[[net_id]]
    
    meta_df_section <- meta_df %>%
      dplyr::filter(section == sid)
    
    clean_one_section(
      net_df = net_df,
      meta_df_section = meta_df_section,
      min_island_size = min_island_size
    )
  })
  
  cleaned_df <- dplyr::bind_rows(cleaned_list)
  
  object[[new_col]] <- cleaned_df$cleaned_label[
    match(colnames(object), cleaned_df$barcode)
  ]
  
  object[[paste0(new_col, "_island_size")]] <- cleaned_df$island_size[
    match(colnames(object), cleaned_df$barcode)
  ]
  
  object
}


#' Generate overlap-aware neighborhood layers around labeled regions.
#'
#' Expands region labels iteratively using spatial neighbors while
#' resolving overlapping assignments between regions.
#'
#' @param object Seurat object.
#' @param column_name Metadata column containing region labels.
#' @param n_layers Number of neighborhood expansion layers.
#' @param exclude_original_labeled_spots_from_other_rings Logical; whether
#' to prevent original labeled spots from being assigned to other regions.
#'
#' @return
#' The input Seurat object with additional metadata columns:
#'
#' \describe{
#'   \item{region_neighbors_n_hits}{
#'   Number of distinct region labels reaching each spot.
#'   }
#'
#'   \item{region_neighbors_is_ambiguous}{
#'   TRUE if a spot was reached by multiple region labels.
#'   }
#'
#'   \item{region_neighbors_n_layer}{
#'   Final non-ambiguous neighborhood layer assignment.
#'   Original labeled spots are not included and remain `NA`.
#'   }
#'
#'   \item{region_neighbors_label}{
#'   Final non-ambiguous region label assignment.
#'   }
#'
#'   \item{region_neighbors_all_hits}{
#'   Raw overlapping region-layer assignments before ambiguity filtering,
#'   stored as semicolon-separated `label:n_layer` pairs.
#'   }
#' }
#'
#' @export
make_region_neighbors_overlap_aware <- function(
    object,
    column_name,
    n_layers = 1,
    exclude_original_labeled_spots_from_other_rings = TRUE
) {
  
  if (missing(column_name)) {
    stop("You must provide 'column_name'")
  }
  
  if (!column_name %in% colnames(object@meta.data)) {
    stop(paste0("Column '", column_name, "' not found in metadata"))
  }
  
  base_labels <- object[[column_name]][, 1]
  labels <- unique(na.omit(base_labels))
  
  ring_mat <- matrix(
    NA_integer_,
    nrow = ncol(object),
    ncol = length(labels),
    dimnames = list(colnames(object), labels)
  )
  
  for (label in labels) {
    
    message("Processing label: ", label)
    
    temp_col <- paste0("tmp_expand_", make.names(label))
    
    object[[temp_col]] <- ifelse(base_labels == label, label, NA)
    
    for (i in seq_len(n_layers)) {
      
      key <- paste0("n", i, "_", make.names(label))
      
      object <- RegionNeighbors(
        object,
        column_name = temp_col,
        column_labels = label,
        mode = "outer",
        column_key = key
      )
      
      md <- object@meta.data
      
      nb_cols <- grep(
        paste0("^", key),
        colnames(md),
        value = TRUE
      )
      
      if (length(nb_cols) == 0) {
        warning("No RegionNeighbors columns found for key: ", key)
        next
      }
      
      is_neighbor <- rowSums(
        md[, nb_cols, drop = FALSE] == as.character(label),
        na.rm = TRUE
      ) > 0
      
      not_already_in_this_expansion <- is.na(object[[temp_col]][, 1])
      
      if (exclude_original_labeled_spots_from_other_rings) {
        allowed_spot <- is.na(base_labels)
      } else {
        allowed_spot <- base_labels != label | is.na(base_labels)
      }
      
      new_ring_spots <- is_neighbor &
        not_already_in_this_expansion &
        allowed_spot
      
      ring_mat[new_ring_spots, label] <- i
      
      object[[temp_col]][new_ring_spots, 1] <- label
    }
  }
  
  region_neighbors_n_hits <- rowSums(!is.na(ring_mat))
  
  region_neighbors_label <- apply(ring_mat, 1, function(x) {
    hits <- names(x)[!is.na(x)]
    if (length(hits) == 1) hits else NA_character_
  })
  
  region_neighbors_n_layer <- apply(ring_mat, 1, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) == 1) vals else NA_integer_
  })
  
  region_neighbors_all_hits <- apply(ring_mat, 1, function(x) {
    valid <- !is.na(x)
    
    if (!any(valid)) {
      return(NA_character_)
    }
    
    paste0(
      names(x)[valid],
      ":n",
      x[valid],
      collapse = ";"
    )
  })
  
  object$region_neighbors_n_hits <- region_neighbors_n_hits
  object$region_neighbors_is_ambiguous <- region_neighbors_n_hits > 1
  
  object$region_neighbors_n_layer <- ifelse(
    region_neighbors_n_hits == 1,
    region_neighbors_n_layer,
    NA_integer_
  )
  
  object$region_neighbors_label <- ifelse(
    region_neighbors_n_hits == 1,
    region_neighbors_label,
    NA_character_
  )
  
  object$region_neighbors_all_hits <- region_neighbors_all_hits
  
  # original_region_spots <- !is.na(base_labels)
  # 
  # object$region_neighbors_n_layer[original_region_spots] <- 0
  # object$region_neighbors_label[original_region_spots] <- base_labels[original_region_spots]
  # object$region_neighbors_is_ambiguous[original_region_spots] <- FALSE
  
  cols_to_remove <- grep(
    "^(tmp_expand_|n[0-9]+_)",
    colnames(object@meta.data),
    value = TRUE
  )
  
  object@meta.data <- object@meta.data[
    ,
    !colnames(object@meta.data) %in% cols_to_remove,
    drop = FALSE
  ]
  
  return(object)
}


#' Plot cumulative region-neighbor layers by label
#'
#' Visualizes spots assigned to region-neighbor layers up to a selected
#' maximum layer. Spots are colored by their final non-ambiguous
#' `region_neighbors_label`.
#'
#' @param object Seurat object containing region-neighbor metadata generated by
#' `make_region_neighbors_overlap_aware()`.
#' @param max_layer Maximum neighborhood layer to include in the plot.
#' For example, `1` plots n1 only, while `3` plots n1+n2+n3.
#' @param pt_size Point size passed to `MapLabels()`.
#' @param section_number Optional section number passed to `MapLabels()`.
#' Defaults to `NULL`, which plots all sections.
#' @param ncol Number of columns for faceted section plots.
#' @param drop_na Logical; whether to drop spots with `NA` labels from the plot.
#' @param image_use A character specifying image type to use.
#' Passed to `MapLabels()`. Defaults to `NULL`.
#'
#' @return
#' A ggplot/patchwork object showing cumulative neighborhood layers colored
#' by region label.
#'
#' @export
plot_region_neighbor_layers <- function(
    object,
    max_layer,
    pt_size = 1,
    section_number = NULL,
    ncol = NULL,
    drop_na = TRUE,
    image_use = NULL
) {
  
  plot_col <- paste0("plot_region_neighbors_1_", max_layer)
  
  object[[plot_col]] <- ifelse(
    object$region_neighbors_n_layer %in% seq_len(max_layer),
    object$region_neighbors_label,
    NA_character_
  )
  
  legend_title <- ifelse(
    max_layer == 1,
    "Region neighbors n1",
    paste0("Region neighbors n1-n", max_layer)
  )
  
  MapLabels(
    object,
    column_name = plot_col,
    drop_na = drop_na,
    label_by = "sample_id",
    section_number = section_number,
    ncol = ncol,
    pt_size = pt_size,
    image_use = image_use,
    mar = c(0, 0, 0, 0)
  ) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 14)
    ) &
    ggplot2::guides(fill = ggplot2::guide_legend(
      title = legend_title,
      override.aes = list(size = 4)
    ))
}