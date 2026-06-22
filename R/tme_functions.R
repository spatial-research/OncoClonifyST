# Helper function to compute a median distance cut-off.
#
# This internal helper returns a median `maxDist` value for a given sample/section
# depending on the chosen method. It extracts pixel coordinates from the `Staffli`
# object, computes the first nearest neighbour distances for all barcodes
# belonging to the specified section, takes their median and multiplies by
# 1.2.  The resulting value can be passed directly to `GetSpatialNetwork()`
# or `RegionNeighbors()` to obtain a median instead of minimum neighbour distance.
#
# @param object A Seurat object containing spatial data and a Staffli object.
# @param sample_id A single character identifying the section (matching the
#        values in `sample_col`).
# @param sample_col The name of the column in metadata that stores sample IDs.
#        Defaults to "sample_id".
# @param method Character string: either "minimum" or "median".  See
#        details above.
# @return A numeric distance threshold when `method = "median"`, otherwise
#        `NULL` if `method = "minimum"`.
# @importFrom semla GetStaffli
# @importFrom dbscan kNN
calculate_median_maxDist <- function(object, sid) {
  st_obj <- semla::GetStaffli(object)
  spots <- colnames(object)[object@tools$Staffli@meta_data$sampleID == sid]
  coords <- st_obj@meta_data[spots, c("pxl_col_in_fullres", "pxl_row_in_fullres")]
  coords <- coords[complete.cases(coords), , drop = FALSE]
  if (nrow(coords) < 2) return(NULL)
  nn_res <- dbscan::kNN(coords, k = 1)
  return(median(nn_res$dist) * 1.2)
}

#' Remove spatial singleton labels
#'
#' Removes labeled spots that do not have any neighboring spot with the same
#' label within each tissue section. A spot is considered a singleton if it has
#' a non-missing label but zero same-label neighbors in the spatial neighbor
#' graph.
#'
#' @param object Seurat object.
#' @param label_col Metadata column containing labels to clean.
#' @param sample_col Metadata column identifying tissue sections or samples.
#' Defaults to `"sample_id"`.
#' @param nNeighbors Number of spatial neighbors to use when constructing the
#' spatial neighbor graph. Passed to `GetSpatialNetwork()`.
#' @param maxDist_method Character string specifying how to compute the
#' distance cut‐off.  Either "minimum" (use semla’s default based on
#' the minimum nearest neighbour distance) or "median" (compute a
#' cut‐off from the median nearest neighbour distance).  Defaults to
#' "minimum".
#' @param new_col Name of the metadata column where cleaned labels will be
#' stored. If `NULL`, defaults to `paste0(label_col, "_nosingletons")`.
#' @param section_map Named character vector mapping sample IDs to spatial
#' network section IDs, for example `c("A1" = "-1", "A2" = "-2")`.
#'
#' @return Seurat object with cleaned labels and same-label neighbor counts.
#'
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom semla GetSpatialNetwork
#' @importFrom magrittr %>%
#'
#' @export
remove_singletons <- function(
    object,
    label_col,
    sample_col = "sample_id",
    nNeighbors = 6,
    maxDist_method = c("minimum","median"),
    new_col = NULL,
    section_map = NULL
) {
  if (missing(label_col)) {
    stop("You must provide 'label_col'")
  }

  maxDist_method <- match.arg(maxDist_method)

  if (is.null(new_col)) {
    new_col <- paste0(label_col, "_nosingletons")
  }

  if (!sample_col %in% colnames(object@meta.data)) {
    stop(sprintf(
      "Column '%s' not found in object metadata. Please supply a valid sample_col.",
      sample_col
    ))
  }

  if (is.null(section_map)) {
    stop("Please provide a section_map, e.g. c('A1'='-1', 'A2'='-2', ...)")
  }

  if (maxDist_method == "median") {
    section_ids <- unique(object@tools$Staffli@meta_data$sampleID)
    spatnet <- setNames(
      lapply(section_ids, function(sid) {
        md <- calculate_median_maxDist(object, sid)
        sub_obj <- subset(object, cells = colnames(object)[object@tools$Staffli@meta_data$sampleID == sid])
        GetSpatialNetwork(sub_obj, nNeighbors = nNeighbors, maxDist = md)[[1]]
      }),
      section_ids
    )
  } else {
    spatnet <- GetSpatialNetwork(object, nNeighbors = nNeighbors)
  }

  meta_df <- object[[]] %>%
    rownames_to_column("barcode") %>%
    select(barcode, all_of(sample_col), all_of(label_col)) %>%
    rename(
      section = all_of(sample_col),
      label = all_of(label_col)
    )

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
        same_label_neighbor = !is.na(from_label) &
          !is.na(to_label) &
          from_label == to_label
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
        n_same_label_neighbors = coalesce(n_same_label_neighbors, 0L),
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

  object[[paste0(new_col, "_n_same_label_neighbors")]] <-
    cleaned_df$n_same_label_neighbors[
      match(colnames(object), cleaned_df$barcode)
    ]

  object
}


#' Remove small spatial islands of labels
#'
#' Removes connected components of same-label spots that are smaller than a
#' specified minimum island size within each tissue section.
#'
#' @param object Seurat object.
#' @param label_col Metadata column containing labels to clean.
#' Defaults to `"clone_na0"`.
#' @param sample_col Metadata column identifying tissue sections or samples.
#' Defaults to `"sample_id"`.
#' @param nNeighbors Number of spatial neighbors to use when constructing the
#' spatial neighbor graph. Passed to `GetSpatialNetwork()`.
#' @param min_island_size Minimum connected island size to retain.
#' @param maxDist_method Character string specifying how to compute the
#' distance cut‐off.  Either "minimum" or "median".  Defaults to
#' "minimum".
#' @param new_col Name of the metadata column where cleaned labels will be
#' stored. If `NULL`, defaults to `paste0(label_col, "_no_small_islands")`.
#' @param section_map Named character vector mapping sample IDs to spatial
#' network section IDs, for example `c("A1" = "-1", "A2" = "-2")`.
#'
#' @return Seurat object with cleaned labels and island sizes.
#'
#' @import dplyr
#' @importFrom tibble rownames_to_column tibble
#' @importFrom semla GetSpatialNetwork
#' @importFrom igraph graph_from_data_frame components
#' @importFrom magrittr %>%
#'
#' @export
remove_small_spatial_islands <- function(
    object,
    label_col = "clone_na0",
    sample_col = "sample_id",
    nNeighbors = 6,
    min_island_size = 2,
    maxDist_method = c("minimum", "median"),
    new_col = NULL,
    section_map = NULL
) {

  maxDist_method <- match.arg(maxDist_method)

  if (is.null(new_col)) {
    new_col <- paste0(label_col, "_no_small_islands")
  }

  if (!sample_col %in% colnames(object@meta.data)) {
    stop(sprintf(
      "Column '%s' not found in object metadata. Please supply a valid sample_col.",
      sample_col
    ))
  }

  if (is.null(section_map)) {
    stop("Please provide a section_map, e.g. c('A1'='-1', 'A2'='-2', ...)")
  }

  if (maxDist_method == "median") {
    section_ids <- unique(object@tools$Staffli@meta_data$sampleID)
    spatnet <- setNames(
      lapply(section_ids, function(sid) {
        md <- calculate_median_maxDist(object, sid)
        sub_obj <- subset(object, cells = colnames(object)[object@tools$Staffli@meta_data$sampleID == sid])
        GetSpatialNetwork(sub_obj, nNeighbors = nNeighbors, maxDist = md)[[1]]
      }),
      section_ids
    )
  } else {
    spatnet <- GetSpatialNetwork(object, nNeighbors = nNeighbors)
  }

  meta_df <- object[[]] %>%
    rownames_to_column("barcode") %>%
    select(
      barcode,
      all_of(sample_col),
      all_of(label_col)
    ) %>%
    rename(
      section = all_of(sample_col),
      label = all_of(label_col)
    )

  clean_one_section <- function(net_df, meta_df_section, min_island_size) {
    same_label_edges <- net_df %>%
      select(from, to) %>%
      left_join(
        meta_df_section %>% select(barcode, from_label = label),
        by = c("from" = "barcode")
      ) %>%
      left_join(
        meta_df_section %>% select(barcode, to_label = label),
        by = c("to" = "barcode")
      ) %>%
      filter(
        !is.na(from_label),
        !is.na(to_label),
        from_label == to_label
      ) %>%
      select(from, to)

    g <- graph_from_data_frame(
      same_label_edges,
      directed = FALSE,
      vertices = meta_df_section %>% select(name = barcode)
    )

    comp <- components(g)

    comp_df <- tibble(
      barcode = names(comp$membership),
      component = as.integer(comp$membership)
    ) %>%
      left_join(
        tibble(
          component = seq_along(comp$csize),
          island_size = as.integer(comp$csize)
        ),
        by = "component"
      )

    meta_df_section %>%
      left_join(comp_df, by = "barcode") %>%
      mutate(
        island_size = case_when(
          is.na(label) ~ NA_integer_,
          TRUE ~ coalesce(island_size, 1L)
        ),
        cleaned_label = case_when(
          !is.na(label) & island_size < min_island_size ~ NA,
          TRUE ~ label
        )
      ) %>%
      select(barcode, cleaned_label, island_size)
  }

  cleaned_list <- lapply(names(spatnet), function(net_id) {
    sid <- names(section_map)[match(paste0("-", net_id), section_map)]

    if (is.na(sid)) {
      stop(paste("No sample_id mapping found for network section:", net_id))
    }

    net_df <- spatnet[[net_id]]
    meta_df_section <- meta_df %>% filter(section == sid)

    clean_one_section(
      net_df = net_df,
      meta_df_section = meta_df_section,
      min_island_size = min_island_size
    )
  })

  cleaned_df <- bind_rows(cleaned_list)

  object[[new_col]] <- cleaned_df$cleaned_label[
    match(colnames(object), cleaned_df$barcode)
  ]

  object[[paste0(new_col, "_island_size")]] <- cleaned_df$island_size[
    match(colnames(object), cleaned_df$barcode)
  ]

  object
}


#' Generate overlap-aware neighborhood layers around labeled regions
#'
#' Expands region labels iteratively using spatial neighbors while resolving
#' overlapping assignments between regions.
#'
#' @param object Seurat object.
#' @param column_name Metadata column containing region labels.
#' @param n_layers Number of neighborhood expansion layers.
#' @param exclude_original_labeled_spots_from_other_rings Logical; whether to
#' prevent original labeled spots from being assigned to other regions.
#' @param maxDist_method Character string specifying how to compute the
#' distance cut‑off passed to `RegionNeighbors()`.  Either "minimum"
#' (default; use semla’s minimum nearest neighbour threshold) or
#' "median" (compute the median nearest neighbour distance × 1.2).
#'
#' @return Seurat object with additional region-neighbor metadata columns.
#'
#' @importFrom semla RegionNeighbors
#' @importFrom stats na.omit
#'
#' @export
make_region_neighbors_overlap_aware <- function(
    object,
    column_name,
    n_layers = 1,
    exclude_original_labeled_spots_from_other_rings = TRUE,
    maxDist_method = c("minimum", "median")
) {
  if (missing(column_name)) {
    stop("You must provide 'column_name'")
  }

  if (!column_name %in% colnames(object@meta.data)) {
    stop(paste0("Column '", column_name, "' not found in metadata"))
  }

  maxDist_method <- match.arg(maxDist_method)

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

    for (sid in unique(object@tools$Staffli@meta_data$sampleID)) {
      cells <- which(object@tools$Staffli@meta_data$sampleID == sid)
      if (length(cells) == 0)
        next

      sub_obj <- subset(object, cells = colnames(object)[cells])

      temp_col_name <- paste0("tmp_expand_", make.names(label))
      temp_col <- rep(NA_character_, length(cells))
      temp_col[base_labels[cells] == label] <- label
      sub_obj[[temp_col_name]] <- temp_col

      md <- NULL
      if (maxDist_method == "median") {
        md <- calculate_median_maxDist(object, sid)
      }

      for (i in seq_len(n_layers)) {
        key <- paste0("n", i, "_", make.names(label))

        sub_obj <- semla::RegionNeighbors(
          sub_obj,
          column_name = temp_col_name,
          column_labels = label,
          mode = "outer",
          column_key = key,
          maxDist = md
        )

        md_data <- sub_obj@meta.data
        nb_cols <- grep(paste0("^", key), colnames(md_data), value = TRUE)
        if (length(nb_cols) == 0) {
          warning("No RegionNeighbors columns found for key: ", key)
          next
        }

        is_neighbor <- rowSums(
          md_data[, nb_cols, drop = FALSE] == as.character(label),
          na.rm = TRUE
        ) > 0

        not_already <- is.na(sub_obj[[temp_col_name]][, 1])

        if (exclude_original_labeled_spots_from_other_rings) {
          allowed <- is.na(base_labels[cells])
        } else {
          allowed <- base_labels[cells] != label | is.na(base_labels[cells])
        }

        new_ring <- is_neighbor & not_already & allowed

        if (any(new_ring))
          ring_mat[cells[new_ring], label] <- i

        tmp <- sub_obj[[temp_col_name]][, 1]
        tmp[new_ring] <- label
        sub_obj[[temp_col_name]] <- tmp
      }
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

  object
}


#' Plot cumulative region-neighbor layers by label
#'
#' Visualizes spots assigned to region-neighbor layers up to a selected maximum
#' layer. Spots are colored by their final non-ambiguous
#' `region_neighbors_label`.
#'
#' @param object Seurat object containing region-neighbor metadata generated by
#' `make_region_neighbors_overlap_aware()`.
#' @param max_layer Maximum neighborhood layer to include in the plot.
#' @param pt_size Point size passed to `MapLabels()`.
#' @param section_number Optional section number passed to `MapLabels()`.
#' @param ncol Number of columns for faceted section plots.
#' @param drop_na Logical; whether to drop spots with `NA` labels from the plot.
#' @param image_use Character specifying image type to use. Passed to
#' `MapLabels()`. Defaults to `NULL`.
#'
#' @return A ggplot/patchwork object showing cumulative neighborhood layers.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom semla MapLabels
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
    plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) &
    guides(fill = guide_legend(
      title = legend_title,
      override.aes = list(size = 4)
    ))
}
