## ============================================================
## BLOCK 0 — Packages + paths + functions
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(grid)
})

base_dir <- "C:/Users/deise/OneDrive/Área de Trabalho/análises R/07_wetlands_fapesp_sp"

imp_sentinel_path  <- file.path(base_dir, "table/varImp/rf_importance_depressaoPeriferica_SENTINEL_v2.csv")
spec_sentinel_path <- file.path(base_dir, "table/spectralLibrary/SENTINEL_AUIS_SPECTRAL_LIBRARY.csv")

imp_embed_path     <- file.path(base_dir, "table/varImp/rf_importance_depressaoPeriferica_EMBEDDINGS_v2.csv")
spec_embed_path    <- file.path(base_dir, "table/spectralLibrary/EMBEDDINGS_AUIS_SPECTRAL_LIBRARY.csv")

plot_top50_importance <- function(df, title = "Top 50") {
  df <- df %>%
    mutate(
      importance = as.numeric(importance),
      variable   = as.character(variable)
    )
  
  top50_vars <- df %>%
    group_by(variable) %>%
    summarise(med_imp = median(importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med_imp)) %>%
    slice_head(n = 50) %>%
    pull(variable)
  
  df_top50 <- df %>%
    filter(variable %in% top50_vars)
  
  ord <- df_top50 %>%
    group_by(variable) %>%
    summarise(med_imp = median(importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(med_imp) %>%
    pull(variable)
  
  df_top50$variable <- factor(df_top50$variable, levels = ord)
  
  ggplot(df_top50, aes(x = variable, y = importance)) +
    geom_boxplot(outlier.shape = NA) +
    coord_flip() +
    xlab(NULL) +
    ylab("Variable Importance") +
    theme_minimal() +
    ggtitle(title)
}

prep_pca <- function(spec_df, ref_col = "reference") {
  drop_cols <- intersect(names(spec_df), c(".geo", "system.index"))
  if (length(drop_cols) > 0) spec_df <- spec_df %>% select(-all_of(drop_cols))
  
  ref <- spec_df[[ref_col]]
  
  X <- spec_df %>%
    select(-all_of(ref_col)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  
  # remove columns that are entirely NA
  non_all_na <- vapply(X, function(v) !all(is.na(v)), logical(1))
  X <- X[, non_all_na, drop = FALSE]
  
  # remove constant / zero-variance columns
  non_const <- vapply(X, function(v) {
    v <- v[is.finite(v)]
    length(v) > 1 && sd(v) > 0
  }, logical(1))
  X <- X[, non_const, drop = FALSE]
  
  # remove rows with NA
  keep <- complete.cases(X)
  X <- X[keep, , drop = FALSE]
  ref <- ref[keep]
  
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  ve  <- (pca$sdev^2) / sum(pca$sdev^2)
  
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$reference <- ref
  
  list(pca = pca, ve = ve, scores = scores, X = X)
}

get_top_loadings <- function(pca, scores, n = 30) {
  loadings <- as.data.frame(pca$rotation[, 1:2, drop = FALSE]) %>%
    rownames_to_column("var") %>%
    mutate(contrib = sqrt(PC1^2 + PC2^2)) %>%
    arrange(desc(contrib)) %>%
    slice_head(n = n)
  
  range_scores <- max(c(abs(scores$PC1), abs(scores$PC2)), na.rm = TRUE)
  range_load   <- max(c(abs(loadings$PC1), abs(loadings$PC2)), na.rm = TRUE)
  mult <- if (!is.finite(range_load) || range_load == 0) 1 else 0.8 * range_scores / range_load
  
  loadings <- loadings %>%
    mutate(PC1 = PC1 * mult, PC2 = PC2 * mult) %>%
    filter(is.finite(PC1), is.finite(PC2)) %>%
    filter(!(PC1 == 0 & PC2 == 0))
  
  loadings
}

recode_reference <- function(x) {
  dplyr::recode(as.character(x),
                "3"  = "Forest",
                "11" = "Wetland",
                "9"  = "GIW",
                "15" = "Pasture",
                "18" = "Agriculture",
                "25" = "Bare",
                "33" = "Water",
                .default = as.character(x)
  )
}

plot_pca_biplot <- function(scores, loadings, ve, title = NULL) {
  scores <- scores %>%
    mutate(reference = recode_reference(reference)) %>%
    filter(is.finite(PC1), is.finite(PC2))
  
  ggplot(scores, aes(x = PC1, y = PC2, color = as.factor(reference))) +
    geom_point(alpha = 0.5, size = 2) +
    geom_segment(
      data = loadings,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.02, "npc"))
    ) +
    ggrepel::geom_text_repel(
      data = loadings,
      aes(x = PC1, y = PC2, label = var),
      inherit.aes = FALSE,
      size = 3,
      max.overlaps = Inf,
      force = 1,
      box.padding = 0.3,
      point.padding = 0.2
    ) +
    xlab(paste0("PC1 (", round(100 * ve[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(100 * ve[2], 1), "%)")) +
    theme_minimal() +
    ggtitle(if (is.null(title)) "" else title) +
    scale_color_manual(values = c(
      "#E974ED", "#db4d4f", "#1f8d49", "#7a5900", "#edde8e", "#2532e4", "#519799"
    ))
}

plot_pca_scores <- function(scores, ve, title = NULL) {
  scores <- scores %>%
    mutate(reference = recode_reference(reference)) %>%
    filter(is.finite(PC1), is.finite(PC2))
  
  ggplot(scores, aes(x = PC1, y = PC2, color = as.factor(reference))) +
    geom_point(alpha = 0.5, size = 2) +
    xlab(paste0("PC1 (", round(100 * ve[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(100 * ve[2], 1), "%)")) +
    theme_minimal() +
    ggtitle(if (is.null(title)) "" else title) +
    scale_color_manual(values = c(
      "#E974ED", "#db4d4f", "#1f8d49", "#7a5900", "#edde8e", "#2532e4", "#519799"
    ))
}



## ============================================================
## BLOCK 1 — PLOT 1: Top 50 Variable Importance (SENTINEL)
## ============================================================

sentinel_imp <- read.csv(imp_sentinel_path)
sentinel_imp <- subset(sentinel_imp, year > 2023)

p1 <- plot_top50_importance(sentinel_imp, "Top 50 in Sentinel")
print(p1)

## ============================================================
## BLOCK 2 — PCA preprocessing (SENTINEL)  
## ============================================================

sentinel_spec <- read.csv(spec_sentinel_path)
sentinel_pca_obj <- prep_pca(sentinel_spec, ref_col = "reference")
sentinel_loadings <- get_top_loadings(sentinel_pca_obj$pca, sentinel_pca_obj$scores, n = 30)

## ============================================================
## BLOCK 3 — PLOT 2: PCA BIPLOT (SENTINEL) 
## ============================================================

p2 <- plot_pca_biplot(
  scores   = sentinel_pca_obj$scores,
  loadings = sentinel_loadings,
  ve       = sentinel_pca_obj$ve,
  title    = "Sentinel PCA biplot (top 30 loadings)"
)
print(p2)


## ============================================================
## BLOCK 4 — PLOT 3: PCA Scores (SENTINEL 
## ============================================================

p3 <- plot_pca_scores(
  scores = sentinel_pca_obj$scores,
  ve     = sentinel_pca_obj$ve,
  title  = "Sentinel PCA scores"
)
print(p3)

## ============================================================
## BLOCK 5 — PLOT 4: Top 50 Variable Importance (EMBEDDINGS)
## ============================================================

embed_imp <- read.csv(imp_embed_path)

p4 <- plot_top50_importance(embed_imp, "Top 50 in Embeddings")
print(p4)

## ============================================================
## BLOCK 6 — PCA preprocessing (EMBEDDINGS)
## ============================================================

embed_spec <- read.csv(spec_embed_path)
embed_pca_obj <- prep_pca(embed_spec, ref_col = "reference")
embed_loadings <- get_top_loadings(embed_pca_obj$pca, embed_pca_obj$scores, n = 30)

## ============================================================
## BLOCK 7 — PLOT 5: PCA BIPLOT (EMBEDDINGS)
## ============================================================

p5 <- plot_pca_biplot(
  scores   = embed_pca_obj$scores,
  loadings = embed_loadings,
  ve       = embed_pca_obj$ve,
  title    = "Embeddings PCA biplot (top 30 loadings)"
)
print(p5)

## ============================================================
## BLOCK 8 — PLOT 6: PCA Scores (EMBEDDINGS)
## ============================================================

p6 <- plot_pca_scores(
  scores = embed_pca_obj$scores,
  ve     = embed_pca_obj$ve,
  title  = "Embeddings PCA scores"
)
print(p6)
