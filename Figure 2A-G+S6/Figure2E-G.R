library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(scales)
library(ggpubr)
library(rstatix)

# ── Data ──────────────────────────────────────────────────────────────────────

ofDDAmv <- cbind(read.csv("csv/E033_DDA_Fusion_Missing_Vals.csv"),   acquisition = "DDA", MS = "Orbitrap Fusion")
ttDDAmv <- cbind(read.csv("csv/E033_DDA_TimsTOF_Missing_Vals.csv"),  acquisition = "DDA", MS = "timsTOF")
ttDIAmv <- cbind(read.csv("csv/E033_DIA_TimsTOF_Missing_Vals.csv"),  acquisition = "DIA", MS = "timsTOF")

datasets <- list(ofDDAmv, ttDDAmv, ttDIAmv)

# ── Constants ─────────────────────────────────────────────────────────────────

ds_levels <- c("DDA Orbitrap Fusion", "DDA timsTOF", "DIA timsTOF")
ds_labels <- c("Orbitrap Fusion\nDDA", "timsTOF\nDDA", "timsTOF\nDIA")

pal_n <- c("0" = "#D9D9D9", "1" = "#92C5DE", "2" = "#2166AC", "3" = "#053061")

n_labels_bar    <- c("0" = "0 (not detected)", "1" = "1", "2" = "2", "3" = "3")
n_labels_violin <- c("0" = "0 (imputed)",      "1" = "1", "2" = "2", "3" = "3")

# ── Helpers ───────────────────────────────────────────────────────────────────

add_missing_cols <- function(df) {
  gr_cols  <- colnames(df)[3:5]
  igg_cols <- colnames(df)[6:8]
  
  df$MissingGR  <- rowSums(is.na(df[, gr_cols]))
  df$MissingIgG <- rowSums(is.na(df[, igg_cols]))
  df$MeanGR     <- ifelse(df$MissingGR  < 3, rowMeans(df[, gr_cols],  na.rm = TRUE), NA)
  df$MeanIgG    <- ifelse(df$MissingIgG < 3, rowMeans(df[, igg_cols], na.rm = TRUE), NA)
  df
}

prep_grid <- function(datasets) {
  df_all   <- bind_rows(lapply(datasets, add_missing_cols))
  prot_col <- colnames(df_all)[1]
  
  df_all$Dataset     <- factor(paste(df_all$acquisition, df_all$MS), levels = ds_levels)
  df_all$significant <- as.logical(df_all$significant)
  df_all$significant[is.na(df_all$significant)] <- FALSE
  
  df_grid <- df_all %>%
    select(all_of(prot_col), Dataset, MeanGR, MeanIgG, MissingGR, MissingIgG, significant) %>%
    complete(.data[[prot_col]], Dataset) %>%
    mutate(
      significant = ifelse(is.na(significant), FALSE, significant),
      MissingGR   = ifelse(is.na(MissingGR),  3L, as.integer(MissingGR)),
      MissingIgG  = ifelse(is.na(MissingIgG), 3L, as.integer(MissingIgG)),
      log_GR      = ifelse(MeanGR  > 0, log10(MeanGR),  NA),
      log_IgG     = ifelse(MeanIgG > 0, log10(MeanIgG), NA),
      N_det_GR    = factor(3L - MissingGR,  levels = 0:3),
      N_det_IgG   = factor(3L - MissingIgG, levels = 0:3)
    )
  
  list(grid = df_grid, prot_col = prot_col)
}

# Median-shift each dataset to the reference (DIA timsTOF).
# GR and IgG are normalised against their own core sets independently,
# since the proteins quantified in all replicates differ between conditions.
normalise_grid <- function(df_grid, prot_col) {
  ref <- ds_levels[3]
  
  core_for <- function(missing_col, log_col) {
    df_grid %>%
      filter(.data[[missing_col]] == 0, !is.na(.data[[log_col]])) %>%
      group_by(.data[[prot_col]]) %>%
      filter(n_distinct(Dataset) == length(ds_levels)) %>%
      pull(.data[[prot_col]]) %>% unique()
  }
  
  compute_shifts <- function(core_prots, log_col, shift_col) {
    df_grid %>%
      filter(.data[[prot_col]] %in% core_prots, !is.na(.data[[log_col]])) %>%
      group_by(Dataset) %>%
      summarise(med = median(.data[[log_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(!!shift_col := med - med[Dataset == ref])
  }
  
  shifts_GR  <- compute_shifts(core_for("MissingGR",  "log_GR"),  "log_GR",  "shift_GR")
  shifts_IgG <- compute_shifts(core_for("MissingIgG", "log_IgG"), "log_IgG", "shift_IgG")
  
  df_grid %>%
    left_join(shifts_GR  %>% select(Dataset, shift_GR),  by = "Dataset") %>%
    left_join(shifts_IgG %>% select(Dataset, shift_IgG), by = "Dataset") %>%
    mutate(
      Norm_GR  = log_GR  - shift_GR,
      Norm_IgG = log_IgG - shift_IgG,
      # Exclude proteins with no Replicate coverage from plots (not imputed)
      Plot_GR  = ifelse(MissingGR  == 3, NA, Norm_GR),
      Plot_IgG = ifelse(MissingIgG == 3, NA, Norm_IgG)
    )
}

# ── Panel A: stacked bar — detection completeness ─────────────────────────────

make_panel_A <- function(datasets) {
  out <- prep_grid(datasets)
  
  out$grid %>%
    pivot_longer(c(N_det_GR, N_det_IgG), names_to = "Group", values_to = "N_detected") %>%
    mutate(Group = factor(ifelse(Group == "N_det_GR", "GR", "IgG"), levels = c("GR", "IgG"))) %>%
    group_by(Dataset, Group, N_detected) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Dataset, Group) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = Dataset, y = prop, fill = N_detected)) +
    geom_col(position = "fill", colour = "grey50", linewidth = 0.2, width = 0.72) +
    facet_wrap(~ Group) +
    scale_fill_manual(values = pal_n, labels = n_labels_bar, name = "Replicate coverage") +
    scale_x_discrete(labels = ds_labels) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(x = NULL, y = "Proportion of proteins") +
    theme_pubr() +
    theme(
      strip.background = element_blank(),
      axis.text.x      = element_text(angle = 90, vjust = 0.7, hjust = 1, size = 10),
      legend.position  = "none"
    )
}

# ── Panel B: detection heatmap ────────────────────────────────────────────────

make_panel_B <- function(datasets) {
  out      <- prep_grid(datasets)
  df_grid  <- out$grid
  prot_col <- out$prot_col
  
  # Sort proteins by detection completeness, DIA first
  prot_order <- df_grid %>%
    select(all_of(prot_col), Dataset, MissingGR) %>%
    pivot_wider(names_from = Dataset, values_from = MissingGR, values_fill = 3) %>%
    arrange(
      .data[[ds_levels[3]]],
      .data[[ds_levels[2]]],
      .data[[ds_levels[1]]]
    ) %>%
    pull(.data[[prot_col]])
  
  df_grid %>%
    pivot_longer(c(MissingGR, MissingIgG), names_to = "Group", values_to = "Missing") %>%
    mutate(
      Group    = factor(ifelse(Group == "MissingGR", "GR", "IgG"), levels = c("GR", "IgG")),
      Detected = factor(3L - Missing, levels = 0:3)
    ) %>%
    mutate(across(all_of(prot_col), ~ factor(., levels = prot_order))) %>%
    ggplot(aes(x = Dataset, y = .data[[prot_col]], fill = Detected)) +
    geom_tile(colour = NA) +
    facet_wrap(~ Group, ncol = 2) +
    scale_fill_manual(values = pal_n, labels = n_labels_bar, name = "Replicate coverage", drop = FALSE) +
    scale_x_discrete(position = "bottom", labels = ds_labels, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    guides(fill = guide_legend(reverse = TRUE, keyheight = unit(0.4, "cm"), keywidth = unit(0.3, "cm"))) +
    labs(x = NULL, y = NULL) +
    theme_pubr() +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.x      = element_text(angle = 90, vjust = 0.7, hjust = 1, size = 10),
      strip.background = element_blank(),
      axis.line        = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.spacing    = unit(1, "lines"),
      legend.position  = "right"
    )
}

# ── Panel C: intensity distributions by replicate count ───────────────────────

make_panel_C <- function(datasets) {
  out      <- prep_grid(datasets)
  df_grid  <- normalise_grid(out$grid, out$prot_col)
  prot_col <- out$prot_col
  
  df_plot <- bind_rows(
    df_grid %>% select(all_of(prot_col), Dataset, N_detected = N_det_GR,  Intensity = Plot_GR)  %>% mutate(Group = "GR"),
    df_grid %>% select(all_of(prot_col), Dataset, N_detected = N_det_IgG, Intensity = Plot_IgG) %>% mutate(Group = "IgG")
  ) %>%
    filter(!is.na(Intensity), !is.nan(Intensity), N_detected != 0) %>%
    mutate(
      Group      = factor(Group,      levels = c("GR", "IgG")),
      N_detected = factor(N_detected, levels = c(1, 2, 3))
    )
  
  # One-way ANOVA + Tukey HSD: N=1, N=2 vs N=3 per platform and condition
  stat_results <- df_plot %>%
    group_by(Group, Dataset) %>%
    filter(n_distinct(N_detected) >= 2) %>%
    group_by(Group, Dataset, N_detected) %>%
    filter(n() >= 3) %>%
    ungroup() %>%
    group_by(Group, Dataset) %>%
    tukey_hsd(Intensity ~ N_detected) %>%
    filter(group2 == "3") %>%
    add_xy_position(x = "Dataset", dodge = 0.85, fun = "max") %>%
    mutate(sig_label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ))
  
  # Violin widths are scaled to count within each replicate group relative to the
  # global maximum, preserving proportionality across panels A and C.
  count_table <- df_plot %>%
    count(Group, Dataset, N_detected, name = "n_rep") %>%
    mutate(
      rel_width = n_rep / max(n_rep),
      offset    = (as.numeric(as.character(N_detected)) - 2) * 0.3
    )
  
  df_plot <- df_plot %>%
    left_join(count_table, by = c("Group", "Dataset", "N_detected")) %>%
    mutate(x_numeric = as.numeric(Dataset) + offset)
  
  ggplot(df_plot, aes(x = x_numeric, y = Intensity, fill = N_detected)) +
    geom_violin(
      aes(group = interaction(Dataset, N_detected), width = rel_width * 0.25),
      scale = "width", alpha = 0.85, trim = TRUE, linewidth = 0.20
    ) +
    stat_pvalue_manual(
      data         = stat_results %>% filter(sig_label != "ns"),
      label        = "sig_label",
      tip.length   = 0.02,
      bracket.size = 0.3,
      size         = 4,
      colour       = "grey20",
      hide.ns      = TRUE
    ) +
    facet_wrap(~ Group, ncol = 1, scales = "free_x") +
    scale_x_continuous(
      breaks = seq_along(ds_labels),
      labels = ds_labels,
      expand = expansion(mult = c(0.1, 0.1))
    ) +
    scale_fill_manual(
      values = c("1" = "#92C5DE", "2" = "#2166AC", "3" = "#053061"),
      labels = c("1" = "1", "2" = "2", "3" = "3"),
      name   = "Replicates\ndetected"
    ) +
    scale_y_continuous(
      name   = expression(Normalised~intensity~(log[10])),
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    labs(x = NULL) +
    theme_pubr() +
    theme(
      strip.background = element_blank(),
      axis.text.x      = element_text(size = 10),
      legend.position  = "none",
      plot.margin      = margin(2, 2, 2, 2, "pt"),
      panel.spacing    = unit(0.5, "lines")
    )
}

# ── Assemble figure ───────────────────────────────────────────────────────────

p_A <- make_panel_A(datasets)
p_B <- make_panel_B(datasets)
p_C <- make_panel_C(datasets)

#svg(filename="Figure2E-G.svg", width=180, height=110, units="mm")
(p_A | p_B) / p_C +
  plot_layout(heights = c(1, 2), widths = c(1, 1.5), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"),
                  legend.position = "bottom",
                  legend.direction = "horizontal")
  ) &
  theme( legend.direction = "horizontal")
#dev.off()

