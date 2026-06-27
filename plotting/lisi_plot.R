source("functions/prepare_env.R")

data_dir <- "../../data/BrainOmicsData/integration/"
res_dir <- check_dir("results/lisi")
fig_dir <- check_dir("figures/lisi")

if (!file.exists(file.path(res_dir, "lisi_results.rds"))) {
  lisi_data <- readRDS(file.path(data_dir, "lisi_data.rds"))
  raw <- compute_lisi(
    X = lisi_data[["umap_raw"]],
    meta_data = lisi_data[["meta_data"]],
    label_colnames = "Dataset"
  )
  harmony_lisi <- compute_lisi(
    X = lisi_data[["umap_harmony"]],
    meta_data = lisi_data[["meta_data"]],
    label_colnames = "Dataset"
  )
  rpca_lisi <- compute_lisi(
    X = lisi_data[["umap_rpca"]],
    meta_data = lisi_data[["meta_data"]],
    label_colnames = "Dataset"
  )
  lisi_results <- cbind(
    raw,
    harmony_lisi,
    rpca_lisi
  )
  names(lisi_results) <- c("Raw", "Harmony", "RPCA")

  saveRDS(lisi_results, file.path(res_dir, "lisi_results.rds"))
} else {
  lisi_results <- readRDS(file.path(res_dir, "lisi_results.rds"))
}

lisi_long <- tidyr::gather(
  lisi_results,
  key = "Method",
  value = "LISI"
)
lisi_long$Method <- factor(
  lisi_long$Method,
  levels = c("Raw", "Harmony", "RPCA")
)

p <- ggplot(lisi_long, aes(x = Method, y = LISI)) +
  geom_boxplot(
    aes(fill = Method),
    width = 0.5,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c("Raw" = "#E41A1C", "Harmony" = "#4DAF4A", "RPCA" = "#377EB8")
  ) +
  labs(
    x = "",
    y = "LISI"
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Harmony", "RPCA"),
      c("Raw", "RPCA"),
      c("Raw", "Harmony")
    ),
    label = "p.signif",
    step.increase = -0.12,
    vjust = -0.22
  ) +
  ylim(0, max(lisi_long$LISI) * 1.2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(
  file.path(fig_dir, "lisi_plot_3.pdf"),
  p,
  width = 2,
  height = 3
)

lisi_long_2 <- lisi_long[lisi_long$Method %in% c("Raw", "RPCA"), ]
lisi_long_2$Method <- factor(
  lisi_long_2$Method,
  levels = c("Raw", "RPCA")
)

p2 <- ggplot(lisi_long_2, aes(x = Method, y = LISI)) +
  geom_boxplot(
    aes(fill = Method),
    width = 0.5,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c("Raw" = "#E41A1C", "RPCA" = "#377EB8")
  ) +
  labs(
    x = "",
    y = "LISI"
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Raw", "RPCA")
    ),
    label = "p.signif",
    step.increase = -0.12,
    vjust = -0.22
  ) +
  ylim(0, max(lisi_long_2$LISI) * 1.2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(
  file.path(fig_dir, "lisi_plot_2.pdf"),
  p2,
  width = 1.5,
  height = 2
)
