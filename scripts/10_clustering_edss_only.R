# ============================================================================
# Script: 10_clustering_edss_only.R
# Descrizione: Clustering dei pazienti usando solo variabili EDSS
#              (numerico, binario, 3-class)
# Autore:  Daria Simonetti
# ============================================================================

source("scripts/01_setup.R")

if (!require(cluster, quietly = TRUE)) install.packages("cluster")
if (!require(dendextend, quietly = TRUE)) install.packages("dendextend")
if (!require(scales, quietly = TRUE)) install.packages("scales")

library(cluster)
library(dendextend)
library(scales)

set.seed(123)

df <- readRDS("outputs/data/data_complete.rds")
cat("Righe:", nrow(df), "| Colonne:", ncol(df), "\n")

if (!"row_id" %in% names(df)) {
  df <- df %>% mutate(row_id = row_number())
}

# Crea edss_bin01 se non esiste
if (!"edss_bin01" %in% names(df)) {
  df <- df %>%
    mutate(
      edss_bin01 = case_when(
        is.na(edss) ~ NA_integer_,
        edss <= 2.0 ~ 0L,
        edss > 2.0 ~ 1L
      )
    )
  cat("Creata variabile edss_bin01 (0 se <=2.0, 1 se >2.0)\n")
}

# Verifica variabili
cat("\nVariabili EDSS disponibili:\n")
cat("  edss (numerico):", sum(! is.na(df$edss)), "valori validi\n")
cat("  edss_bin01 (0/1):", sum(! is.na(df$edss_bin01)), "valori validi\n")
cat("  edss_3class (factor):", sum(!is.na(df$edss_3class)), "valori validi\n")

# Contatori file salvati
saved_csv <- character()
saved_png <- character()

# Helper: silhouette media
calc_avg_silhouette <- function(dist_mat, clusters) {
  if (length(unique(clusters)) < 2) return(NA_real_)
  sil <- silhouette(clusters, dist_mat)
  mean(sil[, "sil_width"])
}

# Helper: dendrogramma "leggibile"
# - Rami colorati per cluster
# - foglie colorate per una label interpretabile (edss_3class o edss_bin)
save_colored_dendrogram <- function(hc, k, file_path, main, sub, labels_vec = NULL) {
  dend <- as.dendrogram(hc)
  
  # Colora rami per cluster
  dend <- dendextend::color_branches(dend, k = k)
  
  # Colora le foglie per label
  if (!is.null(labels_vec)) {
    # labels_vec: vettore nominato names=row_id, values=label
    leaf_vals <- as.character(labels_vec[order.dendrogram(dend)])
    leaf_vals <- factor(leaf_vals)
    
    # palette
    pal_values <- c(project_colors[3], project_colors[2], project_colors[1], project_colors[4])
    pal <- setNames(pal_values[seq_len(nlevels(leaf_vals))], levels(leaf_vals))
    
    dend <- dendextend::set(dend, "labels_col", pal[as.character(leaf_vals)])
    dend <- dendextend::set(dend, "labels_cex", 0.6)
  } else {
    dend <- dendextend::set(dend, "labels_cex", 0.4)
  }
  
  png(file_path, width = 11, height = 7, units = "in", res = 300)
  par(mar = c(5, 4, 4, 2))
  plot(dend, main = main, sub = sub, ylab = "Altezza")
  rect.hclust(hc, k = k, border = "gray50")
  dev.off()
}

# 1. CLUSTERING SU EDSS NUMERICO (1D)

df_edss <- df %>% filter(!is.na(edss))
cat("Osservazioni con EDSS valido:", nrow(df_edss), "\n")

# Matrice per clustering (1 colonna)
X_edss <- as.matrix(df_edss$edss)
rownames(X_edss) <- df_edss$row_id

# Risultati per entrambe le distanze
results_numeric <- tibble(
  distance_type = character(),
  method = character(),
  k = integer(),
  avg_silhouette = numeric(),
  n = integer()
)

best_clusters <- list()

# --- Loop su tipi di distanza ---
distance_types <- c("manhattan", "euclidean")

for (dist_type in distance_types) {
  
  cat("\nDistanza:", dist_type, "\n")
  
  # Calcola distanza
  dist_mat <- dist(X_edss, method = dist_type)
  
  # Clustering gerarchico (average linkage)
  hc <- hclust(dist_mat, method = "average")
  
  # Valuta k in {2, 3, 4}
  k_values <- 2:4
  sil_values <- numeric(length(k_values))
  
  for (i in seq_along(k_values)) {
    k <- k_values[i]
    clusters <- cutree(hc, k = k)
    sil_values[i] <- calc_avg_silhouette(dist_mat, clusters)
    
    results_numeric <- bind_rows(results_numeric, tibble(
      distance_type = dist_type,
      method = "average",
      k = k,
      avg_silhouette = round(sil_values[i], 4),
      n = nrow(df_edss)
    ))
    
    cat("k =", k, "-> silhouette =", round(sil_values[i], 4), "\n")
  }
  
  # Seleziona k migliore
  best_k <- k_values[which.max(sil_values)]
  best_sil <- max(sil_values)
  
  cat("Migliore:  k =", best_k, "(silhouette =", round(best_sil, 4), ")\n")
  
  # Salva assegnazioni cluster per k migliore
  best_clusters[[dist_type]] <- list(
    k = best_k,
    silhouette = best_sil,
    clusters = cutree(hc, k = best_k),
    hc = hc,
    dist_mat = dist_mat
  )
  
  # --- Dendrogramma ---
  leaf_label <- df_edss$edss_3class
  names(leaf_label) <- df_edss$row_id
  
  png_path <- paste0("outputs/figures/dendrogram_edss_numeric_", dist_type, "_colored.png")
  
  save_colored_dendrogram(
    hc = hc,
    k = best_k,
    file_path = png_path,
    main = paste("Dendrogramma EDSS numerico -", str_to_title(dist_type)),
    sub = paste0(
      "Average linkage | k_best=", best_k,
      " | Rami=cluster | Foglie=EDSS 3-class"
    ),
    labels_vec = leaf_label
  )
  
  saved_png <- c(saved_png, png_path)
}

write_csv(results_numeric, "outputs/tables/edss_numeric_hclust_results.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_numeric_hclust_results.csv")
cat("\nSalvato: outputs/tables/edss_numeric_hclust_results.csv\n")

# --- Crea tabella cluster per distanza con miglior silhouette ---
# Scegli la distanza con silhouette media migliore
best_dist <- names(best_clusters)[which.max(sapply(best_clusters, function(x) x$silhouette))]

clusters_best_df <- tibble(
  row_id = df_edss$row_id,
  edss = df_edss$edss,
  edss_3class = as.character(df_edss$edss_3class),
  edss_bin = as.character(df_edss$edss_bin),
  cluster = best_clusters[[best_dist]]$clusters,
  distance_type = best_dist,
  k_best = best_clusters[[best_dist]]$k,
  silhouette_avg = round(best_clusters[[best_dist]]$silhouette, 4)
)

write_csv(clusters_best_df, "outputs/tables/edss_numeric_clusters_best.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_numeric_clusters_best.csv")
cat("outputs/tables/edss_numeric_clusters_best.csv\n")

# Riepilogo cluster (numerico), utile per spiegazione
edss_numeric_summary <- clusters_best_df %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    min = min(edss),
    q1 = quantile(edss, 0.25),
    median = median(edss),
    q3 = quantile(edss, 0.75),
    max = max(edss),
    .groups = "drop"
  )

write_csv(edss_numeric_summary, "outputs/tables/edss_numeric_cluster_summary.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_numeric_cluster_summary.csv")
cat("Salvato: outputs/tables/edss_numeric_cluster_summary.csv\n")

# --- Grafico: strip/jitter plot + boxplot ---

p_edss_clusters <- clusters_best_df %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = cluster, y = edss, fill = cluster)) +
  geom_boxplot(alpha = 0.45, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.65, size = 2) +
  scale_fill_manual(values = project_colors[1:best_clusters[[best_dist]]$k]) +
  labs(
    title = "Clustering su EDSS numerico",
    subtitle = paste0("Distanza: ", best_dist, " | average linkage | k=", best_clusters[[best_dist]]$k,
                      " | silhouette=", round(best_clusters[[best_dist]]$silhouette, 3)),
    x = "Cluster",
    y = "EDSS"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),legend.position = "none")

ggsave("outputs/figures/cluster_edss_numeric.png", plot = p_edss_clusters,
       width = 8, height = 6, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/cluster_edss_numeric.png")
cat("outputs/figures/cluster_edss_numeric.png\n")

# Grafico 2: EDSS ordinato + colore cluster
p_edss_ordered <- clusters_best_df %>%
  arrange(edss) %>%
  mutate(rank = row_number(), cluster = factor(cluster)) %>%
  ggplot(aes(x = rank, y = edss, color = cluster)) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(values = project_colors[1:best_clusters[[best_dist]]$k]) +
  labs(
    title = "EDSS ordinato con assegnazione cluster",
    subtitle = paste0("Distanza: ", best_dist, " | k=", best_clusters[[best_dist]]$k,
                      " | silhouette=", round(best_clusters[[best_dist]]$silhouette, 3)),
    x = "Pazienti (ordinati per EDSS)",
    y = "EDSS"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave("outputs/figures/edss_ordered_by_cluster.png", plot = p_edss_ordered,
       width = 9, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/edss_ordered_by_cluster.png")
cat("Salvato: outputs/figures/edss_ordered_by_cluster.png\n")

# 2. CLUSTERING SU EDSS BINARIO (0/1)

# ----------------------------------------------------------------------------
# NOTA METODOLOGICA: Simple Matching Distance
# ----------------------------------------------------------------------------
# Per variabili binarie, la distanza "Simple Matching" considera:
# - Match 0-0: distanza = 0 (entrambi assenti)
# - Match 1-1: distanza = 0 (entrambi presenti)
# - Mismatch 0-1 o 1-0: distanza = 1
#
# Su una singola variabile binaria, Simple Matching coincide con la distanza
# di Hamming:  d(a,b) = |a - b| = 0 se uguali, 1 se diversi.
#
# A differenza di Jaccard, Simple Matching valorizza anche i match 0-0,
# il che è appropriato quando "non avere la caratteristica" è informativo.
# ----------------------------------------------------------------------------

cat("\nDistanza:  Simple Matching (Hamming su 1 bit)\n")

# Filtra osservazioni valide
df_bin <- df %>% filter(!is.na(edss_bin01))
cat("   Osservazioni con edss_bin01 valido:", nrow(df_bin), "\n")

# Calcola distanza Simple Matching manualmente (su 1 variabile = |a - b|)
# Equivalente a dist() con method="manhattan" su 0/1
X_bin <- as.matrix(df_bin$edss_bin01)
rownames(X_bin) <- df_bin$row_id

dist_bin <- dist(X_bin, method = "manhattan")  # Su 0/1, manhattan = Simple Matching

# Clustering gerarchico
hc_bin <- hclust(dist_bin, method = "average")

# Taglio a k=2 (naturale per variabile binaria)
clusters_bin <- cutree(hc_bin, k = 2)

# Calcola silhouette
sil_bin <- calc_avg_silhouette(dist_bin, clusters_bin)
cat("   Silhouette media (k=2):", round(sil_bin, 4), "\n")

clusters_bin_df <- tibble(
  row_id = df_bin$row_id,
  edss_bin01 = df_bin$edss_bin01,
  edss_bin_label = if_else(df_bin$edss_bin01 == 0L, "<=2.0", ">2.0"),
  cluster = clusters_bin
)

write_csv(clusters_bin_df, "outputs/tables/edss_bin_clusters.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_bin_clusters.csv")
cat("Salvato: outputs/tables/edss_bin_clusters.csv\n")

# Crosstab
crosstab_bin <- clusters_bin_df %>%
  count(cluster, edss_bin_label, name = "n") %>%
  group_by(cluster) %>%
  mutate(perc = round(n / sum(n) * 100, 2)) %>%
  ungroup()

write_csv(crosstab_bin, "outputs/tables/edss_bin_cluster_vs_label.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_bin_cluster_vs_label.csv")
cat("Salvato: outputs/tables/edss_bin_cluster_vs_label.csv\n")

# --- Dendrogramma ---

save_colored_dendrogram(
  hc = hc_bin,
  k = 2,
  file_path = "outputs/figures/dendrogram_edss_bin01_colored.png",
  main = "Dendrogramma EDSS binario (0/1)",
  sub = "Simple Matching | average linkage | k=2 (1 variabile binaria: cluster coincide con classi)",
  labels_vec = setNames(clusters_bin_df$edss_bin_label, clusters_bin_df$row_id)
)
saved_png <- c(saved_png, "outputs/figures/dendrogram_edss_bin01_colored.png")
cat("Salvato: outputs/figures/dendrogram_edss_bin01_colored.png\n")

# Grafico percentuale (stacked 100%) per comprensione immediata
p_bin_clusters <- clusters_bin_df %>%
  mutate(cluster = factor(cluster)) %>%
  count(cluster, edss_bin_label) %>%
  group_by(cluster) %>%
  mutate(perc = n / sum(n) * 100) %>%
  ungroup() %>%
  ggplot(aes(x = cluster, y = perc, fill = edss_bin_label)) +
  geom_col(position = "fill", alpha = 0.85) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(values = c("<=2.0" = project_colors[3], ">2.0" = project_colors[4]),
                    name = "EDSS bin") +
  labs(
    title = "Cluster vs EDSS binario (composizione %)",
    subtitle = paste0("Simple Matching | k=2 | silhouette=", round(sil_bin, 3)),
    x = "Cluster",
    y = "Percentuale nel cluster"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/cluster_edss_bin_composition.png", plot = p_bin_clusters,
       width = 8, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/cluster_edss_bin_composition.png")
cat("Salvato: outputs/figures/cluster_edss_bin_composition.png\n")

# --- Tabella incrociata cluster vs edss_bin01 ---
crosstab_bin <- clusters_bin_df %>%
  count(cluster, edss_bin01, name = "n") %>%
  group_by(cluster) %>%
  mutate(perc = round(n / sum(n) * 100, 2)) %>%
  ungroup() %>%
  rename(label = edss_bin01)

write_csv(crosstab_bin, "outputs/tables/edss_bin_cluster_vs_label.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_bin_cluster_vs_label.csv")
cat("outputs/tables/edss_bin_cluster_vs_label.csv\n")

print(crosstab_bin)

# --- Grafico: barplot counts per cluster ---
cat("\nCreazione grafico cluster EDSS binario...\n")

p_bin_clusters <- clusters_bin_df %>%
  mutate(
    cluster = factor(cluster),
    edss_bin01 = factor(edss_bin01, labels = c("class0 (<=2.0)", "class1 (>2.0)"))
  ) %>%
  count(cluster, edss_bin01) %>%
  ggplot(aes(x = cluster, y = n, fill = edss_bin01)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(project_colors[3], project_colors[4]), name = "EDSS Class") +
  labs(
    title = "Clustering EDSS binario",
    subtitle = paste0("Distanza: Simple Matching | k = 2 | Silhouette = ", round(sil_bin, 3)),
    x = "Cluster",
    y = "Conteggio"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/cluster_edss_bin.png", plot = p_bin_clusters,
       width = 8, height = 6, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/cluster_edss_bin.png")
cat("outputs/figures/cluster_edss_bin. png\n")


# 4. CLUSTERING SU EDSS 3-CLASS (ordinale)

# ----------------------------------------------------------------------------
# NOTA METODOLOGICA: Distanza Manhattan su ordinale
# ----------------------------------------------------------------------------
# Per variabili ordinali (con un ordine naturale), la distanza Manhattan
# è appropriata perché rispetta l'ordine:  la distanza tra "normal" e "severe"
# è maggiore della distanza tra "normal" e "mild".
#
# Mapping: normal=1, mild=2, severe=3
# d(normal, mild) = |1-2| = 1
# d(normal, severe) = |1-3| = 2
# d(mild, severe) = |2-3| = 1
#
# Questo cattura correttamente la "distanza" clinica tra le classi.
# ----------------------------------------------------------------------------

cat("\nDistanza: Manhattan su ordinale (normal=1, mild=2, severe=3)\n")

df_3class <- df %>% filter(!is.na(edss_3class))
cat("  Osservazioni valide:", nrow(df_3class), "\n")
cat("  Livelli edss_3class:", paste(levels(df_3class$edss_3class), collapse = " | "), "\n")

df_3class <- df_3class %>%
  mutate(edss_3class_ord = as.integer(edss_3class))

X_3class <- as.matrix(df_3class$edss_3class_ord)
rownames(X_3class) <- df_3class$row_id

dist_3class <- dist(X_3class, method = "manhattan")
hc_3class <- hclust(dist_3class, method = "average")
clusters_3class <- cutree(hc_3class, k = 3)

sil_3class <- calc_avg_silhouette(dist_3class, clusters_3class)
cat("  Silhouette media (k=3):", round(sil_3class, 4), "\n")

clusters_3class_df <- tibble(
  row_id = df_3class$row_id,
  edss_3class = as.character(df_3class$edss_3class),
  edss_3class_ord = df_3class$edss_3class_ord,
  cluster = clusters_3class
)

write_csv(clusters_3class_df, "outputs/tables/edss_3class_clusters.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_3class_clusters.csv")
cat("Salvato: outputs/tables/edss_3class_clusters.csv\n")

crosstab_3class <- clusters_3class_df %>%
  count(cluster, edss_3class, name = "n") %>%
  group_by(cluster) %>%
  mutate(perc = round(n / sum(n) * 100, 2)) %>%
  ungroup()

write_csv(crosstab_3class, "outputs/tables/edss_3class_cluster_vs_label.csv")
saved_csv <- c(saved_csv, "outputs/tables/edss_3class_cluster_vs_label.csv")
cat("Salvato: outputs/tables/edss_3class_cluster_vs_label.csv\n")

# Dendrogramma colorato con foglie=edss_3class
save_colored_dendrogram(
  hc = hc_3class,
  k = 3,
  file_path = "outputs/figures/dendrogram_edss_3class_colored.png",
  main = "Dendrogramma EDSS 3-class (ordinale)",
  sub = "Manhattan ordinale | average linkage | k=3 (1 variabile discreta: cluster coincide con classi)",
  labels_vec = setNames(clusters_3class_df$edss_3class, clusters_3class_df$row_id)
)
saved_png <- c(saved_png, "outputs/figures/dendrogram_edss_3class_colored.png")
cat("Salvato: outputs/figures/dendrogram_edss_3class_colored.png\n")

# Grafico percentuale stacked 100%
p_3class_clusters <- clusters_3class_df %>%
  mutate(cluster = factor(cluster),
         edss_3class = factor(edss_3class, levels = c("normal (0-2.0)", "mild (2.5-4.0)", "severe (>4.0)"))) %>%
  count(cluster, edss_3class) %>%
  group_by(cluster) %>%
  mutate(perc = n / sum(n) * 100) %>%
  ungroup() %>%
  ggplot(aes(x = cluster, y = perc, fill = edss_3class)) +
  geom_col(position = "fill", alpha = 0.85) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(
    values = c("normal (0-2.0)" = project_colors[3],
               "mild (2.5-4.0)" = project_colors[2],
               "severe (>4.0)" = project_colors[6]),
    name = "EDSS 3-class"
  ) +
  labs(
    title = "Cluster vs EDSS 3-class (composizione %)",
    subtitle = paste0("Manhattan ordinale | k=3 | silhouette=", round(sil_3class, 3)),
    x = "Cluster",
    y = "Percentuale nel cluster"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/cluster_edss_3class_composition.png", plot = p_3class_clusters,
       width = 9, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/cluster_edss_3class_composition.png")
cat("Salvato: outputs/figures/cluster_edss_3class_composition.png\n")


# 4. RIEPILOGO FINALE

cat("\nRisultati clustering:\n")
cat("--------------------------------------------------------------------------------\n")

cat("\n1) EDSS Numerico:\n")
cat("   Distanza migliore:", best_dist, "\n")
cat("   k ottimale:", best_clusters[[best_dist]]$k, "\n")
cat("   Silhouette:", round(best_clusters[[best_dist]]$silhouette, 4), "\n")

cat("\n2) EDSS Binario (0/1):\n")
cat("   Distanza:  Simple Matching\n")
cat("   k: 2\n")
cat("   Silhouette:", round(sil_bin, 4), "\n")

cat("\n3) EDSS 3-class (ordinale):\n")
cat("   Distanza: Manhattan\n")
cat("   k: 3\n")
cat("   Silhouette:", round(sil_3class, 4), "\n")

cat("\nFile CSV salvati (", length(saved_csv), "):\n")
for (f in saved_csv) {
  cat(f, "\n")
}

cat("\nFile PNG salvati (", length(saved_png), "):\n")
for (f in saved_png) {
  cat(f, "\n")
}

cat("
Note metodologiche:
   - EDSS numerico: distanze Manhattan/Euclidea, k selezionato per silhouette.
     In più: grafico EDSS ordinato per visualizzare chiaramente la separazione.
   - EDSS binario: Simple Matching (= Hamming su 1 bit), k=2 naturale
   - EDSS 3-class: Manhattan su ordinale (rispetta l'ordine delle classi)
   - EDSS binario e 3-class: con una sola variabile discreta i cluster coincidono
     (quasi) con le classi; dendrogrammi e silhouette risultano 'perfetti' ma non
     aggiungono struttura oltre alla definizione delle classi.
   - Tutti usano clustering gerarchico con average linkage.
")
