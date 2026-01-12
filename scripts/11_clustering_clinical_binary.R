# ============================================================================
# Script: 11_clustering_clinical_binary. R
# Descrizione: Clustering dei pazienti usando feature cliniche binarie
#              con esplorazione di metriche e algoritmi multipli
# Autore:   Daria Simonetti
# ============================================================================

source("scripts/01_setup.R")

if (!require(tidyverse, quietly = TRUE)) install.packages("tidyverse")
if (!require(vegan, quietly = TRUE)) install.packages("vegan")
if (!require(cluster, quietly = TRUE)) install.packages("cluster")
if (!require(dendextend, quietly = TRUE)) install.packages("dendextend")
if (!require(scales, quietly = TRUE)) install.packages("scales")
if (!require(colorspace, quietly = TRUE)) install.packages("colorspace")

library(tidyverse)
library(vegan)
library(cluster)
library(dendextend)
library(scales)
library(colorspace)

set.seed(123)

df <- readRDS("outputs/data/data_complete.rds")

cat("Righe:", nrow(df), "| Colonne:", ncol(df), "\n")

if (!"row_id" %in% names(df)) {
  df <- df %>% mutate(row_id = row_number())
}

N <- nrow(df)

# Contatori file salvati
saved_csv <- character()
saved_png <- character()


# HELPER
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  format(round(p, 3), nsmall = 3)
}

calc_avg_silhouette <- function(dist_obj, clusters) {
  if (length(unique(clusters)) < 2) return(NA_real_)
  sil <- silhouette(clusters, dist_obj)
  mean(sil[, "sil_width"])
}

save_colored_dendrogram <- function(hc, k, file_path, main, sub) {
  dend <- as.dendrogram(hc)
  dend <- dendextend::color_branches(dend, k = k)
  
  png(file_path, width = 12, height = 8, units = "in", res = 300)
  par(mar = c(5, 4, 4, 2))
  plot(dend, main = main, sub = sub, ylab = "Altezza", xlab = "Pazienti", leaflab = "none")
  rect.hclust(hc, k = k, border = "gray50")
  dev.off()
}


# 1. DEFINIZIONE FEATURE BINARIE


# Binarie cliniche
clinical_binary <- c("pyramidal", "cerebella", "brain_stem", "sensory", 
                     "sphincters", "visual", "mental", "speech", 
                     "motor_system", "sensory_system", "coordination", 
                     "gait", "bowel_bladder", "mobility", "mental_state", 
                     "optic_disc", "nystagmus", "ocular_move", "swallowing")

# Binarie presenting
presenting_binary <- names(df)[grepl("^present_", names(df))]

# Tutte le binarie (intersezione con colonne esistenti)
all_binary <- intersect(c(clinical_binary, presenting_binary), names(df))

cat("Cliniche:", length(intersect(clinical_binary, names(df))), "\n")
cat("Presenting:", length(presenting_binary), "\n")
cat("Totale binarie:", length(all_binary), "\n")

if (length(all_binary) == 0) stop("Nessuna feature binaria trovata: controlla i nomi colonne.")


# 2. CALCOLO FREQUENZE E FILTRO FEATURE


# Conversione robusta a numerico 0/1: factor/character/numeric
df_binary <- df %>%
  select(all_of(all_binary)) %>%
  mutate(across(everything(), ~as.numeric(as.character(.x))))

# Frequenze
freq_df <- df_binary %>%
  pivot_longer(everything(), names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  summarise(
    n_1 = sum(value == 1, na.rm = TRUE),
    n_0 = sum(value == 0, na.rm = TRUE),
    n_na = sum(is.na(value)),
    .groups = "drop"
  ) %>%
  mutate(freq_1 = ifelse((n_1 + n_0) > 0, n_1 / (n_1 + n_0), NA_real_)) %>%
  arrange(desc(freq_1))

write_csv(freq_df, "outputs/tables/binary_freq_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/binary_freq_11.csv")


# --- Filtro feature rare/ubiquitarie ---
min_count <- 5
max_count <- N - 5

feature_selection <- freq_df %>%
  mutate(
    reason = case_when(
      n_1 < min_count ~ "dropped_low_freq",
      n_1 > max_count ~ "dropped_high_freq",
      TRUE ~ "kept"
    )
  )

# Salva lista completa
write_csv(feature_selection %>% select(feature, n_1, n_0, n_na, freq_1, reason),
          "outputs/tables/features_clustering_all_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/features_clustering_all_11.csv")

# Salva lista filtrata
features_filtered <- feature_selection %>% filter(reason == "kept")
write_csv(features_filtered %>% select(feature, n_1, freq_1),
          "outputs/tables/features_clustering_filtered_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/features_clustering_filtered_11.csv")

# Riepilogo
n_kept <- sum(feature_selection$reason == "kept")
n_dropped <- sum(feature_selection$reason != "kept")

cat("Feature kept:", n_kept, "\n")
cat("Feature dropped:", n_dropped, "\n")

if (n_kept < 2) {
  stop("Troppe poche feature dopo filtro (kept < 2). Riduci il filtro o controlla i dati.")
}

# Mostra feature droppate
dropped_features <- feature_selection %>% filter(reason != "kept")
if (nrow(dropped_features) > 0) {
  cat("Feature droppate:\n")
  for (i in 1:nrow(dropped_features)) {
    cat("-", dropped_features$feature[i], 
        "(n_1 =", dropped_features$n_1[i], ",", dropped_features$reason[i], ")\n")
  }
}


# 3. COSTRUZIONE MATRICE X_bin

features_to_use <- features_filtered$feature

X_bin <- df %>%
  select(all_of(features_to_use)) %>%
  mutate(across(everything(), ~as.numeric(as.character(.x)))) %>%
  as.matrix()

# ----------------------------------------------------------------------------
# GESTIONE NA:  sostituiamo con 0
# Motivazione: per distanze binarie, NA viene trattato come "assenza" della
# caratteristica.  Questo è coerente con l'interpretazione clinica dove
# un dato mancante su un indicatore spesso significa "non rilevato/presente". 
# Nel nostro dataset n_na dovrebbe essere 0, ma gestiamo per robustezza.
# ----------------------------------------------------------------------------
rownames(X_bin) <- df$row_id

n_na_total <- sum(is.na(X_bin))
if (n_na_total > 0) {
  cat("Trovati", n_na_total, "NA - sostituiti con 0\n")
  X_bin[is.na(X_bin)] <- 0
} else {
  cat("Nessun NA presente\n")
}

cat("Dimensioni X_bin:", nrow(X_bin), "x", ncol(X_bin), "\n")

# Salva dimensioni
dimensions_df <- tibble(
  N = nrow(X_bin),
  p_features = ncol(X_bin),
  min_count = min_count,
  max_count = max_count
)

write_csv(dimensions_df, "outputs/tables/X_bin_dimensions_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/X_bin_dimensions_11.csv")


# 4. CALCOLO DISTANZE

# ----------------------------------------------------------------------------
# A) DISTANZA JACCARD
# ----------------------------------------------------------------------------
# Jaccard ignora i match 0-0, considerando solo: 
# J = a11 / (a11 + a10 + a01)
# Adatta per dati sparsi dove l'assenza congiunta non è informativa. 
# ----------------------------------------------------------------------------

dist_jaccard <- vegdist(X_bin, method = "jaccard", binary = TRUE)
cat("Calcolata (range:", round(min(dist_jaccard), 3), "-", round(max(dist_jaccard), 3), ")\n")

# ----------------------------------------------------------------------------
# B) DISTANZA SIMPLE MATCHING (Sokal-Michener)
# ----------------------------------------------------------------------------
# Simple Matching considera anche i match 0-0:
# SM = 1 - (a00 + a11) / p
# Appropriata quando l'assenza congiunta di una caratteristica è informativa.
# Implementazione manuale per chiarezza.
# ----------------------------------------------------------------------------

# Funzione per calcolare distanza Simple Matching
calc_simple_matching_dist <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  dm <- matrix(0, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Conta match (0-0 e 1-1)
      matches <- sum(X[i,] == X[j,])
      # Simple Matching Distance = 1 - (matches / p)
      dm[i, j] <- 1 - (matches / p)
      dm[j, i] <- dm[i, j]
    }
  }
  
  rownames(dm) <- rownames(X)
  colnames(dm) <- rownames(X)
  
  return(as.dist(dm))
}

dist_sm <- calc_simple_matching_dist(X_bin)
cat("Simple Matching (range:", round(min(dist_sm), 3), "-", round(max(dist_sm), 3), ")\n")

# Lista distanze per iterazione
distances <- list(
  jaccard = dist_jaccard,
  simple_matching = dist_sm
)

cat("\nRange distanze:\n")
cat("  Jaccard:", round(min(dist_jaccard), 3), "-", round(max(dist_jaccard), 3), "\n")
cat("  Simple Matching:", round(min(dist_sm), 3), "-", round(max(dist_sm), 3), "\n")

# 5. ESPLORAZIONE: ALGORITMI + K (SILHOUETTE GRID)

# Range di k da esplorare
k_range <- 2:6

# Tabella risultati
silhouette_grid <- tibble(
  distance = character(),
  algorithm = character(),
  linkage = character(),
  k = integer(),
  avg_silhouette = numeric(),
  n = integer(),
  p_features = integer()
)

# Salva oggetti clustering per dopo
clustering_objects <- list()

cat("\nEsplorazione configurazioni...\n")

for (dist_name in names(distances)) {
  dist_obj <- distances[[dist_name]]
  cat("\n   Distanza:", dist_name, "\n")
  
  # --- HCLUST con average linkage ---
  cat("      hclust (average)...")
  hc_avg <- hclust(dist_obj, method = "average")
  clustering_objects[[paste0(dist_name, "_hclust_average")]] <- hc_avg
  
  for (k in k_range) {
    clusters <- cutree(hc_avg, k = k)
    silhouette_grid <- bind_rows(silhouette_grid, tibble(
      distance = dist_name,
      algorithm = "hclust",
      linkage = "average",
      k = k,
      avg_silhouette = calc_avg_silhouette(dist_obj, clusters),
      n = N,
      p_features = n_kept
    ))
  }
  
  # --- HCLUST con complete linkage ---
  cat("      hclust (complete)...")
  hc_comp <- hclust(dist_obj, method = "complete")
  clustering_objects[[paste0(dist_name, "_hclust_complete")]] <- hc_comp
  
  for (k in k_range) {
    clusters <- cutree(hc_comp, k = k)
    silhouette_grid <- bind_rows(silhouette_grid, tibble(
      distance = dist_name,
      algorithm = "hclust",
      linkage = "complete",
      k = k,
      avg_silhouette = calc_avg_silhouette(dist_obj, clusters),
      n = N,
      p_features = n_kept
    ))
  }
  
  # --- PAM (k-medoids) ---
  cat("      PAM...")
  
  for (k in k_range) {
    pam_result <- pam(dist_obj, k = k, diss = TRUE)
    
    silhouette_grid <- bind_rows(silhouette_grid, tibble(
      distance = dist_name,
      algorithm = "pam",
      linkage = NA_character_,
      k = k,
      avg_silhouette =  calc_avg_silhouette(dist_obj, pam_result$clustering),
      n = N,
      p_features = n_kept
    ))
    
    # Salva solo per k=2,3,4 per non appesantire
    if (k <= 4) {
      clustering_objects[[paste0(dist_name, "_pam_k", k)]] <- pam_result
    }
  }
}


# 6. SELEZIONE CONFIGURAZIONE MIGLIORE (con vincolo min_cluster_size)

# Prima:
# Trova configurazione con silhouette massima
# in caso di pareggio (entro tolleranza), preferisci Jaccard
# Poi preferisci k più piccolo

# delta <- 0.02  # tolleranza pratica
# 
# max_sil <- max(silhouette_grid$avg_silhouette, na.rm = TRUE)
# 
# candidates <- silhouette_grid %>%
#   filter(!is.na(avg_silhouette)) %>%
#   filter(avg_silhouette >= max_sil - delta)
# 
# # se ci sono più candidati, preferisci jaccard e k piccolo
# best_config <- candidates %>%
#   mutate(distance_pref = ifelse(distance == "jaccard", 0, 1)) %>%
#   arrange(distance_pref, k) %>%
#   slice(1)
# 
# cat("\nMIGLIORE CONFIGURAZIONE:\n")
# cat("   Distanza:", best_config$distance, "\n")
# cat("   Algoritmo:", best_config$algorithm, "\n")
# cat("   Linkage:", ifelse(is.na(best_config$linkage), "N/A", best_config$linkage), "\n")
# cat("   k:", best_config$k, "\n")
# cat("   Silhouette:", round(best_config$avg_silhouette, 4), "\n")
# 
# dist_best <- distances[[best_config$distance]]
# 
# # Estrai cluster dalla configurazione migliore
# if (best_config$algorithm == "hclust") {
#   obj_key <- paste0(best_config$distance, "_hclust_", best_config$linkage)
#   hc_best <- clustering_objects[[obj_key]]
#   clusters_best <- cutree(hc_best, k = best_config$k)
# } else {
#   obj_key <- paste0(best_config$distance, "_pam_k", best_config$k)
#   if (obj_key %in% names(clustering_objects)) {
#     pam_best <- clustering_objects[[obj_key]]
#   } else {
#     # Ricalcola se non salvato
#     pam_best <- pam(distances[[best_config$distance]], k = best_config$k, diss = TRUE)
#   }
#   clusters_best <- pam_best$clustering
# }
# 
# # Salva assegnazioni cluster
# clusters_best_df <- tibble(
#   row_id = df$row_id,
#   cluster = clusters_best,
#   distance = best_config$distance,
#   algorithm = best_config$algorithm,
#   linkage = ifelse(is.na(best_config$linkage), "N/A", best_config$linkage),
#   k_best = best_config$k,
#   avg_silhouette_best = round(best_config$avg_silhouette, 4)
# )
# 
# write_csv(clusters_best_df, "outputs/tables/clinical_clusters_best_11.csv")
# saved_csv <- c(saved_csv, "outputs/tables/clinical_clusters_best_11.csv")
# 
# # Aggiungi cluster al dataframe principale per profiling
# df <- df %>%
#   left_join(clusters_best_df %>% select(row_id, cluster), by = "row_id")


# Dopo:
# imponiamo min_cluster_size e scartiamo soluzioni con micro-cluster

min_cluster_size <- 5
delta <- 0.02

get_clusters_for_config <- function(distance, algorithm, linkage, k) {
  dist_obj <- distances[[distance]]
  
  if (algorithm == "hclust") {
    hc <- clustering_objects[[paste0(distance, "_hclust_", linkage)]]
    cutree(hc, k = k)
  } else {
    # ricalcolo PAM qui: evita dipendenze da oggetti salvati solo per k<=4
    pam(dist_obj, k = k, diss = TRUE)$clustering
  }
}

# Calcola min cluster size per ogni configurazione
silhouette_grid <- silhouette_grid %>%
  rowwise() %>%
  mutate(
    linkage2 = ifelse(is.na(linkage), NA_character_, linkage),
    min_cluster_n = {
      cl <- get_clusters_for_config(distance, algorithm, linkage2, k)
      as.integer(min(table(cl)))
    }
  ) %>%
  ungroup()

# Salva griglia silhouette completa (con min_cluster_n)
silhouette_grid_out <- silhouette_grid %>%
  mutate(avg_silhouette = round(avg_silhouette, 4))

write_csv(silhouette_grid_out, "outputs/tables/silhouette_grid_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/silhouette_grid_11.csv")

# Tieni solo configurazioni "valide" (cluster non troppo piccoli)
valid_grid <- silhouette_grid %>%
  filter(min_cluster_n >= min_cluster_size)

if (nrow(valid_grid) == 0) {
  stop("Nessuna configurazione produce cluster con dimensione minima >= min_cluster_size. Riduci min_cluster_size.")
}

max_sil_valid <- max(valid_grid$avg_silhouette, na.rm = TRUE)

candidates <- valid_grid %>%
  filter(avg_silhouette >= max_sil_valid - delta)

best_config <- candidates %>%
  mutate(distance_pref = ifelse(distance == "jaccard", 0, 1)) %>%
  arrange(desc(avg_silhouette), distance_pref, k) %>%  # silhouette prima; poi preferenza distanza; poi k piccolo
  slice(1)

cat("\nMIGLIORE CONFIGURAZIONE (con vincolo min_cluster_size):\n")
cat("   Distanza:", best_config$distance, "\n")
cat("   Algoritmo:", best_config$algorithm, "\n")
cat("   Linkage:", ifelse(is.na(best_config$linkage), "N/A", best_config$linkage), "\n")
cat("   k:", best_config$k, "\n")
cat("   Silhouette:", round(best_config$avg_silhouette, 4), "\n")
cat("   Min cluster n:", best_config$min_cluster_n, "\n")

dist_best <- distances[[best_config$distance]]

# Estrai cluster dalla configurazione migliore
if (best_config$algorithm == "hclust") {
  obj_key <- paste0(best_config$distance, "_hclust_", best_config$linkage)
  hc_best <- clustering_objects[[obj_key]]
  clusters_best <- cutree(hc_best, k = best_config$k)
} else {
  pam_best <- pam(dist_best, k = best_config$k, diss = TRUE)
  clusters_best <- pam_best$clustering
}

# Salva assegnazioni cluster
clusters_best_df <- tibble(
  row_id = df$row_id,
  cluster = as.integer(clusters_best),
  distance = best_config$distance,
  algorithm = best_config$algorithm,
  linkage = ifelse(is.na(best_config$linkage), "N/A", best_config$linkage),
  k_best = best_config$k,
  avg_silhouette_best = round(best_config$avg_silhouette, 4),
  min_cluster_n = best_config$min_cluster_n
)

write_csv(clusters_best_df, "outputs/tables/clinical_clusters_best_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/clinical_clusters_best_11.csv")

df <- df %>%
  left_join(clusters_best_df %>% select(row_id, cluster), by = "row_id")

cluster_sizes <- table(df$cluster)

# 7. GRAFICI

# INTERPRETAZIONE: frequenza feature per cluster (TOP DIFFERENZE)
# frequenza per feature e cluster
X_long <- as_tibble(X_bin, rownames = "row_id") %>%
  mutate(row_id = as.integer(row_id)) %>%
  left_join(df %>% select(row_id, cluster), by = "row_id") %>%
  pivot_longer(cols = all_of(features_to_use), names_to = "feature", values_to = "value") %>%
  mutate(value = as.integer(value))

feat_freq <- X_long %>%
  group_by(cluster, feature) %>%
  summarise(
    n = n(),
    n_1 = sum(value == 1, na.rm = TRUE),
    freq = n_1 / n,
    .groups = "drop"
  )

# Calcola differenza max-min tra cluster per ogni feature (misura di "discriminazione")
feat_diff <- feat_freq %>%
  group_by(feature) %>%
  summarise(
    freq_min = min(freq),
    freq_max = max(freq),
    delta_freq = freq_max - freq_min,
    .groups = "drop"
  ) %>%
  arrange(desc(delta_freq))

# Salva tabella completa
feat_freq_out <- feat_freq %>%
  left_join(feat_diff, by = "feature") %>%
  mutate(
    freq = round(freq, 3),
    delta_freq = round(delta_freq, 3)
  ) %>%
  arrange(desc(delta_freq), feature, cluster)

write_csv(feat_freq_out, "outputs/tables/cluster_feature_freq_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/cluster_feature_freq_11.csv")

# Grafico: top 10 feature per delta freq
top_n <- min(10, nrow(feat_diff))
top_features <- feat_diff %>%
  arrange(desc(delta_freq)) %>%
  slice_head(n = top_n) %>%
  pull(feature)

p_top_features <- feat_freq %>%
  filter(feature %in% top_features) %>%
  mutate(
    cluster = factor(cluster),
    feature = factor(feature, levels = top_features)
  ) %>%
  ggplot(aes(x = feature, y = freq, fill = cluster)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Top feature che differenziano i cluster",
    subtitle = paste0(
      "Top ", top_n, " per differenza (max-min).  Size cluster: ",
      paste(names(cluster_sizes), cluster_sizes, sep="=", collapse=" | ")
    ),
    x = "Feature (binaria)",
    y = "Freq (1) nel cluster",
    fill = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("outputs/figures/top_features_by_cluster_11.png", plot = p_top_features,
       width = 11, height = 6.5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/top_features_by_cluster_11.png")

# --- 7.1) Silhouette vs k ---

# 7.1 Silhouette vs k (facet per distanza; colore per algoritmo; linetype per linkage)
sil_plot_df <- silhouette_grid %>%
  mutate(
    distance = factor(distance, levels = c("jaccard", "simple_matching")),
    algorithm = factor(algorithm, levels = c("hclust", "pam")),
    linkage = ifelse(is.na(linkage), "none", linkage),
    linkage = factor(linkage, levels = c("none", "average", "complete"))
  )

p_silhouette <- ggplot(sil_plot_df, aes(x = k, y = avg_silhouette, color = algorithm, linetype = linkage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~distance, ncol = 2) +
  geom_vline(xintercept = best_config$k, linetype = "dashed", color = "gray40", alpha = 0.7) +
  scale_x_continuous(breaks = k_range) +
  labs(
    title = "Silhouette media vs numero di cluster (k)",
    subtitle = paste0(
      "Best: ", best_config$distance, " + ", best_config$algorithm,
      ifelse(is.na(best_config$linkage), "", paste0(" + ", best_config$linkage)),
      " | k=", best_config$k,
      " | silhouette=", round(best_config$avg_silhouette, 3),
      " | min cluster n=", best_config$min_cluster_n
    ),
    x = "k",
    y = "Silhouette media",
    color = "Algoritmo",
    linetype = "Linkage"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("outputs/figures/silhouette_by_k_11.png", plot = p_silhouette,
       width = 11, height = 6.5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/silhouette_by_k_11.png")

# 7.2 PCoA (cmdscale) con centroidi cluster
mds <- cmdscale(dist_best, k = 2, eig = TRUE)
eig <- mds$eig
var_expl <- round(eig[1:2] / sum(eig[eig > 0]) * 100, 1)

pcoa_df <- tibble(
  row_id = df$row_id,
  Dim1 = mds$points[, 1],
  Dim2 = mds$points[, 2],
  cluster = factor(df$cluster),
  binary_sum = df$binary_sum,
  edss = df$edss,
  edss_3class = df$edss_3class
)

centroids <- pcoa_df %>% group_by(cluster) %>% summarise(Dim1 = mean(Dim1), Dim2 = mean(Dim2), .groups = "drop")

p_pcoa <- ggplot(pcoa_df, aes(Dim1, Dim2)) +
  geom_point(aes(color = cluster, size = binary_sum), alpha = 0.78) +
  geom_point(data = centroids, aes(Dim1, Dim2, color = cluster),
             shape = 4, stroke = 1.4, size = 4, inherit.aes = FALSE) +
  scale_color_manual(values = project_colors[1:best_config$k], name = "Cluster") +
  scale_size_continuous(range = c(2, 6), name = "Binary sum") +
  labs(
    title = "PCoA su distanza tra pazienti (feature cliniche binarie)",
    subtitle = paste0(best_config$distance, " + ", best_config$algorithm,
                      " | k=", best_config$k,
                      " | Var: ", var_expl[1], "% + ", var_expl[2], "%"),
    x = paste0("Dim 1 (", var_expl[1], "%)"),
    y = paste0("Dim 2 (", var_expl[2], "%)")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/pcoa_best_clusters_11.png", plot = p_pcoa,
       width = 10, height = 7, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/pcoa_best_clusters_11.png")

# 7.3 Dendrogramma best (solo se hclust)
if (best_config$algorithm == "hclust") {
  save_colored_dendrogram(
    hc = hc_best,
    k = best_config$k,
    file_path = "outputs/figures/dendrogram_best_11.png",
    main = paste0("Dendrogramma best (", best_config$distance, " + ", best_config$linkage, ")"),
    sub = paste0("Rami colorati=cluster | k=", best_config$k, " | silhouette=", round(best_config$avg_silhouette, 3))
  )
  saved_png <- c(saved_png, "outputs/figures/dendrogram_best_11.png")
} else {
  cat("\nDendrogramma skipped: best config è PAM.\n")
}

# 7.4 EDSS numerico per cluster
p_edss_cluster <- df %>%
  filter(!is.na(edss), !is.na(cluster)) %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(cluster, edss, fill = cluster)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.55, size = 2) +
  scale_fill_manual(values = project_colors[1:best_config$k]) +
  labs(
    title = "EDSS per cluster (validazione esterna)",
    subtitle = paste0("Clustering su feature binarie | k=", best_config$k),
    x = "Cluster",
    y = "EDSS"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave("outputs/figures/edss_by_cluster_11.png", plot = p_edss_cluster,
       width = 8, height = 6, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/edss_by_cluster_11.png")

# 7.5 Composizione edss_3class (100% stacked)
if ("edss_3class" %in% names(df)) {
  comp_3 <- df %>%
    filter(!is.na(edss_3class), !is.na(cluster)) %>%
    mutate(
      cluster = factor(cluster),
      edss_3class = factor(edss_3class, levels = c("normal (0-2.0)", "mild (2.5-4.0)", "severe (>4.0)"))
    ) %>%
    count(cluster, edss_3class) %>%
    group_by(cluster) %>%
    mutate(perc = n / sum(n) * 100) %>%
    ungroup()
  
  p_comp_3 <- ggplot(comp_3, aes(cluster, perc, fill = edss_3class)) +
    geom_col(position = "fill", alpha = 0.85) +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    scale_fill_manual(
      values = c("normal (0-2.0)" = project_colors[3],
                 "mild (2.5-4.0)" = project_colors[2],
                 "severe (>4.0)" = project_colors[6]),
      name = "EDSS 3-class"
    ) +
    labs(
      title = "Composizione EDSS 3-class per cluster",
      subtitle = "Percentuali all'interno di ciascun cluster (100% stacked)",
      x = "Cluster",
      y = "Percentuale nel cluster"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("outputs/figures/edss3class_composition_by_cluster_11.png", plot = p_comp_3,
         width = 9, height = 5.5, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/edss3class_composition_by_cluster_11.png")
}

# 7.6 Composizione edss_bin (100% stacked)
if ("edss_bin" %in% names(df)) {
  comp_bin <- df %>%
    filter(!is.na(edss_bin), !is.na(cluster)) %>%
    mutate(cluster = factor(cluster)) %>%
    count(cluster, edss_bin) %>%
    group_by(cluster) %>%
    mutate(perc = n / sum(n) * 100) %>%
    ungroup()
  
  p_comp_bin <- ggplot(comp_bin, aes(cluster, perc, fill = edss_bin)) +
    geom_col(position = "fill", alpha = 0.85) +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    scale_fill_manual(
      values = c("class0 (<=2.0)" = project_colors[3],
                 "class1 (>2.0)" = project_colors[4]),
      name = "EDSS bin"
    ) +
    labs(
      title = "Composizione EDSS bin per cluster",
      subtitle = "Percentuali all'interno di ciascun cluster (100% stacked)",
      x = "Cluster",
      y = "Percentuale nel cluster"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("outputs/figures/edssbin_composition_by_cluster_11.png", plot = p_comp_bin,
         width = 8.5, height = 5.5, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/edssbin_composition_by_cluster_11.png")
}

# 7.7 binary_sum per cluster
if ("binary_sum" %in% names(df)) {
  p_binarysum <- df %>%
    filter(!is.na(binary_sum), !is.na(cluster)) %>%
    mutate(cluster = factor(cluster)) %>%
    ggplot(aes(cluster, binary_sum, fill = cluster)) +
    geom_boxplot(alpha = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.55, size = 2) +
    scale_fill_manual(values = project_colors[1:best_config$k]) +
    labs(
      title = "Binary sum per cluster",
      subtitle = "Somma indicatori clinici positivi (profiling cluster)",
      x = "Cluster",
      y = "Binary sum"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  
  ggsave("outputs/figures/binarysum_by_cluster_11.png", plot = p_binarysum,
         width = 8, height = 6, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/binarysum_by_cluster_11.png")
}

# 8. VALIDAZIONE ESTERNA:  TEST STATISTICI
validation_tests <- tibble(
  outcome = character(),
  test = character(),
  statistic = numeric(),
  p_value = numeric(),
  p_value_label = character(),
  note = character()
)

# 8.1 EDSS numerico ~ cluster (Kruskal)
if ("edss" %in% names(df)) {
  df_test <- df %>% filter(!is.na(edss), !is.na(cluster))
  if (length(unique(df_test$cluster)) >= 2) {
    kw <- kruskal.test(edss ~ factor(cluster), data = df_test)
    validation_tests <- bind_rows(validation_tests, tibble(
      outcome = "edss",
      test = "Kruskal-Wallis",
      statistic = unname(kw$statistic),
      p_value = kw$p.value,
      p_value_label = format_p(kw$p.value),
      note = paste0("df=", kw$parameter)
    ))
  }
}

# helper: chi-square vs fisher (simulated)
test_contingency <- function(tbl_mat) {
  # tbl_mat: matrix
  chi <- suppressWarnings(chisq.test(tbl_mat))
  if (any(chi$expected < 5)) {
    ft <- fisher.test(tbl_mat, simulate.p.value = TRUE, B = 10000)
    list(test = "Fisher (simulated)", statistic = NA_real_, p = ft$p.value, note = "expected<5")
  } else {
    list(test = "Chi-squared", statistic = unname(chi$statistic), p = chi$p.value, note = "")
  }
}

# 8.2 cluster ~ edss_bin
if ("edss_bin" %in% names(df)) {
  tab_bin <- df %>%
    filter(!is.na(cluster), !is.na(edss_bin)) %>%
    count(cluster, edss_bin) %>%
    pivot_wider(names_from = edss_bin, values_from = n, values_fill = 0)
  
  write_csv(tab_bin, "outputs/tables/cluster_vs_edss_bin_11.csv")
  saved_csv <- c(saved_csv, "outputs/tables/cluster_vs_edss_bin_11.csv")
  
  mat_bin <- as.matrix(tab_bin %>% select(-cluster))
  res <- test_contingency(mat_bin)
  
  validation_tests <- bind_rows(validation_tests, tibble(
    outcome = "edss_bin",
    test = res$test,
    statistic = res$statistic,
    p_value = res$p,
    p_value_label = format_p(res$p),
    note = res$note
  ))
}

# 8.3 cluster ~ edss_3class
if ("edss_3class" %in% names(df)) {
  tab_3 <- df %>%
    filter(!is.na(cluster), !is.na(edss_3class)) %>%
    count(cluster, edss_3class) %>%
    pivot_wider(names_from = edss_3class, values_from = n, values_fill = 0)
  
  write_csv(tab_3, "outputs/tables/cluster_vs_edss_3class_11.csv")
  saved_csv <- c(saved_csv, "outputs/tables/cluster_vs_edss_3class_11.csv")
  
  mat_3 <- as.matrix(tab_3 %>% select(-cluster))
  res <- test_contingency(mat_3)
  
  validation_tests <- bind_rows(validation_tests, tibble(
    outcome = "edss_3class",
    test = res$test,
    statistic = res$statistic,
    p_value = res$p,
    p_value_label = format_p(res$p),
    note = res$note
  ))
}

# salva tabella test
validation_tests_out <- validation_tests %>%
  mutate(
    statistic = ifelse(is.na(statistic), NA_real_, round(statistic, 3)),
    p_value = signif(p_value, 3)
  )

write_csv(validation_tests_out, "outputs/tables/cluster_validation_tests_11.csv")
saved_csv <- c(saved_csv, "outputs/tables/cluster_validation_tests_11.csv")

cat("\nValidazione esterna (test):\n")
print(validation_tests_out)


# 9. RIEPILOGO FINALE


cat("\nDataset:\n")
cat("   N pazienti:", N, "\n")
cat("   Feature binarie totali:", length(all_binary), "\n")
cat("   Feature dopo filtro:", n_kept, "\n")

cat("\nMIGLIORE CONFIGURAZIONE:\n")
cat("   Distanza:", best_config$distance, "\n")
cat("   Algoritmo:", best_config$algorithm, "\n")
cat("   Linkage:", ifelse(is.na(best_config$linkage), "N/A", best_config$linkage), "\n")
cat("   k:", best_config$k, "\n")
cat("   Silhouette:", best_config$avg_silhouette, "\n")

cat("\nDistribuzione cluster:\n")
print(table(clusters_best))

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
   - Jaccard: ignora match 0-0 (adatta a feature binarie sparse / presenza-assenza)
   - Simple Matching:  considera anche match 0-0, appropriata quando
     l'assenza congiunta è informativa
   - Confrontati hclust (average/complete) e PAM (k-medoids)
   - “k scelto massimizzando silhouette tra le soluzioni con cluster size minima ≥ 5; a prestazioni simili, preferita Jaccard.”
   - Validazione esterna: differenze tra cluster su EDSS (numerico, binario, 3-class) e binary_sum.

")


# # best hclust su simple_matching average
# hc_sm_avg <- clustering_objects[["simple_matching_hclust_average"]]
# 
# table(cutree(hc_sm_avg, k = 2))
# table(cutree(hc_sm_avg, k = 3))
# table(cutree(hc_sm_avg, k = 4))
# 
# # anche per jaccard average (per confronto)
# hc_j_avg <- clustering_objects[["jaccard_hclust_average"]]
# table(cutree(hc_j_avg, k = 2))
# table(cutree(hc_j_avg, k = 3))
# table(cutree(hc_j_avg, k = 4))
# 
# 
# silhouette_grid %>%
#   filter(algorithm == "pam") %>%
#   arrange(desc(avg_silhouette)) %>%
#   head(10)
# pam_best <- pam(dist_sm, k = 4, diss = TRUE)
# table(pam_best$clustering)
# 
# pam_best <- pam(dist_jaccard, k = 4, diss = TRUE)
# table(pam_best$clustering)
