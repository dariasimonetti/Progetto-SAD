# ============================================================================
# Script: 08_eda_multivariate_precluster. R
# Descrizione: EDA Multivariata pre-clustering - Preparazione e visualizzazione
#              struttura dati per clustering Jaccard
# Autore: Daria-Simonetti
# ============================================================================

source("scripts/01_setup.R")

# Carica pacchetti aggiuntivi
if (!require(vegan, quietly = TRUE)) install.packages("vegan")
library(vegan)

# Riproducibilità
set.seed(123)


# Carica dati
df <- readRDS("outputs/data/data_complete.rds")
cat("   Righe:", nrow(df), "| Colonne:", ncol(df), "\n")

# Aggiungi row_id per identificare i soggetti
df <- df %>%
  mutate(row_id = row_number())

N <- nrow(df)

# Contatori file salvati
saved_csv <- character()
saved_png <- character()


# DEFINIZIONE FEATURE SET PER CLUSTERING (JACCARD)

# --- feature binarie ---

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

cat("   Cliniche:", length(intersect(clinical_binary, names(df))), "\n")
cat("   Presenting:", length(presenting_binary), "\n")
cat("   Totale binarie:", length(all_binary), "\n")

# --- Calcolo frequenze ---

freq_df <- df %>%
  select(all_of(all_binary)) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  summarise(
    n_1 = sum(value == 1, na.rm = TRUE),
    n_0 = sum(value == 0, na.rm = TRUE),
    n_na = sum(is.na(value)),
    .groups = "drop"
  ) %>%
  mutate(
    freq_1 = round(n_1 / (n_1 + n_0), 4)
  ) %>%
  arrange(desc(freq_1))

write_csv(freq_df, "outputs/tables/binary_freq_08.csv")
saved_csv <- c(saved_csv, "outputs/tables/binary_freq_08.csv")

# --- Filtro feature rare/ubiquitarie ---

# Soglie:  n_1 >= 5 e n_1 <= (N-5)
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

# lista completa
features_all <- feature_selection %>%
  mutate(set = "all") %>%
  select(feature, set, reason)

write_csv(features_all, "outputs/tables/features_jaccard_all.csv")
saved_csv <- c(saved_csv, "outputs/tables/features_jaccard_all.csv")

# lista filtrata
features_filtered <- feature_selection %>%
  filter(reason == "kept") %>%
  mutate(set = "filtered") %>%
  select(feature, set, reason)

write_csv(features_filtered, "outputs/tables/features_jaccard_filtered.csv")
saved_csv <- c(saved_csv, "outputs/tables/features_jaccard_filtered.csv")


# riepilogo filtro
n_kept <- sum(feature_selection$reason == "kept")
n_dropped_low <- sum(feature_selection$reason == "dropped_low_freq")
n_dropped_high <- sum(feature_selection$reason == "dropped_high_freq")

cat("\nRiepilogo filtro:\n")
cat("   - Kept:", n_kept, "\n")
cat("   - Dropped (low freq, n_1 <", min_count, "):", n_dropped_low, "\n")
cat("   - Dropped (high freq, n_1 >", max_count, "):", n_dropped_high, "\n")

# Elenco feature droppate
dropped_features <- feature_selection %>%
  filter(reason != "kept")

if (nrow(dropped_features) > 0) {
  cat("\n   Feature droppate:\n")
  for (i in 1:nrow(dropped_features)) {
    cat("   -", dropped_features$feature[i], 
        "(n_1 =", dropped_features$n_1[i], ",", dropped_features$reason[i], ")\n")
  }
}

# --- Creazione X_bin (matrice feature filtrate) ---

features_to_use <- features_filtered$feature

X_bin <- df %>%
  select(all_of(features_to_use)) %>%
  mutate(across(everything(), ~as.numeric(as.character(.x))))

# Gestione NA
n_na_total <- sum(is.na(X_bin))
if (n_na_total > 0) {
  cat("Trovati", n_na_total, "NA - sostituiti con 0 per distanza Jaccard\n")
  X_bin[is.na(X_bin)] <- 0
} else {
  cat("Nessun NA presente\n")
}

cat("Dimensioni X_bin:", nrow(X_bin), "x", ncol(X_bin), "\n")


# B) DISTANZA JACCARD TRA PAZIENTI

# Calcola distanza Jaccard con vegan
jaccard_dist <- vegdist(X_bin, method = "jaccard", binary = TRUE)

# matrice
jaccard_mat <- as.matrix(jaccard_dist)

# Aggiungi row_id come nomi
rownames(jaccard_mat) <- df$row_id
colnames(jaccard_mat) <- df$row_id

# Salva matrice wide
jaccard_df <- as.data.frame(jaccard_mat) %>%
  mutate(row_id = rownames(jaccard_mat)) %>%
  select(row_id, everything())

write_csv(jaccard_df, "outputs/tables/jaccard_distance_matrix.csv")
saved_csv <- c(saved_csv, "outputs/tables/jaccard_distance_matrix.csv")

cat("Range distanze:", round(min(jaccard_dist), 3), "-", 
    round(max(jaccard_dist), 3), "\n")
cat("Media distanza:", round(mean(jaccard_dist), 3), "\n")


# C) HEATMAP DISTANZA JACCARD

# Clustering gerarchico per ordinare i soggetti
hc <- hclust(as.dist(jaccard_mat), method = "average")
order_idx <- hc$order

# Matrice ordibata
jaccard_mat_ordered <- jaccard_mat[order_idx, order_idx]

# Trasforma in long per ggplot
jaccard_long <- jaccard_mat_ordered %>%
  as.data.frame() %>%
  mutate(patient_a = factor(rownames(.), levels = rownames(. ))) %>%
  pivot_longer(
    -patient_a, 
    names_to = "patient_b", 
    values_to = "distance"
  ) %>%
  mutate(patient_b = factor(patient_b, levels = colnames(jaccard_mat_ordered)))

# Heatmap
p_heatmap <- ggplot(jaccard_long, aes(x = patient_a, y = patient_b, fill = distance)) +
  geom_tile() +
  scale_fill_gradient(
    low = project_colors[1], 
    high = project_colors[6],
    limits = c(0, 1),
    name = "Distanza\nJaccard"
  ) +
  labs(
    title = "Matrice delle distanze di Jaccard (pazienti)",
    subtitle = "Ordinati per raggruppamento gerarchico (collegamento medio)",
    x = "Paziente",
    y = "Paziente"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("outputs/figures/heatmap_jaccard_patients.png", plot = p_heatmap,
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/heatmap_jaccard_patients.png")


# ORDINATION (MDS/PCoA) SU DISTANZA JACCARD

# PCoA con cmdscale
mds_result <- cmdscale(as.dist(jaccard_mat), k = 2, eig = TRUE)

# Coordinate
mds_coords <- as.data.frame(mds_result$points)
names(mds_coords) <- c("Dim1", "Dim2")

# Row_id e variabili di profiling
mds_coords <- mds_coords %>%
  mutate(
    row_id = df$row_id,
    edss = df$edss,
    binary_sum = df$binary_sum
  ) %>%
  select(row_id, Dim1, Dim2, edss, binary_sum)

# Varianza spiegata
eig_values <- mds_result$eig
var_explained <- round(eig_values[1:2] / sum(eig_values[eig_values > 0]) * 100, 1)

cat("Varianza spiegata Dim1:", var_explained[1], "%\n")
cat("Varianza spiegata Dim2:", var_explained[2], "%\n")

write_csv(mds_coords, "outputs/tables/jaccard_cmdscale_2d.csv")
saved_csv <- c(saved_csv, "outputs/tables/jaccard_cmdscale_2d.csv")

# --- Grafico scatter colorato per EDSS ---

# Gruppi EDSS 
mds_coords <- mds_coords %>%
  mutate(
    edss_group = cut(edss, 
                     breaks = c(-Inf, 2, 4, Inf), 
                     labels = c("Low (0-2)", "Mid (2.5-4)", "High (4.5+)"))
  )

# Scatter con gradiente EDSS
p_mds_edss <- ggplot(mds_coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = edss, size = binary_sum), alpha = 0.7) +
  scale_color_gradient(
    low = project_colors[3], 
    high = project_colors[2],
    name = "EDSS"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "Binary Sum"
  ) +
  labs(
    title = "PCoA sulla distanza di Jaccard (caratteristiche binarie)",
    subtitle = paste0("Dim1: ", var_explained[1], "% | Dim2: ", var_explained[2], "% varianza spiegata"),
    x = paste0("Dimension 1 (", var_explained[1], "%)"),
    y = paste0("Dimension 2 (", var_explained[2], "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("outputs/figures/cmdscale_jaccard_edss.png", plot = p_mds_edss,
       width = 10, height = 7, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/cmdscale_jaccard_edss.png")

# Versione alternativa con gruppi EDSS discreti
p_mds_groups <- ggplot(mds_coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = edss_group, size = binary_sum), alpha = 0.7) +
  scale_color_manual(
    values = c("Low (0-2)" = project_colors[3], 
               "Mid (2.5-4)" = project_colors[2], 
               "High (4.5+)" = project_colors[6]),
    name = "EDSS Group"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "Binary Sum"
  ) +
  labs(
    title = "PCoA sulla distanza di Jaccard (caratteristiche binarie)",
    subtitle = paste0("Colorato per gruppi EDSS | Dim1: ", var_explained[1], 
                      "% | Dim2: ", var_explained[2], "%"),
    x = paste0("Dimension 1 (", var_explained[1], "%)"),
    y = paste0("Dimension 2 (", var_explained[2], "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("outputs/figures/cmdscale_jaccard_edss_groups.png", plot = p_mds_groups,
       width = 10, height = 7, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/cmdscale_jaccard_edss_groups.png")


#PROFILING DESCRITTIVO LEGGERO


# --- Creazione edss_bin per profiling ---


df <- df %>%
  mutate(
    edss_bin = cut(edss, 
                   breaks = c(-Inf, 2, 4, Inf), 
                   labels = c("low", "mid", "high"))
  )

cat("   Distribuzione edss_bin:\n")
print(table(df$edss_bin, useNA = "ifany"))

# --- Boxplot binary_sum ~ edss_bin ---


p_binary_edss <- df %>%
  filter(! is.na(edss_bin) & ! is.na(binary_sum)) %>%
  ggplot(aes(x = edss_bin, y = binary_sum, fill = edss_bin)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(
    values = c("low" = project_colors[3], "mid" = project_colors[2], "high" = project_colors[6])
  ) +
  labs(
    title = "Binary Sum per gruppi EDSS",
    subtitle = "Low (0-2), Mid (2.5-4), High (4.5+)",
    x = "EDSS Group",
    y = "Binary Sum (numero di indicatori positivi)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/figures/binary_sum_by_edss_bin.png", plot = p_binary_edss,
       width = 8, height = 6, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/binary_sum_by_edss_bin.png")

# --- Prevalenza top 5 feature per edss_bin ---

# Top 5 feature più associate a EDSS (dalla fase 6 / letteratura clinica)
top5_features <- c("pyramidal", "motor_system", "coordination", "sphincters", "gait")

# Verifica che esistono nel dataset
top5_features <- intersect(top5_features, names(df))

if (length(top5_features) >= 3) {
  
  # Calcola frequenza per edss_bin
  freq_by_edss <- df %>%
    filter(! is.na(edss_bin)) %>%
    select(edss_bin, all_of(top5_features)) %>%
    pivot_longer(-edss_bin, names_to = "feature", values_to = "value") %>%
    group_by(edss_bin, feature) %>%
    summarise(
      n_1 = sum(value == 1, na.rm = TRUE),
      n_total = n(),
      freq = n_1 / n_total * 100,
      .groups = "drop"
    )
  
  # Barplot
  p_top5 <- ggplot(freq_by_edss, 
                   aes(x = feature, y = freq, fill = edss_bin)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_text(
      aes(label = paste0(round(freq, 0), "%")),
      position = position_dodge(width = 0.8),
      vjust = -0.3,
      size = 3
    ) +
    scale_fill_manual(
      values = c("low" = project_colors[3], "mid" = project_colors[2], "high" = project_colors[6]),
      name = "EDSS Group"
    ) +
    labs(
      title = "Frequenza delle principali caratteristiche cliniche per gruppo EDSS",
      subtitle = "Caratteristiche maggiormente associate alla progressione della disabilità",
      x = "Clinical Feature",
      y = "Frequenza (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylim(0, 100)
  
  ggsave("outputs/figures/top5_features_freq_by_edss_bin.png", plot = p_top5,
         width = 10, height = 6, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/top5_features_freq_by_edss_bin.png")
  
} else {
  cat("Non abbastanza feature top5 nel dataset\n")
}


# OUTPUT CONSOLE FINALE

cat("\nDataset:\n")
cat("   N (pazienti):", N, "\n")
cat("   Feature binarie totali:", length(all_binary), "\n")
cat("   Feature dopo filtro:", n_kept, "\n")

cat("\nFeature droppate:\n")
if (nrow(dropped_features) > 0) {
  cat("   Low freq (n_1 <", min_count, "):\n")
  low_prev <- dropped_features %>% filter(reason == "dropped_low_freq")
  if (nrow(low_prev) > 0) {
    for (f in low_prev$feature) cat("     -", f, "\n")
  } else {
    cat("     (nessuna)\n")
  }
  
  cat("High freq (n_1 >", max_count, "):\n")
  high_prev <- dropped_features %>% filter(reason == "dropped_high_freq")
  if (nrow(high_prev) > 0) {
    for (f in high_prev$feature) cat("     -", f, "\n")
  } else {
    cat("     (nessuna)\n")
  }
} else {
  cat("   (nessuna feature droppata)\n")
}

cat("\nMDS/PCoA:\n")
cat("   Varianza Dim1:", var_explained[1], "%\n")
cat("   Varianza Dim2:", var_explained[2], "%\n")

cat("\nFile CSV salvati (", length(saved_csv), "):\n")
for (f in saved_csv) {
  cat(f, "\n")
}

cat("\nFile PNG salvati (", length(saved_png), "):\n")
for (f in saved_png) {
  cat(f, "\n")
}

cat("
Note: 
   - Preparo i dati per il clustering ma NON eseguo
     un clustering finale con scelta di k ottimale. 
   - La matrice di distanza Jaccard e le coordinate PCoA sono
     pronte per essere usate in un successivo clustering.
   - Le feature rare/ubiquitarie sono state filtrate per evitare
     che dominino o siano irrilevanti nel calcolo delle distanze. 
")