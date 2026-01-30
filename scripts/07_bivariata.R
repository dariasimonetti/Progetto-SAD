# ============================================================================
# Script:  07_eda_bivariate.R
# Descrizione: EDA Bivariata - Correlazioni, test, associazioni binarie
# Autore:  Daria Simonetti
# ============================================================================


if (!require(viridis, quietly = TRUE)) install.packages("viridis")
library(viridis)
source("scripts/01_setup.R")

if (!require(e1071, quietly = TRUE)) install.packages("e1071")
library(e1071)

# Riproducibilità
set.seed(123)

# Carica dati
df <- readRDS("outputs/data/data_complete.rds")

# DEFINIZIONE VARIABILI

# Numeriche
numeric_vars <- c("age", "age_of_onset", "disease_duration", "edss", "binary_sum")

# Categoriche per confronti
categorical_vars <- c("gender", "medicine", "mri_edss_diff", "comorbidity")

# Binarie cliniche
clinical_binary <- c("pyramidal", "cerebella", "brain_stem", "sensory", 
                     "sphincters", "visual", "mental", "speech", 
                     "motor_system", "sensory_system", "coordination", 
                     "gait", "bowel_bladder", "mobility", "mental_state", 
                     "optic_disc", "nystagmus", "ocular_move", "swallowing")

# Binarie presenting
presenting_binary <- names(df)[grepl("^present_", names(df))]

# Tutte le binarie
all_binary <- c(clinical_binary, presenting_binary)
all_binary <- intersect(all_binary, names(df))

cat("\nVariabili:\n")
cat("   Numeriche:", length(numeric_vars), "\n")
cat("   Categoriche:", length(categorical_vars), "\n")
cat("   Binarie totali:", length(all_binary), "\n")

# Contatori file salvati
saved_csv <- character()
saved_png <- character()


# A)-------------------CORRELAZIONI SPEARMAN CON EDSS---------------------------

# Variabili da correlare con EDSS
corr_vars <- c("age", "age_of_onset", "disease_duration", 
               "log_disease_duration", "binary_sum")

# Funzione per calcolare correlazione Spearman
calc_spearman <- function(data, x_var, y_var = "edss") {
  complete_data <- data %>%
    select(x = all_of(x_var), y = all_of(y_var)) %>%
    filter(!is.na(x) & !is.na(y))
  
  n <- nrow(complete_data)
  
  if (n < 3) {
    return(tibble(x = x_var, y = y_var, n = n, rho = NA_real_, p_value = NA_real_))
  }
  
  test <- cor.test(complete_data$x, complete_data$y, method = "spearman", exact = FALSE)
  
  tibble(
    x = x_var,
    y = y_var,
    n = n,
    rho = round(test$estimate, 3),
    p_value = round(test$p.value, 4)
  )
}

# Calcola correlazioni
spearman_corr <- map_dfr(corr_vars, ~calc_spearman(df, .x))

write_csv(spearman_corr, "outputs/tables/spearman_correlations.csv")
saved_csv <- c(saved_csv, "outputs/tables/spearman_correlations.csv")

cat("outputs/tables/spearman_correlations.csv\n")
print(spearman_corr)


# B)-------------------GRAFICI SCATTER EDSS vs NUMERICHE------------------------

# Funzione per creare scatter plot
create_scatter <- function(data, x_var, y_var = "edss") {
  
  title_x <- str_replace_all(x_var, "_", " ") %>% str_to_title()
  
  p <- data %>%
    filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]])) %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(
      position = position_jitter(height = 0.05, width = 0),
      alpha = 0.6,
      color = plasma(6)[4],
      size = 2
    ) +
    geom_smooth(
      method = "loess",
      se = TRUE,
      color = plasma(6)[3],
      fill = plasma(6)[3],
      alpha = 0.2
    ) +
    labs(
      title = paste("EDSS vs", title_x),
      subtitle = "Correlazione Spearman con smoothing LOESS",
      x = title_x,
      y = "EDSS Score"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  return(p)
}

# Crea scatter per ogni variabile
scatter_vars <- c("age", "age_of_onset", "disease_duration", 
                  "log_disease_duration", "binary_sum")

for (var in scatter_vars) {
  p <- create_scatter(df, var)
  filename <- paste0("outputs/figures/scatter_edss_", var, ".png")
  ggsave(filename, plot = p, width = 8, height = 5, dpi = 300)
  saved_png <- c(saved_png, filename)
}


# C)-------------------CONFRONTI EDSS TRA GRUPPI (TEST NON PARAMETRICI)---------

# tabella risultati
bivariate_tests <- tibble(
  predictor = character(),
  test = character(),
  groups = character(),
  n_total = integer(),
  statistic = numeric(),
  p_value = numeric()
)

# --- Wilcoxon per gender ---
cat("\nTest Wilcoxon:  EDSS ~ gender...\n")

df_gender <- df %>% 
  filter(! is.na(gender) & !is.na(edss) & gender %in% c("F", "M"))

if (nrow(df_gender) >= 10) {
  test_gender <- wilcox.test(edss ~ gender, data = df_gender, exact = FALSE)
  bivariate_tests <- bind_rows(bivariate_tests, tibble(
    predictor = "gender",
    test = "Wilcoxon",
    groups = "F vs M",
    n_total = nrow(df_gender),
    statistic = round(test_gender$statistic, 2),
    p_value = round(test_gender$p.value, 4)
  ))
}

# --- Wilcoxon per comorbidity ---
cat("Test Wilcoxon:  EDSS ~ comorbidity...\n")

df_comorbidity <- df %>% 
  filter(! is.na(comorbidity) & !is.na(edss))

if (nrow(df_comorbidity) >= 10) {
  test_comorbidity <- wilcox.test(edss ~ comorbidity, data = df_comorbidity, exact = FALSE)
  bivariate_tests <- bind_rows(bivariate_tests, tibble(
    predictor = "comorbidity",
    test = "Wilcoxon",
    groups = "No vs Yes",
    n_total = nrow(df_comorbidity),
    statistic = round(test_comorbidity$statistic, 2),
    p_value = round(test_comorbidity$p.value, 4)
  ))
}

# --- Wilcoxon per mri_edss_diff ---
cat("Test Wilcoxon:  EDSS ~ mri_edss_diff...\n")

df_mri <- df %>% 
  filter(!is.na(mri_edss_diff) & !is.na(edss))

if (nrow(df_mri) >= 10) {
  test_mri <- wilcox.test(edss ~ mri_edss_diff, data = df_mri, exact = FALSE)
  bivariate_tests <- bind_rows(bivariate_tests, tibble(
    predictor = "mri_edss_diff",
    test = "Wilcoxon",
    groups = "No vs Yes",
    n_total = nrow(df_mri),
    statistic = round(test_mri$statistic, 2),
    p_value = round(test_mri$p.value, 4)
  ))
}

# --- Kruskal-Wallis per medicine ---
cat("Test Kruskal-Wallis: EDSS ~ medicine...\n")

df_medicine <- df %>% 
  filter(!is.na(medicine) & !is.na(edss))

if (nrow(df_medicine) >= 10) {
  test_medicine <- kruskal.test(edss ~ medicine, data = df_medicine)
  medicine_levels <- paste(levels(df_medicine$medicine), collapse = ", ")
  bivariate_tests <- bind_rows(bivariate_tests, tibble(
    predictor = "medicine",
    test = "Kruskal-Wallis",
    groups = medicine_levels,
    n_total = nrow(df_medicine),
    statistic = round(test_medicine$statistic, 2),
    p_value = round(test_medicine$p.value, 4)
  ))
}

# p-value aggiustato (BH)
bivariate_tests <- bivariate_tests %>%
  mutate(p_adj_bh = round(p.adjust(p_value, method = "BH"), 4))

write_csv(bivariate_tests, "outputs/tables/bivariate_tests_edss_vs_categorical.csv")
saved_csv <- c(saved_csv, "outputs/tables/bivariate_tests_edss_vs_categorical.csv")

cat("outputs/tables/bivariate_tests_edss_vs_categorical.csv\n")
print(bivariate_tests)

# --- Boxplot per ogni variabile categorica ---

# Funzione per boxplot + jitter
create_boxplot <- function(data, x_var, y_var = "edss", rotate_x = FALSE) {
  
  title_x <- str_replace_all(x_var, "_", " ") %>% str_to_title()
  
  data <- data %>%
    filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]])) %>%
    mutate(!!x_var := as.factor(.data[[x_var]]))
  
  n_levels <- nlevels(data[[x_var]])
  
  # Palette plasma dinamica
  fill_colors <- plasma(
    n_levels,
    begin = 0.25,
    end   = 0.75
  )
  
  point_colors <- plasma(
    n_levels,
    begin = 0.30,
    end   = 0.70
  )
  
  p <- ggplot(
    data,
    aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      fill = .data[[x_var]],
      color = .data[[x_var]]
    )
  ) +
    geom_boxplot(
      alpha = 0.5,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.2,
      height = 0,
      alpha = 0.7,
      size = 2
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = point_colors) +
    labs(
      title = paste("EDSS per", title_x),
      x = title_x,
      y = "EDSS Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_blank()
    )
  
  if (rotate_x) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}


# Gender
p_gender <- create_boxplot(df %>% filter(gender %in% c("F", "M")), "gender")
ggsave("outputs/figures/box_edss_by_gender.png", plot = p_gender, 
       width = 8, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/box_edss_by_gender.png")

# Comorbidity
p_comorbidity <- create_boxplot(df, "comorbidity")
ggsave("outputs/figures/box_edss_by_comorbidity.png", plot = p_comorbidity, 
       width = 8, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/box_edss_by_comorbidity.png")

# MRI EDSS diff
p_mri <- create_boxplot(df, "mri_edss_diff")
ggsave("outputs/figures/box_edss_by_mri_edss_diff.png", plot = p_mri, 
       width = 8, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/box_edss_by_mri_edss_diff.png")

# Medicine (ordinate per frequenza, ruota etichette)
p_medicine <- df %>%
  filter(!is.na(medicine) & !is.na(edss)) %>%
  mutate(medicine = fct_infreq(medicine)) %>%
  ggplot(
    aes(
      x = medicine,
      y = edss,
      fill = medicine,
      color = medicine
    )
  ) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.7,
    size = 2
  ) +
  scale_fill_manual(
    values = plasma(
      nlevels(df$medicine),
      begin = 0.20,
      end   = 0.80
    )
  ) +
  scale_color_manual(
    values = plasma(
      nlevels(df$medicine),
      begin = 0.20,
      end   = 0.80
    )
  ) +
  labs(
    title = "EDSS by Medicine",
    x = "Medicine (ordered by frequency)",
    y = "EDSS Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


ggsave("outputs/figures/box_edss_by_medicine.png", plot = p_medicine, 
       width = 8, height = 5, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/box_edss_by_medicine.png")


# C-bis)---------------BINARY_SUM VS EDSS CLASSES-------------------------------

# Funzione helper
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  format(round(p, 3), nsmall = 3)
}

# Tabella risultati
binarysum_edssclass_tests <- tibble(
  outcome = character(),
  predictor = character(),
  test = character(),
  groups = character(),
  n_total = integer(),
  statistic = numeric(),
  p_value = numeric(),
  p_value_label = character(),
  note = character()
)

# --- binary_sum ~ edss_bin (Wilcoxon) ---

if ("edss_bin" %in% names(df)) {
  
  df_edss_bin <- df %>%
    filter(! is.na(edss_bin) & ! is.na(binary_sum))
  
  # Conta per gruppo
  counts_bin <- df_edss_bin %>% count(edss_bin)
  
  # Verifica n >= 5 per entrambi i gruppi
  if (all(counts_bin$n >= 5) && nrow(counts_bin) == 2) {
    
    test_bin <- wilcox.test(binary_sum ~ edss_bin, data = df_edss_bin, exact = FALSE)
    
    binarysum_edssclass_tests <- bind_rows(binarysum_edssclass_tests, tibble(
      outcome = "binary_sum",
      predictor = "edss_bin",
      test = "Wilcoxon",
      groups = paste(levels(df_edss_bin$edss_bin), collapse = " vs "),
      n_total = nrow(df_edss_bin),
      statistic = round(test_bin$statistic, 2),
      p_value = test_bin$p.value,
      p_value_label = format_p(test_bin$p.value),
      note = "tested"
    ))
    
    cat("Test eseguito (n =", nrow(df_edss_bin), ")\n")
    
  } else {
    binarysum_edssclass_tests <- bind_rows(binarysum_edssclass_tests, tibble(
      outcome = "binary_sum",
      predictor = "edss_bin",
      test = "Wilcoxon",
      groups = paste(levels(df$edss_bin), collapse = " vs "),
      n_total = nrow(df_edss_bin),
      statistic = NA_real_,
      p_value = NA_real_,
      p_value_label = "NA",
      note = "skipped_low_counts"
    ))
    cat("Skipped:  n < 5 in almeno un gruppo\n")
  }
  
} else {
  cat("Colonna edss_bin non trovata\n")
}

# --- binary_sum ~ edss_3class (Kruskal-Wallis) ---

if ("edss_3class" %in% names(df)) {
  
  df_edss_3class <- df %>%
    filter(!is.na(edss_3class) & !is.na(binary_sum))
  
  # Conta per gruppo
  counts_3class <- df_edss_3class %>% count(edss_3class)
  
  # Verifica almeno 2 classi con n >= 5
  classes_with_enough <- sum(counts_3class$n >= 5)
  
  if (classes_with_enough >= 2) {
    
    test_3class <- kruskal.test(binary_sum ~ edss_3class, data = df_edss_3class)
    
    binarysum_edssclass_tests <- bind_rows(binarysum_edssclass_tests, tibble(
      outcome = "binary_sum",
      predictor = "edss_3class",
      test = "Kruskal-Wallis",
      groups = paste(levels(df_edss_3class$edss_3class), collapse = ", "),
      n_total = nrow(df_edss_3class),
      statistic = round(test_3class$statistic, 2),
      p_value = test_3class$p.value,
      p_value_label = format_p(test_3class$p.value),
      note = "tested"
    ))
    
    cat("Test eseguito (n =", nrow(df_edss_3class), ")\n")
    
  } else {
    binarysum_edssclass_tests <- bind_rows(binarysum_edssclass_tests, tibble(
      outcome = "binary_sum",
      predictor = "edss_3class",
      test = "Kruskal-Wallis",
      groups = paste(levels(df$edss_3class), collapse = ", "),
      n_total = nrow(df_edss_3class),
      statistic = NA_real_,
      p_value = NA_real_,
      p_value_label = "NA",
      note = "skipped_low_counts"
    ))
    cat("Skipped: meno di 2 classi con n >= 5\n")
  }
  
} else {
  cat("Colonna edss_3class non trovata\n")
}

# Salva tabella test
write_csv(binarysum_edssclass_tests, "outputs/tables/binary_sum_vs_edss_classes_tests.csv")
saved_csv <- c(saved_csv, "outputs/tables/binary_sum_vs_edss_classes_tests.csv")

cat("\n otputs/tables/binary_sum_vs_edss_classes_tests.csv\n")
print(binarysum_edssclass_tests)

# --- Grafici boxplot ---

# Boxplot binary_sum ~ edss_bin
if ("edss_bin" %in% names(df)) {
  
  p_bs_edss_bin <- df %>%
    filter(!is.na(edss_bin) & !is.na(binary_sum)) %>%
    mutate(edss_bin = as.factor(edss_bin)) %>%
    ggplot(
      aes(
        x = edss_bin,
        y = binary_sum,
        fill = edss_bin,
        color = edss_bin
      )
    ) +
    geom_boxplot(
      alpha = 0.5,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.2,
      height = 0,
      alpha = 0.7,
      size = 2
    ) +
    scale_fill_manual(
      values = plasma(2, begin = 0.35, end = 0.65)
    ) +
    scale_color_manual(
      values = plasma(2, begin = 0.35, end = 0.65)
    ) +
    labs(
      title = "Binary Sum per EDSS (classificazione binaria)",
      subtitle = "≤ 2.0 vs > 2.0",
      x = "EDSS Classe Binaria",
      y = "Binary Sum"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  
  ggsave("outputs/figures/box_binary_sum_by_edss_bin.png", plot = p_bs_edss_bin,
         width = 8, height = 5, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/box_binary_sum_by_edss_bin.png")
  cat("outputs/figures/box_binary_sum_by_edss_bin.png\n")
}

# Boxplot binary_sum ~ edss_3class
if ("edss_3class" %in% names(df)) {
  
  p_bs_edss_3class <- df %>%
    filter(!is.na(edss_3class) & !is.na(binary_sum)) %>%
    mutate(edss_3class = as.factor(edss_3class)) %>%
    ggplot(
      aes(
        x = edss_3class,
        y = binary_sum,
        fill = edss_3class,
        color = edss_3class
      )
    ) +
    geom_boxplot(
      alpha = 0.5,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.2,
      height = 0,
      alpha = 0.7,
      size = 2
    ) +
    scale_fill_manual(
      values = plasma(3, begin = 0.30, end = 0.70)
    ) +
    scale_color_manual(
      values = plasma(3, begin = 0.30, end = 0.70)
    ) +
    labs(
      title = "Binary Sum per EDSS (3-classi)",
      subtitle = "normal, mild, severe",
      x = "EDSS 3-Classi",
      y = "Binary Sum"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 15, hjust = 1),
      legend.position = "none"
    )
  
  
  ggsave("outputs/figures/box_binary_sum_by_edss_3class.png", plot = p_bs_edss_3class,
         width = 8, height = 5, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/box_binary_sum_by_edss_3class.png")
  cat("outputs/figures/box_binary_sum_by_edss_3class.png\n")
}

# D)-------------------EDSS VS VARIABILI BINARIE (WILCOXON BATCH)---------------

# Funzione per test Wilcoxon su variabile binaria
wilcox_binary <- function(data, feature_name, outcome = "edss") {
  
  # Estrai dati
  feature_vals <- data[[feature_name]]
  outcome_vals <- data[[outcome]]
  
  # Conta
  n_1 <- sum(feature_vals == 1, na.rm = TRUE)
  n_0 <- sum(feature_vals == 0, na.rm = TRUE)
  n_na <- sum(is.na(feature_vals))
  
  # Mediane
  median_1 <- median(outcome_vals[feature_vals == 1], na.rm = TRUE)
  median_0 <- median(outcome_vals[feature_vals == 0], na.rm = TRUE)
  
  # Test solo se n sufficiente
  if (n_1 >= 5 & n_0 >= 5) {
    df_test <- data %>%
      filter(! is.na(.data[[feature_name]]) & !is.na(.data[[outcome]]))
    
    test <- wilcox.test(
      df_test[[outcome]][df_test[[feature_name]] == 1],
      df_test[[outcome]][df_test[[feature_name]] == 0],
      exact = FALSE
    )
    
    return(tibble(
      feature = feature_name,
      n_1 = n_1,
      n_0 = n_0,
      n_na = n_na,
      median_edss_1 = round(median_1, 2),
      median_edss_0 = round(median_0, 2),
      statistic = round(test$statistic, 2),
      p_value = round(test$p.value, 4),
      note = "tested"
    ))
  } else {
    return(tibble(
      feature = feature_name,
      n_1 = n_1,
      n_0 = n_0,
      n_na = n_na,
      median_edss_1 = round(median_1, 2),
      median_edss_0 = round(median_0, 2),
      statistic = NA_real_,
      p_value = NA_real_,
      note = "skipped_low_counts"
    ))
  }
}

# Applica a tutte le binarie
wilcox_results <- map_dfr(all_binary, ~wilcox_binary(df, .x))

# Calcola p_adj_bh solo per i test eseguiti
wilcox_results <- wilcox_results %>%
  mutate(
    p_adj_bh = if_else(
      note == "tested",
      round(p.adjust(p_value[note == "tested"], method = "BH")[cumsum(note == "tested")], 4),
      NA_real_
    )
  )

# Ricalcola p_adj correttamente
tested_pvals <- wilcox_results %>% filter(note == "tested") %>% pull(p_value)
tested_padj <- p.adjust(tested_pvals, method = "BH")

wilcox_results <- wilcox_results %>%
  mutate(row_id = row_number()) %>%
  left_join(
    wilcox_results %>%
      filter(note == "tested") %>%
      mutate(p_adj_bh = round(p.adjust(p_value, method = "BH"), 4)) %>%
      select(feature, p_adj_bh),
    by = "feature",
    suffix = c("", "_new")
  ) %>%
  mutate(p_adj_bh = coalesce(p_adj_bh_new, NA_real_)) %>%
  select(-p_adj_bh_new, -row_id)

# Ordina:  p_value crescente, skipped in fondo
wilcox_results <- wilcox_results %>%
  arrange(note != "tested", p_value)

write_csv(wilcox_results, "outputs/tables/wilcoxon_edss_vs_binary.csv")
saved_csv <- c(saved_csv, "outputs/tables/wilcoxon_edss_vs_binary.csv")


n_tested <- sum(wilcox_results$note == "tested")
n_skipped <- sum(wilcox_results$note == "skipped_low_counts")
cat("   Testate:", n_tested, "| Skipped:", n_skipped, "\n")

# --- Grafico top 10 median diff ---

top10_binary <- wilcox_results %>%
  filter(note == "tested") %>%
  slice_min(p_value, n = 10) %>%
  mutate(
    median_diff = median_edss_1 - median_edss_0,
    feature = fct_reorder(feature, median_diff)
  )

if (nrow(top10_binary) > 0) {
  
  # Palette prisma (gradiente viola-giallo centrale)
  prisma_pal <- plasma(6, begin = 0.25, end = 0.75)
  
  p_top10_grad <- ggplot(top10_binary, aes(x = median_diff, y = feature, color = median_diff)) +
    geom_segment(aes(x = 0, xend = median_diff, y = feature, yend = feature),
                 linewidth = 1) +
    geom_point(size = 4) +
    scale_color_gradientn(
      colors = prisma_pal,
      guide = guide_colorbar(title = "Differenza Mediana")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Top 10 Binary Features: Differenza mediana EDSS",
      subtitle = "Differenza tra mediana EDSS (gruppo=1) e mediana EDSS (gruppo=0)",
      x = "Differenza mediana EDSS (gruppo 1 - gruppo 0)",
      y = "Feature"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  ggsave("outputs/figures/top10_binary_edss_median_diff_prisma_grad.png", plot = p_top10_grad,
         width = 10, height = 4, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/top10_binary_edss_median_diff_prisma_grad.png")
}

# E)-------------------ASSOCIAZIONI TRA VARIABILI BINARIE-----------------------

# --- E1) Selezione feature binarie ---

# Costruisco dataframe binario
bin_df <- df %>%
  select(all_of(all_binary)) %>%
  mutate(across(everything(), as.numeric))

# Nota: NA -> 0 (nel dataset attuale n_na = 0, quindi non cambia)
# se ci fossero NA, questa scelta potrebbe influenzare le associazioni
bin_df[is.na(bin_df)] <- 0

cat("   Dimensioni:", nrow(bin_df), "x", ncol(bin_df), "\n")

# --- E2) Esclusione feature a varianza zero ---

variances <- vapply(bin_df, function(x) stats::var(x, na.rm = TRUE), numeric(1))
zero_var_features <- names(variances)[is.na(variances) | variances == 0]

if (length(zero_var_features) > 0) {
  cat("Feature rimosse (varianza = 0):", paste(zero_var_features, collapse = ", "), "\n")
  bin_df_filtered <- bin_df %>% select(-all_of(zero_var_features))
} else {
  cat("Nessuna feature a varianza zero\n")
  bin_df_filtered <- bin_df
}

cat("Features per matrici:", ncol(bin_df_filtered), "\n")

# --- E3) Matrice PHI (correlazione Pearson su 0/1) ---

phi_mat <- cor(bin_df_filtered, method = "pearson")

# Salva matrice wide
phi_df <- as.data.frame(phi_mat) %>%
  mutate(feature = rownames(phi_mat)) %>%
  select(feature, everything())

write_csv(phi_df, "outputs/tables/phi_matrix.csv")
saved_csv <- c(saved_csv, "outputs/tables/phi_matrix.csv")

# Heatmap PHI
phi_long <- phi_mat %>%
  as.data.frame() %>%
  mutate(feature_a = rownames(. )) %>%
  pivot_longer(-feature_a, names_to = "feature_b", values_to = "phi")

p_phi_heat <- ggplot(phi_long, aes(x = feature_a, y = feature_b, fill = phi)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#807DBA",   # viola soft
    mid = "#F7F7F7",   # bianco
    high = "#FDB863",  # arancio chiaro
    midpoint = 0,
    limits = c(-1, 1),
    name = "PHI"
  )  +
  labs(
    title = "Matrice di correlazione PHI \n (caratteristiche binarie)",
    x = "", y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  coord_fixed()

ggsave("outputs/figures/heatmap_phi_binaries.png", plot = p_phi_heat,
       width = 10, height = 9, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/heatmap_phi_binaries.png")


# --- E5) Top coppie più associate ---

# Funzione per estrarre upper triangle
extract_pairs <- function(mat, value_name) {
  n <- nrow(mat)
  pairs <- tibble()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      pairs <- bind_rows(pairs, tibble(
        feature_a = rownames(mat)[i],
        feature_b = colnames(mat)[j],
        value = mat[i, j]
      ))
    }
  }
  
  names(pairs)[3] <- value_name
  return(pairs)
}

# Top 20 PHI (per |phi|)
phi_pairs <- extract_pairs(phi_mat, "phi") %>%
  filter(! is.na(phi)) %>%
  mutate(abs_phi = abs(phi)) %>%
  slice_max(abs_phi, n = 20) %>%
  select(-abs_phi) %>%
  arrange(desc(abs(phi)))

write_csv(phi_pairs, "outputs/tables/top20_phi_pairs.csv")
saved_csv <- c(saved_csv, "outputs/tables/top20_phi_pairs.csv")


# Top 20 PHI con colori più soft
p_top_phi <- ggplot(phi_pairs, aes(x = abs(phi), 
                                   y = fct_reorder(paste(feature_a, feature_b, sep=" - "), abs(phi)), 
                                   fill = abs(phi))) +
  geom_col() +
  scale_fill_gradient(low = "#BFD3E6", high = "#2C7BB6") +
  labs(title="Top 20 PHI", x="Valore assoluto Phi", y="Coppie di feature") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold"), legend.position="none")

# Combinazione unica con patchwork
figura_unica <- (p_phi_heat) | (p_top_phi)

# Salvataggio
ggsave("outputs/figures/figura_unica_binary_associations_softcolors.png", plot = figura_unica,
       width = 11, height = 4, dpi = 300)


# OUTPUT CONSOLE FINALE

cat("\nDataset:\n")
cat("   Righe:", nrow(df), "\n")
cat("   Colonne:", ncol(df), "\n")

cat("\nNA per variabili principali:\n")
main_vars <- c("edss", "age", "age_of_onset", "disease_duration", 
               "gender", "medicine", "comorbidity")
na_summary <- sapply(df[main_vars], function(x) sum(is.na(x)))
print(na_summary)

cat("\nBinarie testate vs skipped:\n")
cat("Testate:", n_tested, "\n")
cat("Skipped (low counts):", n_skipped, "\n")

cat("\nhead(skewness_df):\n")
print(head(skewness_df))

cat("\nhead(spearman_corr):\n")
print(head(spearman_corr))

cat("\nbivariate_tests_edss_vs_categoricals:\n")
print(bivariate_tests)

cat("\nhead(wilcox_results):\n")
print(head(wilcox_results, 10))

cat("\nFile CSV salvati (", length(saved_csv), "):\n")
for (f in saved_csv) {
  cat(f, "\n")
}

cat("\nFile PNG salvati (", length(saved_png), "):\n")
for (f in saved_png) {
  cat(f, "\n")
}

cat("\nbinary_sum_vs_edss_classes_tests:\n")
print(binarysum_edssclass_tests)
