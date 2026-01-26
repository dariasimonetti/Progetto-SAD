# ============================================================================
# Script: 06_eda_univariate.R
# Descrizione: EDA Univariata - Tabelle e grafici descrittivi
# Autore:  Daria Simonetti
# ============================================================================

if (!require(viridis, quietly = TRUE)) install.packages("viridis")
library(viridis)
source("scripts/01_setup.R")

df <- readRDS("outputs/data/data_complete.rds")

cat("Righe:", nrow(df), "| Colonne:", ncol(df), "\n")

# Verifica dell'assenza di "N" in gender

if ("gender" %in% names(df)) {
  gender_n_count <- sum(df$gender == "N", na.rm = TRUE)
  
  if (gender_n_count > 0) {
    cat("WARNING: Trovate", gender_n_count, "righe con gender = 'N'!\n")
    cat("Distribuzione attuale gender:\n")
    print(table(df$gender, useNA = "ifany"))
  } else {
    cat("Nessun valore 'N' in gender\n")
  }
} else {
  cat("Colonna 'gender' non trovata nel dataset\n")
}

# DEFINIZIONE VARIABILI

# Numeriche
numeric_vars <- c("age", "age_of_onset", "disease_duration", "edss", "binary_sum")

# Categoriche
categorical_vars <- c("gender", "medicine", "mri_edss_diff", "comorbidity", 
                      "edss_bin", "edss_3class")

# Binarie cliniche
clinical_binary <- c("pyramidal", "cerebella", "brain_stem", "sensory", 
                     "sphincters", "visual", "mental", "speech", 
                     "motor_system", "sensory_system", "coordination", 
                     "gait", "bowel_bladder", "mobility", "mental_state", 
                     "optic_disc", "nystagmus", "ocular_move", "swallowing")

# Binarie presenting
presenting_binary <- names(df)[grepl("^present_", names(df))]

cat("\nVariabili identificate:\n")
cat("   Numeriche:", length(numeric_vars), "\n")
cat("   Categoriche:", length(categorical_vars), "\n")
cat("   Binarie cliniche:", length(clinical_binary), "\n")
cat("   Binarie presenting:", length(presenting_binary), "\n")

# Contatori file salvati
saved_csv <- character()
saved_png <- character()


# A)-------------------------TABELLE CSV----------------------------------------

# A1) Summary x numeriche


summary_numeric <- df %>%
  select(all_of(numeric_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    n_valid = sum(! is.na(value)),
    n_na = sum(is.na(value)),
    mean = round(mean(value, na.rm = TRUE), 2), #media
    sd = round(sd(value, na.rm = TRUE), 2), #deviazione standard
    variance = round(var(value, na.rm = TRUE), 2), #varianza
    median = round(median(value, na.rm = TRUE), 2), #mediana
    q1 = round(quantile(value, 0.25, na.rm = TRUE), 2),
    q3 = round(quantile(value, 0.75, na.rm = TRUE), 2),
    iqr = round(IQR(value, na.rm = TRUE), 2),
    min = round(min(value, na.rm = TRUE), 2),
    max = round(max(value, na.rm = TRUE), 2),
    range = round(max(value, na.rm = TRUE) - min(value, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  # Riordina secondo l'ordine originale
  
  mutate(variable = factor(variable, levels = numeric_vars)) %>%
  arrange(variable) %>%
  mutate(variable = as.character(variable))

write_csv(summary_numeric, "outputs/tables/summary_numeric.csv")
saved_csv <- c(saved_csv, "outputs/tables/summary_numeric.csv")


# A2) Frequenze x categoriche


# Funzione per calcolare frequenze
calc_freq <- function(data, var_name, keep_order = FALSE) {
  
  # Calcola NA
  n_na <- sum(is.na(data[[var_name]]))
  
  # Frequenze
  freq_df <- data %>%
    filter(!is.na(.data[[var_name]])) %>%
    count(.data[[var_name]], name = "n") %>%
    rename(level = 1) %>%
    mutate(
      level = as.character(level),
      f = round(n/sum(n), 2),
      perc = round(n / sum(n) * 100, 2)
    )
  
  # Ordina per n decrescente
  if (! keep_order) {
    freq_df <- freq_df %>% arrange(desc(n))
  }
  
  list(freq = freq_df, n_na = n_na)
}

# Salva frequenze per ogni variabile categorica
for (var in categorical_vars) {
  
  # Per edss_group mantieni ordine livelli
  keep_order <- (var == "edss_group")
  
  result <- calc_freq(df, var, keep_order = keep_order)
  
  output_path <- paste0("outputs/tables/freq_", var, ".csv")
  write_csv(result$freq, output_path)
  saved_csv <- c(saved_csv, output_path)
  
}


# A3) Frequenza binarie cliniche


# Filtra solo colonne che esistono nel dataset
clinical_binary_exist <- intersect(clinical_binary, names(df))

freq_clinical <- df %>%
  select(all_of(clinical_binary_exist)) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  summarise(
    n_1 = sum(value == 1, na.rm = TRUE), #freq assoluta 1
    n_0 = sum(value == 0, na.rm = TRUE), #freq assoluta 0
    n_na = sum(is.na(value)),
    .groups = "drop"
  ) %>%
  mutate(
    f_1 = round(n_1/(n_1+n_0), 2),
    f_0 = round(n_0/(n_1+n_0), 2),
    perc_1 = round(n_1 / (n_1 + n_0) * 100, 2) #freq rel %
  ) %>%
  arrange(desc(perc_1))

write_csv(freq_clinical, "outputs/tables/freq_clinical_binary.csv")
saved_csv <- c(saved_csv, "outputs/tables/freq_clinical_binary.csv")


# A4) Freqeunze binarie presenting


if (length(presenting_binary) > 0) {
  
  freq_presenting <- df %>%
    select(all_of(presenting_binary)) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    summarise(
      n_1 = sum(value == 1, na.rm = TRUE),
      n_0 = sum(value == 0, na.rm = TRUE),
      n_na = sum(is.na(value)),
      .groups = "drop"
    ) %>%
    mutate(
      f_1 = round(n_1/(n_1+n_0), 2),
      f_0 = round(n_0/(n_1+n_0), 2),
      perc_1 = round(n_1 / (n_1 + n_0) * 100, 2)
    ) %>%
    arrange(desc(perc_1))
  
  write_csv(freq_presenting, "outputs/tables/freq_presenting_binary.csv")
  saved_csv <- c(saved_csv, "outputs/tables/freq_presenting_binary.csv")
  
} else {
  cat("Nessuna colonna present_* trovata nel dataset\n")
  freq_presenting <- NULL
}


# B)-------------------------GRAFICI PNG----------------------------------------

# B1) Numeriche (istogramma + boxplot + KDP)

plot_numeric_full <- function(data, var_name, bins = 15) {
  
  # Titolo formattato
  title_clean <- str_replace_all(var_name, "_", " ") %>% str_to_title()
  # Larghezza fissa dei bin
  bin_width <- 5
  
  # Creiamo i break dei bin
  breaks <- seq(floor(min(data[[var_name]], na.rm = TRUE) / bin_width) * bin_width,
                ceiling(max(data[[var_name]], na.rm = TRUE) / bin_width) * bin_width,
                by = bin_width)
  
  # Creiamo un nuovo fattore con gli intervalli
  data$bin_interval <- cut(
    data[[var_name]],
    breaks = breaks,
    right = FALSE,   # include l'estremo sinistro e esclude il destro
    include.lowest = TRUE
  )
  
  # Etichette degli intervalli
  levels(data$bin_interval) <- paste(head(breaks, -1), tail(breaks, -1), sep = "-")
  
  # Istogramma con asse x categoriale
  p_hist <- ggplot(data, aes(x = bin_interval)) +
    geom_bar(aes(y = after_stat(count)),
             fill = "#1F9E89FF", color = "white", alpha = 0.7) +
    geom_text(aes(y = after_stat(count), label = after_stat(count)),
              stat = "count",
              vjust = -0.5,
              color = "#482878FF",
              size = 3) +
    labs(
      title = paste("Distribuzione:", title_clean),
      x = title_clean,
      y = "Frequenza"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Boxplot
  p_box <- ggplot(data, aes(y = .data[[var_name]])) +
    geom_boxplot(fill = "#FDE725FF", alpha = 0.7, width = 0.2, outlier.shape = 21 ,outlier.size = 2.5,
                 outlier.color = "#482878FF",
                 outlier.fill = "#1F9E89FF",
                 outlier.alpha = 0.8) +
    labs(y = title_clean, x = "") +
    theme_classic()
  
  p_box <- p_box +
    theme(plot.margin = margin(5, 5, 5, 30))
  
  # Kernel Density Plot
  p_density <- ggplot(data, aes(x = .data[[var_name]])) +
    geom_density(color = "#482878FF", fill = "#1F9E89FF", alpha = 0.3, linewidth = 1) +
    labs(title = paste("Densità di kernel:", title_clean),
         x = title_clean,
         y = "Densità") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12))
  
  # colonna destra più stretta
  right_col <- p_box / p_density + plot_layout(heights = c(1.5, 1))
  p_combined <- p_hist | right_col + plot_layout(widths = c(4, 1.2))
  
  

  
  return(p_combined)
}

# Crea i grafici per le variabili numeriche
for (var in c("age", "age_of_onset", "disease_duration")) {
  
  p <- plot_numeric_full(df, var, bins = 15)
  
  output_path <- paste0("outputs/figures/full_dist_", var, ".png")
  ggsave(output_path, plot = p, width = 10, height = 7, dpi = 300)
  saved_png <- c(saved_png, output_path)
  
}

# EDSS:  boxplot + jitter + barlot frequenze

p_box <- ggplot(df, aes(x = "", y = edss)) +
  geom_boxplot(fill = "#1F9E89FF", alpha = 0.5, width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "#482878FF", alpha = 0.6, size = 2) +
  labs(
    title = "Distribuzione EDSS",
    subtitle = "Boxplot con punti individuali",
    x = "",
    y = "EDSS Score"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

p_bar <- ggplot(df, aes(x = factor(edss), fill = factor(edss))) +
  geom_bar() +
  geom_text(aes(y = after_stat(count), label = after_stat(count)),
            stat = "count",
            vjust = -0.5,
            color = "#482878FF",
            size = 3) +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  labs(
    title = "Frequenze EDSS",
    x = "EDSS Score",
    y = "Frequenza"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")


p_combined <- (p_box | p_bar) + plot_layout(widths = c(2, 2))

ggsave("outputs/figures/edss_box_jitter_bar.png", plot = p_combined, 
       width = 8, height = 6, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/edss_box_jitter_bar.png")

# Binary sum: barplot conteggi

# Crea funzione che genera n colori interpolati
palette <- colorRampPalette(project_colors)
palette_19 <- palette(19)

p_binary_sum <- df %>%
  filter(!is.na(binary_sum)) %>%
  count(binary_sum) %>%
  mutate(
    perc = n / sum(n),                    # calcola la percentuale
    label = paste0(n, " (", round(perc, 2), ")") # crea label: "conteggio (perc%)"
  ) %>%
  ggplot(aes(x = factor(binary_sum), y = n, fill = factor(binary_sum))) +
  geom_col(alpha = 0.8) +
  scale_fill_viridis_d(option =  "C", end = 0.9) +
  geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
  labs(
    title = "Distribuzione Binary Sum",
    subtitle = "Somma indicatori clinici binari",
    x = "Binary Sum",
    y = "Frequenza"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/figures/binary_sum_counts.png", plot = p_binary_sum, 
       width = 8, height = 6, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/binary_sum_counts.png")


# B2) Categoriche: barplot in percentuale

plot_bar_perc <- function(data, var_name, order_by_freq = TRUE,
                          rotate_labels = FALSE, show_n = TRUE,
                          palette = NULL) {
  
  title_clean <- str_replace_all(var_name, "_", " ") %>% str_to_title()
  
  plot_data <- data %>%
    filter(!is.na(.data[[var_name]])) %>%
    count(.data[[var_name]], name = "n") %>%
    mutate(
      perc = n / sum(n),
      level = as.factor(.data[[var_name]])
    )
  
  if (order_by_freq) {
    plot_data <- plot_data %>%
      mutate(level = fct_reorder(level, n, .desc = TRUE))
  }
  
  p <- ggplot(plot_data, aes(x = level, y = perc, fill = level)) +
    geom_col(alpha = 0.85, color = "white", linewidth = 0.3) +
    labs(
      title = paste("Distribuzione:", title_clean),
      x = title_clean,
      y = "Frequenza Relativa",
      fill = title_clean
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  if (is.null(palette)) {
    p <- p + scale_fill_viridis_d(option = "viridis", end = 0.9)
  } else {
    p <- p + scale_fill_manual(values = palette)
  }
  
  if (show_n) {
    p <- p + geom_text(
      aes(label = paste0(round(perc*100, 1), "%\n(n=", n, ")")),
      vjust = -0.25, size = 3
    )
  }
  
  if (rotate_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  }
  
  return(p)
}

# Gender
gender_palette <- c(
  "F" = "#482878FF",
  "M" = "#1F9E89FF"
)

p_gender <- plot_bar_perc(df, "gender", order_by_freq = TRUE, palette = gender_palette)
ggsave("outputs/figures/bar_gender_perc.png", plot = p_gender, 
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_gender_perc.png")

# Medicine
p_medicine <- plot_bar_perc(df, "medicine", order_by_freq = TRUE)
ggsave("outputs/figures/bar_medicine_perc.png", plot = p_medicine, 
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_medicine_perc.png")

# MRI EDSS diff
p_mri <- plot_bar_perc(df, "mri_edss_diff", order_by_freq = FALSE)
ggsave("outputs/figures/bar_mri_edss_diff_perc.png", plot = p_mri, 
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_mri_edss_diff_perc.png")

# Comorbidity
p_comorbidity <- plot_bar_perc(df, "comorbidity", order_by_freq = FALSE)
ggsave("outputs/figures/bar_comorbidity_perc.png", plot = p_comorbidity, 
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_comorbidity_perc.png")

# EDSS bin (mantieni ordine livelli)
p_edss_bin <- plot_bar_perc(df, "edss_bin", order_by_freq = FALSE)
ggsave("outputs/figures/bar_edss_bin_perc.png", plot = p_edss_bin, 
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_edss_bin_perc.png")

# EDSS 3class (mantieni ordine livelli)
p_edss_3class <- plot_bar_perc(df, "edss_3class", order_by_freq = FALSE)
ggsave("outputs/figures/bar_edss_3class_perc.png", plot = p_edss_3class, 
       width = 9, height = 8, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_edss_3class_perc.png")


# B3) Frequenze binarie:  barplot orizzontale


# Cliniche
p_freq_clinical <- freq_clinical %>%
  mutate(feature = fct_reorder(feature, perc_1)) %>%
  ggplot(aes(x = feature, y = perc_1, fill = feature)) +
  geom_col(alpha = 0.85, color = "white", linewidth = 0.3) +
  scale_fill_viridis_d(option = "D", end = 0.95, direction = -1) +
  guides(fill = "none") +
  geom_text(aes(label = paste0(perc_1, "% (n=", n_1, ")")), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Frequenza Indicatori Clinici Binari",
    subtitle = "Percentuale pazienti con valore = 1 per ogni indicatore",
    x = "Indicatore",
    y = "Frequenza (%)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("outputs/figures/bar_frequency_clinical_binary.png", plot = p_freq_clinical, 
       width = 9, height = 7, dpi = 300)
saved_png <- c(saved_png, "outputs/figures/bar_frequency_clinical_binary.png")


# Presenting
if (! is.null(freq_presenting) && nrow(freq_presenting) > 0) {
  
  p_freq_presenting <- freq_presenting %>%
    mutate(feature = fct_reorder(feature, perc_1)) %>%
    ggplot(aes(x = feature, y = perc_1, fill = feature)) +
    geom_col(alpha = 0.85, color = "white", linewidth = 0.3) +
    scale_fill_viridis_d(option = "D", end = 0.95, direction = -1) +
    guides(fill = "none") +
    geom_text(aes(label = paste0(perc_1, "% (n=", n_1, ")")), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Frequenza Presenting Symptoms (Multi-label)",
      subtitle = "Percentuale pazienti con sintomo presente",
      x = "Sintomo",
      y = "Frequenza (%)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("outputs/figures/bar_frequency_presenting.png", plot = p_freq_presenting, 
         width = 9, height = 7, dpi = 300)
  saved_png <- c(saved_png, "outputs/figures/bar_frequency_presenting.png")
  
} else {
  cat("Nessun grafico presenting (colonne non presenti)\n")
}

# D) OUTPUT CONSOLE FINALE



cat("\nFile CSV salvati (", length(saved_csv), "):\n")
for (f in saved_csv) {
  cat(f, "\n")
}

cat("\nFile PNG salvati (", length(saved_png), "):\n")
for (f in saved_png) {
  cat(f, "\n")
}

cat("\nhead(summary_numeric):\n")
cat("--------------------------------------------------------------------------------\n")
print(head(summary_numeric))

