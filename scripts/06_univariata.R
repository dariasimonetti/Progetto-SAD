# ============================================================================
# Script:  06_univariata. R
# Descrizione: Analisi esplorativa univariata - Tabelle e Grafici
# Autore: Daria-Simonetti
# ============================================================================

source("scripts/01_setup.R")


input_path <- "outputs/data/dataset_ready_with_symptoms.rds"

if (!file.exists(input_path)) {
  stop("File non trovato: ", input_path)
}

df <- readRDS(input_path)
cat("   Righe:", nrow(df), "| Colonne:", ncol(df), "\n")


# A) TABELLE DESCRITTIVE

# ----------------------------------------------------------------------------
# A1) Summary numeriche
# ----------------------------------------------------------------------------

numeric_vars <- c("age", "age_of_onset", "disease_duration", "edss", "binary_sum")

summary_numeric <- df %>%
  select(all_of(numeric_vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = round(mean(value, na.rm = TRUE), 2),
    sd = round(sd(value, na.rm = TRUE), 2),
    median = round(median(value, na.rm = TRUE), 2),
    iqr = round(IQR(value, na.rm = TRUE), 2),
    min = round(min(value, na.rm = TRUE), 2),
    max = round(max(value, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  # Mantiene ordine originale
  
  mutate(variable = factor(variable, levels = numeric_vars)) %>%
  arrange(variable) %>%
  mutate(variable = as.character(variable))

write_csv(summary_numeric, "outputs/tables/summary_numeric.csv")

# ----------------------------------------------------------------------------
# A2) Frequenze categoriche
# ----------------------------------------------------------------------------

categorical_vars <- c("gender", "medicine", "mri_edss_diff", "comorbidity", "edss_group")

# Funzione per creare tabella frequenze
create_freq_table <- function(data, var_name, keep_level_order = FALSE) {
  freq_table <- data %>%
    filter(!is.na(!!sym(var_name))) %>%
    count(!!sym(var_name), name = "n") %>%
    mutate(perc = round(n / sum(n) * 100, 2)) %>%
    rename(level = !!sym(var_name)) %>%
    mutate(level = as.character(level))
  
  if (! keep_level_order) {
    freq_table <- freq_table %>% arrange(desc(n))
  }
  
  return(freq_table)
}

# Genera e salva tabelle
for (var in categorical_vars) {
  # Mantiene ordine livelli per gender e edss_group
  keep_order <- var %in% c("gender", "edss_group")
  
  freq_table <- create_freq_table(df, var, keep_level_order = keep_order)
  
  file_path <- paste0("outputs/tables/freq_", var, ". csv")
  write_csv(freq_table, file_path)
}

# ----------------------------------------------------------------------------
# A3) Prevalenza variabili binarie cliniche
# ----------------------------------------------------------------------------

# Identificazione colonne binarie cliniche (escludi present_*)
clinical_binary_cols <- c(
  "pyramidal", "cerebella", "brain_stem", "sensory", "sphincters",
  "visual", "mental", "speech", "motor_system", "sensory_system",
  "coordination", "gait", "bowel_bladder", "mobility", "mental_state",
  "optic_disc", "nystagmus", "ocular_move", "swallowing"
)

# Filtra solo colonne esistenti
clinical_binary_cols <- clinical_binary_cols[clinical_binary_cols %in% names(df)]

prevalence_clinical <- tibble(
  feature = clinical_binary_cols,
  n_1 = sapply(clinical_binary_cols, function(col) sum(df[[col]] == 1, na.rm = TRUE)),
  n_total = sapply(clinical_binary_cols, function(col) sum(! is.na(df[[col]])))
) %>%
  mutate(
    perc_1 = round(n_1 / n_total * 100, 2)
  ) %>%
  select(feature, n_1, perc_1) %>%
  arrange(desc(perc_1))

write_csv(prevalence_clinical, "outputs/tables/prevalence_clinical_binary.csv")

# ----------------------------------------------------------------------------
# A4) Prevalenza presenting symptoms (present_*)
# ----------------------------------------------------------------------------

presenting_cols <- names(df)[str_detect(names(df), "^present_")]

prevalence_presenting <- tibble(
  feature = presenting_cols,
  n_1 = sapply(presenting_cols, function(col) sum(df[[col]] == 1, na.rm = TRUE)),
  n_total = sapply(presenting_cols, function(col) sum(!is.na(df[[col]])))
) %>%
  mutate(
    perc_1 = round(n_1 / n_total * 100, 2)
  ) %>%
  select(feature, n_1, perc_1) %>%
  arrange(desc(perc_1))

write_csv(prevalence_presenting, "outputs/tables/prevalence_presenting_binary.csv")

# B) GRAFICI

# Parametri grafici comuni
fig_width <- 8
fig_height <- 6
fig_dpi <- 300

# Colori
fill_color <- project_colors[1]
fill_color_alt <- project_colors[5]

# ----------------------------------------------------------------------------
# B1) Grafici numeriche
# ----------------------------------------------------------------------------

cat("\n📈 B1) Grafici variabili numeriche.. .\n")

# --- Istogramma + densità:  AGE ---
p_age <- ggplot(df, aes(x = age)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 5, 
                 fill = fill_color, 
                 color = "white",
                 alpha = 0.7) +
  geom_density(color = project_colors[4], linewidth = 1.2) +
  labs(
    title = "Distribuzione dell'età",
    subtitle = paste0("n = ", sum(! is.na(df$age)), " | Media = ", round(mean(df$age, na.rm = TRUE), 1), " anni"),
    x = "Età (anni)",
    y = "Densità"
  )

ggsave("outputs/figures/hist_age.png", p_age, width = fig_width, height = fig_height, dpi = fig_dpi)

# --- Istogramma + densità:  DISEASE_DURATION ---
p_duration <- ggplot(df, aes(x = disease_duration)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 2, 
                 fill = fill_color, 
                 color = "white",
                 alpha = 0.7) +
  geom_density(color = project_colors[4], linewidth = 1.2) +
  labs(
    title = "Distribuzione della durata di malattia",
    subtitle = paste0("n = ", sum(!is.na(df$disease_duration)), 
                      " | Media = ", round(mean(df$disease_duration, na.rm = TRUE), 1), " anni"),
    x = "Durata di malattia (anni)",
    y = "Densità"
  )

ggsave("outputs/figures/hist_disease_duration.png", p_duration, width = fig_width, height = fig_height, dpi = fig_dpi)

# --- Boxplot + jitter: EDSS ---
p_edss <- ggplot(df, aes(x = "", y = edss)) +
  geom_boxplot(fill = fill_color, alpha = 0.5, width = 0.5, outlier.shape = NA) +
  geom_jitter(color = fill_color, alpha = 0.7, width = 0.15, size = 2) +
  labs(
    title = "Distribuzione EDSS",
    subtitle = paste0("n = ", sum(! is.na(df$edss)), 
                      " | Mediana = ", median(df$edss, na.rm = TRUE)),
    x = "",
    y = "EDSS Score"
  ) +
  coord_flip()

ggsave("outputs/figures/box_edss_jitter.png", p_edss, width = fig_width, height = 5, dpi = fig_dpi)

# --- Barplot conteggi: BINARY_SUM ---
binary_sum_counts <- df %>%
  filter(!is.na(binary_sum)) %>%
  count(binary_sum) %>%
  mutate(perc = round(n / sum(n) * 100, 1))

p_binary_sum <- ggplot(binary_sum_counts, aes(x = factor(binary_sum), y = n)) +
  geom_col(fill = fill_color, alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  labs(
    title = "Distribuzione Binary Sum (somma indicatori clinici)",
    subtitle = paste0("n = ", sum(binary_sum_counts$n)),
    x = "Binary Sum",
    y = "Frequenza"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("outputs/figures/bar_binary_sum_counts.png", p_binary_sum, width = fig_width, height = fig_height, dpi = fig_dpi)

# ----------------------------------------------------------------------------
# B2) Grafici categoriche
# ----------------------------------------------------------------------------

# Funzione per creare barplot in percentuale
create_perc_barplot <- function(data, var_name, title, xlab, 
                                order_by_freq = TRUE, horizontal = FALSE) {
  
  plot_data <- data %>%
    filter(!is.na(!!sym(var_name))) %>%
    count(!!sym(var_name), name = "n") %>%
    mutate(
      perc = round(n / sum(n) * 100, 1),
      label = paste0(n, " (", perc, "%)")
    )
  
  if (order_by_freq) {
    plot_data <- plot_data %>%
      mutate(!! sym(var_name) := fct_reorder(!! sym(var_name), n))
  }
  
  p <- ggplot(plot_data, aes(x = !!sym(var_name), y = perc)) +
    geom_col(fill = fill_color, alpha = 0.8) +
    geom_text(aes(label = label), 
              hjust = ifelse(horizontal, -0.1, 0.5),
              vjust = ifelse(horizontal, 0.5, -0.5),
              size = 3.5) +
    labs(
      title = title,
      subtitle = paste0("n = ", sum(plot_data$n)),
      x = xlab,
      y = "Percentuale (%)"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  if (horizontal) {
    p <- p + coord_flip()
  }
  
  return(p)
}

# --- MEDICINE ---
p_medicine <- create_perc_barplot(df, "medicine", 
                                  "Distribuzione Farmaci",
                                  "Tipo di farmaco",
                                  order_by_freq = TRUE)
ggsave("outputs/figures/bar_medicine_perc.png", p_medicine, width = fig_width, height = fig_height, dpi = fig_dpi)

# --- GENDER ---
p_gender <- df %>%
  filter(!is.na(gender)) %>%
  count(gender) %>%
  mutate(
    perc = round(n / sum(n) * 100, 1),
    label = paste0(n, " (", perc, "%)")
  ) %>%
  ggplot(aes(x = gender, y = perc)) +
  geom_col(fill = fill_color, alpha = 0.8) +
  geom_text(aes(label = label), vjust = -0.5, size = 4) +
  labs(
    title = "Distribuzione per Genere",
    subtitle = paste0("n = ", sum(! is.na(df$gender)), 
                      " | NA = ", sum(is.na(df$gender))),
    x = "Genere",
    y = "Percentuale (%)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("outputs/figures/bar_gender_perc.png", p_gender, width = 6, height = fig_height, dpi = fig_dpi)

# --- MRI_EDSS_DIFF ---
p_mri <- create_perc_barplot(df, "mri_edss_diff",
                             "MRI-EDSS Differenza < 2 mesi",
                             "Risposta",
                             order_by_freq = FALSE)
ggsave("outputs/figures/bar_mri_edss_diff_perc.png", p_mri, width = 6, height = fig_height, dpi = fig_dpi)

# --- COMORBIDITY ---
p_comorbidity <- create_perc_barplot(df, "comorbidity",
                                     "Presenza di Comorbidità",
                                     "Comorbidità",
                                     order_by_freq = FALSE)
ggsave("outputs/figures/bar_comorbidity_perc.png", p_comorbidity, width = 6, height = fig_height, dpi = fig_dpi)

# --- EDSS_GROUP ---
p_edss_group <- df %>%
  filter(!is.na(edss_group)) %>%
  count(edss_group) %>%
  mutate(
    perc = round(n / sum(n) * 100, 1),
    label = paste0(n, " (", perc, "%)")
  ) %>%
  ggplot(aes(x = edss_group, y = perc)) +
  geom_col(fill = fill_color, alpha = 0.8) +
  geom_text(aes(label = label), vjust = -0.5, size = 4) +
  labs(
    title = "Distribuzione EDSS per Gruppo di Severità",
    subtitle = paste0("n = ", sum(!is.na(df$edss_group))),
    x = "Gruppo EDSS",
    y = "Percentuale (%)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("outputs/figures/bar_edss_group_perc.png", p_edss_group, width = fig_width, height = fig_height, dpi = fig_dpi)

# ----------------------------------------------------------------------------
# B3) Grafici prevalenze binarie
# ----------------------------------------------------------------------------

# --- Prevalenze cliniche ---
p_prev_clinical <- prevalence_clinical %>%
  mutate(feature = fct_reorder(feature, perc_1)) %>%
  ggplot(aes(x = feature, y = perc_1)) +
  geom_col(fill = fill_color, alpha = 0.8) +
  geom_text(aes(label = paste0(n_1)), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Prevalenza Indicatori Clinici Binari",
    subtitle = paste0("n = ", nrow(df), " pazienti"),
    x = "",
    y = "Prevalenza (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))

ggsave("outputs/figures/bar_prevalence_clinical_binary.png", p_prev_clinical, 
       width = 9, height = 8, dpi = fig_dpi)

# --- Prevalenze presenting symptoms ---
p_prev_presenting <- prevalence_presenting %>%
  mutate(
    feature = str_remove(feature, "^present_"),
    feature = str_to_title(feature),
    feature = fct_reorder(feature, perc_1)
  ) %>%
  ggplot(aes(x = feature, y = perc_1)) +
  geom_col(fill = fill_color_alt, alpha = 0.8) +
  geom_text(aes(label = paste0(n_1, " (", perc_1, "%)")), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "Prevalenza Sintomi di Presentazione",
    subtitle = paste0("n = ", nrow(df), " pazienti | Encoding multi-label"),
    x = "",
    y = "Prevalenza (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.15)))

ggsave("outputs/figures/bar_prevalence_presenting.png", p_prev_presenting, 
       width = 9, height = 6, dpi = fig_dpi)

# --- Fine script ---