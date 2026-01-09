# ============================================================================
# Script: 04_feature_engineering.R
# Descrizione: Feature engineering - creazione nuove variabili
# Autore: Daria Simonetti
# ============================================================================

source("scripts/01_setup.R")


# CARICAMENTO DATI PULITI

df <- readRDS("outputs/data/cleaned.rds")
cat("   Righe:", nrow(df), "| Colonne:", ncol(df), "\n")


# 1. CREAZIONE VARIABILE DISEASE DURATION


df <- df %>%
  mutate(
    disease_duration = age - age_of_onset,
    # Se < 0 (non dovrebbe succedere dopo cleaning), imposta NA
    disease_duration = if_else(disease_duration < 0, NA_real_, disease_duration)
  )

# Statistiche
cat("Disease_duration = age - age_of_onset\n")
cat("   Range:", min(df$disease_duration, na.rm = TRUE), "-", 
    max(df$disease_duration, na.rm = TRUE), "anni\n")
cat("   Media:", round(mean(df$disease_duration, na.rm = TRUE), 2), "anni\n")
cat("   NA:", sum(is.na(df$disease_duration)), "\n")


# 2. IDENTIFICAZIONE COLONNE BINARIE 0/1

# Colonne da escludere dall'identificazione binaria
exclude_cols <- c("id", "age", "age_of_onset", "edss", "disease_duration")

# Funzione per verificare se una colonna è binaria (solo 0 e 1)
is_binary_01 <- function(x) {
  # Deve essere numeric
  if (! is.numeric(x)) return(FALSE)
  
  # Valori unici (escludendo NA)
  unique_vals <- unique(x[!is.na(x)])
  
  
  # Deve contenere solo 0 e/o 1
  all(unique_vals %in% c(0, 1))
}

# Identifica colonne binarie
binary_cols <- names(df)[sapply(df, is_binary_01)]

# Escludi le colonne specificate
binary_cols <- setdiff(binary_cols, exclude_cols)

cat("trovate", length(binary_cols), "colonne binarie:\n")
cat("   ", paste(binary_cols, collapse = ", "), "\n")


# 3. CREAZIONE BINARY_SUM (somma delle colonne binarie)


# Conta NA per riga nelle colonne binarie
df <- df %>%
  mutate(
    binary_na_count = rowSums(is.na(across(all_of(binary_cols))))
  )

# Regola per gestione NA:
# - Se più del 50% delle colonne binarie è NA, binary_sum = NA
# - Altrimenti, calcola somma con na.rm = TRUE
na_threshold <- length(binary_cols) * 0.5

df <- df %>%
  mutate(
    binary_sum = if_else(
      binary_na_count > na_threshold,
      NA_integer_,
      as.integer(rowSums(across(all_of(binary_cols)), na.rm = TRUE))
    )
  ) %>%
  select(-binary_na_count)  # Rimuovi colonna temporanea

cat("Binary_sum = somma di", length(binary_cols), "colonne binarie\n")
cat("Regola NA:  se >50% colonne binarie sono NA, binary_sum = NA\n")
cat("Range:", min(df$binary_sum, na.rm = TRUE), "-", 
    max(df$binary_sum, na.rm = TRUE), "\n")
cat("   Media:", round(mean(df$binary_sum, na.rm = TRUE), 2), "\n")

# Distribuzione
cat("\n   Distribuzione binary_sum:\n")
print(table(df$binary_sum, useNA = "ifany"))


# 4. CREAZIONE VARIABILI EDSS_BIN E EDSS_3CLASS

# Definizione dei cut-off per EDSS_BIN
# Cut-off professore: <=2.0 vs >2.0

df <- df %>%
  mutate(
    edss_bin = case_when(
      is.na(edss) ~ NA_character_,
      edss <= 2.0 ~ "class0 (<=2.0)",
      edss > 2.0 ~ "class1 (>2.0)"
    ),
    edss_bin = factor(edss_bin, levels = c("class0 (<=2.0)", "class1 (>2.0)"))
  )

cat("edss_bin:  class0 (<=2.0) vs class1 (>2.0)\n")
cat("\n   Distribuzione edss_bin:\n")
print(table(df$edss_bin, useNA = "ifany"))

# Definizione dei cut-off per EDSS_3CLASS
# Cut-off professore: normal (0-2.0), mild (2.5-4.0), severe (>4.0)

df <- df %>%
  mutate(
    edss_3class = case_when(
      is.na(edss) ~ NA_character_,
      edss <= 2.0 ~ "normal (0-2.0)",
      edss > 2.0 & edss <= 4.0 ~ "mild (2.5-4.0)",
      edss > 4.0 ~ "severe (>4.0)"
    ),
    edss_3class = factor(edss_3class, levels = c("normal (0-2.0)", "mild (2.5-4.0)", "severe (>4.0)"))
  )

cat("edss_3class: normal (0-2.0), mild (2.5-4.0), severe (>4.0)\n")
cat("\n   Distribuzione edss_3class:\n")
print(table(df$edss_3class, useNA = "ifany"))

# 5. RIEPILOGO NUOVE VARIABILI

new_vars <- c("disease_duration", "binary_sum", "edss_bin", "edss_3class")

for (var in new_vars) {
  cat("\n>>>", var, ":\n")
  
  if (is.factor(df[[var]])) {
    print(table(df[[var]], useNA = "ifany"))
  } else {
    print(summary(df[[var]]))
  }
}


# 6. SALVATAGGIO DATASET


output_path <- "outputs/data/dataset_ready.rds"
saveRDS(df, output_path)


# 7. SALVATAGGIO TABELLE TECNICHE E REPORT

# --- Lista colonne binarie ---
binary_cols_df <- tibble(feature = binary_cols)
write_csv(binary_cols_df, "outputs/tables/binary_cols.csv")

# --- Distribuzione binary_sum ---
binary_sum_counts <- df %>%
  filter(! is.na(binary_sum)) %>%
  count(binary_sum, name = "n") %>%
  mutate(perc = round(n / sum(n) * 100, 2))

write_csv(binary_sum_counts, "outputs/tables/binary_sum_counts.csv")

# --- Distribuzione edss_bin ---
edss_bin_counts <- df %>%
  filter(! is.na(edss_bin)) %>%
  count(edss_bin, name = "n") %>%
  mutate(
    edss_bin = as.character(edss_bin),
    perc = round(n / sum(n) * 100, 2)
  )

write_csv(edss_bin_counts, "outputs/tables/edss_bin_counts.csv")

# --- Distribuzione edss_3class ---
edss_3class_counts <- df %>%
  filter(!is.na(edss_3class)) %>%
  count(edss_3class, name = "n") %>%
  mutate(
    edss_3class = as.character(edss_3class),
    perc = round(n / sum(n) * 100, 2)
  )

write_csv(edss_3class_counts, "outputs/tables/edss_3class_counts.csv")

# ---Mini-report testuale ---

report_path <- "outputs/logs/feature_engineering_report.txt"

# Calcola statistiche per il report
dd_range <- range(df$disease_duration, na.rm = TRUE)
dd_mean <- mean(df$disease_duration, na.rm = TRUE)
dd_na <- sum(is.na(df$disease_duration))

bs_range <- range(df$binary_sum, na.rm = TRUE)
bs_mean <- mean(df$binary_sum, na.rm = TRUE)
bs_na <- sum(is.na(df$binary_sum))

edss_bin_table <- table(df$edss_bin, useNA = "ifany")
edss_3class_table <- table(df$edss_3class, useNA = "ifany")

# Scrivi report
sink(report_path)

cat("================================================================================\n")
cat("                    FEATURE ENGINEERING REPORT\n")
cat("================================================================================\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Dataset:  outputs/data/analysis_ready.rds\n")
cat("Righe:", nrow(df), "| Colonne:", ncol(df), "\n")
cat("================================================================================\n\n")

cat("--- DISEASE_DURATION ---\n")
cat("Range:", dd_range[1], "-", dd_range[2], "anni\n")
cat("Media:", round(dd_mean, 2), "anni\n")
cat("NA:", dd_na, "(", round(dd_na / nrow(df) * 100, 2), "%)\n\n")

cat("--- BINARY_SUM ---\n")
cat("Range:", bs_range[1], "-", bs_range[2], "\n")
cat("Media:", round(bs_mean, 2), "\n")
cat("NA:", bs_na, "(", round(bs_na / nrow(df) * 100, 2), "%)\n\n")

cat("--- EDSS_BIN (Binary Classification) ---\n")
cat("Cut-off: <=2.0 vs >2.0\n")
for (i in seq_along(edss_bin_table)) {
  level_name <- names(edss_bin_table)[i]
  level_n <- edss_bin_table[i]
  level_perc <- round(level_n / sum(edss_bin_table) * 100, 2)
  cat(level_name, ":", level_n, "(", level_perc, "%)\n")
}
cat("\n")

cat("--- EDSS_3CLASS (Multi-class Classification) ---\n")
cat("Cut-off: normal (0-2.0), mild (2.5-4.0), severe (>4.0)\n")
for (i in seq_along(edss_3class_table)) {
  level_name <- names(edss_3class_table)[i]
  level_n <- edss_3class_table[i]
  level_perc <- round(level_n / sum(edss_3class_table) * 100, 2)
  cat(level_name, ":", level_n, "(", level_perc, "%)\n")
}
cat("\n")

cat("--- COLONNE BINARIE IDENTIFICATE ---\n")
cat("Numero:", length(binary_cols), "\n")
cat("Elenco:\n")
cat(paste(" -", binary_cols, collapse = "\n"), "\n")

cat("\n================================================================================\n")
cat("                         FINE REPORT\n")
cat("================================================================================\n")

sink()

# 8. RIEPILOGO FINALE

# Dataset finale
cat("\nStruttura dataset finale:\n")
cat("--------------------------------------------------------------------------------\n")
cat("Righe:", nrow(df), "| Colonne:", ncol(df), "\n\n")

# Glimpse
glimpse(df)

# Head
cat("\n\nPrime 6 righe (variabili chiave):\n")
cat("--------------------------------------------------------------------------------\n")

df %>%
  select(id, age, age_of_onset, disease_duration, edss, edss_bin, edss_3class, binary_sum) %>%
  head() %>%
  print()

# Summary delle nuove variabili
cat("\n\nSummary nuove variabili:\n")
cat("--------------------------------------------------------------------------------\n")

df %>%
  select(all_of(new_vars)) %>%
  summary() %>%
  print()

# Lista colonne binarie identificate
cat("\n\nColonne binarie identificate (", length(binary_cols), "):\n")
cat("--------------------------------------------------------------------------------\n")
cat(paste(binary_cols, collapse = "\n"), "\n")
