# ============================================================================
# Script: 05_symptom_encoding.R
# Descrizione: Encoding multi-label della colonna symptom
# Autore: Daria Simonetti
# ============================================================================

source("scripts/01_setup.R")

input_path <- "outputs/data/dataset_ready.rds"

if (!file.exists(input_path)) {
  stop("File non trovato: ", input_path,)
}

df <- readRDS(input_path)

cat("Caricato:", input_path, "\n")
cat("   Righe:", nrow(df), "| Colonne:", ncol(df), "\n")


# 1. ESPLORAZIONE COLONNA SYMPTOM

if (!"symptom" %in% names(df)) {
  stop("Colonna 'symptom' non trovata nel dataset")
}

# Valori unici
unique_symptoms <- unique(as.character(df$symptom))
cat("Valori unici:", length(unique_symptoms), "\n")
cat("\n Elenco valori:\n")
for (s in sort(unique_symptoms)) {
  cat("   -", s, "\n")
}


# 2. DEFINIZIONE PATTERN PER OGNI SINTOMO

# Pattern regex (case-insensitive gestito dopo con str_to_lower)
# Usiamo word boundaries \\b per evitare match parziali

symptom_patterns <- list(
  present_motor = "\\bmotor[e]?\\b",
  present_sensory = "\\bsensory\\b",
  present_visual = "\\bvisual\\b",
  present_balance = "\\bbalance\\b",
  present_pain = "\\bpain\\b",
  present_behavioural = "\\bbehaviou?r(al)?\\b",
  present_fatigue = "\\bfatigue\\b",
  present_sexual = "\\bsexual\\b"
)

cat("Pattern definiti per", length(symptom_patterns), "sintomi:\n")
for (name in names(symptom_patterns)) {
  cat("   -", name, ":", symptom_patterns[[name]], "\n")
}

# 3. CREAZIONE COLONNE BINARIE MULTI-LABEL

# Conversione symptom in lowercase per matching robusto
df <- df %>%
  mutate(symptom_lower = str_to_lower(as.character(symptom)))

# Crea ogni colonna binaria
for (col_name in names(symptom_patterns)) {
  pattern <- symptom_patterns[[col_name]]
  
  df <- df %>%
    mutate(
      !!col_name := dplyr::case_when(
        is.na(symptom_lower) ~ NA_integer_,               # regola: NA -> NA
        str_detect(symptom_lower, pattern) ~ 1L,
        TRUE ~ 0L
      )
    )
  
  # Conta frequenze (prevalenze)
  n_positive <- sum(df[[col_name]], na.rm = TRUE) #frequenza assoluta
  denom <- sum(!is.na(df[[col_name]]))
  perc <- round(n_positive / nrow(df) * 100, 1)  #frequenza relativa percentuale
  cat(col_name, ":", n_positive, "casi (", perc, "%)\n")
}

# Rimuovi colonna temporanea
df <- df %>% select(-symptom_lower)


# 4. VERIFICA ENCODING

# Colonne create
new_cols <- names(symptom_patterns)

# alcune righe di esempio
cat("\n   Esempi di encoding:\n")
cat("   ------------------------------------------------------------------------\n")

df %>%
  select(id, symptom, all_of(new_cols)) %>%
  head(10) %>%
  print()

# Verifica righe senza alcun match (tutti 0)
df <- df %>%
  mutate(
    symptom_any_match = rowSums(across(all_of(new_cols)))
  )

no_match_rows <- df %>% filter(symptom_any_match == 0)

if (nrow(no_match_rows) > 0) {
  cat("\nRighe senza alcun match (", nrow(no_match_rows), "):\n")
  no_match_rows %>%
    select(id, symptom) %>%
    print()
} else {
  cat("\nTutte le righe hanno almeno un sintomo codificato\n")
}

# Rimuovi colonna temporanea
df <- df %>% select(-symptom_any_match)

# 5. TABELLA FREQUENZE

frequency_table <- tibble(
  symptom = new_cols,
  n = sapply(new_cols, function(col) sum(df[[col]], na.rm = TRUE)),
  total = nrow(df)
) %>%
  mutate(
    perc = round(n / total * 100, 2)
  ) %>%
  arrange(desc(perc)) %>%
  select(symptom, n, perc)

cat("\nFrequenza sintomi:\n")
print(as.data.frame(frequency_table), row.names = FALSE)

# Salva tabella
frequency_path <- "outputs/tables/symptom_frequence.csv"
write_csv(frequency_table, frequency_path)
cat("\nSalvato:", frequency_path, "\n")


# 6. SALVATAGGIO DATASET

output_path <- "outputs/data/data_complete.rds"
saveRDS(df, output_path)

cat("Salvato:", output_path, "(sovrascritto)\n")

# 7. RIEPILOGO FINALE

cat("\nDataset aggiornato:\n")
cat("   Righe:", nrow(df), "| Colonne:", ncol(df), "\n")

cat("\nNuove colonne create:\n")
for (col in new_cols) {
  cat("   -", col, "\n")
}

cat("\nTop 3 sintomi per frequenza:\n")
frequency_table %>%
  head(3) %>%
  mutate(display = paste0(symptom, ": ", n, " (", perc, "%)")) %>%
  pull(display) %>%
  paste("   -", .) %>%
  cat(sep = "\n")

# Glimpse finale
cat("\n\\nStruttura colonne symptom:\n")
cat("--------------------------------------------------------------------------------\n")

df %>%
  select(id, symptom, all_of(new_cols)) %>%
  glimpse()

