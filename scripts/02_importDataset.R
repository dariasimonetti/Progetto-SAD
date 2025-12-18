# ============================================================================
# Script: 02_importDataset.R
# Descrizione: Import dati e controlli qualità (solo audit, no cleaning)
# Autore:  Daria-Simonetti
# ============================================================================

source("scripts/01_setup.R")


# 1. IMPORT DATI

file_path <- "datasets/Patient_info_table1.xlsx"
if (! file.exists(file_path)) {
  stop("ERRORE: File non trovato:  ", file_path)
}

# pulizia nomi colonne
df <- read_excel(file_path, skip = 1) %>%
  clean_names()


# 2. SALVATAGGIO RDS

rds_path <- "outputs/data/raw_import.rds"
saveRDS(df, rds_path)


# 3. REPORT QUALITÀ DATI

log_path <- "outputs/logs/data_quality_report.txt"
sink(log_path, split = TRUE)  # split = TRUE stampa anche in console

cat("
--------------------------------------------------------------------------------
                         DATA QUALITY REPORT
                         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
--------------------------------------------------------------------------------

File sorgente: ", file_path, "

")

# --- 3.1 Dimensioni dataset ---
cat("
--------------------------------------------------------------------------------
1.  DIMENSIONI DATASET
--------------------------------------------------------------------------------
")
cat("Righe (osservazioni):", nrow(df), "\n")
cat("Colonne (variabili): ", ncol(df), "\n")
cat("\nNomi colonne:\n")
cat(paste(" -", names(df), collapse = "\n"), "\n")

# --- 3.2 Skim summary ---
cat("
--------------------------------------------------------------------------------
2. SUMMARY DATASET
--------------------------------------------------------------------------------
")

skim_output <- skim(df)
print(skim_output)

# --- 3.3 Percentuale NA per colonna ---
cat("
--------------------------------------------------------------------------------
3. PERCENTUALE VALORI MANCANTI (NA) PER COLONNA
--------------------------------------------------------------------------------
")

na_summary <- df %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "colonna", values_to = "n_na") %>%
  mutate(
    n_tot = nrow(df),
    perc_na = round(n_na / n_tot * 100, 2)
  ) %>%
  arrange(desc(perc_na))

print(as.data.frame(na_summary), row.names = FALSE)

cat("\nColonne con NA > 0%:\n")
na_with_missing <- na_summary %>% filter(n_na > 0)
if (nrow(na_with_missing) > 0) {
  print(as.data.frame(na_with_missing), row.names = FALSE)
} else {
  cat("   Nessun valore mancante!\n")
}

# --- 3.4 Conteggi livelli variabili categoriali ---
cat("
--------------------------------------------------------------------------------
4. CONTEGGI LIVELLI VARIABILI CATEGORIALI
--------------------------------------------------------------------------------
")

# Funzione helper per stampare conteggi
print_counts <- function(data, var_name) {
  cat("\n>>> ", var_name, ":\n", sep = "")
  
  if (var_name %in% names(data)) {
    counts <- data %>%
      count(!!sym(var_name), name = "n") %>%
      mutate(perc = round(n / sum(n) * 100, 2)) %>%
      arrange(desc(n))
    print(as.data.frame(counts), row.names = FALSE)
  } else {
    cat("Colonna non trovata nel dataset\n")
    # Cerca colonne simili
    similar <- names(data)[grepl(gsub("_", ".*", var_name), names(data), ignore.case = TRUE)]
    if (length(similar) > 0) {
      cat("   Colonne simili trovate:", paste(similar, collapse = ", "), "\n")
    }
  }
}

# Conteggi per genere, medicine e presenting symptoms
print_counts(df, "gender")
print_counts(df, "types_of_medicines")
print_counts(df, "presenting_symptom")

# Colonna MRI-EDSS < 2 months (cerca la versione clean_names)
# clean_names trasforma "MRI-EDSS < 2 months" in qualcosa come "mri_edss_2_months"
cat("\n>>> MRI-EDSS < 2 months:\n")

# Cerca la colonna con pattern mri.*edss
mri_edss_col <- names(df)[grepl("mri.*edss", names(df), ignore.case = TRUE)]

if (length(mri_edss_col) > 0) {
  cat("   (Nome colonna dopo clean_names: ", mri_edss_col[1], ")\n", sep = "")
  counts <- df %>%
    count(!!sym(mri_edss_col[1]), name = "n") %>%
    mutate(perc = round(n / sum(n) * 100, 2)) %>%
    arrange(desc(n))
  print(as.data.frame(counts), row.names = FALSE)
} else {
  cat("   ⚠️ Colonna MRI-EDSS non trovata\n")
  cat("   Colonne disponibili contenenti 'mri' o 'edss':\n")
  cat("   ", paste(names(df)[grepl("mri|edss", names(df), ignore.case = TRUE)], collapse = ", "), "\n")
}

# --- 3.5 Righe con age_of_onset > age ---
cat("
--------------------------------------------------------------------------------
5. ANOMALIE:  RIGHE CON age_of_onset > age
--------------------------------------------------------------------------------
")

# Cerca le colonne corrette (potrebbero avere nomi diversi)
age_col <- names(df)[grepl("^age$|^age_at", names(df), ignore.case = TRUE)][1]
onset_col <- names(df)[grepl("onset", names(df), ignore.case = TRUE)][1]

if (! is.na(age_col) && !is.na(onset_col)) {
  cat("Colonne utilizzate:  age = '", age_col, "', age_of_onset = '", onset_col, "'\n\n", sep = "")
  
  # Trova righe anomale
  anomalies <- df %>%
    mutate(row_id = row_number()) %>%
    filter(!! sym(onset_col) > !!sym(age_col))
  
  if (nrow(anomalies) > 0) {
    cat("⚠️ Trovate", nrow(anomalies), "righe con age_of_onset > age:\n\n")
    
    # Cerca colonna ID (potrebbe essere 'id', 'patient_id', 'cod', etc.)
    id_col <- names(df)[grepl("^id$|patient.*id|^cod", names(df), ignore.case = TRUE)]
    
    if (length(id_col) > 0) {
      anomalies_display <- anomalies %>%
        select(row_id, any_of(id_col), !!sym(age_col), !!sym(onset_col))
    } else {
      anomalies_display <- anomalies %>%
        select(row_id, !!sym(age_col), !!sym(onset_col))
    }
    
    print(as.data.frame(anomalies_display), row.names = FALSE)
  } else {
    cat("✅ Nessuna anomalia trovata (age_of_onset <= age per tutte le righe)\n")
  }
} else {
  cat("⚠️ Colonne 'age' o 'age_of_onset' non trovate\n")
  cat("Colonne disponibili:\n")
  cat(paste(" -", names(df), collapse = "\n"), "\n")
}

# --- 3.6 Conteggio gender == "N" ---
cat("
--------------------------------------------------------------------------------
6. ANOMALIE: RIGHE CON gender == 'N'
--------------------------------------------------------------------------------
")

if ("gender" %in% names(df)) {
  n_gender_n <- sum(df$gender == "N", na.rm = TRUE)
  
  cat("Conteggio righe con gender = 'N':", n_gender_n, "\n")
  
  if (n_gender_n > 0) {
    cat("Percentuale sul totale:", round(n_gender_n / nrow(df) * 100, 2), "%\n")
    
    # Mostra le righe
    gender_n_rows <- df %>%
      mutate(row_id = row_number()) %>%
      filter(gender == "N")
    
    cat("\nRighe interessate:\n")
    
    # Cerca colonna ID
    id_col <- names(df)[grepl("^id$|patient.*id|^cod", names(df), ignore.case = TRUE)]
    
    if (length(id_col) > 0) {
      print(as.data.frame(gender_n_rows %>% select(row_id, any_of(id_col), gender)), 
            row.names = FALSE)
    } else {
      print(as.data.frame(gender_n_rows %>% select(row_id, gender)), 
            row.names = FALSE)
    }
  }
} else {
  cat("⚠️ Colonna 'gender' non trovata\n")
}

# --- Fine report ---

sink()


cat("

Dataset importato: 
   - Righe:", nrow(df), "
   - Colonne:", ncol(df), "
")

summary(df)
