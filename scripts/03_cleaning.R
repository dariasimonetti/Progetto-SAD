# ============================================================================
# Script: 03_cleaning.R
# Descrizione: Cleaning, standardizzazione e rinomina colonne
# Autore:  Daria Simonetti
# ============================================================================

source("scripts/01_setup.R")


# 1. CARICAMENTO DATI

df <- readRDS("outputs/data/raw_import.rds")

cat("Righe:", nrow(df), "| Colonne:", ncol(df), "\n")


# 2. INIZIALIZZAZIONE CLEANING LOG

cleaning_log <- tibble(
  timestamp = as.POSIXct(character()),
  id = integer(),
  field = character(),
  old_value = character(),
  new_value = character(),
  rule = character()
)

# Funzione helper per aggiungere entry al log
add_log <- function(log, id = NA_integer_, field, old_value, new_value, rule) {
  new_entry <- tibble(
    timestamp = Sys.time(),
    id = as.integer(id),
    field = as.character(field),
    old_value = as.character(old_value),
    new_value = as.character(new_value),
    rule = as.character(rule)
  )
  bind_rows(log, new_entry)
}

# 3. REGOLA 1: gender "N" -> NA

# Trova righe con gender == "N" (case-insensitive)
gender_n_rows <- which(str_to_upper(df$gender) == "N")

if (length(gender_n_rows) > 0) {
  
  # Calcolo moda tra "F" e "M"
  gender_valid <- df$gender[str_to_upper(df$gender) %in% c("F", "M")]
  
  if (length(gender_valid)>0){
    # Conta frequenze
    freq_table <- table(str_to_upper(gender_valid))
    # Identifica la moda (più frequente)
    mode_value <- names(freq_table)[which.max(freq_table)]
    
    rule_desc <- paste0("Regola 1: gender 'N' -> mode(", mode_value, ")")
  } else {
    # Fallback: se non ci sono F/M validi, usa "F"
    mode_value <- "F"
    rule_desc <- "Regola 1: gender 'N' -> F (fallback, no F/M)"
    cat("Nessun valore F/M valido trovato - fallback a 'F'\n")
  }
  
  for (row_idx in gender_n_rows) {
    cleaning_log <- add_log(
      cleaning_log,
      id = df$id[row_idx],
      field = "gender",
      old_value = df$gender[row_idx],
      new_value = mode_value,
      rule = rule_desc
    )
  }
  
  # Applica la correzione
  df <- df %>%
    mutate(gender = if_else(str_to_upper(gender) == "N", mode_value, gender))
  
  cat("Convertite", length(gender_n_rows), "righe con gender='N' in '", 
      mode_value, "'\n", sep = "")
} else {
  cat("Nessuna riga con gender='N' trovata\n")
}


# 4. REGOLA 2:  Swap age/age_of_onset se onset > age


# Trova righe anomale
anomaly_rows <- which(df$age_of_onset > df$age)

if (length(anomaly_rows) > 0) {
  for (row_idx in anomaly_rows) {
    old_age <- df$age[row_idx]
    old_onset <- df$age_of_onset[row_idx]
    diff <- old_onset - old_age
    
    if (diff <= 5) {
      # Swap
      df$age[row_idx] <- old_onset
      df$age_of_onset[row_idx] <- old_age
      
      cleaning_log <- add_log(
        cleaning_log,
        id = df$id[row_idx],
        field = "age",
        old_value = as.character(old_age),
        new_value = as.character(old_onset),
        rule = paste0("Regola 2: swap age/onset (diff=", diff, ")")
      )
      
      cleaning_log <- add_log(
        cleaning_log,
        id = df$id[row_idx],
        field = "age_of_onset",
        old_value = as.character(old_onset),
        new_value = as.character(old_age),
        rule = paste0("Regola 2: swap age/onset (diff=", diff, ")")
      )
      
      cat("D", df$id[row_idx], ": swap age=", old_age, "<-> onset=", old_onset, "\n")
      
    } else {
      # Differenza > 5: imposta onset = NA
      df$age_of_onset[row_idx] <- NA_real_
      
      cleaning_log <- add_log(
        cleaning_log,
        id = df$id[row_idx],
        field = "age_of_onset",
        old_value = as.character(old_onset),
        new_value = "NA",
        rule = paste0("Regola 2: onset -> NA (diff=", diff, " > 5)")
      )
      
      cat("ID", df$id[row_idx], ":  onset=", old_onset, "-> NA (diff > 5)\n")
    }
  }
} else {
  cat("Nessuna riga con age_of_onset > age trovata\n")
}


# 5. REGOLA 3: EDSS numerico robusto


# Funzione per convertire EDSS (gestisce virgole e validazione)
clean_edss <- function(x) {
  # Se già numerico, controlla solo range e step
  if (is.numeric(x)) {
    result <- x
  } else {
    # Converti da stringa (gestisci virgola come separatore decimale)
    result <- as.numeric(str_replace(as.character(x), ",", "."))
  }
  
  # deve essere tra 0 e 10 e multiplo di 0.5
  valid <- !is.na(result) & result >= 0 & result <= 10 & (result %% 0.5 == 0)
  result[!valid & ! is.na(result)] <- NA_real_
  
  return(result)
}

# Salva valori originali per confronto
edss_original <- df$edss

# Applica cleaning
df$edss <- clean_edss(df$edss)

# Log delle modifiche
edss_changed <- which(
  (! is.na(edss_original) & is.na(df$edss)) |
    (! is.na(edss_original) & !is.na(df$edss) & as.character(edss_original) != as.character(df$edss))
)

if (length(edss_changed) > 0) {
  for (row_idx in edss_changed) {
    cleaning_log <- add_log(
      cleaning_log,
      id = df$id[row_idx],
      field = "edss",
      old_value = as.character(edss_original[row_idx]),
      new_value = as.character(df$edss[row_idx]),
      rule = "Regola 3: normalizzazione EDSS"
    )
  }
  cat("Normalizzati", length(edss_changed), "valori EDSS\n")
} else {
  cat("Tutti i valori EDSS già validi (numeric, 0-10, step 0.5)\n")
}

# Verifica finale
cat("Range EDSS:", min(df$edss, na.rm = TRUE), "-", max(df$edss, na.rm = TRUE), "\n")


# 6. REGOLA 4:  String cleaning

# str_squish() su tutte le colonne character
char_cols <- names(df)[sapply(df, is.character)]

df <- df %>%
  mutate(across(all_of(char_cols), str_squish))

cat("str_squish() applicato a:", paste(char_cols, collapse = ", "), "\n")

# Normalizza presenting_symptom
# - "Motore" -> "Motor" (case-insensitive)
# - Normalizza separatori (", " " & " " ," "&" ecc.)

if ("presenting_symptom" %in% names(df)) {
  
  # Salva originali per log
  symptom_original <- df$presenting_symptom
  
  df <- df %>%
    mutate(
      presenting_symptom = str_replace_all(
        presenting_symptom,
        regex("\\bMotore\\b", ignore_case = TRUE),
        "Motor"
      ),
      # Normalizza virgole e separatori: ", " o " ," o " & " con o senza spazi -> ", "
      presenting_symptom = str_replace_all(presenting_symptom, "\\s*&\\s*", ", "),
      presenting_symptom = str_replace_all(presenting_symptom, "\\s*,\\s*", ", ")
    )
  
  # Log modifiche
  symptom_changed <- which(symptom_original != df$presenting_symptom)
  
  if (length(symptom_changed) > 0) {
    for (row_idx in symptom_changed) {
      cleaning_log <- add_log(
        cleaning_log,
        id = df$id[row_idx],
        field = "presenting_symptom",
        old_value = symptom_original[row_idx],
        new_value = df$presenting_symptom[row_idx],
        rule = "Regola D: normalizzazione symptom"
      )
    }
    cat("Normalizzati", length(symptom_changed), "valori in presenting_symptom\n")
  }
}

# 7. RINOMINA COLONNE

# mapping rinomina
rename_map <- c(
  mri_edss_diff = "does_the_time_difference_between_mri_acquisition_and_edss_two_months",
  comorbidity = "dose_the_patient_has_co_moroidity",
  medicine = "types_of_medicines",
  symptom = "presenting_symptom",
  bowel_bladder = "bowel_and_bladder_function",
  optic_disc = "optic_discs",
  ocular_move = "ocular_movement"
)

# Log delle rinominazioni
for (new_name in names(rename_map)) {
  old_name <- rename_map[new_name]
  if (old_name %in% names(df)) {
    cleaning_log <- add_log(
      cleaning_log,
      id = NA_integer_,
      field = "column_name",
      old_value = old_name,
      new_value = new_name,
      rule = "Rinomina colonna"
    )
    cat("   ", old_name, " -> ", new_name, "\n")
  }
}

# Applica rinomina
df <- df %>%
  dplyr::rename(dplyr::any_of(rename_map))

# 8. GESTIONE COLONNA fields (NEAR-ZERO VARIANCE)

if ("fields" %in% names(df)) {
  # Verifica se è costante (tutti 0, senza NA)
  fields_values <- unique(df$fields[!is.na(df$fields)])
  fields_na_count <- sum(is.na(df$fields))
  
  is_constant <- length(fields_values) == 1 && fields_na_count == 0
  
  if (is_constant) {
    constant_value <- fields_values[1]
    
    cleaning_log <- add_log(
      cleaning_log,
      id = NA_integer_,
      field = "fields",
      old_value = paste0("constant (all ", constant_value, ")"),
      new_value = "REMOVED",
      rule = "Near-zero variance:  colonna costante rimossa"
    )
    
    # Rimuovi la colonna
    df <- df %>% select(-fields)
    
    cat("Colonna 'fields' rimossa:  valore costante =", constant_value, "\n")
    cat("Registrato nel cleaning log\n")
  } else {
    cat("Colonna 'fields' non è costante, mantenuta nel dataset\n")
  }
} else {
  cat("Colonna 'fields' non presente nel dataset\n")
}


# 9. CONVERSIONE A FACTOR


# gender:  factor con livelli c("F", "M")
df <- df %>%
  mutate(
    gender = factor(gender, levels = c("F", "M"))
    # Valori non in c("F", "M") diventano automaticamente NA
  )

# mri_edss_diff: factor con livelli c("No", "Yes")
df <- df %>%
  mutate(
    # Normalizza case prima di convertire
    mri_edss_diff = str_to_title(mri_edss_diff),
    mri_edss_diff = factor(mri_edss_diff, levels = c("No", "Yes"))
  )


# medicine: factor ordinato per frequenza
df <- df %>%
  mutate(
    medicine = fct_infreq(factor(medicine))
  )


# 9.4 symptom: factor (ordine alfabetico)
df <- df %>%
  mutate(
    symptom = factor(symptom)
  )


# comorbidity: factor con livelli c("No", "Yes")
df <- df %>%
  mutate(
    comorbidity = str_to_title(comorbidity),
    comorbidity = factor(comorbidity, levels = c("No", "Yes"))
  )


# 10. SALVATAGGIO CLEANING LOG

log_output_path <- "outputs/logs/cleaning_log.csv"
write_csv(cleaning_log, log_output_path)

cat("Salvato:", log_output_path, "\n")
cat("Totale operazioni registrate:", nrow(cleaning_log), "\n")


# 11. SALVATAGGIO DATASET PULITO

output_path <- "outputs/data/cleaned.rds"
saveRDS(df, output_path)


# 12. CONTROLLI FINALI E RIEPILOGO


# Glimpse del dataset
cat("\nStruttura dataset finale:\n")
cat("--------------------------------------------------------------------------------\n")
glimpse(df)

# Conteggi livelli factor
cat("\n\nConteggi livelli variabili factor:\n")
cat("--------------------------------------------------------------------------------\n")

cat("\n>>> gender:\n")
print(table(df$gender, useNA = "ifany"))

cat("\n>>> medicine:\n")
print(table(df$medicine, useNA = "ifany"))

cat("\n>>> mri_edss_diff:\n")
print(table(df$mri_edss_diff, useNA = "ifany"))

cat("\n>>> comorbidity:\n")
print(table(df$comorbidity, useNA = "ifany"))

cat("\n>>> symptom (primi 10 livelli):\n")
symptom_counts <- sort(table(df$symptom), decreasing = TRUE)
print(head(symptom_counts, 10))

# Verifica anomalia age_of_onset > age
cat("\n\nVerifica finale:  righe con age_of_onset > age:\n")
cat("--------------------------------------------------------------------------------\n")

anomalies_post <- df %>% filter(age_of_onset > age)

if (nrow(anomalies_post) == 0) {
  cat("NESSUNA anomalia trovata.  Tutte le righe hanno age_of_onset <= age\n")
} else {
  cat("ATTENZIONE: Trovate", nrow(anomalies_post), "righe anomale:\n")
  print(anomalies_post %>% select(id, age, age_of_onset))
}

# Riepilogo cleaning log
cat("\n\nRiepilogo Cleaning Log:\n")
cat("--------------------------------------------------------------------------------\n")

log_summary <- cleaning_log %>%
  count(rule, name = "n_operazioni") %>%
  arrange(desc(n_operazioni))

print(as.data.frame(log_summary), row.names = FALSE)
