# ============================================================================
# Script: 09_data_dictionary_and_tables.R
# Descrizione: Data dictionary automatico + Table 1 descrittiva
# Autore:  Daria-Simonetti
# ============================================================================

source("scripts/01_setup.R")

# Carica dati
df <- readRDS("outputs/data/data_complete.rds")

cat("   Righe:", nrow(df), "| Colonne:", ncol(df), "\n")

# Contatori file salvati
saved_csv <- character()


# DATA DICTIONARY AUTOMATICO

# Funzione per ottenere la classe principale
get_main_class <- function(x) {
  cls <- class(x)
  # Restituisci la prima classe (più specifica)
  cls[1]
}

# Funzione per generare example_values
get_example_values <- function(x, max_examples = 5) {
  
  # Rimuovi NA per gli esempi
  x_clean <- x[! is.na(x)]
  
  if (length(x_clean) == 0) {
    return("(all NA)")
  }
  
  if (is.numeric(x)) {
    # Per numeriche: valori unici arrotondati
    unique_vals <- unique(round(x_clean, 2))
    unique_vals <- sort(unique_vals)
    
    if (length(unique_vals) > max_examples) {
      # Mostra range + alcuni esempi
      examples <- c(
        head(unique_vals, 2),
        "...",
        tail(unique_vals, 2)
      )
    } else {
      examples <- unique_vals
    }
    
    return(paste(examples, collapse = "; "))
    
  } else if (is.factor(x) || is.character(x)) {
    # Per factor/character: livelli più comuni
    freq_table <- sort(table(x_clean), decreasing = TRUE)
    top_values <- names(head(freq_table, max_examples))
    
    if (length(freq_table) > max_examples) {
      return(paste(c(top_values, "..."), collapse = "; "))
    } else {
      return(paste(top_values, collapse = "; "))
    }
    
  } else if (is.logical(x)) {
    return(paste(unique(x_clean), collapse = "; "))
    
  } else {
    return(paste(head(unique(as.character(x_clean)), max_examples), collapse = "; "))
  }
}

# Costruisci data dictionary
data_dictionary <- tibble(
  variable = names(df)
) %>%
  rowwise() %>%
  mutate(
    class = get_main_class(df[[variable]]),
    n_missing = sum(is.na(df[[variable]])),
    pct_missing = round(n_missing / nrow(df) * 100, 1),
    n_unique = n_distinct(df[[variable]], na.rm = TRUE),
    example_values = get_example_values(df[[variable]])
  ) %>%
  ungroup() %>%
  arrange(variable)

# Salva
write_csv(data_dictionary, "outputs/tables/data_dictionary.csv")
saved_csv <- c(saved_csv, "outputs/tables/data_dictionary.csv")

cat("   Variabili totali:", nrow(data_dictionary), "\n")

# Riepilogo classi
cat("\nRiepilogo per classe:\n")
class_summary <- data_dictionary %>%
  count(class, name = "n_vars") %>%
  arrange(desc(n_vars))
print(as.data.frame(class_summary), row.names = FALSE)


# TABLE 1 DESCRITTIVA DEL CAMPIONE

# Inizializza tabella
table1 <- tibble(
  variable = character(),
  level = character(),
  value = character()
)

# --- Helper funct ---

# Funzione per formattare mediana [IQR]
format_median_iqr <- function(x, digits = 1) {
  x_clean <- x[!is.na(x)]
  
  if (length(x_clean) == 0) {
    return("NA")
  }
  
  med <- round(median(x_clean), digits)
  q1 <- round(quantile(x_clean, 0.25), digits)
  q3 <- round(quantile(x_clean, 0.75), digits)
  
  paste0(med, " [", q1, "–", q3, "]")
}

# Funzione per formattare conteggi e percentuali
format_n_pct <- function(n, total, digits = 1) {
  pct <- round(n / total * 100, digits)
  paste0(n, " (", pct, "%)")
}

# Funzione per aggiungere variabile categorica a table1
add_categorical <- function(table1, data, var_name, var_label = NULL) {
  
  if (is.null(var_label)) var_label <- var_name
  
  # Conteggi per livello
  counts <- data %>%
    count(!! sym(var_name), name = "n") %>%
    rename(level = 1) %>%
    mutate(level = as.character(level))
  
  # Gestisci NA come livello separato
  n_na <- sum(is.na(data[[var_name]]))
  
  if (n_na > 0) {
    counts <- counts %>%
      filter(! is.na(level)) %>%
      bind_rows(tibble(level = "Missing", n = n_na))
  }
  
  # Calcola percentuali
  total <- nrow(data)
  
  new_rows <- counts %>%
    mutate(
      variable = var_label,
      value = format_n_pct(n, total)
    ) %>%
    select(variable, level, value)
  
  bind_rows(table1, new_rows)
}

# --- Costruisci Table 1 ---

N <- nrow(df)

# 1. N totale
table1 <- bind_rows(table1, tibble(
  variable = "N",
  level = "",
  value = as.character(N)
))

# 2. Age
table1 <- bind_rows(table1, tibble(
  variable = "Age, years",
  level = "",
  value = format_median_iqr(df$age)
))

# 2. Age of onset
table1 <- bind_rows(table1, tibble(
  variable = "Age of onset, years",
  level = "",
  value = format_median_iqr(df$age_of_onset)
))

# Aggiungi missing se presenti
n_missing_age <- sum(is.na(df$age))
if (n_missing_age > 0) {
  table1 <- bind_rows(table1, tibble(
    variable = "Age, years",
    level = "Missing",
    value = format_n_pct(n_missing_age, N)
  ))
}

# 3. Disease duration
table1 <- bind_rows(table1, tibble(
  variable = "Disease duration, years",
  level = "",
  value = format_median_iqr(df$disease_duration)
))

n_missing_dd <- sum(is.na(df$disease_duration))
if (n_missing_dd > 0) {
  table1 <- bind_rows(table1, tibble(
    variable = "Disease duration, years",
    level = "Missing",
    value = format_n_pct(n_missing_dd, N)
  ))
}

# 4. EDSS
table1 <- bind_rows(table1, tibble(
  variable = "EDSS",
  level = "",
  value = format_median_iqr(df$edss)
))

n_missing_edss <- sum(is.na(df$edss))
if (n_missing_edss > 0) {
  table1 <- bind_rows(table1, tibble(
    variable = "EDSS",
    level = "Missing",
    value = format_n_pct(n_missing_edss, N)
  ))
}

# 5. Gender
table1 <- add_categorical(table1, df, "gender", "Gender")

# 6. Medicine
table1 <- add_categorical(table1, df, "medicine", "Medicine")

# 7. MRI EDSS diff
if ("mri_edss_diff" %in% names(df)) {
  table1 <- add_categorical(table1, df, "mri_edss_diff", "MRI-EDSS < 2 months")
}

# 8. Comorbidity
table1 <- add_categorical(table1, df, "comorbidity", "Comorbidity")

# Salva
write_csv(table1, "outputs/tables/table1_demographics.csv")
saved_csv <- c(saved_csv, "outputs/tables/table1_demographics.csv")

cat("   Righe in Table 1:", nrow(table1), "\n")


# C) OUTPUT CONSOLE FINALE


cat("\nFile salvati:\n")
for (f in saved_csv) {
  cat(f, "\n")
}

cat("\nData Dictionary (head):\n")
cat("--------------------------------------------------------------------------------\n")
print(head(data_dictionary, 10))

cat("\n\nTable 1 - Caratteristiche del campione:\n")
cat("--------------------------------------------------------------------------------\n")
print(as.data.frame(table1), row.names = FALSE)

cat("

Note: 
   - Il data dictionary contiene metadati per tutte le variabili
   - Table 1 è formattata per inclusione diretta nella relazione
   - Le variabili numeriche sono riportate come mediana [IQR]
   - Le variabili categoriche mostrano n (%)
   - I missing sono esplicitamente indicati dove presenti

")
