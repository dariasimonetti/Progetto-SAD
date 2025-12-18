# ============================================================================
# Script: 01_setup.R
# Descrizione: Setup iniziale del progetto
# Autore: Daria Simonetti
# ============================================================================

cat("\n Working directory:", getwd(), "\n\n")

# --- Creazione struttura cartelle ---
directories <- c(
  "datasets",
  "scripts",
  "outputs/data",
  "outputs/figures",
  "outputs/tables",
  "outputs/logs"
)

for (dir in directories) {
  if (! dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Creata:", dir, "\n")
  }
}

# --- Caricamento pacchetti ---
packages <- c(
  "tidyverse", "readxl", "janitor", "skimr", "naniar",
  "stringr", "lubridate", "ggplot2", "patchwork", 
  "psych", "broom", "ggcorrplot"
)

load_pkg <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(packages, load_pkg))
cat("\nPacchetti caricati\n")

# --- Tema ggplot2 globale ---
theme_set(
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
)

# Palette colori
project_colors <- c("#a2d2ff", "#fec89a", "#ccd5ae", "#eac4d5", 
                    "#cfbaf0", "#ff928b", "#809bce", "#d5bdaf")

# --- Opzioni globali ---
options(scipen = 999)

cat("\nSetup completato!\n")