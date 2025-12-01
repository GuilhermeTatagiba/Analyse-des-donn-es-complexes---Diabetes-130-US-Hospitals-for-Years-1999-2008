# Project: Diabetes Hospital Dataset Analysis
# Authors: Rim El Fatihi & Yuetong Lu & Guilherme Peres Tatagiba

rm(list=ls())
#setwd("C:/Users/guipe/Modelisation statistique/Analyse-des-donn-es-complexes---Diabetes-130-US-Hospitals-for-Years-1999-2008")
#getwd()

# Chargement des bibliothèques
library(dplyr)
library(ggplot2)
library(VIM)         # kNN
library(mice)        # mice (EM)
library(forecast)    # Forecasting
library(caret)       # Préprocessing

# 1) Chargement des données
diabetes <- read.csv("dataset/diabetic_data.csv", na.strings = c("", "NA", "?"))
cat("Dimensions initiales :", dim(diabetes), "\n")
gc()

# 2) Prétraitement initial
cat("\n== Prétraitement initial ==\n")

# 2.1 Supprimer 'weight' (trop de NA)
if("weight" %in% names(diabetes)) {
  diabetes$weight <- NULL
  cat("-> colonne 'weight' supprimée\n")
}

# 2.2 Convertir chaînes en facteurs et nettoyer entrées vides
for(v in names(diabetes)) {
  if(is.character(diabetes[[v]])) {
    diabetes[[v]] <- trimws(diabetes[[v]])
    diabetes[[v]][diabetes[[v]] == ""] <- NA
    diabetes[[v]] <- as.factor(diabetes[[v]])
  }
}

# 2.3 Remplacer NA massifs par catégories explicites (utile pour kNN/mice)
if("payer_code" %in% names(diabetes)) {
  diabetes$payer_code <- addNA(as.factor(diabetes$payer_code))
  levels(diabetes$payer_code)[is.na(levels(diabetes$payer_code))] <- "Unknown_payer"
  cat("-> payer_code : NA => 'Unknown_payer'\n")
}
if("medical_specialty" %in% names(diabetes)) {
  diabetes$medical_specialty <- addNA(as.factor(diabetes$medical_specialty))
  levels(diabetes$medical_specialty)[is.na(levels(diabetes$medical_specialty))] <- "Unknown_specialty"
  cat("-> medical_specialty : NA => 'Unknown_specialty'\n")
}

gc()

# 3) Regroupement des codes diag_1/diag_2/diag_3 (ICD-9) en macro-catégories
cat("\n== Regroupement ICD-9 pour diag_1/2/3 ==\n")
group_diag_icd9 <- function(code) {
  if(is.na(code) || code == "") return(NA)
  s <- as.character(code)
  s <- trimws(s)
  # Codes commençant par 'V' ou 'E'
  if(grepl("^[Vv]", s)) return("V_codes")
  if(grepl("^[Ee]", s)) return("E_codes")
  # essayer extraire partie numérique (premiers 3 chars)
  num <- suppressWarnings(as.numeric(substr(s, 1, 3)))
  if(is.na(num)) return(NA)
  if (num >= 1   & num <= 139) return("infectious")
  if (num >= 140 & num <= 239) return("neoplasms")
  if (num >= 240 & num <= 279) return("endocrine")
  if (num >= 280 & num <= 289) return("blood")
  if (num >= 290 & num <= 319) return("mental")
  if (num >= 320 & num <= 389) return("nervous")
  if (num >= 390 & num <= 459) return("circulatory")
  if (num >= 460 & num <= 519) return("respiratory")
  if (num >= 520 & num <= 579) return("digestive")
  if (num >= 580 & num <= 629) return("genitourinary")
  if (num >= 630 & num <= 679) return("pregnancy")
  if (num >= 680 & num <= 709) return("skin")
  if (num >= 710 & num <= 739) return("musculoskeletal")
  if (num >= 740 & num <= 759) return("congenital")
  if (num >= 760 & num <= 779) return("perinatal")
  if (num >= 780 & num <= 799) return("symptoms")
  if (num >= 800 & num <= 999) return("injury")
  return("other")
}

for(col in c("diag_1","diag_2","diag_3")) {
  if(col %in% names(diabetes)) {
    new_col <- paste0(col, "_group")
    diabetes[[new_col]] <- factor(sapply(as.character(diabetes[[col]]), group_diag_icd9))
    cat("-> créé :", new_col, "avec", length(levels(diabetes[[new_col]])), "niveaux\n")
  }
}
gc()

# 4) Vérifier variables avec NA après preprocessing
missing_summary <- sapply(diabetes, function(x) sum(is.na(x)))
missing_vars <- names(missing_summary[missing_summary > 0])
cat("\nVariables avec NA (après preprocessing) :\n")
print(missing_summary[missing_summary > 0])
gc()

# 5) IMPUTATION 1 : kNN optimisé (variable-par-variable)
cat("\n== IMPUTATION kNN (variable par variable, imp_var=FALSE) ==\n")

diab_knn <- diabetes

#Suppression : 

diab_knn$diag_1 <- NULL
diab_knn$diag_2 <- NULL
diab_knn$diag_3 <- NULL

# Choisir variables numériques robustes pour distance
numeric_candidates <- intersect(names(diab_knn)[sapply(diab_knn, is.numeric)],
                                c("time_in_hospital","num_lab_procedures","num_procedures","num_medications","number_diagnoses"))
if(length(numeric_candidates) < 1) {
  numeric_candidates <- names(diab_knn)[sapply(diab_knn, is.numeric)][1:min(3, sum(sapply(diab_knn, is.numeric)))]
}
dist_vars <- numeric_candidates
cat("Variables utilisées pour la distance kNN :", paste(dist_vars, collapse = ", "), "\n")

# recalculer liste NA car on a créé des diag_group
missing_summary <- sapply(diab_knn, function(x) sum(is.na(x)))
vars_to_impute <- names(missing_summary[missing_summary > 0])
cat("Nombre de variables à imputer (kNN loop):", length(vars_to_impute), "\n")
print(vars_to_impute)

for(var in vars_to_impute) {
  cat("kNN -> imputer :", var, " (NA:", sum(is.na(diab_knn[[var]])), ")\n")
  # sauter si déjà plein
  if(sum(is.na(diab_knn[[var]])) == 0) next
  # préparer subset minimal
  cols_for_kNN <- unique(c(var, dist_vars))
  subset_df <- diab_knn[, cols_for_kNN, drop = FALSE]
  tryCatch({
    out <- kNN(subset_df, variable = var, k = 5, imp_var = FALSE, dist_var = dist_vars, weightDist = FALSE)
    # remplacer colonne entière (s'il s'agit d'un factor, kNN garde le type)
    diab_knn[[var]] <- out[[var]]
    cat("  -> kNN succès\n")
    gc()
  }, error = function(e) {
    cat("  -> kNN erreur sur", var, ":", conditionMessage(e), " -> fallback median/mode\n")
    if(is.numeric(diab_knn[[var]])) {
      diab_knn[[var]][is.na(diab_knn[[var]])] <- median(diab_knn[[var]], na.rm = TRUE)
    } else {
      moda <- names(sort(table(diab_knn[[var]]), decreasing = TRUE))[1]
      diab_knn[[var]][is.na(diab_knn[[var]])] <- moda
      diab_knn[[var]] <- factor(diab_knn[[var]])
    }
    gc()
  })
}
cat("NA restants après kNN :", sum(is.na(diab_knn)), "\n")
saveRDS(diab_knn, "diabetes_after_knn.rds")
gc()

#############################################
# 6) IMPUTATION 2 : EM avec mice (norm + polyreg)
#############################################

# Copier le dataset kNN
diab_em <- diabetes

# 🚨 Supprimer les anciennes variables diag_* (trop de catégories)
diab_em$diag_1 <- NULL
diab_em$diag_2 <- NULL
diab_em$diag_3 <- NULL

# S'assurer que toutes les variables qualitatives sont des facteurs
for (v in names(diab_em)) {
  if (is.character(diab_em[[v]])) {
    diab_em[[v]] <- factor(diab_em[[v]])
  }
}

# Méthodes automatiques
methods_em <- make.method(diab_em)

# Attribution manuelle : numérique → norm, factor → polyreg
for (v in names(diab_em)) {
  if (any(is.na(diab_em[[v]]))) {
    if (is.numeric(diab_em[[v]])) {
      methods_em[v] <- "norm"
    } else {
      methods_em[v] <- "cart"
    }
  } else {
    methods_em[v] <- ""
  }
}

# Matrice de prédicteurs : tout utiliser sauf ID
predM_em <- make.predictorMatrix(diab_em)
predM_em[, "encounter_id"] <- 0
predM_em["encounter_id", ] <- 0

# Imputation EM
imp_em <- mice(diab_em, method = methods_em, predictorMatrix = predM_em,
               m = 1, maxit = 6, print = TRUE)

diab_em_final <- complete(imp_em)


#############################################
# 7) IMPUTATION 3 : Random Forest (mice, method = "rf")
#############################################

diab_rf <- diabetes

# 🚨 Supprimer diag_1/2/3 (trop de catégories)
diab_rf$diag_1 <- NULL
diab_rf$diag_2 <- NULL
diab_rf$diag_3 <- NULL

# Conversion des variables qualitatives en facteurs
for (v in names(diab_rf)) {
  if (is.character(diab_rf[[v]])) {
    diab_rf[[v]] <- factor(diab_rf[[v]])
  }
}

# Méthodes RF
methods_rf <- make.method(diab_rf)

for (v in names(diab_rf)) {
  if (any(is.na(diab_rf[[v]]))) {
    methods_rf[v] <- "rf"
  } else {
    methods_rf[v] <- ""
  }
}

# Matrice de prédicteurs
predM_rf <- make.predictorMatrix(diab_rf)
predM_rf[, "encounter_id"] <- 0
predM_rf["encounter_id", ] <- 0

# Imputation Random Forest
imp_rf <- mice(diab_rf, method = methods_rf, predictorMatrix = predM_rf,
               m = 1, maxit = 5, print = TRUE)

diab_rf_final <- complete(imp_rf)


# Vérifier NA résiduels et appliquer ranger si nécessaire (sur les trois versions : knn/em)
res_na_knn <- colSums(is.na(diab_knn))
res_na_em  <- colSums(is.na(diab_em))

cat("NA résiduels (kNN) :\n"); print(res_na_knn[res_na_knn > 0])
cat("NA résiduels (EM)  :\n"); print(res_na_em[res_na_em > 0])


# Imputer résiduels kNN
vars_knn_res <- names(res_na_knn[res_na_knn > 0])
if(length(vars_knn_res) > 0) {
  cat("Imputation ranger pour variables restantes dans kNN...\n")
  pred_set <- intersect(names(diab_knn)[sapply(diab_knn, is.numeric)],
                        c("time_in_hospital","num_lab_procedures","num_procedures","num_medications","number_diagnoses"))
  if(length(pred_set) < 1) pred_set <- names(diab_knn)[sapply(diab_knn, is.numeric)][1:min(3, sum(sapply(diab_knn, is.numeric)))]
  for(v in vars_knn_res) {
    cat(" -> ranger imputer", v, "\n")
    try({
      diab_knn <- impute_with_ranger(diab_knn, target = v, pred_vars = pred_set, sample_train = 15000, ntree = 50)
      cat("    NA restants pour", v, ":", sum(is.na(diab_knn[[v]])), "\n")
    }, silent = TRUE)
    gc()
  }
}

# Imputer résiduels EM
vars_em_res <- names(res_na_em[res_na_em > 0])
if(length(vars_em_res) > 0) {
  cat("Imputation ranger pour variables restantes dans EM...\n")
  pred_set2 <- intersect(names(diab_em)[sapply(diab_em, is.numeric)],
                         c("time_in_hospital","num_lab_procedures","num_procedures","num_medications","number_diagnoses"))
  if(length(pred_set2) < 1) pred_set2 <- names(diab_em)[sapply(diab_em, is.numeric)][1:min(3, sum(sapply(diab_em, is.numeric)))]
  for(v in vars_em_res) {
    cat(" -> ranger imputer", v, "\n")
    try({
      diab_em <- impute_with_ranger(diab_em, target = v, pred_vars = pred_set2, sample_train = 15000, ntree = 50)
      cat("    NA restants pour", v, ":", sum(is.na(diab_em[[v]])), "\n")
    }, silent = TRUE)
    gc()
  }
}

# 9) Fallback final : médiane / mode pour garantir 0 NA
final_fix <- function(df) {
  for(nm in names(df)) {
    if(any(is.na(df[[nm]]))) {
      if(is.numeric(df[[nm]])) {
        df[[nm]][is.na(df[[nm]])] <- median(df[[nm]], na.rm = TRUE)
      } else {
        moda <- names(which.max(table(df[[nm]])))
        df[[nm]][is.na(df[[nm]])] <- moda
        df[[nm]] <- factor(df[[nm]])
      }
    }
  }
  return(df)
}

diab_knn_final <- final_fix(diab_knn)
diab_em_final  <- final_fix(diab_em)
# pour MissForest, on n'a imputé que l'échantillon ; ici on n'applique pas missForest sur tout le dataset (trop lourd)
# tu peux choisir entre diab_knn_final / diab_em_final pour les étapes suivantes.

cat("\nNA finaux (kNN/EM) :", sum(is.na(diab_knn_final)), "/", sum(is.na(diab_em_final)), "\n")

# 10) Sauvegardes
saveRDS(diab_knn_final, "diabetes_imputed_knn_final.rds")
saveRDS(diab_em_final,  "diabetes_imputed_em_final.rds")
write.csv(diab_knn_final, "diabetes_imputed_knn_final.csv", row.names = FALSE)
write.csv(diab_em_final,  "diabetes_imputed_em_final.csv",  row.names = FALSE)
saveRDS(diab_rf_final, "diabetes_missforest_sample_imputed.rds")
write.csv(diab_rf_final, "diabetes_missforest_sample_imputed.csv", row.names = FALSE)

cat("\n== Imputation terminée: fichiers écrits ==\n")
cat("Résumé NA finaux (kNN / EM):", sum(is.na(diab_knn_final)), "/", sum(is.na(diab_em_final)), "\n")


