### ============================================
### 1. IMPORTATION & NETTOYAGE DES DONNÉES
### ============================================

# Chargement des librairies
library(dplyr)
library(stringr)

# ------------------------------
# IMPORTATION DE LA BASE
# ------------------------------
df <- read.csv("diabetic_data.csv", na.strings = c("?", "", "NA"))

# Aperçu général
dim(df)                  # nombre de lignes / colonnes
summary(df$readmitted)   # type de la variable cible brute


# ===================================================
# 1) SUPPRESSION DES VARIABLES INUTILES OU TROP VIDES
# ===================================================

# Liste des colonnes à retirer :
# - IDs (encounter_id, patient_nbr)
# - Variables quasi vides (weight, payer_code, medical_specialty)
# - Les 20 colonnes de médicaments (trop complexes + trop de NA)
cols_to_remove <- c(
  "encounter_id", "patient_nbr",
  "weight", "payer_code", "medical_specialty",
  "metformin", "repaglinide", "nateglinide", "chlorpropamide",
  "glimepiride", "acetohexamide", "glipizide", "glyburide",
  "tolbutamide", "pioglitazone", "rosiglitazone", "acarbose",
  "miglitol", "troglitazone", "tolazamide", "insulin",
  "glyburide.metformin", "glipizide.metformin", "metformin.rosiglitazone",
  "metformin.pioglitazone"
)

# Suppression
df <- df %>% select(-one_of(cols_to_remove))

# Vérification
dim(df)


# ======================================================
# 2) TRANSFORMATION DE LA VARIABLE CIBLE EN BINAIRE
# ======================================================

# readmitted prend 3 valeurs : "<30", ">30", "NO"
# Nous créons une version binaire : 1 = réadmission < 30 jours, 0 sinon
df$readmit_binary <- ifelse(df$readmitted == "NO", 0, 1)

# Vérification
table(df$readmit_binary)
prop.table(table(df$readmit_binary))


# =======================================
# 3) GESTION DES VALEURS MANQUANTES (NA)
# =======================================

# Fonction : renvoie la modalité la plus fréquente (mode)
mode_impute <- function(x){
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

# Application de l'imputation :
# - pour les variables texte/factor : remplacer NA par la valeur la plus fréquente
# - pour les variables numériques : remplacer NA par la médiane
df <- df %>%
  mutate(across(where(is.character), ~replace(., is.na(.), mode_impute(.)))) %>%
  mutate(across(where(is.factor), ~replace(., is.na(.), mode_impute(.)))) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), median(., na.rm = TRUE))))

# ===============================
# 4) ENCODAGE VARIABLES MÉDICALES
# ===============================

# Nettoyage de A1Cresult : renommer les catégories pour les rendre lisibles
df$A1Cresult <- factor(df$A1Cresult,
                       levels = c("None", "Norm", ">7", ">8"),
                       labels = c("None", "Normal", "High", "VeryHigh"))

# Nettoyage de max_glu_serum
df$max_glu_serum <- factor(df$max_glu_serum,
                           levels = c("None", "Norm", ">200", ">300"),
                           labels = c("None", "Normal", "High", "VeryHigh"))


# ==============================================
# 5) REGROUPEMENT DIAGNOSTIQUES 
#    diag1, diag2, diag3 -> diag_1_group, diag_2_group, diag_3_group
# ==============================================

# Fonction qui transforme un code médical ICD en grande catégorie
regroup_diag <- function(code){
  if(is.na(code)) return("Other")
  first <- substr(code, 1, 1)   # on lit juste le premier chiffre du code
  
  if(first == "2") return("Neoplasms")      # Cancers
  if(first == "3") return("Blood")          # Sang / immunité
  if(first == "4") return("Endocrine")      # Diabète, thyroïde...
  if(first == "5") return("Mental")         # Psychiatrie
  if(first == "6") return("Nervous")        # Neurologie
  if(first == "7") return("Circulatory")    # Cœur, vascular
  if(first == "8") return("Respiratory")    # Poumons
  if(first == "9") return("Digestive")      # Tube digestif, reins
  
  return("Other")
}

# Conversion en caractere
df$diag_1 <- as.character(df$diag_1)
df$diag_2 <- as.character(df$diag_2)
df$diag_3 <- as.character(df$diag_3)

# Application du regroupement
df$diag_1_group <- factor(sapply(df$diag_1, regroup_diag))
df$diag_2_group <- factor(sapply(df$diag_2, regroup_diag))
df$diag_3_group <- factor(sapply(df$diag_3, regroup_diag))

# Vérification
table(df$diag_1_group)
table(df$diag_2_group)
table(df$diag_3_group)



# ======================
# 6) TRAITEMENT DE L'ÂGE
# ======================

# L'âge est déjà en classes ordonnées, on le convertit proprement
df$age <- factor(df$age,
                 levels = c("[0-10)", "[10-20)", "[20-30)", "[30-40)", "[40-50)",
                            "[50-60)", "[60-70)", "[70-80)", "[80-90)", "[90-100)"),
                 ordered = TRUE)

df <- df %>% select(-diag_1, -diag_2, -diag_3)


View(df)

# ============================================
# RÉCAPITULATIF DES VARIABLES GARDÉES / SUPPRIMÉES
# ============================================

# Variables SUPPRIMÉES :
# - encounter_id, patient_nbr (identifiants)
# - weight, payer_code, medical_specialty (trop de NA)
# - 20 colonnes de médicaments (trop complexes, trop de NA, faible valeur prédictive)

# Variables CONSERVÉES :
# - readmitted (brut) + readmit_binary (cible binaire)
# - A1Cresult (renommée)
# - max_glu_serum (renommée)
# - age
# - time_in_hospital
# - num_lab_procedures
# - num_procedures
# - number_emergency
# - number_inpatient
# - diag_1, diag_2, diag_3 (brutes)
# - diag_1_group, diag_2_group, diag_3_group (regroupées proprement)
# - race, gender, admission_type, discharge_disposition, etc. (toutes conservées)

### ============================================
### 2. ANALYSE EXPLORATOIRE DU DATAFRAME COMPLET
### ============================================

library(dplyr)
library(ggplot2)

### --------------------------------------------
### 2.0 APERÇU GLOBAL DU DATAFRAME
### --------------------------------------------

# Structure de df : type de chaque variable, quelques exemples de valeurs
str(df)

# Résumé statistique global (numériques + qualitatives)
summary(df)

df <- df %>% select(-examide, -citoglipton, -glimepiride.pioglitazone)


### --------------------------------------------
### 2.1 ANALYSE DE LA VARIABLE CIBLE
### --------------------------------------------

## 2.1.1 Variable readmitted (3 modalités : "<30", ">30", "NO")
## -> description brute, juste pour comprendre la répartition des types de réadmission

# Effectifs
table(df$readmitted)

# Proportions
prop.table(table(df$readmitted))

# Petit barplot de la version 3 classes (brute)
ggplot(df, aes(x = readmitted)) +
  geom_bar() +
  labs(title = "Distribution de la variable readmitted (3 modalités)",
       x = "readmitted",
       y = "Nombre de patients") +
  theme_minimal()


## 2.1.2 Variable binaire readmit_binary (0/1 : 1 = réadmission, 0 = pas de réadmission)
## -> C'EST LA VRAIE VARIABLE CIBLE DU PROJET

# Effectifs
table(df$readmit_binary)

# Proportions
prop.table(table(df$readmit_binary))

# Barplot de la cible binaire
ggplot(df, aes(x = factor(readmit_binary))) +
  geom_bar() +
  labs(x = "Réadmission (0 = Non, 1 = Oui)",
       y = "Nombre de patients",
       title = "Distribution de la variable cible binaire (réadmis ou non)") +
  theme_minimal()

### --------------------------------------------
### 2.2 DESCRIPTION DES VARIABLES NUMÉRIQUES
### --------------------------------------------

# Sélection de toutes les variables numériques
vars_num <- df %>%
  select(where(is.numeric))

# Vérifier quelles sont les variables numériques
colnames(vars_num)

# Résumé statistique (min, max, moyenne, médiane, quartiles, etc.)
summary(vars_num)

# Histogrammes pour toutes les variables numériques
for (v in colnames(vars_num)) {
  print(
    ggplot(df, aes_string(x = v)) +
      geom_histogram(bins = 30) +
      labs(title = paste("Histogramme de", v),
           x = v,
           y = "Effectif") +
      theme_minimal()
  )
}

### --------------------------------------------
### 2.3 DESCRIPTION DES VARIABLES QUALITATIVES
### --------------------------------------------

# Sélection de toutes les variables qualitatives (factor ou character)
vars_cat <- df %>%
  select(where(~ is.factor(.) || is.character(.)))

# Vérifier quelles sont les variables qualitatives
colnames(vars_cat)

# Pour chaque variable qualitative : tableau d'effectifs + proportions
for (v in colnames(vars_cat)) {
  cat("\n\n==============================\n")
  cat("Variable :", v, "\n")
  cat("==============================\n")
  print(table(vars_cat[[v]]))
  cat("\nProportions :\n")
  print(prop.table(table(vars_cat[[v]])))
}

# Quelques barplots pour les variables qualitatives CLÉS

# Âge
ggplot(df, aes(x = age)) +
  geom_bar() +
  labs(title = "Distribution de l'âge",
       x = "Classe d'âge",
       y = "Nombre de patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sexe
ggplot(df, aes(x = gender)) +
  geom_bar() +
  labs(title = "Répartition par sexe",
       x = "Sexe",
       y = "Nombre de patients") +
  theme_minimal()

# Résultat A1C
ggplot(df, aes(x = A1Cresult)) +
  geom_bar() +
  labs(title = "Distribution de A1Cresult",
       x = "Résultat A1C",
       y = "Nombre de patients") +
  theme_minimal()

# Diagnostic principal regroupé
ggplot(df, aes(x = diag_1_group)) +
  geom_bar() +
  labs(title = "Répartition des diagnostics principaux (diag_1_group)",
       x = "Groupe diagnostique",
       y = "Nombre de patients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### ---------------------------------------------------------
### 2.4 PREMIÈRES RELATIONS VISUELLES AVEC LA VARIABLE CIBLE
###      (DESCRIPTIF SEULEMENT, SANS TESTS STATISTIQUES)
### ---------------------------------------------------------

## 2.4.1 Variables numériques vs readmit_binary : boxplots

# On exclut la cible binaire des variables numériques explicatives
vars_num_explicatives <- vars_num %>% select(-readmit_binary)

for (v in colnames(vars_num_explicatives)) {
  print(
    ggplot(df, aes(x = factor(readmit_binary), y = .data[[v]])) +
      geom_boxplot() +
      labs(title = paste("Boxplot de", v, "selon le statut de réadmission"),
           x = "Réadmission (0 = Non, 1 = Oui)",
           y = v) +
      theme_minimal()
  )
}


## 2.4.2 Variables qualitatives vs readmit_binary : barplots empilés

# Pour chaque variable qualitative, on affiche la proportion de réadmis/non-réadmis
for (v in colnames(vars_cat)) {
  print(
    ggplot(df, aes(x = .data[[v]], fill = factor(readmit_binary))) +
      geom_bar(position = "fill") +
      labs(title = paste("Proportion de réadmission selon", v),
           x = v,
           y = "Proportion de patients",
           fill = "Réadmission (0 = Non, 1 = Oui)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}

prop.table(table(df$readmit_binary))
df <- df %>% select(-readmitted)

###########################################################
### PARTIE 3 — Sélection univariée des variables
### Objectif : tester si chaque variable explicative
### est associée à la réadmission (readmit_binary)
### en utilisant uniquement les méthodes vues en cours :
###  - Test du Chi-deux (variables qualitatives)
###  - Test de Wilcoxon (variables numériques)
###  - Correction de Benjamini–Hochberg (contrôle du FDR)
###########################################################

library(dplyr)

###########################################################
### 1. Séparer les variables selon leur type
###########################################################

# Variables numériques (on retire la cible)
vars_num <- df %>% 
  select(where(is.numeric)) %>% 
  select(-readmit_binary)

# Variables qualitatives (character ou factor)
vars_cat <- df %>%
  select(where(~ is.factor(.) || is.character(.)))

###########################################################
### 2. TESTS POUR LES VARIABLES QUALITATIVES
### Méthode : test du Chi-deux d’indépendance
### H0 : pas d’association avec la réadmission
### H1 : association significative
###########################################################

pvals_cat <- c()   # vecteur pour stocker les p-values

for (v in colnames(vars_cat)) {
  
  # Table de contingence : modalités de la variable vs réadmission
  tab <- table(vars_cat[[v]], df$readmit_binary)
  
  # Test du Chi-deux (warnings supprimés pour les petites fréquences)
  test <- suppressWarnings(chisq.test(tab))
  
  # Stockage de la p-value brute
  pvals_cat[v] <- test$p.value
}

# Correction du risque de faux positifs (FDR)
pvals_cat_BH <- p.adjust(pvals_cat, method = "BH")

# Tableau récapitulatif propre
cat_results <- data.frame(
  variable = names(pvals_cat),
  p_value = pvals_cat,
  p_value_BH = pvals_cat_BH
) %>% arrange(p_value_BH)
cat_results
###########################################################
### 3. TESTS POUR LES VARIABLES NUMÉRIQUES
### Méthode : test de Wilcoxon (non paramétrique)
### Justification : distributions asymétriques dans le dataset
###########################################################

pvals_num <- c()  # vecteur pour stocker les p-values

for (v in colnames(vars_num)) {
  
  # Séparer la variable selon les groupes (0 = non réadmis, 1 = réadmis)
  x0 <- vars_num[df$readmit_binary == 0, v]
  x1 <- vars_num[df$readmit_binary == 1, v]
  
  # Test de Wilcoxon pour comparer les distributions
  test <- wilcox.test(x0, x1)
  
  # Stockage de la p-value brute
  pvals_num[v] <- test$p.value
}

# Correction BH sur les variables numériques
pvals_num_BH <- p.adjust(pvals_num, method = "BH")

# Tableau propre
num_results <- data.frame(
  variable = names(pvals_num),
  p_value = pvals_num,
  p_value_BH = pvals_num_BH
) %>% arrange(p_value_BH)

###########################################################
### 4. AFFICHAGE DES VARIABLES SIGNIFICATIVES
### (p-value BH < 0.05)
### Ce sont les variables retenues pour la suite
### → elles seront utilisées dans la modélisation multivariée
###########################################################


cat("Variables qualitatives significatives (Chi-deux + BH): ")
print(filter(cat_results, p_value_BH < 0.05))


cat("Variables numériques significatives (Wilcoxon + BH): ")
print(filter(num_results, p_value_BH < 0.05))



############################################################
### PARTIE 4 — Sélection multivariée par régularisation
### Objectifs :
###  - Identifier les variables réellement prédictives
###  - Réduire la dimension
###  - Gérer la multicolinéarité
###  - Éviter le sur-apprentissage (overfitting)
###
### Méthodes utilisées :
###  - Ridge Regression (pénalisation L2)
###  - LASSO (pénalisation L1)
###  - Elastic-Net (compromis L1/L2)
############################################################

library(glmnet)
library(caret)
library(dplyr)

############################################################
### 1. Préparation des données
###    → glmnet exige des variables numériques uniquement
###    → On encode donc toutes les variables catégorielles
###      en "dummy variables" via model.matrix
############################################################
df_model <- df
# Création de la matrice X (numérique uniquement)
# model.matrix crée automatiquement :
#   - variables dummy pour les facteurs
#   - un intercept (qu’on retire avec [, -1])
X <- model.matrix(readmit_binary ~ ., data = df_model)[, -1]

# Variable cible
y <- df_model$readmit_binary

############################################################
### 2. Séparation en échantillon d'entraînement et test (70/30)
###    → permet d'évaluer les modèles ensuite
############################################################

set.seed(123)  # reproductibilité
train_index <- createDataPartition(y, p = 0.7, list = FALSE)

X_train <- X[train_index, ]
X_test  <- X[-train_index, ]

y_train <- y[train_index]
y_test  <- y[-train_index]

############################################################
### 3. MODELE RIDGE (alpha = 0)
###    → Ne sélectionne pas les variables
###    → Mais stabilise les coefficients
###    → Utile pour voir "l'importance globale"
############################################################

set.seed(123)

ridge_cv <- cv.glmnet(
  X_train, y_train,
  alpha = 0,              # alpha = 0 → Ridge
  family = "binomial"     # modèle logistique pénalisé
)

# Lambda optimal (celui qui minimise l’erreur)
ridge_lambda <- ridge_cv$lambda.min

cat("===== RIDGE : meilleur lambda =====")
print(ridge_lambda)

# Extraction des coefficients Ridge
ridge_coef <- coef(ridge_cv, s = "lambda.min")

cat("===== Coefficients Ridge =====")
print(ridge_coef)

############################################################
### 4. MODELE LASSO (alpha = 1)
###    → Sélection automatique des variables
###    → Les coefficients inutiles deviennent exactement 0
###    → Modèle "sparse" = épuré, très utile pour l'interprétation
############################################################

set.seed(123)

lasso_cv <- cv.glmnet(
  X_train, y_train,
  alpha = 1,              # alpha = 1 → LASSO
  family = "binomial"
)

lasso_lambda <- lasso_cv$lambda.min

cat("===== LASSO : meilleur lambda =====")
print(lasso_lambda)

# Coefficients LASSO
lasso_coef <- coef(lasso_cv, s = "lambda.min")

cat("===== Coefficients LASSO =====")
print(lasso_coef)

# Liste des variables sélectionnées (coeff ≠ 0)
selected_lasso <- rownames(lasso_coef)[lasso_coef[, 1] != 0]

cat("===== Variables sélectionnées par LASSO =====")
print(selected_lasso)

############################################################
### 5. MODELE ELASTIC NET (alpha = 0.5)
###    → Combine L1 et L2
###    → Utile quand plusieurs variables sont corrélées
###    → Plus stable que LASSO seul
############################################################

set.seed(123)

enet_cv <- cv.glmnet(
  X_train, y_train,
  alpha = 0.5,            # compromis Lasso / Ridge
  family = "binomial"
)

enet_lambda <- enet_cv$lambda.min

cat("===== Elastic Net : meilleur lambda =====")
print(enet_lambda)

# Coefficients Elastic-Net
enet_coef <- coef(enet_cv, s = "lambda.min")

cat("===== Coefficients Elastic Net =====")
print(enet_coef)

############################################################
### 6. GRAPHIQUES DE VALIDATION
###    → montrent comment l’erreur varie selon lambda
###    → utilisés en cours pour choisir le lambda optimal
############################################################

plot(ridge_cv)
title("Validation croisée — Ridge Regression")

plot(lasso_cv)
title("Validation croisée — LASSO")

plot(enet_cv)
title("Validation croisée — Elastic Net")

