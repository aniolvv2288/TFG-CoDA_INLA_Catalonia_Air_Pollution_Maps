## Aniol Vilamala Vidal

################################################################################
#  ----- IMPUTACIÓ DELS MISSINGS ----- 
################################################################################

dades <- readRDS("dades_preprocessades.rds")

contaminants <- c("NO2", "PM10", "PM2.5", "O3", "CO", "SO2")

X <- dades %>%
  dplyr::select(all_of(contaminants))

## ----- KNN -----

imp_knn <- robCompositions::impKNNa(X, k = 5)

imp_knn

## ----- REGRESSIÓ ROBUSTA -----

imp_lts <- robCompositions::impCoda(X, method = "ltsReg", k = 5)

imp_lts


################################################################################
# ----- ANÀLISI DE SENSIBILITAT -----
################################################################################

impKNN <- data.frame(acomp(imp_knn$xImp))
impLTS <- imp_lts$xImp

## ----- BOXPLOT COMPARATIU PER CONTAMINANT -----

impKNN_long <- impKNN %>%
  mutate(Mètode = "KNN") %>%
  melt(id.vars = "Mètode", variable.name = "Contaminant", value.name = "Valor")

impLTS_long <- impLTS %>%
  mutate(Mètode = "LTS") %>%
  melt(id.vars = "Mètode", variable.name = "Contaminant", value.name = "Valor")

na_matrix <- is.na(as.matrix(X))
tipus_valor <- ifelse(as.vector(na_matrix), "Imputada", "Original")

df_comp <- data.frame(
  KNN = as.vector(as.matrix(impKNN)),
  LTS = as.vector(as.matrix(impLTS)),
  Contaminant = rep(colnames(impKNN), each = nrow(impKNN)),
  Tipus = tipus_valor,
  TIPUS_ESTACIO = rep(dades$TIPUS_ESTACIO, times = length(contaminants)),
  AREA_URBANA  = rep(dades$AREA_URBANA,  times = length(contaminants))
)

df_comp_long <- melt(df_comp,
                     id.vars = c("Contaminant", "Tipus", "TIPUS_ESTACIO", "AREA_URBANA"),
                     variable.name = "Mètode",
                     value.name = "Valor")

## ----- BOXPLOT COMPARATIU -----

ggplot(df_comp_long, aes(x = Contaminant, y = Valor, fill = Mètode)) +
  geom_boxplot(alpha = 0.7,
               outlier.color = "red",
               position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 10) +
  labs(title = "Comparació KNN vs LTS per contaminant",
       x = "Contaminant",
       y = "Valor imputat") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## ----- BOXPLOT COMPARATIU ESTRATIFICAT PER TIPUS ESTACIÓ I ÀREA URBANA -----

ggplot(df_comp_long, aes(x = TIPUS_ESTACIO, y = Valor, fill = Mètode)) +
  geom_boxplot(alpha = 0.7,
               outlier.color = "grey30",
               position = position_dodge(width = 0.8)) +
  facet_grid(AREA_URBANA ~ Contaminant, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 10) +
  labs(title = "Comparació de mètodes d’imputació per contaminant",
       subtitle = "Estratificat per tipus d’estació i àrea urbana",
       x = "Tipus d’estació",
       y = "Valor imputat",
       fill = "Mètode") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Boxplot segons els 2 contaminants més diferents segons la imputació (CO i O3):

ggplot(df_comp_long %>% dplyr::filter(Contaminant %in% c("CO", "O3")), 
       aes(x = TIPUS_ESTACIO, y = Valor, fill = Mètode)) +
  geom_boxplot(alpha = 0.7,
               outlier.color = "red",
               position = position_dodge(width = 0.8)) +
  facet_grid(AREA_URBANA ~ Contaminant, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 10) +
    labs(
      title = expression("Comparació imputacions per CO i " * O[3]),
      subtitle = "Estratificat per tipus d’estació i àrea urbana",
      x = "Tipus d’estació",
      y = "Valor imputat",
      fill = "Mètode"
    ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1))

# Boxplots de les dades originals pels 2 contaminants

ggplot(
  df_comp_long %>%
    filter(Contaminant %in% c("CO", "O3"),
           Tipus == "Original"),
  aes(x = TIPUS_ESTACIO, y = Valor, fill = Mètode)
) +
  geom_boxplot(alpha = 0.7,
               outlier.color = "red",
               position = position_dodge(width = 0.8)) +
  facet_grid(AREA_URBANA ~ Contaminant, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 10) +
  labs(
    title = expression("Distribucions originals de CO i " * O[3]),
    subtitle = "Estratificat per tipus d’estació i àrea urbana",
    x = "Tipus d’estació",
    y = "Valor observat",
    fill = "Mètode"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

## ----- SCATTERPLOT COMPARATIU -----

ggplot(df_comp, aes(x = KNN, y = LTS,
                    color = Contaminant,
                    shape = Tipus)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1,
              color = "black", linetype = "dashed") +
  scale_shape_manual(values = c("Original" = 16, "Imputada" = 8)) +
  theme_minimal(base_size = 10) +
  labs(title = "Comparació imputacions KNN vs LTS",
       x = "Imputació KNN",
       y = "Imputació LTS",
       color = "Contaminant",
       shape = "Tipus de dada") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## ----- DIFERÈNCIES PCA ---

pca_knn <- princomp(clr(acomp(impKNN)))
pca_lts <- princomp(clr(acomp(impLTS)))

scores_knn <- as.data.frame(pca_knn$scores); scores_knn$Mètode <- "KNN"
scores_lts <- as.data.frame(pca_lts$scores); scores_lts$Mètode <- "LTS"
scores_all <- rbind(scores_knn, scores_lts)

ggplot(scores_all, aes(x = Comp.1, y = Comp.2, color = Mètode)) +
  geom_point(alpha = 0.6) +
  theme_minimal(base_size = 13) +
  labs(title = "PCA KNN vs LTS",
       x = "Component 1",
       y = "Component 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Variabilitat total explicada segons el mètode:

var_knn <- pca_knn$sdev^2
var_lts <- pca_lts$sdev^2
pct_knn <- round(100 * var_knn / sum(var_knn), 2); pct_knn # 80.11
pct_lts <- round(100 * var_lts / sum(var_lts), 2); pct_lts # 79.25

## ----- CONCORDANÇA GLOBAL -----

cor_total <- cor(as.vector(as.matrix(impKNN)),
                 as.vector(as.matrix(impLTS)))

rmse_total <- sqrt(mean((as.vector(as.matrix(impKNN)) -
                           as.vector(as.matrix(impLTS)))^2))

resultats <- data.frame(
  Mètrica = c("Correlació global", "RMSE global"),
  Valor = c(round(cor_total, 3), round(rmse_total, 4))
)

titol_html <- "<span style='color:#1b1b1b; font-weight:bold; font-size:15px;'>Comparació global entre KNN i LTS</span>"

kable(resultats, format = "html", align = "c", caption = titol_html) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"),
                full_width = FALSE, position = "center") %>%
  row_spec(0, bold = TRUE, color = "white", background = "#2e7d32") %>%
  column_spec(1, width = "10em") %>%  
  column_spec(2, width = "6em", color = "black", background = "#f2f2f2", bold = TRUE)


################################################################################
# ----- EXPORTACIÓ CONJUNT IMPUTAT TRIAT (KNN) -----
################################################################################

dades_imputades <- dades
dades_imputades[, contaminants] <- imp_knn$xImp

saveRDS(dades_imputades, file = "dades_imputades.rds")

