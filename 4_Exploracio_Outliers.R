## Aniol Vilamala Vidal

################################################################################
# ----- DETECCIÓ D'OUTLIERS MULTIVARIANTS -----
################################################################################

dades <- readRDS("dades_preprocessades.rds")

dades_imputades <- readRDS("dades_imputades.rds")

contaminants <- c("NO2", "PM10", "PM2.5", "O3", "CO", "SO2")

X <- dades_imputades %>%
  dplyr::select(all_of(contaminants)) 

# set.seed(1714)
# out <- OutlierClassifier1(
#   acomp(X),
#   alpha = 0.05,              
#   type = "all",            
#   corrected = TRUE,          
#   RedCorrected = FALSE,      
#   robust = TRUE              
# )

outlier <- ifelse(out == "ok", "No outlier", "Outlier")

outlier_quin <- sapply(out, function(class) {
  class <- as.character(class) 
  if (class == "ok") {
    return("ok")
  } else if (class == "000000") {
    return("?")
  } else {
    bits <- strsplit(class, "")[[1]]
    idx <- which(bits == "1")
    paste(contaminants[idx], collapse = "_")
  }
})

dades_outliers <- dades_imputades %>%
  mutate(OUTLIER = outlier,
         OUTLIER_BIN = ifelse(outlier == "Outlier", 1, 0),
         OUTLIER_QUIN = outlier_quin)

casos_outlier <- dades_outliers %>%
  filter(OUTLIER == "Outlier")
dim(casos_outlier)

################################################################################
# ----- OUTLIERS SEGONS TIPUS D'ESTACIÓ I ÀREA URBANA -----
################################################################################

table(casos_outlier$TIPUS_ESTACIO, casos_outlier$AREA_URBANA)

resum_outliers <- dades_outliers %>%
  group_by(AREA_URBANA, TIPUS_ESTACIO, OUTLIER) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(AREA_URBANA, TIPUS_ESTACIO) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(AREA_URBANA = factor(AREA_URBANA, levels = c("urban", "suburban", "rural")))

ggplot(resum_outliers, aes(x = TIPUS_ESTACIO, y = prop, fill = OUTLIER)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ AREA_URBANA) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No outlier" = "steelblue", "Outlier" = "red")) +
  labs(title = "% d’outliers segons tipus d'estació i àrea urbana",
       x = "Tipus d’estació", y = "Proporció", fill = "Categoria") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
# ----- DINS DELS OUTLIERS, QUINS SÓN ORIGINALS I QUINS IMPUTATS -----
################################################################################

dades <- dades %>%
  mutate(ID = row_number())

dades_imputades <- dades_imputades %>%
  mutate(ID = row_number())

dades_outliers <- dades_outliers %>%
  mutate(ID = row_number())

info_originals <- dades %>%
  dplyr::select(ID, all_of(contaminants)) %>%
  pivot_longer(cols = all_of(contaminants),
               names_to = "contaminant", 
               values_to = "valor") %>%
  group_by(ID) %>%
  summarise(n_originals = sum(!is.na(valor)),
            .groups = "drop")

dades_outliers <- dades_outliers %>%
  left_join(info_originals, by = "ID")

casos_outlier <- dades_outliers %>%
  filter(OUTLIER == "Outlier")

table(casos_outlier$n_originals)

ggplot(casos_outlier, aes(x = factor(n_originals))) +
  geom_bar(fill = "steelblue", color = "white", alpha = 0.8) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.4, size = 4) +
  labs(title = "Outliers segons el nombre de valors originals",
       x = "Nombre de contaminants no imputats",
       y = "Nombre d'outliers") +
  coord_cartesian(ylim = c(0, 50)) + 
  theme_minimal(base_size = 13)

ggplot(casos_outlier, aes(x = factor(n_originals), fill = TIPUS_ESTACIO)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  labs(title = "Outliers segons el nombre de valors originals",
       subtitle = "Estratificat segons el tipus d'estació",
       x = "Nombre de contaminants no imputats",
       y = "Freqüència",
       fill = "Tipus d'estació") +
  theme_minimal(base_size = 13)

ggplot(casos_outlier, aes(x = factor(n_originals), fill = AREA_URBANA)) +
  geom_bar(position = "dodge", alpha = 0.8) +
  labs(title = "Outliers segons el nombre de valors originals",
       subtitle = "Estratificat segons l'àrea urbana",
       x = "Nombre de contaminants no imputats",
       y = "Freqüència",
       fill = "Tipus d'estació") +
  theme_minimal(base_size = 13)


################################################################################
# ----- ELIMINACIÓ OUTLIERS AMB < 3 CONTAMINANTS ORIGINALS -----
################################################################################

exclosos <- dades_outliers %>%
  filter(n_originals <= 3 & OUTLIER == "Outlier")

dades_netes <- dades_outliers %>%
  filter(!(ID %in% exclosos$ID))


################################################################################
## ----- VISUALITZACIÓ OUTLIERS A CATALUNYA-----
################################################################################

cat <- gisco_get_nuts(
  year = 2021,
  resolution = "3",
  nuts_level = 2,
  country = "ES"
) %>%
  filter(NUTS_NAME %in% c("Cataluña", "Catalonia")) %>%
  st_transform(4326)

prov_cat <- gisco_get_nuts(
  year = 2021,
  resolution = "3",
  nuts_level = 3,
  country = "ES"
) %>%
  filter(NUTS_NAME %in% c("Barcelona", "Girona", "Lleida", "Tarragona")) %>%
  st_transform(4326)

outlier_per_estacio <- dades_netes %>%
  group_by(`NOM ESTACIO`, LATITUD, LONGITUD, TIPUS_ESTACIO, AREA_URBANA) %>%
  summarise(
    pct_outliers = 100 * mean(OUTLIER_BIN == 1, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  geom_sf(data = cat, fill = "grey95", color = "black", size = 0.3) +
  geom_sf(data = prov_cat, fill = NA, color = "grey40", linetype = "dashed") +
  geom_point(
    data = outlier_per_estacio,
    aes(x = LONGITUD, y = LATITUD, color = pct_outliers),
    size = 3
  ) +
  scale_color_gradient(
    name = "% d'outliers",
    low = "steelblue",
    high = "red",
    limits = c(0, 100)
  ) +
  labs(
    title = "Localització i % d'observacions outlier per estació",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


################################################################################
# ----- EXPORTACIÓ CONJUNT DE DADES FINAL -----
################################################################################

saveRDS(dades_netes, file = "dades_finals.rds")

