## Aniol Vilamala Vidal

################################################################################
# INCORPORACIÓ DE LES DADES
################################################################################
 
XVPCA_2020_2024 <- read_csv("XVPCA_2020_2024.csv")

# Les dades originals (crues) pesen > 100 MB i no es poden pujar a GitHub

################################################################################
# PROCESSAMENT INICIAL
################################################################################

## SELECCIÓ VARIABLES RELLEVANTS: 

dades <- XVPCA_2020_2024 %>%
  select(-c(`CODI EOI`, `CODI INE`, `CODI COMARCA`, `NOM COMARCA`, 
            MAGNITUD, UNITATS, MUNICIPI, Georeferència)
  ) %>%
  mutate(`TIPUS ESTACIO` = as_factor(`TIPUS ESTACIO`),
         `AREA URBANA` = as_factor(`AREA URBANA`))

## SELECCIÓ DELS CONTAMINANTS: 

contaminants <- c("NO2", "PM10", "PM2.5", "O3", "CO", "SO2") 
dades <- dades %>%
  filter(CONTAMINANT %in% contaminants)

# Perdem dades de les estacions "Cubelles (Poliesportiu)" -08074005- i 
# "Amposta" -43014001- perquè no disposem de quin contaminant s'avalua a cada entrada

## RECODIFICACIÓ VARIABLE ANY I MES:

dades <- dades %>%
  mutate(
    DATA = dmy(DATA),
    ANY = year(DATA),
    MES = month(DATA)
  )

## RECODIFICACIÓ VARIABLES LATITUD i LONGITUD:

dades <- dades %>%
  mutate(LONGITUD = as.numeric(gsub(",", ".", LONGITUD)),
         LATITUD = as.numeric(
           sub("^([0-9]{2})([0-9]+)$", "\\1.\\2", as.character(LATITUD))
         ))

## MITJANA DIÀRIA PER ESTACIÓ DE CADA CONTAMINANT + FORA "NaN"

dades <- dades %>%
  mutate(
    CONT_DIARIA = rowMeans(select(., `01h`:`24h`), na.rm = TRUE),
    DATA = as.Date(DATA)  # assegura't que tens una columna de data
  ) %>%
  filter(!is.nan(CONT_DIARIA))

## MITJANA MENSUAL PER ESTACIÓ DE CADA CONTAMINANT

dades <- dades %>%
  group_by(`NOM ESTACIO`, CONTAMINANT, ANY, MES) %>%
  summarise(
    MITJANA_MENSUAL = mean(CONT_DIARIA, na.rm = TRUE),
    TIPUS_ESTACIO   = first(`TIPUS ESTACIO`),
    AREA_URBANA     = first(`AREA URBANA`),
    ALTITUD         = first(ALTITUD),
    LATITUD         = first(LATITUD),
    LONGITUD        = first(LONGITUD),
    .groups = "drop"
  ) %>%
  arrange(`NOM ESTACIO`, CONTAMINANT, ANY, MES)

## DADES EN FORMAT AMPLE (1 FILA PER ESTACIÓ, ANY i MES) + CONVERSIÓ UNITATS CO

dades <- dades %>%
  pivot_wider(names_from = CONTAMINANT, values_from = MITJANA_MENSUAL) %>%
  mutate(CO = CO * 1000)

################################################################################
# DETECCIÓ CASOS EXTRANYS I TRACTAMENT
################################################################################

## ESTACIONS AMB MÉS D'UN "TIPUS D'ESTACIÓ"

dades %>%
  count(`NOM ESTACIO`, ANY, MES) %>%
  filter(n > 1)

# Tant a finals del 2021 com per tot el 2022 a Sant Andreu de la Barca tenim que la única estació 
# és mesurada com a suburban i urban al mateix temps. Tractem aquest cas fent que 
# l'estació només pugui ser urbana i els valors de contaminants presents als casos
# suburban es combinin amb els d'urban per obtenir casos més complets.
# Com que només tenim informació per les PM10 canviem aquests casos.

SAB_suburban <- dades %>%
  filter(`NOM ESTACIO` == "Sant Andreu de la Barca",
         ANY %in% c(2021, 2022),
         AREA_URBANA == "suburban") %>%
  select(ANY, MES, PM10)

for(i in seq_len(nrow(SAB_suburban))) {
  any_i <- SAB_suburban$ANY[i]
  mes_i <- SAB_suburban$MES[i]
  val_i <- SAB_suburban$PM10[i]
  
  dades$PM10[
    dades$`NOM ESTACIO` == "Sant Andreu de la Barca" &
      dades$ANY == any_i &
      dades$MES == mes_i &
      dades$AREA_URBANA == "urban"
  ] <- val_i
}

dades <- dades %>%
  filter(!( `NOM ESTACIO` == "Sant Andreu de la Barca" &
              ANY %in% c(2021, 2022) &
              AREA_URBANA == "suburban"))

length(unique(dades$`NOM ESTACIO`))
dim(dades)

################################################################################
# TRACTAMENT DELS MISSINGS - SELECCIÓ CASOS AMB MÍNIM 2 OBSERVACIONS
################################################################################

dades <- dades %>%
  rowwise() %>%
  mutate(
    n_contaminants_valids = sum(!is.na(c_across(c(NO2, PM10, PM2.5, O3, SO2, CO))))
  ) %>%
  ungroup() %>%
  dplyr::filter(n_contaminants_valids >= 2) %>%
  dplyr::select(-n_contaminants_valids)

length(unique(dades$`NOM ESTACIO`))
dim(dades)

saveRDS(dades, file = "dades_2_preprocessades.rds")

################################################################################
# RESUM PRE TRACTAMENT DE MISSINGS
################################################################################

dades <- readRDS("dades_preprocessades.rds")

## BASE DE DADES:

DT::datatable(dades)

## RESUM DELS CONTAMINANTS:

resum_cont <- dades %>%
  select(-c(ALTITUD, ANY, LATITUD, LONGITUD, MES)) %>%
  select(where(is.numeric)) %>%
  summarise(
    across(
      everything(),
      list(
        n              = ~sum(!is.na(.)),
        n_missing      = ~sum(is.na(.)),
        pct_missing    = ~round(100 * mean(is.na(.)), 2),
        mean           = ~mean(., na.rm = TRUE),
        median         = ~median(., na.rm = TRUE),
        sd             = ~sd(., na.rm = TRUE),
        min            = ~min(., na.rm = TRUE),
        max            = ~max(., na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    )
  ) %>%
  pivot_longer(everything(),
               names_to = c("variable", "estat"),
               names_sep = "__",
               values_to = "valor") %>%
  pivot_wider(names_from = estat, values_from = valor) %>%
  arrange(desc(n))

print(resum_cont)

