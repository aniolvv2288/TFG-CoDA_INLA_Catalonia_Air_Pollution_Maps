# Aniol Vilamala Vidal

################################################################################
## ----- ALTITUD -----
################################################################################

# r <- geodata::elevation_30s(country = "ESP", path = tempdir())
# r <- terra::project(r, "EPSG:25831")
# r_cat <- terra::crop(r, vect(m))
# r_cat <- terra::mask(r_cat, vect(m))

altituds_catalunya <- read_csv("altituds_catalunya.csv")
dim(altituds_catalunya)

cat <- gisco_get_nuts(
  year = 2021,
  resolution = "3",
  nuts_level = 2,
  country = "ES"
) %>%
  filter(NUTS_NAME %in% c("Cataluña", "Catalonia")) %>%
  st_transform(4326)
cat_utm <- st_transform(cat, 25831)
ggplot(cat_utm) + 
  geom_sf(fill = "lightgreen", color = "black") + 
  theme_bw() +
  coord_sf(datum = st_crs(cat_utm))
m <- cat_utm


p <- st_as_sf(data.frame(long = altituds_catalunya$`Coordenades UTM (metres). UTMx`,
                         lat = altituds_catalunya$`Coordenades UTM (metres). UTMy`),
              coords = c("long", "lat"),
              crs = 25831)
altituds_catalunya[, c("x", "y")] <- st_coordinates(p)

## PLOT MAPA + ESTACIONS
ind <- st_intersects(m, p, sparse = FALSE)
altituds_catalunya <- altituds_catalunya[ind[[1]], ]

ggplot(m) + geom_sf() + coord_sf(datum = st_crs(m)) +
  geom_point(data = altituds_catalunya, aes(x = x, y = y)) + theme_bw()


## DATA FRAME DEPURAT
altituds_catalunya <- altituds_catalunya %>%
  rename("Altitud"=`Altitud (metres)`) %>%
  select(-`Coordenades UTM (metres). UTMx`, -`Coordenades UTM (metres). UTMy`)


################################################################################
## ----- ÀREA URBANA -----
################################################################################

area_urbana <- read_csv("area_urbana.csv")
area_urbana <- area_urbana %>%
  filter(Any %in% 2024) %>%
  mutate(
    AREA_URBANA = case_when(
      Degurba == 1 ~ "urban",
      Degurba == 2 ~ "suburban",
      Degurba == 3 ~ "rural"
    )
  ) %>%
  select(-Any, -Degurba)


## CREACIÓ DATA FRAME COVARIABLES

anys <- 2020:2024
mesos <- 1:12
temps <- expand.grid(Any = anys, Mes = mesos)

covariables <- left_join(altituds_catalunya, area_urbana, by="Codi")
covariables <- covariables %>% 
  dplyr::select(Codi, Nom, Altitud, AREA_URBANA, x, y)

covariables <- covariables %>%
  tidyr::crossing(temps) %>%
  dplyr::select(Codi, Nom, Any, Mes, Altitud, AREA_URBANA, x, y)
  

################################################################################
## ----- TEMPERATURES MÀXIMES -----
################################################################################

temp_max <- read_csv("temperatures_maximes.csv")

temp_max <- temp_max %>%
  select(-`Media mensual. Media anual`) %>%
  
  pivot_longer(
    cols = starts_with("Media mensual"),
    names_to = "Mes",
    values_to = "Temp_max"
  ) %>%
  
  mutate(
    Mes = str_remove(Mes, "Media mensual\\. "),
    Mes = str_trim(Mes)
  ) %>%
  
  mutate(
    Mes = case_match(
      Mes,
      "Enero" ~ 1,
      "Febrero" ~ 2,
      "Marzo" ~ 3,
      "Abril" ~ 4,
      "Mayo" ~ 5,
      "Junio" ~ 6,
      "Julio" ~ 7,
      "Agosto" ~ 8,
      "Septiembre" ~ 9,
      "Octubre" ~ 10,
      "Noviembre" ~ 11,
      "Diciembre" ~ 12
    )
  )

## ASSIGNACIÓ DE POBLES A COMARCA:

municipis_comarca <- read_csv("mpiscatalunya.csv")

### Passar els municipis del Lluçanès a la comarca anterior a la seva formació

municipis_comarca <- municipis_comarca %>%
  mutate(
    Comarca = if_else(Comarca == "Lluçanès", "Osona", Comarca)
  ) %>%
  select(-`Codi comarca`)

###

temp_max <- temp_max %>%
  left_join(
    municipis_comarca %>% select(Nom, Comarca),
    by = "Comarca"
  )

## AFEGIR A DATA FRAME COVARIABLES:

covariables <- covariables %>%
  left_join(
    temp_max,
    by = c("Nom", "Any", "Mes")
  )

covariables <- covariables %>%
  dplyr::select(Codi, Nom, Any, Mes, Altitud, AREA_URBANA, Temp_max, 
                x, y)


################################################################################
## ----- TEMPERATURES MÍNIMES -----
################################################################################

temp_min <- read_csv("temperatures_minimes.csv")

temp_min <- temp_min %>%
  select(-`Mitjana mensual. Mitjana anual`) %>%
  
  pivot_longer(
    cols = starts_with("Mitjana mensual"),
    names_to = "Mes",
    values_to = "Temp_min"
  ) %>%
  
  mutate(
    Mes = str_remove(Mes, "Mitjana mensual\\. "),
    Mes = str_trim(Mes)
  ) %>%
  
  mutate(
    Mes = case_match(
      Mes,
      "Gener" ~ 1,
      "Febrer" ~ 2,
      "Març" ~ 3,
      "Abril" ~ 4,
      "Maig" ~ 5,
      "Juny" ~ 6,
      "Juliol" ~ 7,
      "Agost" ~ 8,
      "Setembre" ~ 9,
      "Octubre" ~ 10,
      "Novembre" ~ 11,
      "Desembre" ~ 12
    )
  )

## ASSIGNACIÓ DE POBLES A COMARCA:

municipis_comarca <- read_csv("mpiscatalunya.csv")

temp_min <- temp_min %>%
  left_join(
    municipis_comarca %>% select(Nom, Comarca),
    by = "Comarca"
  )

## AFEGIR A DATA FRAME COVARIABLES:

covariables <- covariables %>%
  left_join(
    temp_min,
    by = c("Nom", "Any", "Mes")
  )

covariables <- covariables %>%
  dplyr::select(Codi, Nom, Any, Mes, Altitud, AREA_URBANA, Temp_max, Temp_min, 
                x, y)


################################################################################
## ----- PLUVIOMETRIA -----
################################################################################

pluviometria <- read_csv("pluviometria.csv")

pluviometria <- pluviometria %>%
  dplyr::select(-Total) %>%
  
  pivot_longer(
    cols = c(Gener, Febrer, Març, Abril, Maig, Juny, Juliol, Agost,
             Setembre, Octubre, Novembre, Desembre),
    names_to = "Mes",
    values_to = "Pluja"
  ) %>%
  
  mutate(
    Mes = case_match(
      Mes,
      "Gener" ~ 1,
      "Febrer" ~ 2,
      "Març" ~ 3,
      "Abril" ~ 4,
      "Maig" ~ 5,
      "Juny" ~ 6,
      "Juliol" ~ 7,
      "Agost" ~ 8,
      "Setembre" ~ 9,
      "Octubre" ~ 10,
      "Novembre" ~ 11,
      "Desembre" ~ 12
    )
  )

## ASSIGNACIÓ DE POBLES A COMARCA

pluviometria <- pluviometria %>%
  left_join(
    municipis_comarca %>% select(Nom, Comarca),
    by = "Comarca"
  )

## AFEGIR A DATA FRAME COVARIABLES:

covariables <- covariables %>%
  left_join(
    pluviometria,
    by = c("Nom", "Any", "Mes")
  )

covariables <- covariables %>%
  dplyr::select(Codi, Nom, Any, Mes, Altitud, AREA_URBANA, Temp_max, Temp_min, 
                Pluja, x, y) %>%
  mutate(
    SEASON = case_when(
      Mes %in% c(12, 1, 2) ~ "Hivern",
      Mes %in% c(3, 4, 5)  ~ "Primavera",
      Mes %in% c(6, 7, 8)  ~ "Estiu",
      Mes %in% c(9, 10, 11) ~ "Tardor"
    ),
    SEASON = factor(
      SEASON,
      levels = c("Hivern", "Primavera", "Estiu", "Tardor"),
      ordered = TRUE
    )
  ) %>%
  relocate(SEASON, .after = Mes)

# saveRDS(covariables, "covariables.rds")


################################################################################
## ----- ACTUALITZACIÓ DADES FINALS AMB LES NOVES VARIABLES -----
################################################################################

## AFEGIR MUNICIPI I COMARCA DE CADA ESTACIÓ

dades_finals <- readRDS("dades_finals.rds")

XVPCA_2020_2024 <- read_csv("XVPCA_2020_2024.csv")

info_estacions <- XVPCA_2020_2024 %>%
  select(`NOM ESTACIO`, MUNICIPI, Comarca = `NOM COMARCA`) %>%
  distinct()

dades_finals <- dades_finals %>%
  left_join(info_estacions, by = "NOM ESTACIO") %>%
  mutate(
    MUNICIPI = ifelse(
      MUNICIPI == "Vandellòs i l'Hospitalet de l'",
      "Vandellòs i l'Hospitalet de l'Infant",
      MUNICIPI
    )
  ) %>%
  dplyr::select(
    ID, MUNICIPI, Comarca, `NOM ESTACIO`, ANY, MES, TIPUS_ESTACIO,
    AREA_URBANA, ALTITUD, LONGITUD, LATITUD,
    O3, NO2, CO, SO2, PM10, PM2.5
  )

## AFEGIR COVARIABLES

dades_finals <- dades_finals %>%
  left_join(
    covariables %>% 
      select(Nom, Any, Mes, Temp_max, Temp_min, Pluja),
    by = c("MUNICIPI" = "Nom", "ANY" = "Any", "MES" = "Mes")
  ) %>%
  relocate(Temp_max, Temp_min, Pluja, .after = ALTITUD) %>%
  relocate(SEASON, .after = MES)

# saveRDS(dades_finals, "dades_finals.rds")

