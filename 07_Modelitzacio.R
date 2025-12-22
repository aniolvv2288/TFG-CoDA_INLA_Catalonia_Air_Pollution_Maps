# Aniol Vilamala Vidal

################################################################################
## ----- MAPA -----
################################################################################

cat <- gisco_get_nuts(
  year = 2021,
  resolution = "3",
  nuts_level = 2,
  country = "ES"
) %>%
  filter(NUTS_NAME %in% c("Cataluña", "Catalonia")) %>%
  st_transform(4326)

# Plot simple en WGS84
ggplot(cat) + 
  geom_sf(fill = "lightblue", color = "black") + 
  theme_bw()

# Transformar a UTM 31N (EPSG:25831)
cat_utm <- st_transform(cat, 25831)

ggplot(cat_utm) + 
  geom_sf(fill = "lightgreen", color = "black") + 
  theme_bw() +
  coord_sf(datum = st_crs(cat_utm))

m <- cat_utm

################################################################################
## ----- INCORPORACIÓ DADES I CREACIÓ DE LA VARIABLE RESPOSTA -----
################################################################################

dades_finals <- readRDS("dades_finals.rds")
d <- dades_finals
contaminants <- c("NO2", "PM10", "PM2.5", "O3", "CO", "SO2")

## VARIABLE RESPOSTA

X <- dades_finals %>%
  dplyr::select(all_of(contaminants))

D <- ncol(X)
Total_basal <- (1 / sqrt(D)) * rowSums(log(X))

d$Total_basal <- Total_basal

## LOCALITZACIÓ ESTACIONS

p <- st_as_sf(data.frame(long = d$LONGITUD, lat = d$LATITUD),
              coords = c("long", "lat"))
st_crs(p) <- st_crs(4326)
p <- p %>% st_transform(25831)
d[, c("x", "y")] <- st_coordinates(p)

################################################################################
## ----- VARIABLE TEMPERATURA -----
################################################################################

cor.test(d$Temp_max, d$Temp_min) # Correlació molt alta: 0.8903017

d$Amplitud_termica <- d$Temp_max - d$Temp_min
d$Temp_mitj <- (d$Temp_max + d$Temp_min)/2


################################################################################
## ----- CONSTRUCCIÓ DE LA MALLA -----
################################################################################

coo <- cbind(d$x, d$y)
bnd <- fmesher::fm_nonconvex_hull(st_coordinates(m)[, 1:2])

mesh_specs <- list(
  list(max.edge = c(10000, 35000), cutoff = 4000),   
  list(max.edge = c(15000, 45000), cutoff = 4500),  
  list(max.edge = c(19000, 50000), cutoff = 5000),  
  list(max.edge = c(20000, 60000), cutoff = 5000),   
  list(max.edge = c(25000, 70000), cutoff = 5000),  
  list(max.edge = c(30000, 70000), cutoff = 7000)    
)

meshes <- lapply(mesh_specs, function(spec) {
  fmesher::fm_mesh_2d_inla(
    loc = coo,
    boundary = bnd,
    max.edge = spec$max.edge,
    cutoff = spec$cutoff
  )
})

par(mfrow = c(2, 3), mar = c(2, 2, 2, 2))

for (i in seq_along(meshes)) {
  plot(meshes[[i]], main = paste("Mesh", i))
  #points(coo, col = "red", pch = 20)
}

sapply(meshes, function(m) m$n)

mesh <- fmesher::fm_mesh_2d_inla(loc = coo, boundary = bnd,
                                 max.edge = c(30000, 70000),
                                 cutoff = 7000)
par(mfrow = c(1,1))
plot(mesh, main = "Mesh")
points(coo, col = "red", pch = 20)

################################################################################
## ----- MODEL SPDE SOBRE LA MALLA + CAMP ESPACIAL + MATRIU DE PROJECCIÓ -----
################################################################################

# spde <- inla.spde2.pcmatern(
#   mesh = mesh, alpha = 2, constr = TRUE,
#   prior.range = c(10000, 0.01), # P(range < 10000) = 0.01
#   prior.sigma = c(3, 0.01)      # P(sigma > 3) = 0.01
# )

spde <- inla.spde2.matern(mesh, alpha=2)

indexs <- inla.spde.make.index("s", n.spde = spde$n.spde)

A <- inla.spde.make.A(mesh = mesh, loc = coo)

################################################################################
## ----- STACK -----
################################################################################

stk <- inla.stack(
  tag = "dat",
  data=list(y = d$Total_basal), 
  A=list(A,1),
  effects=list(c(list(b0 = 1),
                 indexs),
               list(ALTITUD = d$ALTITUD,
                    AREA_URBANA = d$AREA_URBANA,
                    ANY = as.integer(d$ANY),
                    MES = as.integer(d$MES),
                    SEASON = as.integer(d$SEASON),
                    MUNICIPI = as.factor(d$MUNICIPI),
                    TEMP_MAX = d$Temp_max,
                    TEMP_MIN = d$Temp_min,
                    PLUJA = d$Pluja,
                    TEMP_MITJ = d$Temp_mitj,
                    AMPL_TERM = d$Amplitud_termica
                    )
               )
  )


################################################################################
## ----- MODELS -----
################################################################################

## EFECTES FIXES + CAMP ESPACIAL

formula1 <- y ~ -1 + b0 + f(s, model = spde)

formula2 <- y ~ -1 + b0 + ALTITUD + f(s, model = spde)

formula3 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + f(s, model = spde)

formula4 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + f(s, model = spde)

formula5 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN + 
  f(s, model = spde)

formula6 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN + PLUJA +
  f(s, model = spde)

## INCORPORACIÓ EFECTES ALEATORIS D'ANY I ESTACIÓ

formula7 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN + PLUJA +
  f(s, model = spde) + f(MUNICIPI, model = "iid")

formula8 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN + PLUJA +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1")

formula9 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN + PLUJA +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "seasonal", season.length = 4)

formula10 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN + PLUJA +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "rw1", cyclic = TRUE)

formula11 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA + TEMP_MAX + TEMP_MIN +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "seasonal", season.length = 4)

formula12 <- y ~ -1 + b0 + ALTITUD + AREA_URBANA +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "seasonal", season.length = 4) + 
  f(TEMP_MAX, model = "rw1") + f(TEMP_MIN, model = "rw1")

formula13 <- y ~ -1 + b0 + ALTITUD +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "seasonal", season.length = 4) + 
  f(TEMP_MAX, model = "rw1") + f(TEMP_MIN, model = "rw1")

formula14 <- y ~ -1 + b0 + AREA_URBANA +
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "seasonal", season.length = 4) + 
  f(TEMP_MAX, model = "rw1") + f(TEMP_MIN, model = "rw1")

formula15 <- y ~ -1 + b0 + 
  f(s, model = spde) + f(MUNICIPI, model = "iid") + f(ANY, model = "ar1") +
  f(SEASON, model = "seasonal", season.length = 4) + 
  f(TEMP_MAX, model = "rw1") + f(TEMP_MIN, model = "rw1")

## INCORPORACIÓ EFECTES ALEATORIS AMB MES


################################################################################
## ----- RESULTATS DELS MODELS -----
################################################################################

models <- list()

for (i in 1:15) {
  f <- get(paste0("formula", i))
  m <- inla(
    f,
    data = inla.stack.data(stk, spde = spde),
    control.predictor = list(compute = TRUE, A = inla.stack.A(stk)),
    control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  models[[paste0("model", i)]] <- m
}

models <- lapply(models, summary)

extract_metrics <- function(m, name) {
  data.frame(
    model = name,
    DIC   = m$dic$dic,
    pD    = m$dic$p.eff,
    WAIC  = m$waic$waic,
    pWAIC = m$waic$p.eff,
    loglik = m$mlik[1,1],
    stringsAsFactors = FALSE
  )
}

metrics <- do.call(rbind, lapply(seq_along(models), function(i) {
  extract_metrics(models[[i]], paste0("model", i))
}))

metrics <- metrics[order(metrics$WAIC), ]
print(metrics)

