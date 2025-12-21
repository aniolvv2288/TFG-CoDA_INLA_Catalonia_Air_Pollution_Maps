## Aniol Vilamala Vidal

################################################################################
# ----- MAPA -----
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
# ----- INCORPORACIÓ DADES -----
################################################################################

## BASE DE DADES

dades_finals <- readRDS("dades_finals.rds")

dades_finals <- dades_finals %>%
  dplyr::select(-OUTLIER, -OUTLIER_BIN, -OUTLIER_QUIN, -n_originals)

d <- dades_finals %>%
  dplyr::select(`NOM ESTACIO`, ANY, MES, ALTITUD, AREA_URBANA, TIPUS_ESTACIO,
                LATITUD, LONGITUD, NO2) %>%
  dplyr::group_by(`NOM ESTACIO`, ANY, ALTITUD, AREA_URBANA, TIPUS_ESTACIO,
                  LATITUD, LONGITUD) %>%
  dplyr::summarise(NO2 = mean(NO2, na.rm = TRUE), .groups = "drop")

r <- geodata::elevation_30s(country = "ESP", path = tempdir())
r <- terra::project(r, "EPSG:25831")
r_cat <- terra::crop(r, vect(m))
r_cat <- terra::mask(r_cat, vect(m))

## LOCALITZACIÓ ESTACIONS

p <- st_as_sf(data.frame(long = d$LONGITUD, lat = d$LATITUD),
              coords = c("long", "lat"))
st_crs(p) <- st_crs(4326)
p <- p %>% st_transform(25831)
d[, c("x", "y")] <- st_coordinates(p)

## PLOT MAPA + ESTACIONS

ind <- st_intersects(m, p)
d <- d[ind[[1]], ]

ggplot(m) + geom_sf() + coord_sf(datum = st_crs(m)) +
  geom_point(data = d, aes(x = x, y = y)) + theme_bw()

################################################################################
# ----- CONSTRUCCIÓ DE LA MALLA -----
################################################################################

library(INLA)

coo <- cbind(d$x, d$y)
bnd <- inla.nonconvex.hull(st_coordinates(m)[, 1:2])

mesh <- inla.mesh.2d(
  loc = coo, boundary = bnd,
  max.edge = c(25000, 100000), cutoff = 5500
)

mesh$n

plot(mesh)
points(coo, col = "red")


################################################################################
# ----- MODEL SPDE SOBRE LA MALLA + CAMP ESPACIAL + MATRIU DE PROJECCIÓ -----
################################################################################

spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2, constr = TRUE,
  prior.range = c(10000, 0.01), # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01
)

timesn <- length(unique(d$ANY))

indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn)

group <- d$ANY - min(d$ANY) + 1
A <- inla.spde.make.A(mesh = mesh, loc = coo, group = group)


################################################################################
# ----- DADES DE PREDICCIÓ -----
################################################################################

bb <- st_bbox(m)
x <- seq(bb$xmin - 1, bb$xmax + 1, length.out = 35)
y <- seq(bb$ymin - 1, bb$ymax + 1, length.out = 35)
dp <- as.matrix(expand.grid(x, y))
plot(dp, asp = 1)

p <- st_as_sf(data.frame(x = dp[, 1], y = dp[, 2]),
              coords = c("x", "y")
)
st_crs(p) <- st_crs(25831)
ind <- st_intersects(m, p)
dp <- dp[ind[[1]], ]
plot(dp, asp = 1)

dp <- rbind(cbind(dp, 1), cbind(dp, 2), cbind(dp, 3), cbind(dp, 4), cbind(dp, 5))
dim(dp)

coop <- dp[, 1:2]
groupp <- dp[, 3]
Ap <- inla.spde.make.A(mesh = mesh, loc = coop, group = groupp)

# dp <- as.data.frame(dp)
# dp$altitud <- extract(r_cat, coop)[, 1]


################################################################################
# ----- STACKS DE PREDICCIÓ -----
################################################################################

stk.e <- inla.stack(
  tag = "est",
  data = list(y = d$NO2),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(d))), s = indexs)
)

stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(dp))), s = indexs)
)

stk.full <- inla.stack(stk.e, stk.p)



################################################################################
# ----- FÓRMULA DEL MODEL -----
################################################################################

rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))

formula <- y ~ 0 + b0 + f(s,
                          model = spde, group = s.group,
                          control.group = list(model = "ar1", hyper = rprior)
)

res <- inla(formula,
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE,
              A = inla.stack.A(stk.full)
            )
)

################################################################################
# ----- RESULTATS -----
################################################################################

summary(res)


################################################################################
# ----- MAPA -----
################################################################################

index <- inla.stack.index(stack = stk.full, tag = "pred")$data

dp <- data.frame(dp)
names(dp) <- c("x", "y", "time")

dp$pred_mean <- res$summary.fitted.values[index, "mean"]
dp$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
dp$pred_ul <- res$summary.fitted.values[index, "0.975quant"]

library(reshape2)
dpm <- melt(dp,
            id.vars = c("x", "y", "time"),
            measure.vars = c("pred_mean", "pred_ll", "pred_ul")
)

dpm <- dpm %>%
  mutate(year = case_when(
    time == 1 ~ 2020,
    time == 2 ~ 2021,
    time == 3 ~ 2022,
    time == 4 ~ 2023,
    time == 5 ~ 2024
  ))

ggplot(m) + 
  geom_sf() + 
  coord_sf(datum = NA) +
  geom_tile(data = dpm, aes(x = x, y = y, fill = value)) +
  labs(x = "", y = "") +
  facet_grid(year ~ variable) + 
  scale_fill_gradientn(
    name = "NO2",
    colours = c("green", "yellow", "red")
  ) +
  theme_bw()

