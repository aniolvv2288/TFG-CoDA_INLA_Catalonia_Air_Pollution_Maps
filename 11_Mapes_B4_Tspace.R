# Aniol Vilamala Vidal

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

ggplot(cat) + 
  geom_sf(fill = "lightblue", color = "black") + 
  theme_bw()

cat_utm <- st_transform(cat, 25831)

ggplot(cat_utm) + 
  geom_sf(fill = "lightgreen", color = "black") + 
  theme_bw() +
  coord_sf(datum = st_crs(cat_utm))

m <- cat_utm


################################################################################
# ----- INCORPORACIÓ DADES -----
################################################################################

dades_finals <- readRDS("dades_finals.rds")
d <- dades_finals
contaminants <- c("O3", "NO2", "CO", "SO2", "PM10", "PM2.5")

## TOTAL T-SPACE

X <- dades_finals %>%
  dplyr::select(all_of(contaminants))

D <- ncol(X)
Total_basal <- (1 / sqrt(D)) * rowSums(log(X))

d$Total_basal <- Total_basal

## BALANCES ILR

V <- matrix(0, nrow = 6, ncol = 5)
rownames(V) <- contaminants
colnames(V) <- paste0("B", 1:5)

V["SO2", "B1"] <-  1
V[c("NO2","PM10","PM2.5","O3","CO"), "B1"] <- -1

V[c("NO2","O3","CO"), "B2"] <-  1
V[c("PM10","PM2.5"), "B2"] <- -1

V[c("NO2","O3"), "B3"] <-  1
V["CO", "B3"] <- -1

V["NO2", "B4"] <-  1
V["O3",  "B4"] <- -1

V["PM2.5", "B5"] <-  1
V["PM10",  "B5"] <- -1

V

X <- d %>%
  dplyr::select(all_of(contaminants)) %>%
  acomp

bal_base <- gsi.buildilrBase(V)

crossprod(bal_base)

ilr_bal <- ilr(X, V=bal_base)

ilr_df <- as.data.frame(ilr_bal)

d$id <- seq_len(nrow(d))
ilr_df$id <- seq_len(nrow(ilr_df))

d <- d %>% 
  left_join(ilr_df, by = "id") %>%
  dplyr::select(-all_of(contaminants), -id) %>%
  rename(Total = Total_basal)

## RASTER ALTITUD

r <- geodata::elevation_30s(country = "ESP", path = tempdir())
r <- terra::project(r, "EPSG:25831")
r_cat <- terra::crop(r, vect(m))
r_cat <- terra::mask(r_cat, vect(m))

## COVARIABLES

covariables <- readRDS("covariables.rds")
covariables <- covariables %>%
  dplyr::select(-Altitud)

## LOCALITZACIÓ ESTACIONS

p <- st_as_sf(data.frame(long = d$LONGITUD, lat = d$LATITUD),
              coords = c("long", "lat"))
st_crs(p) <- st_crs(4326)
p <- p %>% st_transform(25831)
d[, c("x", "y")] <- st_coordinates(p)


################################################################################
# ----- CONSTRUCCIÓ DE LA MALLA -----
################################################################################

coo <- cbind(d$x, d$y)
bnd <- fmesher::fm_nonconvex_hull(st_coordinates(m)[, 1:2])

mesh <- fmesher::fm_mesh_2d_inla(loc = coo, boundary = bnd,
                                 max.edge = c(30000, 70000), cutoff = 1500)

mesh$n           
par(mfrow=c(1,1))
plot(mesh, main = "Mesh")
points(coo, col = "red", pch = 20)


################################################################################
# ----- MODEL SPDE SOBRE LA MALLA + CAMP ESPACIAL + MATRIU DE PROJECCIÓ -----
################################################################################

spde <- inla.spde2.matern(mesh, alpha=2)

d$YM <- d$ANY * 12 + d$MES

timesn <- length(unique(d$YM))

indexs <- inla.spde.make.index(
  name = "s",
  n.spde = spde$n.spde,
  n.group = timesn
)

group <- d$YM - min(d$YM) + 1

A <- inla.spde.make.A(
  mesh = mesh,
  loc = coo,
  group = group
)


################################################################################
# ----- DADES DE PREDICCIÓ (dp) AMB ESPAI–TEMPS COMPLET (ANY–MES) -------------
################################################################################

bb <- st_bbox(m)
x <- seq(bb$xmin - 1, bb$xmax + 1, length.out = 100)
y <- seq(bb$ymin - 1, bb$ymax + 1, length.out = 100)
dp <- as.matrix(expand.grid(x, y))

p <- st_as_sf(data.frame(x = dp[, 1], y = dp[, 2]),
              coords = c("x", "y"), crs = 25831)
ind <- st_intersects(m, p)
dp <- dp[ind[[1]], ]
dp <- as.data.frame(dp)
names(dp) <- c("Var1", "Var2")

p_pred <- st_as_sf(dp, coords = c("Var1", "Var2"), crs = 25831)
p_pred_v <- terra::vect(p_pred)
alt_raster <- terra::extract(r_cat, p_pred_v)
coords <- sf::st_coordinates(p_pred)

df_alt <- data.frame(
  Var1 = coords[, 1],
  Var2 = coords[, 2],
  Altitud = alt_raster[, 2]
)

dp <- dp %>%
  left_join(df_alt, by = c("Var1", "Var2")) %>%
  mutate(Altitud = pmax(Altitud, 0)) %>%
  filter(!is.na(Altitud)) %>%
  distinct(Var1, Var2, .keep_all = TRUE)

dp <- tidyr::crossing(
  dp,
  ANY = 2020:2024,
  MES = 1:12
)

cov_pts <- covariables %>%
  distinct(Codi, x, y) %>%
  st_as_sf(coords = c("x", "y"), crs = 25831)

p_pred <- st_as_sf(dp, coords = c("Var1", "Var2"), crs = 25831)
nearest <- st_nearest_feature(p_pred, cov_pts)
dp$Codi <- cov_pts$Codi[nearest]

dp <- dp %>%
  mutate(YM = ANY * 12 + MES)

covariables2 <- covariables %>%
  mutate(YM = Any * 12 + Mes)

dp <- dp %>%
  left_join(
    covariables2,
    by = c("Codi", "YM")
  ) %>%
  select(
    Var1, Var2, ANY, MES, YM, Codi,
    Nom, Altitud, AREA_URBANA, x, y,
    Temp_max, Temp_min, Pluja
  ) %>%
  group_by(
    Var1, Var2, ANY, MES, YM, Codi,
    Nom, Altitud, AREA_URBANA, x, y
  ) %>%
  summarise(
    Temp_max = mean(Temp_max, na.rm = TRUE),
    Temp_min = mean(Temp_min, na.rm = TRUE),
    Pluja    = mean(Pluja, na.rm = TRUE),
    .groups = "drop"
  )

coop   <- as.matrix(dp[, c("Var1", "Var2")])
groupp <- dp$YM - min(dp$YM) + 1

Ap <- inla.spde.make.A(
  mesh = mesh,
  loc  = coop,
  group = groupp
)


################################################################################
# ----- STACKS DE PREDICCIÓ -----
################################################################################

stk.e <- inla.stack(
  tag = "est",
  data = list(y = d$B4),
  A = list(1, A),
  effects = list(
    data.frame(b0 = rep(1, nrow(d)), 
               ALTITUD = d$ALTITUD,
               MUNICIPI = as.factor(d$MUNICIPI),
               AREA_URBANA = d$AREA_URBANA,
               TEMP_MAX = d$Temp_max,
               TEMP_MIN = d$Temp_min),
    s = indexs
  )
)

stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(
    data.frame(b0 = rep(1, nrow(dp)), 
               ALTITUD = dp$Altitud,                  
               MUNICIPI = as.factor(dp$Nom),
               AREA_URBANA = dp$AREA_URBANA,
               TEMP_MAX = dp$Temp_max,
               TEMP_MIN = dp$Temp_min),
    s = indexs
  )
)

stk.full <- inla.stack(stk.e, stk.p)


################################################################################
# ----- FÓRMULA DEL MODEL -----
################################################################################

rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.85)))

formula <- y ~ -1 + b0 + ALTITUD + TEMP_MAX + TEMP_MIN + 
  f(s, model = spde, group = s.group,
    control.group = list(model = "ar1", hyper = rprior)) + 
  f(MUNICIPI, model = "iid")

if (!file.exists("model_mes_any_inla_b4.rds")) {
  resultat <- inla(
    formula,
    data = inla.stack.data(stk.full),
    control.predictor = list(
      compute = TRUE,
      A = inla.stack.A(stk.full)
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      config = TRUE
    )
  )
  
  saveRDS(resultat, "model_mes_any_inla_b4.rds")
  
} else {
  resultat <- readRDS("model_mes_any_inla_b4.rds")
}


################################################################################
# ----- RESULTATS -----
################################################################################

summary(resultat)

## DISTRIBUCIÓ DE LA MARGINAL DE LA POSTERIOR PELS EFECTES FIXES

par(mfrow=c(2,2))
plot(resultat$marginals.fixed[[1]], ty = "l", xlab = expression(beta[0]), 
     ylab = "Density") 
plot(resultat$marginals.fixed[[2]], ty = "l", xlab = expression(beta[Altitud]), 
     ylab = "Density") 
plot(resultat$marginals.fixed[[3]], ty = "l", xlab = expression(beta[Temp_max]), 
     ylab = "Density") 
plot(resultat$marginals.fixed[[4]], ty = "l", xlab = expression(beta[Temp_min]), 
     ylab = "Density") 

## DISTRIBUCIÓ DE LA MARGINAL DE LA POSTERIOR PELS HIPERPARÀMETRES

resultat.res<-inla.spde2.result(resultat, 's', spde, do.transf=TRUE)

par(mfrow=c(1,3))
plot(resultat.res$marginals.variance.nominal[[1]], ty = "l", 
     xlab = expression(sigma[randomfield]^2), ylab = "Density")  # Contra rapidesa a la que cau correlació camp espacial
plot(resultat.res$marginals.kappa[[1]], type = "l", 
     xlab = expression(kappa), ylab = "Density")                 # Variància del camp espacial
plot(resultat.res$marginals.range.nominal[[1]], type = "l", 
     xlab = "range nominal", ylab = "Density")                   # Distància on les observacions ja no estan correlacionades
par(mfrow=c(1,1))

################################################################################
# ----- VALIDACIÓ DEL MODEL -----
################################################################################

dades <- d %>%
  mutate(YM = ANY * 12 + MES) %>%
  select(MUNICIPI, YM, ANY, MES, ALTITUD, Temp_max, Temp_min, AREA_URBANA, Pluja,
         x, y, B4)

coords <- cbind(dades$x, dades$y)

## DADES EN TRAIN-TEST

smp_size <- floor(0.75 * nrow(dades))
set.seed(1714)
train_ind <- sample(seq_len(nrow(dades)), size = smp_size, replace = FALSE)

train <- dades[train_ind, ]
test <- dades[-train_ind, ]
test$B4 <- NA

train_coords <- coords[train_ind,]
test_coords <- coords[-train_ind,]

## MATRIUS I STACKS

group_train <- train$YM - min(train$YM) + 1
group_test  <- test$YM  - min(test$YM) + 1

Ae <- inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(train_coords),
  group = group_train
)

App <- inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(test_coords),
  group = group_test
)

stk.e <- inla.stack(
  tag = "est",
  data = list(y = train$B4),
  A = list(1, Ae),
  effects = list(
    data.frame(b0 = rep(1, nrow(train)), 
               ALTITUD = train$ALTITUD,
               MUNICIPI = as.factor(train$MUNICIPI),
               AREA_URBANA = train$AREA_URBANA,
               TEMP_MAX = train$Temp_max,
               TEMP_MIN = train$Temp_min),
    s = indexs
  )
)

stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, App),
  effects = list(
    data.frame(b0 = rep(1, nrow(test)), 
               ALTITUD = test$ALTITUD,                        
               MUNICIPI = as.factor(test$MUNICIPI),
               AREA_URBANA = test$AREA_URBANA,
               TEMP_MAX = test$Temp_max,
               TEMP_MIN = test$Temp_min),
    s = indexs
  )
)

stk.f <- inla.stack(stk.e, stk.p)

## FÓRMULA I MODEL

prior_any <- list(theta = list(prior = "pccor1", param = c(0, 0.85)))

formula <- y ~ -1 + b0 + ALTITUD + TEMP_MAX + TEMP_MIN +
  f(s, model = spde, group = s.group,
    control.group = list(model = "ar1", hyper = prior_any)) + 
  f(MUNICIPI, model = "iid")

if (!file.exists("val_model_mes_any_inla_b4.rds")) {
  p.res <- inla(
    formula,
    data = inla.stack.data(stk.f),
    control.predictor = list(
      compute = TRUE,
      A = inla.stack.A(stk.f)
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      config = TRUE
    )
  )
  
  saveRDS(p.res, "val_model_mes_any_inla_b4.rds")
  
} else {
  p.res <- readRDS("val_model_mes_any_inla_b4.rds")
}


## PREDICCIONS

index.pred <- inla.stack.index(stk.f, "pred")$data
pred <- p.res$summary.linear.predictor[index.pred, "mean"]
obs <- dades[-train_ind, ]$B4

corr <- cor(obs, pred, use = "complete.obs")
rmse <- sqrt(mean((pred - obs)^2, na.rm = TRUE))
mae  <- mean(abs(pred - obs), na.rm = TRUE)

df_plot <- data.frame(
  obs  = obs,
  pred = pred
)

metrics_text <- sprintf("Correlació = %.3f\nRMSE = %.3f\nMAE = %.3f",
                        corr, rmse, mae)

ggplot(df_plot, aes(x = obs, y = pred)) +
  geom_point(size = 2, alpha = 0.8, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  annotate("label",
           x = min(df_plot$obs, na.rm = TRUE),
           y = max(df_plot$pred, na.rm = TRUE),
           label = metrics_text,
           hjust = 0, vjust = 1,
           size = 4.5,
           fill = "white",        
           alpha = 0.9) +         
  labs(
    x = "Observats",
    y = "Predits",
    title = "Comparació entre valors observats i predits"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    plot.title = element_text(hjust = 0.5)
  )


################################################################################
# ----- MAPA -----
################################################################################

if (!file.exists("prediccio_model_inla_mes_any_b4.rds")) {
  
  index <- inla.stack.index(stack = stk.full, tag = "pred")$data
  
  dp <- data.frame(dp)
  names(dp) <- c("x", "y", "any", "mes", "time")
  
  dp$pred_mean <- resultat$summary.fitted.values[index, "mean"]
  dp$pred_ll   <- resultat$summary.fitted.values[index, "0.025quant"]
  dp$pred_ul   <- resultat$summary.fitted.values[index, "0.975quant"]
  
  dpm <- reshape2::melt(
    dp,
    id.vars = c("x", "y", "time"),
    measure.vars = c("pred_mean", "pred_ll", "pred_ul")
  )
  
  t0 <- min(dpm$time)
  
  dpm2 <- dpm %>%
    mutate(
      mes_index = time - t0,
      date = as.Date("2020-01-01") %m+% months(mes_index)
    )
  
  saveRDS(dpm2, "prediccio_model_inla_mes_any_b4.rds")
  
} else {
  
  dpm2 <- readRDS("prediccio_model_inla_mes_any_b4.rds")
}

dpm_febrer <- dpm2 %>%
  filter(month(date) == 2) %>%
  mutate(
    Any = year(date),
    variable = factor(variable, levels = c("pred_mean", "pred_ll", "pred_ul"))
  )

ggplot(m) +
  geom_sf() +
  coord_sf(datum = NA) +
  geom_tile(
    data = dpm_febrer,
    aes(x = x, y = y, fill = value)
  ) +
  facet_grid(
    rows = vars(Any),
    cols = vars(variable),
    labeller = labeller(
      variable = c(
        pred_mean = "Mitjana\nbalança 4",
        pred_ll   = "Límit inferior\n(IC 95%)",
        pred_ul   = "Límit superior\n(IC 95%)"
      )
    )
  ) +
  scale_fill_gradientn(
    name = "Coordenades \nB4 T-space",
    colours = c("green", "yellow", "red"),
    limits = c(min(dpm2$value), max(dpm2$value)),
    breaks = scales::pretty_breaks(n = 5),
    oob = scales::squish
  ) +
  labs(
    title = "Component B4 del T-space al febrer (2020–2024)",
    x = "", y = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0, "lines"),
    plot.title = element_text(
      size = 15
    )
  ) +
  guides(fill = guide_colourbar(title.position = "top"))


dpm_setembre <- dpm2 %>%
  filter(month(date) == 9) %>%
  mutate(
    Any = year(date),
    variable = factor(variable, levels = c("pred_mean", "pred_ll", "pred_ul"))
  )

ggplot(m) +
  geom_sf() +
  coord_sf(datum = NA) +
  geom_tile(
    data = dpm_setembre,
    aes(x = x, y = y, fill = value)
  ) +
  facet_grid(
    rows = vars(Any),
    cols = vars(variable),
    labeller = labeller(
      variable = c(
        pred_mean = "Mitjana\nbalança 4",
        pred_ll   = "Límit inferior\n(IC 95%)",
        pred_ul   = "Límit superior\n(IC 95%)"
      )
    )
  ) +
  scale_fill_gradientn(
    name = "Coordenades \nB4 T-space",
    colours = c("green", "yellow", "red"),
    limits = c(min(dpm2$value), max(dpm2$value)),
    breaks = scales::pretty_breaks(n = 5),
    oob = scales::squish
  ) +
  labs(
    title = "Component B4 del T-space al setembre (2020–2024)",
    x = "", y = ""
  ) +
  theme_bw() + 
  theme(
    legend.position = "right",
    legend.box = "vertical",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0, "lines"),
    plot.title = element_text(
      size = 15
    )
  ) +
  guides(fill = guide_colourbar(title.position = "top"))

