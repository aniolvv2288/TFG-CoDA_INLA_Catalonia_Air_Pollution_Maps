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

dades_finals <- readRDS("dades_finals.rds")
d <- dades_finals
contaminants <- c("NO2", "PM10", "PM2.5", "O3", "CO", "SO2")

## VARIABLE RESPOSTA

X <- dades_finals %>%
  dplyr::select(all_of(contaminants))

D <- ncol(X)
Total_basal <- (1 / sqrt(D)) * rowSums(log(X))

d$Total_basal <- Total_basal

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
                                 max.edge = c(30000, 70000),
                                 cutoff = 7000)

par(mfrow=c(1,1))
plot(mesh, main = "Mesh")
points(coo, col = "red", pch = 20)


################################################################################
# ----- MODEL SPDE SOBRE LA MALLA + CAMP ESPACIAL + MATRIU DE PROJECCIÓ -----
################################################################################

spde <- inla.spde2.matern(mesh, alpha=2)

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
x <- seq(bb$xmin - 1, bb$xmax + 1, length.out = 100)
y <- seq(bb$ymin - 1, bb$ymax + 1, length.out = 100)
dp <- as.matrix(expand.grid(x, y))
plot(dp, asp = 1)

p <- st_as_sf(data.frame(x = dp[, 1], y = dp[, 2]),
              coords = c("x", "y")
)
st_crs(p) <- st_crs(25831)
ind <- st_intersects(m, p)
dp <- dp[ind[[1]], ]
plot(dp, asp = 1)

## ALTITUD AMB RASTER

dp <- as.data.frame(dp)
p_pred <- st_as_sf(dp, coords = c("Var1", "Var2"), crs = 25831)
p_pred_v <- terra::vect(p_pred)
alt_raster <- terra::extract(r_cat, p_pred_v)
coords <- sf::st_coordinates(p_pred)

df_alt <- data.frame(
  x = coords[, 1],
  y = coords[, 2],
  Altitud = alt_raster[, 2]
)

df_alt <- df_alt %>%
  rename(Var1 = x,
         Var2 = y,
         Altitud_raster = Altitud)

# dp <- dp %>%
#   left_join(df_alt, by = c("Var1","Var2")) %>%
#   filter(!is.na(Altitud_raster)) %>%
#   filter(Altitud_raster >= 0) %>%
#   distinct(Var1, Var2, .keep_all = TRUE) %>%
#   rename(Altitud = Altitud_raster)

dp <- dp %>%
  left_join(df_alt, by = c("Var1", "Var2")) %>%
  mutate(
    Altitud_raster = ifelse(is.na(Altitud_raster), NA, pmax(Altitud_raster, 0))
  ) %>%
  filter(!is.na(Altitud_raster)) %>%
  distinct(Var1, Var2, .keep_all = TRUE) %>%
  rename(Altitud = Altitud_raster)

dp <- bind_rows(
  dp %>% mutate(V3 = 1),
  dp %>% mutate(V3 = 2),
  dp %>% mutate(V3 = 3),
  dp %>% mutate(V3 = 4),
  dp %>% mutate(V3 = 5)
)

dim(dp)[1]/5

dp <- dp %>% 
  relocate(V3, .after=Var2)

## UNIÓ COVARIABLES I DADES DE PREDICCIÓ
cov_pts <- covariables %>%
  distinct(Codi, x, y) %>%
  st_as_sf(coords = c("x","y"), crs = 25831)

dp <- as.data.frame(dp)
p_pred <- st_as_sf(dp, coords = c("Var1", "Var2"), crs = 25831)

nearest <- st_nearest_feature(p_pred, cov_pts)
dp$Codi <- cov_pts$Codi[nearest]

dp <- dp %>%
  mutate(
    Any = case_when(
      V3 == 1 ~ 2020,
      V3 == 2 ~ 2021,
      V3 == 3 ~ 2022,
      V3 == 4 ~ 2023,
      V3 == 5 ~ 2024,
    )
  )


################################################################################
# ----- BUCLE PER MES AMB CADA MODEL, VALIDACIÓ I GENERACIÓ DE MAPA -----
################################################################################

dp_original <- dp

MESOS <- 1:12

for (mes in MESOS) {
  
  cat("\n=============================\n")
  cat("   PROCESSANT MES:", mes, "\n")
  cat("=============================\n")
  
  dp_mes <- dp_original
  dp_mes$MES <- mes
  
  # Join amb covariables del mes
  dp_mes <- dp_mes %>%
    left_join(
      covariables %>% filter(Mes == mes),
      by = c("Codi", "Any")
    ) %>%
    select(
      Var1, Var2, V3, Codi, Any, MES,
      Nom, Altitud, AREA_URBANA, x, y,
      Temp_max, Temp_min, Pluja
    ) %>%
    group_by(
      Var1, Var2, V3, Codi, Any, MES,
      Nom, Altitud, AREA_URBANA, x, y
    ) %>%
    summarise(
      Temp_max = mean(Temp_max, na.rm = TRUE),
      Temp_min = mean(Temp_min, na.rm = TRUE),
      Pluja    = mean(Pluja, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Matriu A de predicció
  coop   <- as.matrix(dp_mes[, c("Var1", "Var2")])
  groupp <- dp_mes$V3
  
  Ap <- inla.spde.make.A(
    mesh = mesh,
    loc  = coop,
    group = groupp
  )
  
  ##############################################################################
  # STACK DE DADES
  ##############################################################################
  
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = d$Total_basal),
    A = list(1, A),
    effects = list(
      data.frame(
        b0          = 1,
        ALTITUD     = d$ALTITUD,     
        MUNICIPI    = as.factor(d$MUNICIPI),
        TEMP_MAX    = d$Temp_max,
        TEMP_MIN    = d$Temp_min
      ),
      s = indexs
    )
  )
  
  ##############################################################################
  # STACK DE PREDICCIÓ
  ##############################################################################
  
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA),
    A = list(1, Ap),
    effects = list(
      data.frame(
        b0          = 1,
        ALTITUD     = dp_mes$Altitud,
        MUNICIPI    = factor(dp_mes$Nom),
        TEMP_MAX    = dp_mes$Temp_max,
        TEMP_MIN    = dp_mes$Temp_min
      ),
      s = indexs
    )
  )
  
  stk.full <- inla.stack(stk.e, stk.p)
  
  ##############################################################################
  # MODEL
  ##############################################################################
  
  rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.85)))
  
  formula <- y ~ -1 + b0 + ALTITUD + TEMP_MAX + TEMP_MIN +
    f(s, model = spde, group = s.group,
      control.group = list(model = "ar1", hyper = rprior)) +
    f(MUNICIPI, model = "iid") 
  
  resultat <- inla(
    formula,
    data = inla.stack.data(stk.full),
    control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)),
    control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  
  ##############################################################################
  # RESULTATS MODEL
  ##############################################################################
  
  cat("   Resultats del model: \n")
  
  print(summary(resultat))
  
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
  
  
  ##############################################################################
  # VALIDACIÓ DEL MODEL
  ##############################################################################
  
  dades <- d %>%
    filter(MES == mes) %>%
    select(MUNICIPI, ANY, MES, ALTITUD, Temp_max, Temp_min, AREA_URBANA, Pluja,
           x, y, Total_basal)
  
  coords <- cbind(dades$x, dades$y)
  
  ## DADES EN TRAIN-TEST
  
  smp_size <- floor(0.75 * nrow(dades))
  set.seed(1714)
  train_ind <- sample(seq_len(nrow(dades)), size = smp_size, replace = FALSE)
  
  train <- dades[train_ind, ]
  test <- dades[-train_ind, ]
  test$Total_basal <- NA
  
  train_coords <- coords[train_ind,]
  test_coords <- coords[-train_ind,]
  
  ## MATRIUS I STACKS
  
  group_train <- train$ANY - min(d$ANY) + 1
  group_test  <- test$ANY  - min(d$ANY) + 1
  
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
    data = list(y = train$Total_basal),
    A = list(1, Ae),
    effects = list(
      data.frame(b0 = rep(1, nrow(train)), 
                 ALTITUD = train$ALTITUD,
                 MUNICIPI = as.factor(train$MUNICIPI),
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
  
  p.res <- inla(
    formula,
    data = inla.stack.data(stk.f),
    control.predictor = list(compute = TRUE, A = inla.stack.A(stk.f)),
    control.compute   = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  
  ## PREDICCIONS 
  
  index.pred <- inla.stack.index(stk.f, "pred")$data
  pred <- p.res$summary.linear.predictor[index.pred, "mean"]
  obs <- dades[-train_ind, ]$Total_basal
  
  corr <- cor(obs, pred, use = "complete.obs")
  rmse <- sqrt(mean((pred - obs)^2, na.rm = TRUE))
  mae  <- mean(abs(pred - obs), na.rm = TRUE)
  
  df_plot <- data.frame(
    obs  = obs,
    pred = pred
  )
  
  metrics_text <- sprintf("Correlació = %.3f\nRMSE = %.3f\nMAE = %.3f",
                          corr, rmse, mae)
  print(
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
        title = paste("Valors observats i predits del mes -", mes),
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5)
      )
  )
  
  
  ##############################################################################
  # PREDICCIONS I MAPA
  ##############################################################################
  
  index <- inla.stack.index(stack = stk.full, tag = "pred")$data
  
  dp_plot <- dp_mes
  names(dp_plot)[1:3] <- c("x", "y", "time")
  
  dp_plot$p_mitj <- resultat$summary.fitted.values[index, "mean"]
  dp_plot$p_inf   <- resultat$summary.fitted.values[index, "0.025quant"]
  dp_plot$p_sup   <- resultat$summary.fitted.values[index, "0.975quant"]
  
  # MAPA
  dpm <- reshape2::melt(
    dp_plot,
    id.vars = c("x", "y", "time"),
    measure.vars = c("p_mitj", "p_inf", "p_sup")
  )
  
  dpm <- dpm %>%
    mutate(year = case_when(
      time == 1 ~ 2020,
      time == 2 ~ 2021,
      time == 3 ~ 2022,
      time == 4 ~ 2023,
      time == 5 ~ 2024
    ))
  
  print(max(dpm$value))
  print(min(dpm$value))
  
  print(
    ggplot(m) +
      geom_sf() +
      coord_sf(datum = NA) +
      geom_tile(data = dpm, aes(x = x, y = y, fill = value)) +
      labs(
        title = paste("Contaminació del mes -", mes),
        x = "", y = ""
      ) +
      facet_grid(year ~ variable) +
      scale_fill_gradientn(
        name = "Contaminació",
        colours = c("green", "yellow", "red"),
        limits = c(-0.2, 12),
        oob = scales::squish
      ) +
      theme_bw()
  )
}
