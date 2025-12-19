## Aniol Vilamala Vidal

################################################################################
# ----- MATRIU DE VARIACIÓ DELS CONTAMINANTS -----
################################################################################

dades_finals <- readRDS("dades_finals.rds")

contaminants <- c("NO2", "PM10", "PM2.5", "O3", "CO", "SO2")
expressions <- c("NO₂", "PM₁₀", "PM₂.₅", "O₃", "CO", "SO₂")

X <- dades_finals %>%
  dplyr::select(all_of(contaminants))

matriu_var <- compositions::variation(compositions::acomp(X))
matriu_var <- round(matriu_var, 2)
diag(matriu_var) <- NA

df_var <- as.data.frame(matriu_var)
rownames(df_var) <- expressions
colnames(df_var) <- expressions

df_var_fmt <- df_var %>%
  tibble::rownames_to_column("Contaminants") %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "",
                                       cell_spec(formatC(.x, format = "f", digits = 2),
                                                 format = "html"))))

kable(df_var_fmt, format = "html", escape = FALSE, align = "c") %>%
  kable_styling(full_width = FALSE,
                bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "lightblue") %>%
  column_spec(1, bold = TRUE, color = "black", background = "lightblue") %>%
  column_spec(2:ncol(df_var_fmt), width = "7em") 


################################################################################
# ----- CODA BIPLOTS -----
################################################################################

pca <- X %>%
  clr %>% 
  princomp

fviz_pca_biplot(pca,
                geom.ind = "point",                
                palette = "jco",                  
                addEllipses = FALSE,                
                col.var = "red", 
                col.ind = "lightblue",
                repel = TRUE) +
  ggplot2::ggtitle("CoDA Biplot: casos i contaminants")    

grup_combinat <- interaction(dades_finals$TIPUS_ESTACIO, dades_finals$AREA_URBANA)

levels(grup_combinat) <- c(
  "Background-Suburban",
  "Traffic-Suburban",
  "Industrial-Suburban",
  "Background-Rural",
  "Traffic-Rural",
  "Industrial-Rural",
  "Background-Urban",
  "Traffic-Urban",
  "Industrial-Urban"
)


fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = grup_combinat,     
             palette = "lancet",
             addEllipses = FALSE,
             legend.title = "Tipus d'estació - \nÀrea urbana",
             pointshape = 16) +
  ggplot2::ggtitle("CoDA PCA: casos segons tipus d’estació i àrea urbana")        

################################################################################
# ----- TERNARY DIAGRAMS -----
################################################################################

X_2 <- X %>% 
  dplyr::select(PM10, O3, NO2)


plot(acomp(X_2), pca=TRUE, col.pca="red", axes=TRUE, col = as.numeric(dades$AREA_URBANA))


wind.cf1 = X %>% acomp %>% ClusterFinder1

sort(wind.cf1$typeTbl, decreasing = T)

################################################################################
# ----- NORMALITAT DE LES DADES -----
################################################################################

X %>% acomp %>% qqnorm

################################################################################
# ----- DESCRIPTIVA DE LA TRASNFORMACIÓ ILR -----
################################################################################

## ----- DESCRPIPTIVA BÀSICA ----- 

X_ilr <- compositions::ilr(as.matrix(X))

sum_mat <- summary(X_ilr)
summary_ilr <- as.data.frame(t(sum_mat))  
summary_ilr$Variance <- apply(X_ilr, 2, var, na.rm = TRUE)

summary_ilr$Coordinate <- rownames(summary_ilr)

summary_ilr <- summary_ilr %>%
  dplyr::select(Coordinate, Min., `1st Qu.`, Median, Mean, `3rd Qu.`, Max., Variance) %>%
  dplyr::rename(
    Coordenada = Coordinate,
    Mín. = Min.,
    Q1 = `1st Qu.`,
    Mediana = Median,
    Mitjana = Mean,
    Q3 = `3rd Qu.`,
    Màx. = Max.,
    Variància = Variance
  )

summary_ilr_fmt <- summary_ilr %>%
  mutate(across(where(is.numeric),
                ~ cell_spec(formatC(.x, format = "f", digits = 2),
                            format = "html"))) %>%
  dplyr::select(-Coordenada)

kable(summary_ilr_fmt, format = "html", escape = FALSE, align = "c") %>%
  kable_styling(full_width = FALSE,
                bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  row_spec(0, bold = TRUE, color = "black", background = "lightblue") %>%
  column_spec(1:ncol(summary_ilr_fmt), width = "4em") %>%
  column_spec(1, bold = TRUE, color = "black", background = "lightblue")

## ----- NORMALITAT DE LA TRANSFORMACIÓ -----

pairs(acomp(X))

pairs(
  acomp(X),
  panel = vp.lrdensityplot,      
  upper.panel = vp.lrdensityplot,
  lower.panel = vp.lrboxplot,    
  diag.panel = NULL,             
  col = 4,                      
  alpha = 0.05,                  
  main = "Pairs plot composicional dels contaminants atmosfèrics"
)

pairs(
  acomp(X),
  panel = vp.lrdensityplot,
  upper.panel = vp.lrboxplot,
  labels = contaminants,
  cex.labels = 1.2,
  font.labels = 2,
  gap = 0.5,
  row1attop = TRUE,
  main = "Relacions log-ratio entre contaminants (CoDA)"
)



