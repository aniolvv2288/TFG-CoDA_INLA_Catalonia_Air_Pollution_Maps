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
# ----- T-SPACE -----
################################################################################

## BALANCES ILR

d <- dades_finals

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

## TOTAL T-SPACE

X <- dades_finals %>%
  dplyr::select(all_of(contaminants))

D <- ncol(X)
Total_basal <- (1 / sqrt(D)) * rowSums(log(X))

d$Total_basal <- Total_basal


################################################################################
# ----- RESUM T-SPACE -----
################################################################################

ilr_df <- as.data.frame(ilr_bal)

d$id <- seq_len(nrow(d))
ilr_df$id <- seq_len(nrow(ilr_df))

d <- d %>% 
  left_join(ilr_df, by = "id") %>%
  dplyr::select(-all_of(contaminants), -id) %>%
  rename(Total = Total_basal)

taula_stats <- d %>%
  select(B1, B2, B3, B4, B5, Total) %>%
  summarise(
    across(
      everything(),
      list(
        Min = ~min(.x, na.rm = TRUE),
        Q1 = ~quantile(.x, 0.25, na.rm = TRUE),
        Median = ~median(.x, na.rm = TRUE),
        Mean = ~mean(.x, na.rm = TRUE),
        Q3 = ~quantile(.x, 0.75, na.rm = TRUE),
        Max = ~max(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = c("Variable", "Estadistic"),
    names_sep = "_",
    values_to = "Valor"
  ) %>%
  tidyr::pivot_wider(
    names_from = Estadistic,
    values_from = Valor
  )

taula_stats


################################################################################
# ----- EVOLUCIÓ COORDENADES DEL T-SPACE SEGONS EL MES I ANY -----
################################################################################

d_plot <- d %>%
  select(ANY, MES, Total, B1, B2, B3, B4, B5) %>%
  mutate(
    MES = factor(MES, levels = 1:12,
                 labels = c("Gen","Feb","Mar","Abr","Mai","Jun",
                            "Jul","Ago","Set","Oct","Nov","Des"))
  ) %>%
  rename(
    "Balança 1" = B1,
    "Balança 2" = B2,
    "Balança 3" = B3,
    "Balança 4" = B4,
    "Balança 5" = B5
  ) %>%
  pivot_longer(
    cols = c(Total, `Balança 1`, `Balança 2`, `Balança 3`, `Balança 4`, `Balança 5`),
    names_to = "Variable",
    values_to = "Valor"
  ) %>%
  group_by(Variable, ANY, MES) %>%       
  summarise(Valor = mean(Valor), .groups = "drop")

ggplot(d_plot, aes(x = MES, y = Valor, color = factor(ANY), group = ANY)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  scale_color_viridis_d(option = "C", end = 0.85) +
  
  labs(
    title = "Evolució mensual de les coordenades del T-space segons l'any",
    x = "Mes",
    y = "Valor mitjà",
    color = "Any"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",     
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
