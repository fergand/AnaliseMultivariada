rm(list = ls())

pak::pak("pgmm")
data("coffee", package = "pgmm")


str(coffee)
?coffee


pak::pak("MVar.pt")
data(DataMix, package = "MVar.pt")
data(Data_Cafes, package = "MVar.pt")

str(DataMix)
str(Data_Cafes)


pak::pak("metan")
data(data_ge2, package = "metan")




# Pacotes iniciais --------------------------------------------------------

pacman::p_load(tidyverse, readxl, janitor, corrplot, MVN, candisc, factoextra, MASS)

# Dados Café: Profa. Dra Sonia Maria Costa Celestino ----------------------

df_solo_bruto <- read_tsv("dados/Dados do solo.tsv",locale = locale(decimal_mark = ","))
df_met_bruto <- read_csv("dados/Dados metereológicos.csv",locale = locale(decimal_mark = ","))
df_cafe_bruto <- read_xlsx("dados/Dados qualidade café.xlsx")

str(df_solo_bruto);str(df_met_bruto);str(df_cafe_bruto)
summary(df_solo_bruto)
summary(df_met_bruto)
summary(df_cafe_bruto)


# Limpeza inicial dos dados -----------------------------------------------

df_cafe <- df_cafe_bruto %>% 
  clean_names() %>% 
  rename(
    amostra = codigo_da_amostra,
    umidade = teor_de_agua_g_100g,
    brix = solidos_soluveis_totais_graus_brix,
    cloreto = cloreto_percent,
    cafeina = cafeina_percent_b_s,
    ac_clorogenico = ac_clorogenico_percent_b_s,
    pen16 = peneira_16_percent,
    pen17 = peneira_17_percent,
    pen18 = peneira_18_percent,
    pen19 = peneira_19_percent
  ) %>% 
  mutate(
    amostra = as.factor(amostra),
    pen18_19 = pen18+pen19,
    class_grao = factor(ifelse(pen17 + pen18_19 >= median(pen17 + pen18_19), "Graudo", "Medio_Miudo"))
  )


summary(df_cafe)  

boxplot(df_cafe)

library(car)
df_cafe %>% 
  dplyr::select(-amostra,-class_grao) %>% 
  cor() %>% 
  corrplot(method = "ellipse", type = "upper", 
           addCoef.col = "black", tl.col = "black", number.cex = 0.8,
           title = "Matriz de Correlação das Respostas", mar = c(0,0,2,0))


vars_y <- c("umidade", "brix", "cloreto", "cafeina", "ac_clorogenico", 
            "pen16", "pen17", "pen18_19")

corrplot(cor(as.matrix(df_cafe[, vars_y])), method = "ellipse", type = "upper", 
         addCoef.col = "black", tl.col = "black", number.cex = 0.75)


Y_mat <- as.matrix(df_cafe[, vars_y])
acp_cafe <- prcomp(Y_mat, center = TRUE, scale. = TRUE)
fviz_pca_biplot(acp_cafe, 
                habillage = df_cafe$class_grao,
                addEllipses = TRUE, palette = "Dark2",
                repel = TRUE, title = "Biplot ACP - Qualidade e Granulometria do Café")

bloco_fisico  <- scale(df_cafe[, c("pen16", "pen17", "pen18_19")])
bloco_quimico <- scale(df_cafe[, c("umidade", "brix", "cloreto", "cafeina", "ac_clorogenico")])
acc_res <- candisc::cancor(bloco_fisico, bloco_quimico)
summary(acc_res)

# Extrair os escores canônicos do 1º par (Xcan1 e Ycan1)
df_escores_acc <- data.frame(
  Amostra     = df_cafe$amostra,
  Classe_Grao = df_cafe$class_grao, # ou df_cafe$class_grao
  Xcan1       = acc_res$scores$X[, 1],
  Ycan1       = acc_res$scores$Y[, 1]
)

# Gráfico de dispersão dos escores canônicos com ggplot2
ggplot(df_escores_acc, aes(x = Xcan1, y = Ycan1, color = Classe_Grao)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "gray30", linetype = "dashed") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "1º Par Canônico: Granulometria (Xcan1) vs. Química (Ycan1)",
    subtitle = "Correlação Canônica r = 0,527 (R² = 27,8%, p = 0,064)",
    x = "Variável Canônica Física (Xcan1 - Peneiras)",
    y = "Variável Canônica Química (Ycan1 - Compostos)",
    color = "Padrão de Grão"
  ) +
  ggthemes::theme_clean()

# Extrair as correlações (cargas canônicas) das variáveis originais com o 1º par
cargas_x <- data.frame(
  Variavel = rownames(acc_res$structure$X.xscores),
  Bloco    = "Bloco X (Físico - Peneiras)",
  Carga    = acc_res$structure$X.xscores[, 1]
)

cargas_y <- data.frame(
  Variavel = rownames(acc_res$structure$Y.yscores),
  Bloco    = "Bloco Y (Químico)",
  Carga    = acc_res$structure$Y.yscores[, 1]
)

df_cargas <- bind_rows(cargas_x, cargas_y) |>
  mutate(Variavel = fct_reorder(Variavel, Carga))

ggplot(df_cargas, aes(x = Carga, y = Variavel, fill = Carga > 0)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_vline(xintercept = c(-0.3, 0, 0.3)) +
  facet_wrap(~ Bloco, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
  xlim(-1, 1) +
  labs(
    title = "Cargas Canônicas (Correlações Estruturais) no 1º Par",
    subtitle = "Correlação direta de cada variável original com sua respectiva variável canônica",
    x = "Correlação com a 1ª Variável Canônica (r)",
    y = NULL
  )

# Plota Xcan1 vs Ycan1 junto com os vetores de peso das variáveis
plot(acc_res, which = 1, smooth = TRUE)
