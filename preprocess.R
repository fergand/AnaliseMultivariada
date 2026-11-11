rm(list = ls())

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