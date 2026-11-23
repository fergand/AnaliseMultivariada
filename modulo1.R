# carregando os dados no ambiente de trabalho
load("dados/cafe.RData")


pacman::p_load(tidyverse, janitor, MVN, heplots, Hotelling, magrittr)
set_theme(ggthemes::theme_clean())

# Módulo 1: efetuar MANOVA e/ou T2 de Hotelling para comparar cafés de 
# Peneira Alta vs. Peneira Baixa quanto ao vetor químico

df_cafe %>% glimpse()

df_cafe %<>%
  mutate(media_pen_alta = mean(pen17 + pen18_19), .by = amostra) %>%
  mutate(
    grupo_peneira = factor(
      ifelse(media_pen_alta >= median(media_pen_alta), "Peneira Alta", "Peneira Baixa"),
      levels = c("Peneira Alta", "Peneira Baixa")
    )
  )

vars_quim <- c("umidade", "brix", "cloreto", "cafeina", "ac_clorogenico")

print(table(df_cafe$grupo_peneira))

df_cafe %>% 
  group_by(grupo_peneira) %>% 
  dplyr::summarise(
    n = dplyr::n(),
    across(all_of(vars_quim), list(media = mean, dp = sd), .names = "{.col}_{.fn}")
  )


# Verificação dos pressupostos multivariados ------------------------------

Y_alta <- as.matrix(df_cafe %>% filter(grupo_peneira == "Peneira Alta") %>% select(all_of(vars_quim)))
Y_baixa <- as.matrix(df_cafe %>% filter(grupo_peneira == "Peneira Baixa") %>% select(all_of(vars_quim)))

mvnTest::HZ.test(Y_alta, qqplot = T)
MVN::mvn(Y_alta, mvn_test = "mardia")$multivariate_normality
mvnTest::HZ.test(Y_baixa, qqplot = T)
MVN::mvn(Y_baixa, mvn_test = "mardia")$multivariate_normality

# Homogeneidade (pressuposto do T2 e MANOVA)
# H0: Sigma_Alta = Sigma_Baixa
heplots::boxM(df_cafe[, vars_quim], df_cafe$grupo_peneira)
# valor-p > 0.05 -> não rejeita-se H0. portanto, segue homogeneidades das matrizes
# de variância-covariância

n1 <- nrow(Y_alta); n2 <- nrow(Y_baixa); p <- ncol(Y_alta)
m1 <- colMeans(Y_alta); m2 <- colMeans(Y_baixa)
S1 <- cov(Y_alta); S2 <- cov(Y_baixa)

Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)

diff_medias <- m1 - m2
T2_stat <- as.numeric((n1 * n2 / (n1 + n2)) * t(diff_medias) %*% solve(Sp) %*% diff_medias)
F_stat <- ((n1 + n2 - p - 1) / (p * (n1 + n2 - 2))) * T2_stat
p_valor_T2 <- pf(F_stat, df1 = p, df2 = n1 + n2 - p - 1, lower.tail = FALSE)
wilks <- 1 / (1 + T2_stat / (n1 + n2 - 2))

cat(sprintf("T² = %.4f | F(%d, %d) = %.4f | p-valor = %.5f | Wilks eq. = %.5f\n",
            T2_stat, p, n1 + n2 - p - 1, F_stat, p_valor_T2, wilks))


fit_manova <- manova(as.matrix(df_cafe[, vars_quim]) ~ grupo_peneira, data = df_cafe)
print(summary(fit_manova, test = "Wilks"))
print(summary(fit_manova, test = "Pillai"))


alpha <- 0.05
se_diff <- sqrt(diag(Sp) * (1/n1 + 1/n2))

# Constante crítica de Hotelling (Roy-Bose) e de Bonferroni (m = 5 comparações)
c_hotelling  <- sqrt((p * (n1 + n2 - 2) / (n1 + n2 - p - 1)) * qf(1 - alpha, p, n1 + n2 - p - 1))
c_bonferroni <- qt(1 - alpha / (2 * p), df = n1 + n2 - 2)

# Vetor discriminante padronizado (mostra qual variável mais separa os grupos)
a_bruto <- as.numeric(solve(Sp) %*% diff_medias)
a_padronizado <- a_bruto * sqrt(diag(Sp))

df_posthoc <- tibble(
  Variavel     = vars_quim,
  Media_Alta   = m1,
  Media_Baixa  = m2,
  Diferenca    = diff_medias,
  EP           = se_diff,
  IC_Bonf_Inf  = diff_medias - c_bonferroni * se_diff,
  IC_Bonf_Sup  = diff_medias + c_bonferroni * se_diff,
  IC_T2_Inf    = diff_medias - c_hotelling * se_diff,
  IC_T2_Sup    = diff_medias + c_hotelling * se_diff
)

cat("\n ICs SIMULTÂNEOS 95% E PESO DISCRIMINANTE)\n")
print(df_posthoc, width = Inf)


df_cafe %>% 
  mutate(across(all_of(vars_quim), ~ as.numeric(scale(.)))) %>% 
  pivot_longer(all_of(vars_quim), names_to = "Variavel", values_to = "Z_score") %>% 
  group_by(grupo_peneira, Variavel) %>% 
  dplyr::summarise(
    Media_Z = mean(Z_score),
    EP_Z    = sd(Z_score) / sqrt(n()),
    .groups = "drop"
  ) %>% 
  mutate(Variavel = factor(Variavel, levels = vars_quim,
                           labels = c("Umidade", "°Brix", "Cloreto", "Cafeína", "Ác. Clorogênico"))) %>% 
  ggplot(aes(x = Variavel, y = Media_Z, color = grupo_peneira, group = grupo_peneira)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Media_Z - 1.96 * EP_Z, ymax = Media_Z + 1.96 * EP_Z), width = 0.12) +
  scale_color_manual(values = c("Peneira Alta" = "#1b9e77", "Peneira Baixa" = "#d95f02")) +
  labs(
    x = NULL, y = "Média Padronizada (Z-score)", color = "Classificação"
  )


df_posthoc %>% 
  mutate(
    Variavel = factor(Variavel, levels = rev(vars_quim),
                      labels = rev(c("Umidade", "°Brix", "Cloreto", "Cafeína", "Ác. Clorogênico"))),
    Diff_d     = Diferenca /sqrt(diag(Sp)),
    Bonf_Inf_d = IC_Bonf_Inf /sqrt(diag(Sp)),
    Bonf_Sup_d = IC_Bonf_Sup /sqrt(diag(Sp)),
    T2_Inf_d   = IC_T2_Inf / sqrt(diag(Sp)),
    T2_Sup_d   = IC_T2_Sup / sqrt(diag(Sp))
  ) %>% 
  ggplot(aes(y = Variavel, x = Diff_d)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick") +
  geom_errorbar(aes(xmin = T2_Inf_d, xmax = T2_Sup_d), width = 0.25, color = "gray65", linewidth = 1) +
  geom_errorbar(aes(xmin = Bonf_Inf_d, xmax = Bonf_Sup_d), width = 0.12, color = "#1b9e77", linewidth = 1.5) +
  geom_point(size = 3, color = "black") +
  labs(
    x = "Tamanho de Efeito Padronizado (d de Cohen)", y = NULL
  )

