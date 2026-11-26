library(shiny)
library(bslib)
library(plotly)
library(tidyverse)
library(kableExtra)
library(ggthemes); set_theme(theme_clean())
library(rgl)
library(MVA)

# pacotes adicionais para o modulo 2

library(candisc)
library(FactoMineR)
library(factoextra)
library(MASS)
library(caret)
library(pROC)
library(cluster)

## para uma primeira versão, vamos começar com o df do exercício, para depois 
## generalizar

load("dados/cafe.RData") |> as_tibble()
dados <- df_cafe %<>%
  mutate(media_pen_alta = mean(pen17 + pen18 + pen19), .by = amostra) %>%
  mutate(
    grupo_peneira = factor(
      ifelse(media_pen_alta >= median(media_pen_alta), "Peneira Alta", "Peneira Baixa"),
      levels = c("Peneira Alta", "Peneira Baixa")
    )
  )
vars_quim <- c("umidade", "brix", "cloreto", "cafeina", "ac_clorogenico")
vars_fis <- c("pen16", "pen17", "pen18", "pen19")

rotulos_quim <- c(
  "Umidade" = "umidade",
  "°Brix" = "brix",
  "Cloreto" = "cloreto",
  "Cafeína" = "cafeina",
  "Ác. Clorogênico" = "ac_clorogenico"
)

## ui
ui <- page_fluid(
  withMathJax(),
  layout_columns(
    col_widths = c(3, 9),
    
    # ==========================================================================
    # CARD DA ESQUERDA (muda sozinho conforme a aba escolhida)
    # ==========================================================================
    card(
      # Controles do Módulo 1 
      conditionalPanel(
        condition = "input.tabset != 'acc' && input.tabset != 'acp' && input.tabset != 'disc' && input.tabset != 'cluster'",
        card(
          helpText("VETOR DE DIFERENÇA NULO (\\(\\boldsymbol{\\delta}_0 = \\boldsymbol{\\mu}_1 - \\boldsymbol{\\mu}_2\\))"),
          numericInput("mu01", "Δ Umidade",         value = 0, step = 0.1),
          numericInput("mu02", "Δ °Brix",           value = 0, step = 0.5),
          numericInput("mu03", "Δ Cloreto",         value = 0, step = 0.005),
          numericInput("mu04", "Δ Cafeína",         value = 0, step = 0.05),
          numericInput("mu05", "Δ Ác. Clorogênico", value = 0, step = 0.1)
        ),
        br(),
        card(
          "PARÂMETROS",
          sliderInput(
            "alpha",
            "Nível de significância \\(\\alpha\\):",
            value = 0.05, min = 0.01, max = 0.10, step = 0.01
          ),
          fluidRow(
            column(6, actionButton("b", "Bonferroni", width = "100%")),
            column(6, actionButton("h", "Hotelling",  width = "100%"))
          )
        )
      ),
      
      # Controles da ACC
      conditionalPanel(
        condition = "input.tabset == 'acc'",
        strong("PARÂMETROS DA ACC"),
        hr(),
        selectInput(
          "acc_par", "Par canônico para visualização:",
          choices = c("1º Par (u1, v1)" = 1, "2º Par (u2, v2)" = 2),
          selected = 1
        ),
        sliderInput(
          "acc_alpha", "Nível de significância (α):",
          min = 0.01, max = 0.10, value = 0.05, step = 0.01
        ),
        radioButtons(
          "acc_tipo_carga", "Tipo de carga nos gráficos:",
          choices = c("Diretas (X com u / Y com v)" = "direta",
                      "Cruzadas (X com v / Y com u)" = "cruzada"),
          selected = "direta"
        )
      ),
      
      # Controles da ACP
      conditionalPanel(
        condition = "input.tabset == 'acp'",
        strong("CONFIGURAÇÃO DA ACP"),
        hr(),
        radioButtons(
          "acp_modo", "Conjunto de variáveis:",
          choices = c(
            "Apenas Variáveis Químicas" = "quim",
            "Todas as Variáveis (Química + Peneiras)" = "full"
          ),
          selected = "quim"
        ),
        selectInput("acp_eixo1", "Eixo horizontal:", choices = 1:4, selected = 1),
        selectInput("acp_eixo2", "Eixo vertical:",   choices = 1:4, selected = 2),
        checkboxInput("acp_elipses", "Mostrar elipses de confiança (95%)", value = TRUE)
      ),
      
      # Controles da Análise Discriminante
      conditionalPanel(
        condition = "input.tabset == 'disc'",
        strong("MODELO DISCRIMINANTE"),
        hr(),
        radioButtons(
          "disc_metodo", "Função discriminante:",
          choices = c("Linear (LDA)" = "lda", "Quadrática (QDA)" = "qda"),
          selected = "lda"
        ),
        radioButtons(
          "disc_eval", "Método de avaliação:",
          choices = c("Validação Cruzada (LOO-CV)" = "cv",
                      "Ressubstituição (aparente)" = "aparente"),
          selected = "cv"
        ),
        checkboxGroupInput(
          "disc_vars", "Variáveis químicas no modelo:",
          choices = rotulos_quim, selected = vars_quim
        ),
        hr(),
        strong("CLASSIFICAR NOVA AMOSTRA"),
        numericInput("pred_umidade", "Umidade (g/100g)", value = 11.20, step = 0.1),
        numericInput("pred_brix",    "°Brix",            value = 32.00, step = 0.5),
        numericInput("pred_cloreto", "Cloreto (%)",      value = 0.022, step = 0.002),
        numericInput("pred_cafeina", "Cafeína (% b.s.)", value = 0.92,  step = 0.05),
        numericInput("pred_cloro",   "Ác. Clorogênico",  value = 4.20,  step = 0.1),
        uiOutput("disc_pred_nova")
      ),
      
      # Controles do Cluster
      conditionalPanel(
        condition = "input.tabset == 'cluster'",
        strong("PARÂMETROS DO CLUSTER"),
        helpText("Agrupamento realizado sobre as 5 variáveis químicas padronizadas."),
        hr(),
        sliderInput("cl_k", "Número de clusters (k):", min = 2, max = 5, value = 2, step = 1),
        selectInput(
          "cl_metodo_hc", "Método de ligação (Hierárquico):",
          choices = c(
            "Ward.D2" = "ward.D2",
            "Média (UPGMA)" = "average",
            "Completa" = "complete",
            "Simples" = "single"
          ),
          selected = "ward.D2"
        ),
        hr(),
        strong("Correlação Cofenética:"),
        verbatimTextOutput("cl_cofen")
      )
    ),
    
    # ==========================================================================
    # CARD DA DIREITA 
    # ==========================================================================
    navset_card_tab(
      id = "tabset",
      
      # MENU DO MÓDULO 1
      nav_menu(
        "Módulo 1",
        nav_panel(
          "Dados",
          br(),
          "Pré-visualização dos dados carregados (variáveis químicas e classe de peneira)",
          br(), br(),
          tableOutput("table1")
        ),
        nav_panel(
          "Exploratória",
          br(),
          fluidRow(
            column(
              6,
              selectInput(
                "grupo_diag",
                "Selecione o grupo para diagnóstico de distribuição e normalidade:",
                choices = c("Peneira Alta", "Peneira Baixa"),
                selected = "Peneira Alta"
              )
            )
          ),
          fluidRow(
            column(4, plotOutput("hist1", height = "240px")),
            column(4, plotOutput("hist2", height = "240px")),
            column(4, plotOutput("hist3", height = "240px"))
          ),
          fluidRow(
            column(4, plotOutput("hist4", height = "240px")),
            column(4, plotOutput("hist5", height = "240px"))
          ),
          hr(),
          strong("Boxplots Comparativos das Variáveis Químicas (Peneira Alta vs. Peneira Baixa)"),
          br(), br(),
          fluidRow(
            column(4, checkboxInput("box_pontos", "Mostrar pontos individuais", value = TRUE)),
            column(4, checkboxInput("box_medias", "Destacar médias (losango)", value = TRUE))
          ),
          plotOutput("boxplots_exp", height = "420px"),
          hr(),
          fluidRow(
            column(
              6,
              verbatimTextOutput("hz_output"),
              verbatimTextOutput("boxm_output")
            ),
            column(6, plotOutput("qqplot"))
          )
        ),
        nav_panel(
          "Teste T2",
          br(),
          withMathJax(helpText(
            "Teste \\(T^2\\) de Hotelling (2 amostras independentes): ",
            "\\(H_0\\colon \\boldsymbol{\\mu}_1 - \\boldsymbol{\\mu}_2 = \\boldsymbol{\\delta}_0\\) vs. ",
            "\\(H_1\\colon \\boldsymbol{\\mu}_1 - \\boldsymbol{\\mu}_2 \\neq \\boldsymbol{\\delta}_0\\)"
          )),
          br(),
          fluidRow(
            column(3, strong("Estatística \\(T^2\\)"), textOutput("stats_t2")),
            column(3, strong("\\(F\\) calculado"), textOutput("Fcalc")),
            column(3, strong("\\(F\\) crítico"), textOutput("Ftab")),
            column(3, strong("Valor-p"), textOutput("valor_p"))
          ),
          br(),
          fluidRow(
            column(4, strong("Lambda de Wilks (\\(\\Lambda\\))"), textOutput("wilks_val"))
          ),
          hr(),
          fluidRow(
            column(12, uiOutput("t2_conclusao"))
          )
        ),
        nav_panel(
          "ICS e Elipsóide",
          br(),
          uiOutput("titulo_ics"),
          tableOutput("tabics"),
          hr(),
          strong("Projeção 3D do Elipsoide de Confiança (\\(\\boldsymbol{\\mu}_1 - \\boldsymbol{\\mu}_2\\))"),
          fluidRow(
            column(4, selectInput("eixo_x", "Eixo X:", choices = rotulos_quim, selected = "umidade")),
            column(4, selectInput("eixo_y", "Eixo Y:", choices = rotulos_quim, selected = "brix")),
            column(4, selectInput("eixo_z", "Eixo Z:", choices = rotulos_quim, selected = "cafeina"))
          ),
          plotlyOutput("elipse", height = "550px")
        )
      ),
      
      # MENU DO MÓDULO 2
      nav_menu(
        "Módulo 2",
        nav_panel(
          "Correlação Canônica (ACC)",
          value = "acc",
          br(),
          strong("Correlações Canônicas e Teste de Bartlett (Qui-Quadrado)"),
          tableOutput("acc_tabela"),
          uiOutput("acc_conclusao"),
          hr(),
          layout_columns(
            col_widths = c(6, 6),
            plotOutput("acc_cargas_plot", height = "340px"),
            plotlyOutput("acc_dispersao", height = "340px")
          )
        ),
        nav_panel(
          "Componentes Principais (ACP)",
          value = "acp",
          br(),
          layout_columns(
            col_widths = c(7, 5),
            plotOutput("acp_plot_principal", height = "420px"),
            plotOutput("acp_scree", height = "420px")
          ),
          hr(),
          layout_columns(
            col_widths = c(6, 6),
            div(
              strong("Autovalores e Variância Explicada"),
              tableOutput("acp_autovalores")
            ),
            div(
              strong("Cargas (Correlações com os Eixos) e cos² no Plano"),
              tableOutput("acp_cargas")
            )
          )
        ),
        nav_panel(
          "Análise Discriminante",
          value = "disc",
          br(),
          layout_columns(
            col_widths = c(6, 6),
            div(
              strong("Matriz de Confusão e Métricas"),
              verbatimTextOutput("disc_cm")
            ),
            plotOutput("disc_roc", height = "340px")
          ),
          hr(),
          layout_columns(
            col_widths = c(6, 6),
            div(
              strong("Importância das Variáveis (Coef. Padronizados e Cargas LD1)"),
              tableOutput("disc_importancia")
            ),
            plotOutput("disc_densidade", height = "300px")
          )
        ),
        nav_panel(
          "Agrupamento (Cluster)",
          value = "cluster",
          br(),
          layout_columns(
            col_widths = c(6, 6),
            plotOutput("cl_kmeans_plot", height = "340px"),
            plotOutput("cl_dendro_plot", height = "340px")
          ),
          hr(),
          layout_columns(
            col_widths = c(5, 7),
            div(
              strong("K-Means vs. Tamanho de Grão (Teste de Fisher)"),
              verbatimTextOutput("cl_fisher")
            ),
            div(
              strong("Perfil Médio (Químico e Físico) por Cluster K-Means"),
              tableOutput("cl_perfis")
            )
          )
        )
      )
    )
  )
)

## functions para o server
qqmahala <- function(data) {
  xp <- scale(data, scale = FALSE)
  n <- nrow(data)
  p <- ncol(data)
  Sa <- cov(xp)
  d <- diag(xp %*% solve(Sa) %*% t(xp))
  r <- rank(d)
  chi2q <- qchisq((r - 0.5)/n, p)
  plot(
    d, chi2q, pch = 19, 
    main = "Chi-Square Q-Q Plot", 
    xlab = "Squared Mahalanobis Distance",
    ylab = "Chi-Square Quantile"
  )
  abline(0,1, lwd = 2, col = "firebrick2")
}

t2test <- function(delta0, data, vars, alpha) {
  Y1 <- as.matrix(data %>% dplyr::filter(grupo_peneira == "Peneira Alta")  %>% dplyr::select(all_of(vars)))
  Y2 <- as.matrix(data %>% dplyr::filter(grupo_peneira == "Peneira Baixa") %>% dplyr::select(all_of(vars)))
  
  n1 <- nrow(Y1); n2 <- nrow(Y2); p <- ncol(Y1)
  m1 <- colMeans(Y1); m2 <- colMeans(Y2)
  S1 <- cov(Y1);      S2 <- cov(Y2)
  
  Sp <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  diff <- (m1 - m2) - delta0
  
  T2 <- as.numeric((n1 * n2 / (n1 + n2)) * t(diff) %*% solve(Sp) %*% diff)
  Fcalc <- ((n1 + n2 - p - 1) / (p * (n1 + n2 - 2))) * T2
  Ftab <- qf(1 - alpha, p, n1 + n2 - p - 1)
  pvalue <- pf(Fcalc, p, n1 + n2 - p - 1, lower.tail = FALSE)
  wilks <- 1 / (1 + T2 / (n1 + n2 - 2))
  
  data.frame(
    t2 = T2, Fcalc = Fcalc, Ftab = Ftab, valorp = pvalue, wilks = wilks
  )
}

ICS <- function(data, vars, type = "b", alpha = 0.05) {
  Y1 <- as.matrix(data %>% dplyr::filter(grupo_peneira == "Peneira Alta")  %>% dplyr::select(all_of(vars)))
  Y2 <- as.matrix(data %>% dplyr::filter(grupo_peneira == "Peneira Baixa") %>% dplyr::select(all_of(vars)))
  
  n1 <- nrow(Y1); n2 <- nrow(Y2); p <- ncol(Y1)
  m1 <- colMeans(Y1); m2 <- colMeans(Y2)
  Sp <- ((n1 - 1) * cov(Y1) + (n2 - 1) * cov(Y2)) / (n1 + n2 - 2)
  
  if (!(type %in% c("h", "b")))
    stop("Forneça valor para type: 'b' para Bonferroni ou 'h' para Hotelling")
  
  diff_m <- m1 - m2
  
  resultado <- data.frame(
    Media_Alta = numeric(p),
    Media_Baixa = numeric(p),
    Diferenca = numeric(p),
    EP = numeric(p),
    LI = numeric(p),
    LS = numeric(p),
    a = character(p),
    row.names   = vars
  )
  
  for (i in 1:p) {
    a <- matrix(0, p, 1); a[i] <- 1
    var_diff <- as.numeric((1 / n1 + 1 / n2) * t(a) %*% Sp %*% a)
    ep_i <- sqrt(var_diff)
    
    if (type == "h") {
      crit <- sqrt((p * (n1 + n2 - 2) / (n1 + n2 - p - 1)) * qf(1 - alpha, p, n1 + n2 - p - 1))
    } else {
      crit <- qt(1 - alpha / (2 * p), df = n1 + n2 - 2)
    }
    
    erro <- crit * ep_i
    LI <- diff_m[i] - erro
    LS <- diff_m[i] + erro
    
    resultado[i, ] <- list(m1[i], m2[i], diff_m[i], ep_i, LI, LS, paste(a[1:p], collapse = " "))
  }
  resultado
}

elips <- function(delta0, data, vars_all, vars_3d, alpha) {
  Y1 <- as.matrix(data %>% dplyr::filter(grupo_peneira == "Peneira Alta")  %>% dplyr::select(all_of(vars_all)))
  Y2 <- as.matrix(data %>% dplyr::filter(grupo_peneira == "Peneira Baixa") %>% dplyr::select(all_of(vars_all)))
  
  n1 <- nrow(Y1); n2 <- nrow(Y2); p <- ncol(Y1)
  diff_m <- colMeans(Y1) - colMeans(Y2)
  Sp <- ((n1 - 1) * cov(Y1) + (n2 - 1) * cov(Y2)) / (n1 + n2 - 2)
  
  # Índices das 3 variáveis escolhidas para projeção
  idx <- match(vars_3d, vars_all)
  S_diff_3d <- (1 / n1 + 1 / n2) * Sp[idx, idx]
  centro_3d <- diff_m[idx]
  delta_3d <- delta0[idx]
  
  # Raio crítico da projeção exata do elipsoide 5D (T² de Hotelling)
  t_crit <- sqrt((p * (n1 + n2 - 2) / (n1 + n2 - p - 1)) * qf(1 - alpha, p, n1 + n2 - p - 1))
  
  elip <- rgl::ellipse3d(
    as.matrix(S_diff_3d),
    centre = centro_3d,
    t = t_crit
  )
  
  verts <- t(elip$vb[1:3, ])
  faces <- t(elip$ib) - 1
  
  nomes_eixos <- names(rotulos_quim)[match(vars_3d, rotulos_quim)]
  
  plot_ly() %>%
    add_trace(
      type = "mesh3d",
      x = verts[, 1], y = verts[, 2], z = verts[, 3],
      i = faces[, 1], j = faces[, 2], k = faces[, 3],
      color = I("firebrick"),
      opacity = 0.20,
      contour = list(show = TRUE, color = "black", width = 2),
      name = "Elipsoide (1 - α)",
      hoverinfo = "skip"
    ) %>%
    add_markers(
      x = centro_3d[1], y = centro_3d[2], z = centro_3d[3],
      marker = list(color = "blue", size = 5),
      name = "Diferença Amostral"
    ) %>%
    add_markers(
      x = delta_3d[1], y = delta_3d[2], z = delta_3d[3],
      marker = list(color = "green", size = 7),
      name = "Vetor Nulo (δ₀)"
    ) %>%
    layout(scene = list(
      xaxis = list(title = paste("Δ", nomes_eixos[1]), showspikes = FALSE),
      yaxis = list(title = paste("Δ", nomes_eixos[2]), showspikes = FALSE),
      zaxis = list(title = paste("Δ", nomes_eixos[3]), showspikes = FALSE)
    ))
}

## server

server <- function(input, output) {
  
  # Reativos principais
  df <- reactive(dados)
  
  selected_grupo <- reactive({
    df() %>%
      dplyr::filter(grupo_peneira == input$grupo_diag) %>%
      dplyr::select(all_of(vars_quim))
  })
  
  delta0 <- reactive(c(input$mu01, input$mu02, input$mu03, input$mu04, input$mu05))
  alpha <- reactive(input$alpha)
  
  hz_result <- reactive({
    mvnTest::HZ.test(selected_grupo(), qqplot = FALSE)
  })
  
  boxm_result <- reactive({
    heplots::boxM(df()[, vars_quim], df()$grupo_peneira)
  })
  
  t2_result <- reactive({
    t2test(delta0(), df(), vars_quim, alpha())
  })
  
  # Tipo de intervalo simultâneo (b = Bonferroni, h = Hotelling)
  v <- reactiveValues(value = "b")
  observeEvent(input$b, { v$value <- "b" })
  observeEvent(input$h, { v$value <- "h" })
  
  # Diagnósticos de pressupostos
  output$hz_output <- renderPrint({
    res <- hz_result()
    a <- alpha()
    
    cat("Teste de Henze-Zirkler — Grupo:", input$grupo_diag, "\n\n")
    cat("  HZ      :", round(res@HZ, 4), "\n")
    cat("  p-value :", format.pval(res@p.value, digits = 4), "\n\n")
    cat("  Result  :",
        if (res@p.value > a) {
          paste0("Não se rejeita Normalidade Multivariada (α = ", a, ")")
        } else {
          paste0("Rejeita-se Normalidade Multivariada (α = ", a, ")")
        })
  })
  
  output$boxm_output <- renderPrint({
    bm <- boxm_result()
    a <- alpha()
    
    cat("Teste M de Box (Homogeneidade de Covariâncias: Σ₁ = Σ₂)\n\n")
    cat("  Qui-Quadrado :", round(bm$statistic, 4), "| gl =", bm$parameter, "\n")
    cat("  p-value      :", round(bm$p.value, 4), "\n\n")
    cat("  Result       :",
        if (bm$p.value > a) {
          paste0("Não se rejeita homogeneidade de covariâncias (α = ", a, ")")
        } else {
          paste0("Rejeita-se homogeneidade de covariâncias (α = ", a, ")")
        })
  })
  
  # Saídas do Teste T² de Hotelling
  output$stats_t2 <- renderText(round(t2_result()$t2, 4))
  output$Fcalc <- renderText(round(t2_result()$Fcalc, 4))
  output$Ftab <- renderText(round(t2_result()$Ftab, 4))
  output$valor_p <- renderText(round(t2_result()$valorp, 4))
  output$wilks_val <- renderText(round(t2_result()$wilks, 4))
  
  output$t2_conclusao <- renderUI({
    p_value <- t2_result()$valorp
    delta0_txt <- paste(round(delta0(), 3), collapse = ", ")
    
    if (p_value < alpha()) {
      withMathJax(HTML(paste0(
        "Ao nível de significância de <b>", alpha() * 100, "%</b>, <b>rejeita-se \\(H_0\\)</b>: ",
        "há evidências estatísticas de que o vetor de diferenças médias entre Peneira Alta e Peneira Baixa ",
        "difere de \\(\\boldsymbol{\\delta}_0 = (", delta0_txt, ")'\\)."
      )))
    } else {
      withMathJax(HTML(paste0(
        "Ao nível de significância de <b>", alpha() * 100, "%</b>, <b>não se rejeita \\(H_0\\)</b>: ",
        "não há evidências estatísticas suficientes para afirmar que o vetor de diferenças médias entre ",
        "Peneira Alta e Peneira Baixa difira de \\(\\boldsymbol{\\delta}_0 = (", delta0_txt, ")'\\)."
      )))
    }
  })
  
  # Tabelas
  output$table1 <- renderTable({
    df() %>% dplyr::select(amostra, grupo_peneira, all_of(vars_quim))
  }, digits = 4)
  
  output$titulo_ics <- renderUI({
    metodo <- if (v$value == "b") "Bonferroni" else "Hotelling (Roy-Bose)"
    strong(paste0("Intervalos de Confiança Simultâneos (", (1 - alpha()) * 100, "%) — Método: ", metodo))
  })
  
  output$tabics <- renderTable({
    ICS(data = df(), vars = vars_quim, type = v$value, alpha = alpha())
  }, rownames = TRUE, align = "c", digits = 4)
  
  # Gráficos da aba Exploratória
  output$hist1 <- renderPlot({
    hist(selected_grupo()$umidade, main = "Umidade", xlab = "g/100g", col = "gray85", border = "white")
  })
  output$hist2 <- renderPlot({
    hist(selected_grupo()$brix, main = "°Brix", xlab = "°Brix", col = "gray85", border = "white")
  })
  output$hist3 <- renderPlot({
    hist(selected_grupo()$cloreto, main = "Cloreto", xlab = "%", col = "gray85", border = "white")
  })
  output$hist4 <- renderPlot({
    hist(selected_grupo()$cafeina, main = "Cafeína", xlab = "% b.s.", col = "gray85", border = "white")
  })
  output$hist5 <- renderPlot({
    hist(selected_grupo()$ac_clorogenico, main = "Ác. Clorogênico", xlab = "% b.s.", col = "gray85", border = "white")
  })
  
  rotulos_box_quim <- c(
    umidade        = "Umidade (g/100g)",
    brix           = "°Brix",
    cloreto        = "Cloreto (%)",
    cafeina        = "Cafeína (% b.s.)",
    ac_clorogenico = "Ác. Clorogênico (% b.s.)"
  )
  
  output$boxplots_exp <- renderPlot({
    df_long <- df() %>%
      dplyr::select(grupo_peneira, all_of(vars_quim)) %>%
      tidyr::pivot_longer(
        cols = all_of(vars_quim),
        names_to = "variavel",
        values_to = "valor"
      ) %>%
      mutate(variavel = factor(rotulos_box_quim[variavel], levels = rotulos_box_quim[vars_quim]))
    
    p <- ggplot(df_long, aes(x = grupo_peneira, y = valor, fill = grupo_peneira)) +
      geom_boxplot(alpha = 0.65, width = 0.55, outlier.shape = if (input$box_pontos) NA else 19) +
      facet_wrap(~ variavel, scales = "free_y") +
      scale_fill_manual(values = c("Peneira Alta" = "#FC4E07", "Peneira Baixa" = "#00AFBB")) +
      labs(x = NULL, y = NULL, fill = "Grupo") +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "bold")
      )
    
    if (input$box_pontos) {
      p <- p + geom_jitter(width = 0.14, alpha = 0.6, size = 1.8, color = "gray25")
    }
    
    if (input$box_medias) {
      p <- p + stat_summary(
        fun = mean, geom = "point", shape = 23, size = 3,
        fill = "white", color = "black", stroke = 1
      )
    }
    
    p
  })
  
  output$qqplot <- renderPlot({
    qqmahala(selected_grupo())
  }, res = 96)
  
  # Elipsoide 3D
  output$elipse <- renderPlotly({
    vars_3d <- c(input$eixo_x, input$eixo_y, input$eixo_z)
    validate(
      need(length(unique(vars_3d)) == 3, "Selecione três variáveis distintas para visualizar o elipsoide 3D.")
    )
    elips(delta0(), df(), vars_quim, vars_3d, alpha())
  })
  
  # MODULO 2 - SERVER
  
  # --- 2.1 ACC ---
  acc_obj <- reactive({
    Y1n <- scale(as.matrix(df() %>% dplyr::select(all_of(vars_quim))))
    Y2n <- scale(as.matrix(df() %>% dplyr::select(all_of(vars_fis))))
    candisc::cancor(Y1n, Y2n)
  })
  
  acc_bartlett <- reactive({
    cc <- acc_obj()
    n <- nrow(df()); p <- length(vars_quim); q <- length(vars_fis)
    lambda <- cc$cancor^2
    k_max <- min(p, q)
    
    estat <- numeric(k_max); gl <- numeric(k_max); pval <- numeric(k_max)
    for (i in 1:k_max) {
      estat[i] <- -(n - (p + q + 3) / 2) * sum(log(1 - lambda[i:k_max]))
      gl[i]    <- (p - i + 1) * (q - i + 1)
      pval[i]  <- pchisq(estat[i], df = gl[i], lower.tail = FALSE)
    }
    
    eigen_val <- lambda / (1 - lambda)
    tibble(
      Par           = paste0("Par ", 1:k_max),
      Cor_Canonica  = cc$cancor,
      R2_Canonico   = lambda,
      Perc_Assoc    = (eigen_val / sum(eigen_val)) * 100,
      Acum_Perc     = cumsum(Perc_Assoc),
      Qui2_Bartlett = estat,
      gl            = gl,
      p_valor       = pval
    )
  })
  
  output$acc_tabela <- renderTable({
    acc_bartlett()
  }, digits = 4, align = "c")
  
  output$acc_conclusao <- renderUI({
    tab <- acc_bartlett()
    a <- input$acc_alpha
    p1 <- tab$p_valor[1]
    r1 <- round(tab$Cor_Canonica[1], 4)
    
    if (p1 < a) {
      HTML(paste0(
        "<p>Ao nível de <b>", a * 100, "%</b>, a primeira correlação canônica (", r1,
        ") é <b>estatisticamente significativa</b> (p = ", round(p1, 4), ").</p>"
      ))
    } else {
      HTML(paste0(
        "<p>Ao nível de <b>", a * 100, "%</b>, <b>nenhuma</b> correlação canônica atingiu significância ",
        "(1º par: r = ", r1, ", p = ", round(p1, 4), "). ",
        if (p1 < 0.10) "Observa-se tendência marginal (p < 0,10)." else "", "</p>"
      ))
    }
  })
  
  output$acc_cargas_plot <- renderPlot({
    cc <- acc_obj()
    k  <- as.integer(input$acc_par)
    
    if (input$acc_tipo_carga == "direta") {
      c_x <- cc$structure$X.xscores[, k]
      c_y <- cc$structure$Y.yscores[, k]
      sub_txt <- paste0("Cargas Diretas — Par ", k)
    } else {
      c_x <- cc$structure$X.yscores[, k]
      c_y <- cc$structure$Y.xscores[, k]
      sub_txt <- paste0("Cargas Cruzadas — Par ", k)
    }
    
    df_cargas <- bind_rows(
      tibble(Variavel = names(c_x), Carga = as.numeric(c_x), Bloco = paste0("Química (u", k, ")")),
      tibble(Variavel = names(c_y), Carga = as.numeric(c_y), Bloco = paste0("Física (v", k, ")"))
    )
    
    ggplot(df_cargas, aes(x = reorder(Variavel, abs(Carga)), y = Carga, fill = Carga > 0)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_flip() +
      facet_wrap(~ Bloco, scales = "free_y") +
      scale_fill_manual(values = c("TRUE" = "#00AFBB", "FALSE" = "#FC4E07")) +
      ylim(-1, 1) +
      labs(title = sub_txt, x = NULL, y = "Correlação") +
      theme_minimal()
  })
  
  output$acc_dispersao <- renderPlotly({
    cc <- acc_obj()
    k  <- as.integer(input$acc_par)
    
    df_sc <- tibble(
      u = cc$scores$X[, k],
      v = cc$scores$Y[, k],
      Amostra = df()$amostra,
      Grupo = df()$grupo_peneira
    )
    
    p <- ggplot(df_sc, aes(x = u, y = v, color = Grupo)) +
      geom_point(aes(text = paste("Amostra:", Amostra)), size = 2.5, alpha = 0.85) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "gray30", linetype = "dashed") +
      scale_color_manual(values = c("Peneira Alta" = "#FC4E07", "Peneira Baixa" = "#00AFBB")) +
      labs(
        title = paste0("Escores Canônicos (r = ", round(cc$cancor[k], 3), ")"),
        x = paste0("u", k, " (Química)"),
        y = paste0("v", k, " (Física)")
      ) +
      theme_minimal()
    
    ggplotly(p, tooltip = c("text", "color", "x", "y"))
  })
  
  # --- 2.2 ACP ---
  acp_obj <- reactive({
    if (input$acp_modo == "quim") {
      prcomp(df() %>% dplyr::select(all_of(vars_quim)), scale. = TRUE)
    } else if (input$acp_modo == "full") {
      prcomp(df() %>% dplyr::select(all_of(c(vars_quim, vars_fis))), scale. = TRUE)
    } else {
      FactoMineR::PCA(
        df() %>% dplyr::select(all_of(c(vars_quim, vars_fis))),
        quanti.sup = 6:9, scale.unit = TRUE, graph = FALSE
      )
    }
  })
  
  output$acp_scree <- renderPlot({
    factoextra::fviz_eig(acp_obj(), addlabels = TRUE, ylim = c(0, 60),
                         main = "Scree Plot (% Variância Explicada)")
  })
  
  output$acp_plot_principal <- renderPlot({
    eixos <- c(as.integer(input$acp_eixo1), as.integer(input$acp_eixo2))
    validate(need(eixos[1] != eixos[2], "Selecione dois eixos distintos."))
    
    if (input$acp_modo == "sup") {
      factoextra::fviz_pca_var(
        acp_obj(),
        axes = eixos,
        col.var = "black",
        col.quanti.sup = "steelblue",
        repel = TRUE,
        title = "Círculo de Correlações (Química Ativa + Peneiras Suplementares)"
      )
    } else {
      factoextra::fviz_pca_biplot(
        acp_obj(),
        axes = eixos,
        repel = TRUE,
        col.var = "black",
        col.ind = df()$grupo_peneira,
        palette = c("#FC4E07", "#00AFBB"),
        addEllipses = input$acp_elipses,
        ellipse.type = "confidence",
        legend.title = "Grupo",
        title = "Biplot ACP"
      )
    }
  })
  
  output$acp_autovalores <- renderTable({
    factoextra::get_eigenvalue(acp_obj()) %>% 
      as.data.frame()%>% 
      rownames_to_column("Componente") %>%
      dplyr::rename(
        Autovalor  = eigenvalue,
        `Var (%)`  = variance.percent,
        `Acum (%)` = cumulative.variance.percent
      )
  }, digits = 3)
  
  output$acp_cargas <- renderTable({
    eixos <- c(as.integer(input$acp_eixo1), as.integer(input$acp_eixo2))
    validate(need(eixos[1] != eixos[2], ""))
    
    v_acp <- factoextra::get_pca_var(acp_obj())
    tab <- data.frame(
      Variavel   = rownames(v_acp$coord),
      Eixo_X     = v_acp$coord[, eixos[1]],
      Eixo_Y     = v_acp$coord[, eixos[2]],
      cos2_Plano = rowSums(v_acp$cos2[, eixos, drop = FALSE])
    )
    
    if (input$acp_modo == "sup") {
      sup <- acp_obj()$quanti.sup
      tab_sup <- data.frame(
        Variavel   = paste0(rownames(sup$cor), " (sup)"),
        Eixo_X     = sup$cor[, eixos[1]],
        Eixo_Y     = sup$cor[, eixos[2]],
        cos2_Plano = rowSums(sup$cos2[, eixos, drop = FALSE])
      )
      tab <- bind_rows(tab, tab_sup)
    }
    tab
  }, digits = 3)
  
  # --- 2.3 ANÁLISE DISCRIMINANTE ---
  disc_formula <- reactive({
    validate(need(length(input$disc_vars) >= 2, "Selecione ao menos 2 variáveis químicas."))
    as.formula(paste("grupo_peneira ~", paste(input$disc_vars, collapse = " + ")))
  })
  
  disc_fit_full <- reactive({
    if (input$disc_metodo == "lda") {
      MASS::lda(disc_formula(), data = df(), prior = c(0.5, 0.5))
    } else {
      MASS::qda(disc_formula(), data = df(), prior = c(0.5, 0.5))
    }
  })
  
  disc_eval_res <- reactive({
    if (input$disc_eval == "cv") {
      if (input$disc_metodo == "lda") {
        MASS::lda(disc_formula(), data = df(), prior = c(0.5, 0.5), CV = TRUE)
      } else {
        MASS::qda(disc_formula(), data = df(), prior = c(0.5, 0.5), CV = TRUE)
      }
    } else {
      predict(disc_fit_full(), df())
    }
  })
  
  output$disc_cm <- renderPrint({
    res <- disc_eval_res()
    caret::confusionMatrix(data = res$class, reference = df()$grupo_peneira)
  })
  
  output$disc_roc <- renderPlot({
    res <- disc_eval_res()
    roc_obj <- pROC::roc(df()$grupo_peneira, res$posterior[, "Peneira Alta"], quiet = TRUE)
    pROC::plot.roc(
      roc_obj, print.auc = TRUE, auc.polygon = TRUE,
      auc.polygon.col = "#00AFBB33", max.auc.polygon = TRUE,
      main = paste0("Curva ROC (", toupper(input$disc_metodo), " - ", toupper(input$disc_eval), ")")
    )
  })
  
  output$disc_importancia <- renderTable({
    validate(need(input$disc_metodo == "lda",
                  "Coeficientes padronizados e cargas LD1 aplicam-se apenas à LDA."))
    X <- df() %>% dplyr::select(all_of(input$disc_vars))
    fit <- disc_fit_full()
    pred <- predict(fit, df())
    
    g1 <- X[df()$grupo_peneira == "Peneira Alta", , drop = FALSE]
    g2 <- X[df()$grupo_peneira == "Peneira Baixa", , drop = FALSE]
    Sp <- ((nrow(g1) - 1) * cov(g1) + (nrow(g2) - 1) * cov(g2)) / (nrow(df()) - 2)
    
    tibble(
      Variavel    = colnames(X),
      Coef_Bruto  = as.numeric(coef(fit)),
      Coef_Padrao = as.numeric(coef(fit)) * sqrt(diag(Sp)),
      Carga_LD1   = as.numeric(cor(X, pred$x))
    ) %>%
      arrange(desc(abs(Carga_LD1)))
  }, digits = 3)
  
  output$disc_densidade <- renderPlot({
    validate(need(input$disc_metodo == "lda",
                  "O gráfico de densidade dos escores LD1 é exclusivo da LDA."))
    pred <- predict(disc_fit_full(), df())
    df_p <- tibble(LD1 = as.numeric(pred$x), Grupo = df()$grupo_peneira)
    
    ggplot(df_p, aes(x = LD1, fill = Grupo)) +
      geom_density(alpha = 0.45) +
      geom_rug(aes(color = Grupo)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("Peneira Alta" = "#FC4E07", "Peneira Baixa" = "#00AFBB")) +
      scale_color_manual(values = c("Peneira Alta" = "#FC4E07", "Peneira Baixa" = "#00AFBB")) +
      labs(title = "Distribuição dos Escores Discriminantes (LD1)", x = "LD1", y = "Densidade") +
      theme_minimal()
  })
  
  output$disc_pred_nova <- renderUI({
    nova <- data.frame(
      umidade        = input$pred_umidade,
      brix           = input$pred_brix,
      cloreto        = input$pred_cloreto,
      cafeina        = input$pred_cafeina,
      ac_clorogenico = input$pred_cloro
    )
    pred_n  <- predict(disc_fit_full(), newdata = nova)
    p_alta  <- round(pred_n$posterior[1, "Peneira Alta"] * 100, 1)
    p_baixa <- round(pred_n$posterior[1, "Peneira Baixa"] * 100, 1)
    
    HTML(paste0(
      "<b>Classe prevista:</b> ", as.character(pred_n$class), "<br>",
      "P(Peneira Alta) = <b>", p_alta, "%</b> | P(Peneira Baixa) = <b>", p_baixa, "%</b>"
    ))
  })
  
  # --- 2.4 AGRUPAMENTO (CLUSTER) ---
  X_cl_pad <- reactive({
    mat <- scale(df() %>% dplyr::select(all_of(vars_quim)))
    rownames(mat) <- paste0(seq_len(nrow(df())), "_", df()$amostra)
    mat
  })
  
  cl_dist <- reactive(dist(X_cl_pad(), method = "euclidean"))
  cl_hc   <- reactive(hclust(cl_dist(), method = input$cl_metodo_hc))
  
  cl_km <- reactive({
    set.seed(123)
    kmeans(X_cl_pad(), centers = input$cl_k, nstart = 25)
  })
  
  output$cl_cofen <- renderPrint({
    r_cof <- cor(cl_dist(), cophenetic(cl_hc()))
    cat("r cofenético =", round(r_cof, 4))
  })
  
  output$cl_kmeans_plot <- renderPlot({
    factoextra::fviz_cluster(
      cl_km(), data = X_cl_pad(),
      palette = "jco", ellipse.type = "convex",
      repel = TRUE, ggtheme = theme_minimal(),
      main = paste0("K-Means (k = ", input$cl_k, ") sobre Perfil Químico")
    )
  })
  
  output$cl_dendro_plot <- renderPlot({
    factoextra::fviz_dend(
      cl_hc(), k = input$cl_k, cex = 0.5, palette = "jco",
      rect = TRUE, rect_fill = TRUE,
      main = paste0("Dendrograma (", input$cl_metodo_hc, ")"),
      ylab = "Distância Euclidiana"
    )
  })
  
  output$cl_fisher <- renderPrint({
    tab <- table(Cluster_KM = cl_km()$cluster, Grupo = df()$grupo_peneira)
    print(tab)
    cat("\nProporções por Cluster (%):\n")
    print(round(prop.table(tab, margin = 1) * 100, 1))
    cat("\n")
    print(fisher.test(tab, simulate.p.value = (input$cl_k > 2)))
  })
  
  output$cl_perfis <- renderTable({
    df() %>%
      mutate(Cluster = factor(cl_km()$cluster)) %>%
      group_by(Cluster) %>%
      summarise(
        n = n(),
        across(all_of(c(vars_quim, vars_fis)), ~ round(mean(.x), 3)),
        .groups = "drop"
      )
  }, digits = 3)
}

## run app
shinyApp(ui = ui, server = server)