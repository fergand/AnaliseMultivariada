library(shiny)
library(bslib)
library(plotly)
library(tidyverse)
library(kableExtra)
library(ggthemes); set_theme(theme_clean())
library(rgl)
library(MVA)

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
    card(
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
          value = 0.05,
          min = 0.01,
          max = 0.10,
          step = 0.01
        ),
        fluidRow(
          column(6, actionButton("b", "Bonferroni", width = "100%")),
          column(6, actionButton("h", "Hotelling",  width = "100%"))
        )
      )
    ),
    card(
      tabsetPanel(
        id = "tabset",
        tabPanel(
          "Dados",
          br(),
          "Pré-visualização dos dados carregados (variáveis químicas e classe de peneira)",
          br(), br(),
          tableOutput("table1")
        ),
        tabPanel(
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
          fluidRow(
            column(
              6,
              verbatimTextOutput("hz_output"),
              verbatimTextOutput("boxm_output")
            ),
            column(6, plotOutput("qqplot"))
          )
        ),
        tabPanel(
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
        tabPanel(
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
      )
    ),
    col_widths = c(3, 9)
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
}

## run app
shinyApp(ui = ui, server = server)