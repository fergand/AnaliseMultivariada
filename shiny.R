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

dados <- read_csv("Dados/dados.csv") |> as_tibble()

## ui
ui <- fluidPage(
  page_fluid(
    layout_columns(
      card(
        fluidRow(
          card(
            withMathJax(helpText("VETOR DE MÉDIAS NULO (\\(\\boldsymbol{\\mu}_0\\))")),
            numericInput("mu01", "x1", value = 32),
            numericInput("mu02", "x2", value = 85),
            numericInput("mu03", "x3", value = 18)
          ),
          br(), # Espaço extra
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
              column(6,actionButton("b", "Bonferroni")),
              column(6,actionButton("h", "Hotelling"))
            )
          )
        )
      ),
      card("",
        tabsetPanel(
          id = "tabset",
          tabPanel(
            "Dados", 
            "Pré-visualização dos dados carregados",
            br("\n"),
            tableOutput("table1")
          ),
          tabPanel(
            "Exploratória", 
            "Distribuição individual dos dados e diagnóstico de normalidade",
            fluidRow(
              column(4, plotOutput("hist1")),
              column(4, plotOutput("hist2")),
              column(4, plotOutput("hist3"))
            ),
            fluidRow(
              column(6, verbatimTextOutput("hz_output")),
              column(6, plotOutput("qqplot"))
            )
          ),
          # Falta escrever as hipóteses e resultado do teste
          # Testando a hipótese \\(H_0\colon \boldsymbol{\mu} = \boldsymbol{\mu}_0\\)
          tabPanel(
            "Teste T2",
            withMathJax(helpText("Teste \\(T^2\\) de Hotelling: testando a hipótese \\(H_0\\colon \\boldsymbol{\\mu} = \\boldsymbol{\\mu}_0\\) vs \\(H_1\\colon \\boldsymbol{\\mu} \\ne \\boldsymbol{\\mu}_0\\)")),
            br("\n"),
            fluidRow(
              column(4, "Estatística \\(T^2\\)", textOutput("stats_t2")),
              column(4, "F crítico", textOutput("Fcalc")),
              column(4, "Valor-p", textOutput("valor_p"))
            ),
            br("\n"),
            fluidRow(
              column(12, uiOutput("t2_conclusao"))
            )
          ),
          tabPanel(
            "ICS e Elipsóide",
            br("\n"),
            "Intervalos de Confiança Simultâneos",
            tableOutput("tabics"),
            br("\n"),
            "Elipse de confiança",
            plotlyOutput("elipse", height = "600px")
            
          )
        )
      ),
      col_widths = c(3,9)
    )
  )
)

## functions para o server
qqmahala <- function(data) {
  xp <- scale(data, scale = FALSE)
  n <- dim(data)[1]
  p <- dim(data)[2]
  Sa <- cov(xp)
  D <- xp %*% solve(Sa) %*% t(xp)
  d <- diag(D)
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

t2test <- function(mu0, data, alpha) {
  mu0 <- mu0
  yb <- colMeans(data %>% select(-1))
  invS <- solve(cov(data %>% select(-1)))
  
  n <- nrow(data)
  p <- ncol(data %>% select(-1))
  
  ## T2
  T2 <- n*t(yb - mu0) %*% invS %*% (yb-mu0)
  
  ## transformação direta de T2 para F
  Fcalc <- ((n-p)/((n-1)*p))*n*t(yb - mu0) %*% invS %*% (yb-mu0);
  
  ## Ftab
  Ftab <- qf(1 - alpha, p, n-p)
  pvalue <- pf(Fcalc, p, n-p, lower.tail = FALSE)
  
  resultado <-  data.frame(
    t2 = T2, Fcalc = Fcalc, Ftab = Ftab, valorp = pvalue
  )
  resultado
}

ICS <- function(data, type = "b", alpha = 0.05) {
  n <- data |> nrow()
  p <- data |> dplyr::select(-1) |> ncol()
  yb <- data |> dplyr::select(-1) |> colMeans()
  S <- data |> dplyr::select(-1) |> cov()
  
  if (!(type %in% c("h", "b")))
    stop("Forneça valor para type: 'b' para Bonferroni ou 'h' para Hotelling")
  
  resultado <- data.frame(
    media = numeric(p), Var = numeric(p), LI = numeric(p), LS = numeric(p),
    a = character(p),row.names = names(yb)
  )
  
  for (i in 1:p) {
    a <- matrix(0, p, 1); a[i] <- 1
    
    if (type == "h") {
      erro <- sqrt(p*(n-1)/(n*(n-p)) * qf(1-alpha, p, n-p) * t(a) %*% S %*% a)
    } else {
      erro <- sqrt(qt(1 - alpha/(2*p), n-1)^2 * 1/n * t(a) %*% S %*% a)
    }
    
    LI <- t(a) %*% yb - erro
    LS <- t(a) %*% yb + erro
    
    resultado[i,] <- list(yb[i],S[i,i],LI,LS,paste(a[1:p],collapse = " "))
  }
  resultado
}

elips <- function(mu0, data, alpha) {
  
  yb <- colMeans(data %>% dplyr::select(-1))
  S  <- cov(data %>% dplyr::select(-1))
  n <- nrow(data)
  p <- ncol(data %>% dplyr::select(-1))
  
  elip <- rgl::ellipse3d(
    as.matrix(S),
    centre = yb,
    t = sqrt((n-1)*p/(n*(n-p))*qf(1 - alpha, p, n-p))
  )
  
  verts <- t(elip$vb[1:3,])
  faces <- t(elip$ib) - 1
  
  plot_ly() %>%
    add_trace(
      type = "mesh3d",
      x = verts[,1], y = verts[,2], z = verts[,3],
      i = faces[,1], j = faces[,2], k = faces[,3],
      color = "red",
      opacity = 0.15,
      contour = list(show = TRUE, color = "black", width = 2),
      hoverinfo = "skip"
    ) %>%
    add_markers(
      x = mu0[1], y = mu0[2], z = mu0[3],
      marker = list(color = "green", size = 6)
    ) %>%
    layout(scene = list(
      xaxis = list(title = "x1", showspikes = F),
      yaxis = list(title = "x2", showspikes = F),
      zaxis = list(title = "x3", showspikes = F)
    ))
}

## server

server <- function(input, output) {
  # reactives
  df <- reactive(dados)
  selected <- reactive(dados %>% dplyr::select(-1))
  
  hz_result <- reactive({
    mvnTest::HZ.test(selected(), qqplot = F)
  })
  
  t2_result <- reactive({
    
  })
  
  mu0 <- reactive(c(input$mu01, input$mu02, input$mu03))
  alpha <- reactive(input$alpha)
  
  #tint -> tipo de intervalo
  v <- reactiveValues(value = "b")
  observeEvent(input$b, {
    v$value <- "b"
  })
  observeEvent(input$h, {
    v$value <- "h"
  })
  
  # Teste Normalidade
  output$hz_output <- renderPrint({
    res <- hz_result()
    alpha <- input$alpha
    
    cat("Teste de Henze-Zirkler para Normal Multivariada\n\n")
    cat("  HZ      :", res@HZ, "\n")
    cat("  p-value :", res@p.value, "\n\n")
    cat("  Result  :", 
        if (res@p.value > alpha) {
          paste0("Os dados seguem a dist. Normal Multivariada (sig.level = ", alpha, ")")
        } else {
          paste0("Os dados NÃO seguem a dist. Normal Multivariada (sig.level = ", alpha, ")")
        })
  })
  
  # T2 de Hotelling
  output$stats_t2 <- renderText(
    t2test(mu0(), df(), alpha()) %>% dplyr::select(1) %>% pull() %>% round(2)
  )
  output$Fcalc <- renderText(
    t2test(mu0(), df(), alpha()) %>% dplyr::select(2) %>% pull() %>% round(3)
  )
  output$valor_p <- renderText(
    t2test(mu0(), df(), alpha()) %>% dplyr::select(4) %>% pull() %>% round(3)
  )
  output$t2_conclusao <- renderUI({
    p_value <- t2test(mu0(), df(), alpha()) %>% dplyr::select(4) %>% pull() %>% round(3)
    mu0_txt <- paste(round(mu0(), 2), collapse = ", ")
    
    if (p_value < alpha()) {
      withMathJax(HTML(paste0(
        "Há evidências a favor de que a verdadeira média das observações, ao nível de ",
        alpha() * 100, "%, NÃO seja o vetor \\(\\boldsymbol{\\mu_0} = (", mu0_txt, ")\\)"
      )))
    } else {
      withMathJax(HTML(paste0(
        "Há evidências a favor de que a verdadeira média das observações, ao nível de ",
        alpha() * 100, "%, seja o vetor \\(\\boldsymbol{\\mu_0} = (", mu0_txt, ").\\)"
      )))
    }
  })

  
  # tabelas
  output$table1 <- renderTable(df(), auto = T)
  output$tabics <- renderTable({
    ICS(data = df(), type = v$value, alpha = alpha())
    
  }, rownames = T, align = "c")
  
  # graficos
  output$hist1 <- renderPlot({
    selected() %>% select(1) %>% pull() %>% hist(main = "Histograma de X1")
  })
  output$hist2 <- renderPlot({
    selected() %>% select(2) %>% pull() %>% hist(main = "Histograma de X2")
  })
  output$hist3 <- renderPlot({
    selected() %>% select(3) %>% pull %>% hist(main = "Histograma de X3")
  })
  output$qqplot <- renderPlot(
    selected() %>% qqmahala(),
    res = 96
  )
  output$elipse <- renderPlotly({
    elips(mu0(), df(), alpha())
  })
  
  
}


## run app
shinyApp(ui = ui, server = server)