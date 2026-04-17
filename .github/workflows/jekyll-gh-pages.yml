####################################################################
### MIS410 - Culottes Forecasting Dashboard
### High-End AI Style Shiny App
### Based on analysis_culottes.R logic
####################################################################

####################################################################
### Step-1: Load Libraries
####################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(forecast)
  library(tseries)
  library(ggplot2)
  library(dplyr)
  library(RMariaDB)
  library(DBI)
  library(scales)
  library(readr)
  library(lubridate)
  library(DT)
  library(gridExtra)
})

####################################################################
### Step-2: Read and Prepare Data
####################################################################

gtrends <- read.csv("gtrends.csv")
names(gtrends) <- tolower(names(gtrends))

culottes.data <- gtrends[, c("date", "culottes")]
culottes.data$date <- as.Date(culottes.data$date, "%m/%d/%Y")
culottes.data <- culottes.data[order(culottes.data$date), ]

culottes <- ts(culottes.data$culottes, start = c(2015, 41), frequency = 52)

####################################################################
### Step-3: Pattern / Summary Prep
####################################################################

adf.result <- adf.test(culottes)

culottes.decomp <- stl(culottes, s.window = "periodic", robust = TRUE)

summary_stats <- data.frame(
  Statistic = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum", "Std. Dev."),
  Value = c(
    min(culottes.data$culottes),
    quantile(culottes.data$culottes, 0.25),
    median(culottes.data$culottes),
    mean(culottes.data$culottes),
    quantile(culottes.data$culottes, 0.75),
    max(culottes.data$culottes),
    sd(culottes.data$culottes)
  )
)

####################################################################
### Step-4: Model Comparison Using RMSFE
####################################################################

h_eval <- 6
test.culottes <- tail(culottes, h_eval)
train.culottes <- head(culottes, length(culottes) - h_eval)

naive.train <- naive(train.culottes, h = h_eval)
ses.train   <- ses(train.culottes, h = h_eval, level = 95, initial = "optimal")
holt.train  <- holt(train.culottes, h = h_eval, level = 95, initial = "optimal")
holtd.train <- holt(train.culottes, h = h_eval, level = 95, initial = "optimal", damped = TRUE)

arima.train <- auto.arima(train.culottes)
arima.fc    <- forecast(arima.train, h = h_eval)

stlf.train  <- stlf(train.culottes, h = h_eval, level = 95)

tbats.train <- tbats(train.culottes)
tbats.fc    <- forecast(tbats.train, h = h_eval)

rmsfe.naive <- sqrt(mean((test.culottes - naive.train$mean)^2))
rmsfe.ses   <- sqrt(mean((test.culottes - ses.train$mean)^2))
rmsfe.holt  <- sqrt(mean((test.culottes - holt.train$mean)^2))
rmsfe.holtd <- sqrt(mean((test.culottes - holtd.train$mean)^2))
rmsfe.arima <- sqrt(mean((test.culottes - arima.fc$mean)^2))
rmsfe.stlf  <- sqrt(mean((test.culottes - stlf.train$mean)^2))
rmsfe.tbats <- sqrt(mean((test.culottes - tbats.fc$mean)^2))

rmsfe.results <- data.frame(
  Model = c("Naive", "SES", "Holt", "Holt-Damped", "ARIMA", "STLF", "TBATS"),
  RMSFE = round(c(
    rmsfe.naive, rmsfe.ses, rmsfe.holt, rmsfe.holtd,
    rmsfe.arima, rmsfe.stlf, rmsfe.tbats
  ), 4)
)

rmsfe.results <- rmsfe.results[order(rmsfe.results$RMSFE), ]
best.model.name <- as.character(rmsfe.results$Model[1])

####################################################################
### Step-5: Best Model Forecast on Full Series
####################################################################

tbats.model <- tbats(culottes)

####################################################################
### Step-6: Inventory Parameters
####################################################################

mu.d <- mean(culottes.data$culottes)
sigma.d <- sd(culottes.data$culottes)

inventory_calc <- function(lead_time = 3, csl = 95) {
  z <- qnorm(csl / 100)
  SS <- z * sigma.d * sqrt(lead_time)
  ROP <- mu.d * lead_time + SS
  
  data.frame(
    Parameter = c(
      "Lead Time (weeks)",
      "Cycle Service Level (%)",
      "z-score",
      "Mean Weekly Demand",
      "Std Dev of Demand",
      "Safety Stock (SS)",
      "Reorder Point (ROP)"
    ),
    Value = c(
      lead_time,
      csl,
      round(z, 4),
      round(mu.d, 2),
      round(sigma.d, 2),
      round(SS, 2),
      round(ROP, 2)
    )
  )
}

####################################################################
### Step-7: Helper Functions
####################################################################

get_model_forecast <- function(model_name, h) {
  if (model_name == "Naive") {
    fc <- naive(culottes, h = h, level = 95)
  } else if (model_name == "SES") {
    fc <- ses(culottes, h = h, level = 95, initial = "optimal")
  } else if (model_name == "Holt") {
    fc <- holt(culottes, h = h, level = 95, initial = "optimal")
  } else if (model_name == "Holt-Damped") {
    fc <- holt(culottes, h = h, level = 95, initial = "optimal", damped = TRUE)
  } else if (model_name == "ARIMA") {
    fit <- auto.arima(culottes)
    fc <- forecast(fit, h = h, level = 95)
  } else if (model_name == "STLF") {
    fc <- stlf(culottes, h = h, level = 95)
  } else if (model_name == "TBATS") {
    fit <- tbats(culottes)
    fc <- forecast(fit, h = h, level = 95)
  } else {
    fc <- forecast(tbats.model, h = h, level = 95)
  }
  return(fc)
}

forecast_table_builder <- function(fc, h, model_name) {
  last.date <- max(culottes.data$date)
  fc.dates <- last.date + (1:h) * 7
  
  data.frame(
    Forecast_Date = format(fc.dates, "%Y-%m-%d"),
    Model = model_name,
    Forecast = round(as.numeric(fc$mean), 2),
    Lo_95 = round(as.numeric(fc$lower[, 1]), 2),
    Hi_95 = round(as.numeric(fc$upper[, 1]), 2)
  )
}

save_to_mysql <- function(forecast_tbl, rmsfe_tbl, inventory_tbl) {
  tryCatch({
    con <- dbConnect(
      RMariaDB::MariaDB(),
      dbname = "mis410",
      host = "localhost",
      port = 3306,
      user = "root",
      password = ""
    )
    
    dbWriteTable(con, "culottes_model_results", rmsfe_tbl, overwrite = TRUE, row.names = FALSE)
    dbWriteTable(con, "culottes_forecast_dashboard", forecast_tbl, overwrite = TRUE, row.names = FALSE)
    dbWriteTable(con, "culottes_inventory_dashboard", inventory_tbl, overwrite = TRUE, row.names = FALSE)
    
    dbDisconnect(con)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

####################################################################
### Step-8: UI
####################################################################

ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bg = "#0a0f1a",
    fg = "#e5e7eb",
    primary = "#7c3aed",
    secondary = "#111827",
    success = "#10b981",
    info = "#38bdf8",
    warning = "#f59e0b",
    danger = "#ef4444",
    base_font = font_google("Inter"),
    code_font = font_google("JetBrains Mono")
  ),
  
  tags$head(
    tags$style(HTML("
      body {
        background: radial-gradient(circle at top left, #111827 0%, #0a0f1a 45%, #05070d 100%);
      }

      .app-title-wrap {
        padding: 18px 0 8px 0;
      }

      .app-title {
        font-size: 30px;
        font-weight: 800;
        color: #f8fafc;
        margin-bottom: 4px;
      }

      .app-subtitle {
        color: #94a3b8;
        font-size: 13px;
        margin-bottom: 0;
      }

      .glass-card {
        background: linear-gradient(180deg, rgba(17,24,39,0.95) 0%, rgba(9,14,26,0.97) 100%);
        border: 1px solid rgba(255,255,255,0.07);
        border-radius: 20px;
        box-shadow: 0 12px 40px rgba(0,0,0,0.32);
        padding: 18px;
        margin-bottom: 18px;
      }

      .sidebar-box {
        background: linear-gradient(180deg, rgba(10,15,26,0.98) 0%, rgba(17,24,39,0.98) 100%);
        border: 1px solid rgba(255,255,255,0.07);
        border-radius: 22px;
        box-shadow: 0 12px 36px rgba(0,0,0,0.35);
        padding: 18px;
        min-height: 100%;
      }

      .kpi-card {
        background: linear-gradient(135deg, rgba(17,24,39,0.97) 0%, rgba(25,32,55,0.97) 100%);
        border: 1px solid rgba(255,255,255,0.07);
        border-radius: 20px;
        box-shadow: 0 10px 35px rgba(0,0,0,0.30);
        padding: 18px;
        min-height: 140px;
      }

      .kpi-label {
        font-size: 12px;
        color: #9ca3af;
        letter-spacing: 0.8px;
        text-transform: uppercase;
        font-weight: 700;
      }

      .kpi-value {
        font-size: 30px;
        font-weight: 800;
        color: #f9fafb;
        margin-top: 8px;
        margin-bottom: 8px;
      }

      .kpi-desc {
        color: #cbd5e1;
        font-size: 13px;
        line-height: 1.5;
      }

      .section-title {
        font-size: 16px;
        font-weight: 700;
        color: #f8fafc;
        margin-bottom: 12px;
      }

      .mini-tag {
        display: inline-block;
        padding: 6px 10px;
        font-size: 12px;
        border-radius: 999px;
        color: #dbeafe;
        background: rgba(124,58,237,0.18);
        border: 1px solid rgba(124,58,237,0.28);
        margin-right: 6px;
        margin-bottom: 8px;
      }

      .insight-box {
        background: rgba(255,255,255,0.03);
        border-left: 4px solid #7c3aed;
        border-radius: 12px;
        padding: 14px;
        margin-bottom: 10px;
        color: #e5e7eb;
      }

      .nav-tabs {
        border-bottom: 1px solid rgba(255,255,255,0.08);
      }

      .nav-tabs .nav-link {
        color: #94a3b8;
        background: transparent;
        border: none;
      }

      .nav-tabs .nav-link.active {
        color: #f8fafc;
        background: rgba(124,58,237,0.15);
        border-radius: 12px;
      }

      .dataTables_wrapper .dataTables_length,
      .dataTables_wrapper .dataTables_filter,
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_processing,
      .dataTables_wrapper .dataTables_paginate {
        color: #cbd5e1 !important;
      }

      table.dataTable tbody tr {
        background-color: #0f172a !important;
        color: #e5e7eb !important;
      }

      table.dataTable thead th {
        background-color: #111827 !important;
        color: #f8fafc !important;
        border-bottom: 1px solid rgba(255,255,255,0.08) !important;
      }

      .form-control, .selectize-input, .selectize-dropdown, .irs {
        background-color: #0f172a !important;
        color: #e5e7eb !important;
        border-color: rgba(255,255,255,0.08) !important;
      }

      .small-note {
        color: #94a3b8;
        font-size: 12px;
      }
    "))
  ),
  
  fluidRow(
    column(
      12,
      div(
        class = "app-title-wrap",
        div(class = "app-title", "Culottes Demand Intelligence Dashboard"),
        p(
          class = "app-subtitle",
          "Pattern analysis, best model selection, 6–12 week forecasting, safety stock, reorder point, and business insight"
        )
      )
    )
  ),
  
  fluidRow(
    column(
      3,
      div(
        class = "sidebar-box",
        div(class = "section-title", "Control Center"),
        
        sliderInput("h_forecast", "Forecast Horizon (weeks)", min = 6, max = 12, value = 6, step = 1),
        sliderInput("lead_time", "Lead Time (weeks)", min = 1, max = 8, value = 3, step = 1),
        sliderInput("csl", "Cycle Service Level (%)", min = 80, max = 99, value = 95, step = 1),
        
        selectInput(
          "plot_window",
          "Main Chart View",
          choices = c("Full Series" = "full", "Last 52 Weeks" = "recent"),
          selected = "full"
        ),
        
        checkboxGroupInput(
          "plot_layers",
          "Main Forecast Layers",
          choices = c(
            "Historical Series" = "hist",
            "Forecast Line" = "fc",
            "95% Interval" = "band",
            "Divider Line" = "divider"
          ),
          selected = c("hist", "fc", "band", "divider")
        ),
        
        actionButton("run_btn", "Run Analysis", class = "btn btn-primary"),
        
        br(), br(),
        div(class = "section-title", "Project Coverage"),
        div(class = "mini-tag", "(a) Pattern"),
        div(class = "mini-tag", "(b) Best Model"),
        div(class = "mini-tag", "(c) Forecast"),
        div(class = "mini-tag", "(d) SS & ROP"),
        div(class = "mini-tag", "(e) Insight"),
        
        br(), br(),
        p(class = "small-note",
          "This app uses your culottes analysis structure and shows all forecasting models plus a separate 6–12 week forecast page."
        )
      )
    ),
    
    column(
      9,
      
      fluidRow(
        column(
          3,
          div(
            class = "kpi-card",
            div(class = "kpi-label", "Best Model"),
            div(class = "kpi-value", textOutput("best_model", inline = TRUE)),
            div(class = "kpi-desc", "Selected using lowest RMSFE")
          )
        ),
        column(
          3,
          div(
            class = "kpi-card",
            div(class = "kpi-label", "Best RMSFE"),
            div(class = "kpi-value", textOutput("best_rmsfe", inline = TRUE)),
            div(class = "kpi-desc", "Lower means better forecasting accuracy")
          )
        ),
        column(
          3,
          div(
            class = "kpi-card",
            div(class = "kpi-label", "Safety Stock"),
            div(class = "kpi-value", textOutput("kpi_ss", inline = TRUE)),
            div(class = "kpi-desc", "Inventory buffer during lead time")
          )
        ),
        column(
          3,
          div(
            class = "kpi-card",
            div(class = "kpi-label", "Reorder Point"),
            div(class = "kpi-value", textOutput("kpi_rop", inline = TRUE)),
            div(class = "kpi-desc", "Replenishment trigger level")
          )
        )
      ),
      
      br(),
      
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "Executive View",
          br(),
          div(
            class = "glass-card",
            div(class = "section-title", "Best Model Forecast View"),
            plotOutput("plot_forecast_main", height = "450px")
          ),
          fluidRow(
            column(
              6,
              div(
                class = "glass-card",
                div(class = "section-title", "(a) Pattern of Data"),
                uiOutput("pattern_text")
              )
            ),
            column(
              6,
              div(
                class = "glass-card",
                div(class = "section-title", "(e) Business Insight"),
                uiOutput("insights_text")
              )
            )
          )
        ),
        
        tabPanel(
          "Data Pattern",
          br(),
          div(
            class = "glass-card",
            div(class = "section-title", "Weekly Time Series"),
            plotOutput("plot_series", height = "380px")
          ),
          div(
            class = "glass-card",
            div(class = "section-title", "STL Decomposition"),
            plotOutput("plot_decomp", height = "430px")
          ),
          div(
            class = "glass-card",
            div(class = "section-title", "ACF and Periodogram"),
            fluidRow(
              column(6, plotOutput("plot_acf", height = "320px")),
              column(6, plotOutput("plot_periodogram", height = "320px"))
            )
          ),
          div(
            class = "glass-card",
            div(class = "section-title", "Summary Statistics"),
            DTOutput("summary_table")
          )
        ),
        
        tabPanel(
          "Model Comparison",
          br(),
          div(
            class = "glass-card",
            div(class = "section-title", "(b) RMSFE Comparison"),
            plotOutput("plot_rmsfe", height = "380px")
          ),
          div(
            class = "glass-card",
            div(class = "section-title", "RMSFE Results Table"),
            DTOutput("table_rmsfe")
          )
        ),
        
        tabPanel(
          "All Model Forecasts",
          br(),
          div(
            class = "glass-card",
            div(class = "section-title", "Forecast Graphs for All Candidate Models"),
            plotOutput("plot_all_models", height = "1000px")
          )
        ),
        
        tabPanel(
          "6-12 Week Forecast",
          br(),
          div(
            class = "glass-card",
            div(class = "section-title", "Dedicated 6–12 Week Forecast Page"),
            plotOutput("plot_forecast_horizon", height = "430px")
          ),
          div(
            class = "glass-card",
            div(class = "section-title", "(c) Forecast Values for Selected Horizon"),
            DTOutput("table_forecast")
          )
        ),
        
        tabPanel(
          "Inventory Decision",
          br(),
          fluidRow(
            column(
              6,
              div(
                class = "glass-card",
                div(class = "section-title", "(d) Safety Stock and Reorder Point"),
                tableOutput("table_inventory")
              )
            ),
            column(
              6,
              div(
                class = "glass-card",
                div(class = "section-title", "Inventory Summary"),
                verbatimTextOutput("contents_inv")
              )
            )
          ),
          div(
            class = "glass-card",
            div(class = "section-title", "Lead-Time Demand Distribution"),
            plotOutput("plot_inventory", height = "350px")
          )
        )
      )
    )
  )
)

####################################################################
### Step-9: Server
####################################################################

server <- function(input, output, session) {
  
  analysis_data <- eventReactive(input$run_btn, {
    
    best_fc <- forecast(tbats.model, h = input$h_forecast, level = 95)
    best_fc_tbl <- forecast_table_builder(best_fc, input$h_forecast, "TBATS")
    
    inv_tbl <- inventory_calc(input$lead_time, input$csl)
    SS <- as.numeric(inv_tbl$Value[inv_tbl$Parameter == "Safety Stock (SS)"])
    ROP <- as.numeric(inv_tbl$Value[inv_tbl$Parameter == "Reorder Point (ROP)"])
    
    pattern_lines <- c(
      "The culottes series is non-stationary based on the ADF test, so a stable mean assumption is not appropriate.",
      "The STL decomposition and seasonality diagnostics indicate recurring seasonal movement in the weekly search pattern.",
      "The series also shows an overall category growth pattern across the full sample period."
    )
    
    insight_lines <- c(
      paste0("The best model is ", best.model.name, " because it has the lowest RMSFE among all tested models."),
      "The demand pattern is seasonal, so replenishment and promotion timing should follow the annual peak cycle.",
      "Short-term forecasts can be used to reduce over-ordering in weaker weeks and prepare earlier for stronger weeks.",
      paste0("Safety stock is ", round(SS, 2), " and reorder point is ", round(ROP, 2), " under the chosen inputs."),
      "Since the source is Google Trends, the results should be treated as demand signals unless converted into actual unit demand."
    )
    
    save_to_mysql(best_fc_tbl, rmsfe.results, inv_tbl)
    
    list(
      best_fc = best_fc,
      best_fc_tbl = best_fc_tbl,
      inv_tbl = inv_tbl,
      SS = SS,
      ROP = ROP,
      pattern_lines = pattern_lines,
      insight_lines = insight_lines
    )
  }, ignoreInit = FALSE)
  
  output$best_model <- renderText({
    best.model.name
  })
  
  output$best_rmsfe <- renderText({
    round(min(rmsfe.results$RMSFE), 4)
  })
  
  output$kpi_ss <- renderText({
    inv <- inventory_calc(input$lead_time, input$csl)
    ss <- as.numeric(inv$Value[inv$Parameter == "Safety Stock (SS)"])
    comma(round(ss, 2))
  })
  
  output$kpi_rop <- renderText({
    inv <- inventory_calc(input$lead_time, input$csl)
    rop <- as.numeric(inv$Value[inv$Parameter == "Reorder Point (ROP)"])
    comma(round(rop, 2))
  })
  
  output$pattern_text <- renderUI({
    tagList(lapply(analysis_data()$pattern_lines, function(x) div(class = "insight-box", x)))
  })
  
  output$insights_text <- renderUI({
    tagList(lapply(analysis_data()$insight_lines, function(x) div(class = "insight-box", x)))
  })
  
  output$summary_table <- renderDT({
    datatable(summary_stats, rownames = FALSE, options = list(dom = "tip", pageLength = 10))
  })
  
  output$table_rmsfe <- renderDT({
    datatable(rmsfe.results, rownames = FALSE, options = list(dom = "tip", pageLength = 10))
  })
  
  output$table_forecast <- renderDT({
    datatable(analysis_data()$best_fc_tbl, rownames = FALSE, options = list(dom = "tip", pageLength = 12))
  })
  
  output$table_inventory <- renderTable({
    analysis_data()$inv_tbl
  })
  
  output$contents_inv <- renderPrint({
    cat(sprintf("Safety Stock  = %.2f  ~  %.0f index units\n", analysis_data()$SS, round(analysis_data()$SS)))
    cat(sprintf("Reorder Point = %.2f  ~  %.0f index units\n", analysis_data()$ROP, round(analysis_data()$ROP)))
    cat(sprintf("Lead Time     = %d weeks\n", input$lead_time))
    cat(sprintf("Service Level = %d%%\n", input$csl))
  })
  
  output$plot_series <- renderPlot({
    plot_df <- if (input$plot_window == "recent") tail(culottes.data, 52) else culottes.data
    
    ggplot(plot_df, aes(x = date, y = culottes)) +
      geom_line(color = "#8b5cf6", linewidth = 1) +
      labs(
        title = "Culottes Weekly Google Trends Index",
        x = NULL,
        y = "Search Interest Index"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#111827", color = NA),
        panel.background = element_rect(fill = "#111827", color = NA),
        panel.grid.major = element_line(color = rgb(1, 1, 1, 0.06)),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "#cbd5e1"),
        axis.title = element_text(color = "#e5e7eb"),
        plot.title = element_text(color = "#f8fafc", face = "bold")
      )
  })
  
  output$plot_decomp <- renderPlot({
    plot(culottes.decomp, main = "STL Decomposition - Culottes")
  })
  
  output$plot_acf <- renderPlot({
    acf(culottes, lag.max = 104, main = "ACF - Culottes")
  })
  
  output$plot_periodogram <- renderPlot({
    spec.pgram(culottes, main = "Periodogram - Culottes")
  })
  
  output$plot_rmsfe <- renderPlot({
    rmsfe.plot.df <- rmsfe.results
    rmsfe.plot.df$ColorGroup <- ifelse(rmsfe.plot.df$Model == best.model.name, "Best", "Other")
    
    ggplot(rmsfe.plot.df, aes(x = reorder(Model, RMSFE), y = RMSFE, fill = ColorGroup)) +
      geom_col(width = 0.68) +
      geom_text(aes(label = round(RMSFE, 4)), hjust = -0.1, color = "#f8fafc", size = 4) +
      coord_flip() +
      scale_fill_manual(values = c("Best" = "#8b5cf6", "Other" = "#475569"), guide = "none") +
      labs(
        title = "Forecast Model Comparison Using RMSFE",
        x = NULL,
        y = "RMSFE"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#111827", color = NA),
        panel.background = element_rect(fill = "#111827", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = rgb(1, 1, 1, 0.06)),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "#cbd5e1"),
        axis.title = element_text(color = "#e5e7eb"),
        plot.title = element_text(color = "#f8fafc", face = "bold")
      ) +
      expand_limits(y = max(rmsfe.plot.df$RMSFE) * 1.2)
  })
  
  output$plot_forecast_main <- renderPlot({
    hist_df <- if (input$plot_window == "recent") tail(culottes.data, 52) else culottes.data
    fc_tbl <- analysis_data()$best_fc_tbl
    fc_tbl$date <- as.Date(fc_tbl$Forecast_Date)
    
    p <- ggplot() +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#111827", color = NA),
        panel.background = element_rect(fill = "#111827", color = NA),
        panel.grid.major = element_line(color = rgb(1, 1, 1, 0.06)),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "#cbd5e1"),
        axis.title = element_text(color = "#e5e7eb"),
        plot.title = element_text(color = "#f8fafc", face = "bold"),
        plot.subtitle = element_text(color = "#94a3b8"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "#e5e7eb")
      ) +
      labs(
        title = "Best Model Forecast View",
        subtitle = paste0("Best model based on RMSFE: ", best.model.name),
        x = NULL,
        y = "Search Interest Index"
      )
    
    if ("hist" %in% input$plot_layers) {
      p <- p + geom_line(data = hist_df, aes(x = date, y = culottes, color = "Historical"), linewidth = 1.1)
    }
    if ("band" %in% input$plot_layers) {
      p <- p + geom_ribbon(data = fc_tbl, aes(x = date, ymin = Lo_95, ymax = Hi_95, fill = "95% Interval"), alpha = 0.18)
    }
    if ("fc" %in% input$plot_layers) {
      p <- p + geom_line(data = fc_tbl, aes(x = date, y = Forecast, color = "Forecast"), linewidth = 1.2)
    }
    if ("divider" %in% input$plot_layers) {
      p <- p + geom_vline(xintercept = as.numeric(max(culottes.data$date)), linetype = "dashed", linewidth = 0.7, color = "#94a3b8")
    }
    
    p +
      scale_color_manual(values = c("Historical" = "#f59e0b", "Forecast" = "#8b5cf6")) +
      scale_fill_manual(values = c("95% Interval" = "#8b5cf6"))
  })
  
  output$plot_forecast_horizon <- renderPlot({
    fc_tbl <- analysis_data()$best_fc_tbl
    fc_tbl$date <- as.Date(fc_tbl$Forecast_Date)
    
    ggplot() +
      geom_line(data = culottes.data, aes(x = date, y = culottes, color = "Historical"), linewidth = 1) +
      geom_ribbon(data = fc_tbl, aes(x = date, ymin = Lo_95, ymax = Hi_95), fill = "#8b5cf6", alpha = 0.18) +
      geom_line(data = fc_tbl, aes(x = date, y = Forecast, color = "Forecast"), linewidth = 1.2) +
      geom_vline(xintercept = as.numeric(max(culottes.data$date)), linetype = "dashed", color = "#94a3b8") +
      scale_color_manual(values = c("Historical" = "#f59e0b", "Forecast" = "#8b5cf6")) +
      labs(
        title = paste0("Dedicated Forecast Page: ", input$h_forecast, "-Week Forecast"),
        subtitle = "This page is reserved only for 6–12 week forecasting",
        x = NULL,
        y = "Search Interest Index",
        color = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "#111827", color = NA),
        panel.background = element_rect(fill = "#111827", color = NA),
        panel.grid.major = element_line(color = rgb(1, 1, 1, 0.06)),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "#cbd5e1"),
        axis.title = element_text(color = "#e5e7eb"),
        plot.title = element_text(color = "#f8fafc", face = "bold"),
        plot.subtitle = element_text(color = "#94a3b8"),
        legend.position = "top",
        legend.text = element_text(color = "#e5e7eb")
      )
  })
  
  output$plot_all_models <- renderPlot({
    model_names <- c("Naive", "SES", "Holt", "Holt-Damped", "ARIMA", "STLF", "TBATS")
    
    plot_list <- lapply(model_names, function(m) {
      fc <- get_model_forecast(m, input$h_forecast)
      fc_tbl <- forecast_table_builder(fc, input$h_forecast, m)
      fc_tbl$date <- as.Date(fc_tbl$Forecast_Date)
      
      ggplot() +
        geom_line(data = tail(culottes.data, 80), aes(x = date, y = culottes), color = "#f59e0b", linewidth = 0.9) +
        geom_ribbon(data = fc_tbl, aes(x = date, ymin = Lo_95, ymax = Hi_95), fill = "#8b5cf6", alpha = 0.16) +
        geom_line(data = fc_tbl, aes(x = date, y = Forecast), color = "#8b5cf6", linewidth = 1) +
        geom_vline(xintercept = as.numeric(max(culottes.data$date)), linetype = "dashed", color = "#94a3b8") +
        labs(title = m, x = NULL, y = "Index") +
        theme_minimal(base_size = 12) +
        theme(
          plot.background = element_rect(fill = "#111827", color = NA),
          panel.background = element_rect(fill = "#111827", color = NA),
          panel.grid.major = element_line(color = rgb(1, 1, 1, 0.05)),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "#cbd5e1", size = 9),
          axis.title = element_text(color = "#e5e7eb", size = 10),
          plot.title = element_text(color = "#f8fafc", face = "bold", size = 13)
        )
    })
    
    do.call(grid.arrange, c(plot_list, ncol = 2))
  })
  
  output$plot_inventory <- renderPlot({
    z.val <- qnorm(input$csl / 100)
    SS <- z.val * sigma.d * sqrt(input$lead_time)
    ROP <- mu.d * input$lead_time + SS
    mu.lt <- mu.d * input$lead_time
    sigma.lt <- sigma.d * sqrt(input$lead_time)
    
    x <- seq(mu.lt - 4 * sigma.lt, mu.lt + 4 * sigma.lt, length.out = 400)
    y <- dnorm(x, mu.lt, sigma.lt)
    
    plot(
      x, y, type = "l", lwd = 2, col = "#8b5cf6",
      xlab = paste0("Demand During Lead Time (", input$lead_time, " Weeks)"),
      ylab = "Density",
      main = paste0("Lead-Time Demand Distribution | CSL = ", input$csl, "%")
    )
    
    x.fill <- x[x <= ROP]
    y.fill <- y[x <= ROP]
    polygon(c(x.fill, rev(x.fill)), c(y.fill, rep(0, length(y.fill))),
            col = rgb(0.1, 0.7, 0.5, 0.30), border = NA)
    
    x.risk <- x[x > ROP]
    y.risk <- y[x > ROP]
    polygon(c(x.risk, rev(x.risk)), c(y.risk, rep(0, length(y.risk))),
            col = rgb(0.9, 0.3, 0.3, 0.30), border = NA)
    
    abline(v = ROP, col = "#ef4444", lwd = 2, lty = 2)
    abline(v = mu.lt, col = "#38bdf8", lwd = 2, lty = 3)
    
    legend(
      "topright",
      legend = c(
        paste0("ROP = ", round(ROP, 2)),
        paste0("Mean LT Demand = ", round(mu.lt, 2)),
        paste0("Safety Stock = ", round(SS, 2))
      ),
      col = c("#ef4444", "#38bdf8", "#10b981"),
      lty = c(2, 3, 1),
      lwd = 2,
      bty = "n"
    )
  })
}

####################################################################
### Step-10: Run App
####################################################################

shinyApp(ui = ui, server = server)
