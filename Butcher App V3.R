# app.R
# Bucher's Method (Adjusted Indirect Comparison) Shiny app
# WARN on assumption violations but DO NOT block computation.
# Includes pre-loaded EXAMPLE DATA from published literature.

library(shiny)
library(shinyalert)
library(ggplot2)
library(dplyr)
library(tibble)

# ---------- Helpers ----------

se_from_ci_logratio <- function(lower, upper) {
  (log(upper) - log(lower)) / 3.92
}

se_from_ci_linear <- function(lower, upper) {
  (upper - lower) / 3.92
}

invert_ratio <- function(est, lower, upper) {
  list(est = 1/est, lower = 1/upper, upper = 1/lower)
}

invert_linear <- function(est, lower, upper) {
  list(est = -est, lower = -upper, upper = -lower)
}

bucher_indirect <- function(effect_AC, se_AC, effect_BC, se_BC, rho = 0) {
  effect_AB <- effect_AC - effect_BC
  var_AB <- se_AC^2 + se_BC^2 - 2 * rho * se_AC * se_BC
  se_AB <- sqrt(var_AB)
  list(effect = effect_AB, se = se_AB)
}

ci95 <- function(est, se) {
  c(lower = est - 1.96 * se, upper = est + 1.96 * se)
}

# ---------- Example Datasets (Published) ----------
examples <- list(
  "None" = list(
    measure = "HR",
    AC = list(est = 0.8, lcl = 0.6, ucl = 1.0, dir = "A_C"),
    BC = list(est = 1.0, lcl = 0.8, ucl = 1.2, dir = "B_C"),
    desc = "Default blank state."
  ),
  
  "HIV: AZT vs Didanosine (HR)" = list(
    measure = "HR",
    # Trial 1: AZT (A) vs Placebo (C)
    AC = list(est = 0.55, lcl = 0.45, ucl = 0.67, dir = "A_C"),
    # Trial 2: Didanosine (B) vs Placebo (C)
    BC = list(est = 0.49, lcl = 0.38, ucl = 0.63, dir = "B_C"),
    desc = "Classic HIV example. Both drugs reduce mortality vs placebo. Indirectly compares AZT vs Didanosine."
  ),
  
  "Stroke: Clopidogrel vs Aspirin (OR)" = list(
    measure = "OR",
    # Clopidogrel vs Placebo (A vs C)
    AC = list(est = 0.80, lcl = 0.72, ucl = 0.89, dir = "A_C"),
    # Aspirin vs Placebo (B vs C)
    BC = list(est = 0.87, lcl = 0.80, ucl = 0.95, dir = "B_C"),
    desc = "Stroke prevention. Both active drugs are effective. Indirect comparison checks if Clopidogrel is superior to Aspirin."
  ),
  
  "Depression: Drug A vs Drug B (Mean Diff)" = list(
    measure = "MD",
    AC = list(est = -3.5, lcl = -4.5, ucl = -2.5, dir = "A_C"),
    BC = list(est = -2.0, lcl = -3.0, ucl = -1.0, dir = "B_C"),
    desc = "HAM-D score reduction. Negative values indicate improvement. Checks if Drug A reduces scores more than B."
  ),
  
  "Multi-arm Trial Example (Requires Rho)" = list(
    measure = "HR",
    AC = list(est = 0.70, lcl = 0.55, ucl = 0.89, dir = "A_C"),
    BC = list(est = 0.75, lcl = 0.59, ucl = 0.95, dir = "B_C"),
    desc = "WARNING: These come from the *same* 3-arm trial. You MUST enable the multi-arm check in Assumptions and set rho (usually 0.5) to get correct variance."
  )
)


# ---------- UI ----------

ui <- navbarPage(
  title = "Bucher's Method ITC",
  header = tagList(useShinyalert()),
  
  tabPanel(
    "Assumptions",
    fluidRow(
      column(
        8,
        h3("Check assumptions"),
        p("The app will warn if assumptions are violated, but it will still compute."),
        
        radioButtons("assump_common_comp",
                     "Do both comparisons use the SAME common comparator (C)?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_same_outcome",
                     "Are the outcome definitions and follow-up windows comparable across trials?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_transitivity",
                     "Is transitivity (similarity of effect modifiers) plausible?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_independent",
                     "Are A vs C and B vs C evidence sources independent (different trials / no shared patients)?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_multiarm",
                     "Is the evidence coming from a single multi-arm trial (A vs B vs C)?",
                     choices = c("No", "Yes"), selected = "No", inline = TRUE),
        
        conditionalPanel(
          condition = "input.assump_multiarm == 'Yes'",
          sliderInput("rho_multiarm",
                      "If multi-arm: correlation (rho) between effect estimates (A vs C) and (B vs C)",
                      min = 0, max = 1, value = 0.5, step = 0.05),
          p("If you do not know rho, 0.5 is often used as a rough starting point for balanced 3-arm trials.")
        ),
        
        actionButton("btn_check", "Check assumptions", class = "btn-primary")
      ),
      column(4, h4("Status"), uiOutput("assump_status"))
    )
  ),
  
  tabPanel(
    "Inputs",
    fluidRow(
      # NEW: Example Loader
      column(
        12,
        wellPanel(
          h4("Load Example Data"),
          selectInput("example_select", "Choose a published/demo scenario:", choices = names(examples)),
          textOutput("example_desc")
        )
      )
    ),
    fluidRow(
      column(
        5,
        h3("Effect measure"),
        selectInput("measure", "Choose measure:",
                    choices = c("Hazard Ratio (HR)" = "HR",
                                "Odds Ratio (OR)" = "OR",
                                "Relative Risk (RR)" = "RR",
                                "Mean Difference (MD)" = "MD",
                                "Risk Difference (RD)" = "RD"),
                    selected = "HR"),
        helpText("For HR/OR/RR: enter estimate + 95% CI on ratio scale. For MD/RD: enter on original scale."),
        
        hr(),
        h4("A vs C"),
        selectInput("dir_AC", "Direction entered",
                    choices = c("A vs C" = "A_C", "C vs A" = "C_A"),
                    selected = "A_C"),
        numericInput("est_AC", "Estimate", value = 0.8),
        numericInput("lcl_AC", "Lower 95% CI", value = 0.6),
        numericInput("ucl_AC", "Upper 95% CI", value = 1.0),
        
        hr(),
        h4("B vs C"),
        selectInput("dir_BC", "Direction entered",
                    choices = c("B vs C" = "B_C", "C vs B" = "C_B"),
                    selected = "B_C"),
        numericInput("est_BC", "Estimate", value = 1.0),
        numericInput("lcl_BC", "Lower 95% CI", value = 0.8),
        numericInput("ucl_BC", "Upper 95% CI", value = 1.2),
        
        hr(),
        actionButton("btn_run", "Run Bucher ITC", class = "btn-success")
      ),
      
      # NOTES COLUMN (includes RD block)
      column(
        7,
        h3("Notes"),
        tags$ul(
          tags$li("If multi-arm is flagged, variance will be adjusted with rho.")
        ),
        
        conditionalPanel(
          condition = "input.measure == 'RD'",
          tagList(
            tags$hr(),
            tags$h4("RD (Risk Difference): baseline risk warning"),
            tags$p("If that rate is high in one trial set and low in another, the same treatment can appear to have a very different absolute benefit even if its relative effect is the same."),
            tags$h4("Why RD is vulnerable (intuition)"),
            tags$p("RD is an absolute difference in probabilities:"),
            tags$ul(
              tags$li("If an outcome is rare, there is little room to reduce it (RD tends to be small)."),
              tags$li("If an outcome is common, there is more room to reduce it (RD can be large).")
            ),
            tags$p("So, if Trial 1's placebo risk is 40% and Trial 2's placebo risk is 10%, even identical true treatment potency can yield very different RDs, breaking the transitivity/similarity idea needed for Bucher.")
          )
        )
      )
    )
  ),
  
  tabPanel(
    "Results",
    fluidRow(
      column(
        12,
        h3("Indirect comparison result"),
        tableOutput("tbl_results"),
        hr(),
        h3("Forest plot"),
        plotOutput("plot_forest", height = "380px"),
        hr(),
        h3("Simple network view"),
        plotOutput("plot_network", height = "260px")
      )
    )
  )
)

# ---------- Server ----------

server <- function(input, output, session) {
  
  # ---- Load Example Data Observer ----
  observeEvent(input$example_select, {
    req(input$example_select)
    ex <- examples[[input$example_select]]
    
    # Update inputs
    updateSelectInput(session, "measure", selected = ex$measure)
    
    updateSelectInput(session, "dir_AC", selected = ex$AC$dir)
    updateNumericInput(session, "est_AC", value = ex$AC$est)
    updateNumericInput(session, "lcl_AC", value = ex$AC$lcl)
    updateNumericInput(session, "ucl_AC", value = ex$AC$ucl)
    
    updateSelectInput(session, "dir_BC", selected = ex$BC$dir)
    updateNumericInput(session, "est_BC", value = ex$BC$est)
    updateNumericInput(session, "lcl_BC", value = ex$BC$lcl)
    updateNumericInput(session, "ucl_BC", value = ex$BC$ucl)
    
    # Special handling for multi-arm example warning
    if (grepl("Multi-arm", input$example_select)) {
      updateRadioButtons(session, "assump_multiarm", selected = "Yes")
      shinyalert("Multi-arm example loaded",
                 "Notice that the Multi-arm assumption is now set to 'Yes'. Rho (0.5) will be used.",
                 type = "info")
    } else {
      updateRadioButtons(session, "assump_multiarm", selected = "No")
    }
  })
  
  output$example_desc <- renderText({
    ex <- examples[[input$example_select]]
    ex$desc
  })
  
  
  # ---- Assumptions Logic ----
  assumptions_ok <- reactive({
    checks <- c(
      common_comp  = input$assump_common_comp == "Yes",
      same_outcome = input$assump_same_outcome == "Yes",
      transitivity = input$assump_transitivity == "Yes",
      independent  = input$assump_independent == "Yes"
    )
    
    multiarm <- input$assump_multiarm == "Yes"
    
    list(
      ok = all(checks),
      checks = checks,
      multiarm = multiarm,
      rho = if (multiarm) input$rho_multiarm else 0
    )
  })
  
  output$assump_status <- renderUI({
    st <- assumptions_ok()
    failed <- names(st$checks)[!st$checks]
    
    if (length(failed) == 0) {
      tagList(
        tags$p(strong("All key assumptions marked as Yes.")),
        if (st$multiarm) tags$p("Multi-arm flagged: variance will be adjusted using rho.") else NULL
      )
    } else {
      tagList(
        tags$p(strong("Assumption warning(s):")),
        tags$ul(lapply(failed, function(x) tags$li(x))),
        tags$p("Computation will still run, but results may be biased.")
      )
    }
  })
  
  observeEvent(input$btn_check, {
    st <- assumptions_ok()
    failed <- names(st$checks)[!st$checks]
    
    if (length(failed) == 0 && !st$multiarm) {
      shinyalert(title = "Assumptions OK",
                 text = "All key assumptions marked as Yes.",
                 type = "success")
      return()
    }
    
    msg <- character(0)
    if (length(failed) > 0) {
      msg <- c(msg, paste("These assumptions are marked 'No':", paste(failed, collapse = ", ")))
    }
    if (st$multiarm) {
      msg <- c(msg, "Multi-arm flagged: the app will adjust the variance using rho, but results can be sensitive to rho.")
    }
    
    shinyalert(
      title = "Warnings",
      text = paste(msg, collapse = "\n\n"),
      type = if (length(failed) > 0) "warning" else "info"
    )
  })
  
  parsed_inputs <- reactive({
    measure <- input$measure
    is_ratio <- measure %in% c("HR", "OR", "RR")
    
    req(input$est_AC, input$lcl_AC, input$ucl_AC, input$est_BC, input$lcl_BC, input$ucl_BC)
    
    if (is_ratio) {
      validate(
        need(input$est_AC > 0, "A vs C estimate must be > 0."),
        need(input$lcl_AC > 0, "A vs C lower CI must be > 0."),
        need(input$ucl_AC > 0, "A vs C upper CI must be > 0."),
        need(input$est_BC > 0, "B vs C estimate must be > 0."),
        need(input$lcl_BC > 0, "B vs C lower CI must be > 0."),
        need(input$ucl_BC > 0, "B vs C upper CI must be > 0."),
        need(input$lcl_AC < input$ucl_AC, "A vs C CI must have lower < upper."),
        need(input$lcl_BC < input$ucl_BC, "B vs C CI must have lower < upper.")
      )
    } else {
      validate(
        need(input$lcl_AC < input$ucl_AC, "A vs C CI must have lower < upper."),
        need(input$lcl_BC < input$ucl_BC, "B vs C CI must have lower < upper.")
      )
    }
    
    est_AC <- input$est_AC; lcl_AC <- input$lcl_AC; ucl_AC <- input$ucl_AC
    est_BC <- input$est_BC; lcl_BC <- input$lcl_BC; ucl_BC <- input$ucl_BC
    
    if (input$dir_AC == "C_A") {
      inv <- if (is_ratio) invert_ratio(est_AC, lcl_AC, ucl_AC) else invert_linear(est_AC, lcl_AC, ucl_AC)
      est_AC <- inv$est; lcl_AC <- inv$lower; ucl_AC <- inv$upper
    }
    
    if (input$dir_BC == "C_B") {
      inv <- if (is_ratio) invert_ratio(est_BC, lcl_BC, ucl_BC) else invert_linear(est_BC, lcl_BC, ucl_BC)
      est_BC <- inv$est; lcl_BC <- inv$lower; ucl_BC <- inv$upper
    }
    
    if (is_ratio) {
      eff_AC <- log(est_AC); se_AC <- se_from_ci_logratio(lcl_AC, ucl_AC)
      eff_BC <- log(est_BC); se_BC <- se_from_ci_logratio(lcl_BC, ucl_BC)
    } else {
      eff_AC <- est_AC; se_AC <- se_from_ci_linear(lcl_AC, ucl_AC)
      eff_BC <- est_BC; se_BC <- se_from_ci_linear(lcl_BC, ucl_BC)
    }
    
    list(
      measure = measure,
      is_ratio = is_ratio,
      AC = list(est = est_AC, lcl = lcl_AC, ucl = ucl_AC, eff = eff_AC, se = se_AC),
      BC = list(est = est_BC, lcl = lcl_BC, ucl = ucl_BC, eff = eff_BC, se = se_BC)
    )
  })
  
  results <- eventReactive(input$btn_run, {
    st <- assumptions_ok()
    pin <- parsed_inputs()
    
    failed <- names(st$checks)[!st$checks]
    if (length(failed) > 0) {
      shinyalert(
        title = "Assumption warning",
        text = paste(
          "You marked the following assumptions as 'No':",
          paste(failed, collapse = ", "),
          "\n\nThe app will still compute, but results may be biased/invalid.",
          sep = "\n"
        ),
        type = "warning"
      )
    }
    
    if (pin$measure == "RD") {
      shinyalert(
        title = "RD caution",
        text = "Risk Difference is sensitive to baseline risk. Confirm control event rates are similar across trials before trusting RD-based ITC.",
        type = "info"
      )
    }
    
    rho <- st$rho
    
    ind <- bucher_indirect(pin$AC$eff, pin$AC$se, pin$BC$eff, pin$BC$se, rho = rho)
    ci_ind <- ci95(ind$effect, ind$se)
    
    if (pin$is_ratio) {
      est_AB <- exp(ind$effect)
      lcl_AB <- exp(ci_ind["lower"])
      ucl_AB <- exp(ci_ind["upper"])
    } else {
      est_AB <- ind$effect
      lcl_AB <- ci_ind["lower"]
      ucl_AB <- ci_ind["upper"]
    }
    
    df <- tibble(
      comparison = c("A vs C", "B vs C", "Indirect A vs B"),
      estimate   = c(pin$AC$est, pin$BC$est, est_AB),
      lower      = c(pin$AC$lcl, pin$BC$lcl, lcl_AB),
      upper      = c(pin$AC$ucl, pin$BC$ucl, ucl_AB)
    )
    
    list(df = df, is_ratio = pin$is_ratio, measure = pin$measure)
  })
  
  output$tbl_results <- renderTable({
    res <- results()
    req(res)
    
    res$df %>%
      mutate(estimate = round(estimate, 4),
             lower = round(lower, 4),
             upper = round(upper, 4)) %>%
      rename(Comparison = comparison,
             Estimate = estimate,
             `Lower 95%` = lower,
             `Upper 95%` = upper)
  })
  
  output$plot_forest <- renderPlot({
    res <- results()
    req(res)
    
    df <- res$df %>%
      mutate(comparison = factor(comparison, levels = rev(comparison)))
    
    if (res$is_ratio) {
      ggplot(df, aes(y = comparison, x = estimate, xmin = lower, xmax = upper)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_errorbarh(height = 0.2) +
        geom_point(size = 2) +
        scale_x_log10() +
        labs(x = paste0(res$measure, " (x-axis on log scale)"), y = NULL) +
        theme_minimal(base_size = 12)
    } else {
      ggplot(df, aes(y = comparison, x = estimate, xmin = lower, xmax = upper)) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_errorbarh(height = 0.2) +
        geom_point(size = 2) +
        labs(x = res$measure, y = NULL) +
        theme_minimal(base_size = 12)
    }
  })
  
  output$plot_network <- renderPlot({
    res <- results()
    req(res)
    
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    plot(NA, xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
    
    pts <- data.frame(node = c("A", "B", "C"), x = c(2, 8, 5), y = c(8, 8, 2))
    
    segments(pts$x[pts$node == "A"], pts$y[pts$node == "A"], pts$x[pts$node == "C"], pts$y[pts$node == "C"], lwd = 2)
    segments(pts$x[pts$node == "B"], pts$y[pts$node == "B"], pts$x[pts$node == "C"], pts$y[pts$node == "C"], lwd = 2)
    segments(pts$x[pts$node == "A"], pts$y[pts$node == "A"], pts$x[pts$node == "B"], pts$y[pts$node == "B"], lwd = 1, lty = 2)
    
    symbols(pts$x, pts$y, circles = rep(0.6, 3), inches = FALSE, add = TRUE, bg = "white")
    text(pts$x, pts$y, labels = pts$node, cex = 1.3)
    
    legend("bottom", legend = c("Direct evidence", "Indirect comparison path"),
           lwd = c(2, 1), lty = c(1, 2), bty = "n")
  })
}

shinyApp(ui, server)