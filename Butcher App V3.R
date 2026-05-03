# app.R — Bucher's Method Adjusted Indirect Treatment Comparison (ITC)
# Frequentist ITC when no direct head-to-head RCT exists between A and B,
# but both have been compared to a common comparator C.

library(shiny)
library(shinyalert)
library(ggplot2)
library(dplyr)
library(tibble)

# ---------- Helpers ----------

# SE from symmetric 95% CI on log-ratio scale (HR, OR, RR).
# Assumes log-normality; less reliable for very small samples (t-distributed CIs).
se_from_ci_logratio <- function(lower, upper) {
  (log(upper) - log(lower)) / 3.92
}

# SE from symmetric 95% CI on linear scale (MD, RD)
se_from_ci_linear <- function(lower, upper) {
  (upper - lower) / 3.92
}

# Flip C-vs-A ratio direction to A-vs-C: reciprocal, CI bounds swap
invert_ratio <- function(est, lower, upper) {
  list(est = 1/est, lower = 1/upper, upper = 1/lower)
}

# Flip C-vs-A linear direction to A-vs-C: negate all, bounds swap
invert_linear <- function(est, lower, upper) {
  list(est = -est, lower = -upper, upper = -lower)
}

# Bucher indirect comparison.
# effect_* must be on log scale for ratios, raw scale for MD/RD.
# rho corrects variance for shared-control multi-arm trials (typically 0.5).
bucher_indirect <- function(effect_AC, se_AC, effect_BC, se_BC, rho = 0) {
  effect_AB <- effect_AC - effect_BC
  var_AB    <- se_AC^2 + se_BC^2 - 2 * rho * se_AC * se_BC
  se_AB     <- sqrt(var_AB)
  list(effect = effect_AB, se = se_AB)
}

ci95 <- function(est, se) {
  c(lower = est - 1.96 * se, upper = est + 1.96 * se)
}

# Two-sided Z-test p-value
pval_z <- function(effect, se) {
  2 * (1 - pnorm(abs(effect / se)))
}

format_pval <- function(p) {
  if (!is.finite(p)) return("N/A")
  if (p < 0.001)     return("< 0.001")
  formatC(p, digits = 3, format = "f")
}

# ---------- Example datasets ----------
examples <- list(

  "Custom (enter your own data)" = list(
    measure   = "HR",
    trt_names = list(A = "Treatment A", B = "Treatment B", C = "Comparator C"),
    AC        = list(est = NA, lcl = NA, ucl = NA, dir = "A_C"),
    BC        = list(est = NA, lcl = NA, ucl = NA, dir = "B_C"),
    desc      = "Enter your own values in the fields below."
  ),

  "HIV: AZT vs Didanosine (HR)" = list(
    measure   = "HR",
    trt_names = list(A = "AZT", B = "Didanosine", C = "Placebo"),
    AC        = list(est = 0.55, lcl = 0.45, ucl = 0.67, dir = "A_C"),
    BC        = list(est = 0.49, lcl = 0.38, ucl = 0.63, dir = "B_C"),
    desc      = "Classic HIV example (Bucher et al., 1997). Both drugs reduce mortality vs placebo. Indirectly compares AZT vs Didanosine."
  ),

  "Stroke: Clopidogrel vs Aspirin (OR)" = list(
    measure   = "OR",
    trt_names = list(A = "Clopidogrel", B = "Aspirin", C = "Placebo"),
    AC        = list(est = 0.80, lcl = 0.72, ucl = 0.89, dir = "A_C"),
    BC        = list(est = 0.87, lcl = 0.80, ucl = 0.95, dir = "B_C"),
    desc      = "Stroke prevention. Both active drugs effective vs placebo. Checks if Clopidogrel is superior to Aspirin."
  ),

  "Depression: Drug A vs Drug B (Mean Diff)" = list(
    measure   = "MD",
    trt_names = list(A = "Drug A", B = "Drug B", C = "Placebo"),
    AC        = list(est = -3.5, lcl = -4.5, ucl = -2.5, dir = "A_C"),
    BC        = list(est = -2.0, lcl = -3.0, ucl = -1.0, dir = "B_C"),
    desc      = "HAM-D score reduction. Negative values indicate improvement. Checks if Drug A reduces scores more than Drug B."
  ),

  "Multi-arm Trial Example (Requires Rho)" = list(
    measure   = "HR",
    trt_names = list(A = "Arm A", B = "Arm B", C = "Control"),
    AC        = list(est = 0.70, lcl = 0.55, ucl = 0.89, dir = "A_C"),
    BC        = list(est = 0.75, lcl = 0.59, ucl = 0.95, dir = "B_C"),
    desc      = "NOTICE: A vs C and B vs C come from the same 3-arm trial. Enable 'Multi-arm' in Assumptions and set rho (typically 0.5) for correct variance."
  )
)

# ---------- UI ----------
ui <- navbarPage(
  title  = "Bucher's Method ITC",
  header = tagList(useShinyalert()),

  # ---- Assumptions Tab ----
  tabPanel(
    "Assumptions",
    fluidRow(
      column(
        8,
        h3("Check assumptions before running"),
        p("The app will warn if assumptions are violated but will still compute."),

        radioButtons("assump_common_comp",
                     "1. Do both comparisons use the SAME common comparator (C)?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_same_outcome",
                     "2. Are outcome definitions and follow-up windows comparable across trials?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_transitivity",
                     "3. Is transitivity (similarity of effect modifiers) plausible?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_independent",
                     "4. Are A vs C and B vs C from independent trials (no shared patients)?",
                     choices = c("Yes", "No"), selected = "Yes", inline = TRUE),
        radioButtons("assump_multiarm",
                     "5. Does the evidence come from a single multi-arm trial (A vs B vs C in one study)?",
                     choices = c("No", "Yes"), selected = "No", inline = TRUE),

        conditionalPanel(
          condition = "input.assump_multiarm == 'Yes'",
          sliderInput("rho_multiarm",
                      "Correlation (rho) between effect estimates sharing the same control arm",
                      min = 0, max = 1, value = 0.5, step = 0.05),
          helpText("For a balanced 3-arm trial, rho = 0.5 is the standard starting point (Senn, 1995). Rho > 0 for all shared-control designs.")
        ),

        br(),
        actionButton("btn_check", "Check assumptions", class = "btn-primary")
      ),
      column(
        4,
        h4("Status"),
        uiOutput("assump_status")
      )
    )
  ),

  # ---- Inputs Tab ----
  tabPanel(
    "Inputs",

    fluidRow(
      column(
        12,
        wellPanel(
          fluidRow(
            column(
              6,
              h4("Load example data"),
              selectInput("example_select", label = NULL,
                          choices = names(examples), width = "100%"),
              textOutput("example_desc")
            ),
            column(
              6,
              h4("Assumption status"),
              uiOutput("assump_status_inline")
            )
          )
        )
      )
    ),

    fluidRow(
      column(
        12,
        wellPanel(
          h4("Treatment / Comparator labels"),
          helpText("These labels appear in the results table, forest plot, and network diagram."),
          fluidRow(
            column(4, textInput("name_A", "Treatment A name", value = "Treatment A")),
            column(4, textInput("name_B", "Treatment B name", value = "Treatment B")),
            column(4, textInput("name_C", "Comparator C name", value = "Comparator C"))
          )
        )
      )
    ),

    fluidRow(
      column(
        5,
        h3("Effect measure"),
        selectInput("measure", "Choose measure:",
                    choices = c("Hazard Ratio (HR)" = "HR",
                                "Odds Ratio (OR)"   = "OR",
                                "Relative Risk (RR)" = "RR",
                                "Mean Difference (MD)" = "MD",
                                "Risk Difference (RD)" = "RD"),
                    selected = "HR"),
        helpText("HR/OR/RR: enter estimate and 95% CI on the ratio scale (all values > 0). MD/RD: enter on the original scale."),

        hr(),
        h4(uiOutput("label_AC")),
        selectInput("dir_AC", "Direction entered",
                    choices = c("A vs C" = "A_C", "C vs A" = "C_A"),
                    selected = "A_C"),
        numericInput("est_AC", "Estimate", value = NA),
        numericInput("lcl_AC", "Lower 95% CI", value = NA),
        numericInput("ucl_AC", "Upper 95% CI", value = NA),

        hr(),
        h4(uiOutput("label_BC")),
        selectInput("dir_BC", "Direction entered",
                    choices = c("B vs C" = "B_C", "C vs B" = "C_B"),
                    selected = "B_C"),
        numericInput("est_BC", "Estimate", value = NA),
        numericInput("lcl_BC", "Lower 95% CI", value = NA),
        numericInput("ucl_BC", "Upper 95% CI", value = NA),

        hr(),
        actionButton("btn_run", "Run Bucher ITC", class = "btn-success btn-lg")
      ),

      column(
        7,
        h3("Notes"),
        tags$ul(
          tags$li("All fields must be filled before running."),
          tags$li("If multi-arm is flagged in Assumptions, variance is adjusted using rho."),
          tags$li("CIs must satisfy lower < upper. For ratio measures, all values must be > 0.")
        ),

        conditionalPanel(
          condition = "input.measure == 'RD'",
          tagList(
            tags$hr(),
            tags$h4("RD (Risk Difference): baseline risk warning"),
            tags$p(HTML("Risk Difference is an <em>absolute</em> measure sensitive to baseline event rates. If control-arm rates differ between the two trials, the same treatment can appear to have very different absolute effects even when its relative effect is identical.")),
            tags$ul(
              tags$li("Low baseline risk → small room for absolute reduction → small RD."),
              tags$li("High baseline risk → more room for absolute reduction → larger RD.")
            ),
            tags$p("Confirm that control-arm event rates are similar across trials before interpreting RD-based ITC results.")
          )
        )
      )
    )
  ),

  # ---- Results Tab ----
  tabPanel(
    "Results",
    fluidRow(
      column(
        12,
        h3("Indirect comparison result"),
        tableOutput("tbl_results"),
        fluidRow(
          column(3, downloadButton("download_csv",  "Download table (.csv)",        class = "btn-default")),
          column(3, downloadButton("download_plot", "Download forest plot (.png)",  class = "btn-default"))
        ),
        hr(),
        h3("Forest plot"),
        plotOutput("plot_forest", height = "380px"),
        hr(),
        h3("Evidence network"),
        plotOutput("plot_network", height = "300px")
      )
    )
  ),

  # ---- About Tab ----
  tabPanel(
    "About",
    fluidRow(
      column(
        8,
        h3("About this app"),
        p("This application implements ", strong("Bucher's method"), " for adjusted indirect treatment comparison (ITC) — a frequentist approach for estimating the relative efficacy of treatments A and B when no head-to-head RCT exists, but both have been compared to a common comparator C."),

        h4("Statistical method"),
        tags$ul(
          tags$li("Indirect estimate (log scale for ratios): log(HRₐᴮ) = log(HRₐᶜ) − log(HRᴮᶜ)"),
          tags$li("Variance: Var(AB) = Var(AC) + Var(BC) − 2·ρ·SE(AC)·SE(BC). Under independence (ρ = 0): Var(AB) = Var(AC) + Var(BC)."),
          tags$li("SEs are derived from 95% CIs assuming a symmetric normal distribution: SE = (upper − lower) / 3.92. Less reliable for very small samples where CIs are t-derived."),
          tags$li("95% CI for the indirect estimate: estimate ± 1.96 × SE."),
          tags$li("P-value: two-sided Z-test, p = 2 × (1 − Φ(|Z|)), where Z = estimate / SE.")
        ),

        h4("Key assumptions (Bucher, 1997)"),
        tags$ol(
          tags$li(strong("Common comparator:"), " Both trials must share the same comparator C."),
          tags$li(strong("Outcome homogeneity:"), " Outcome definitions and follow-up must be comparable."),
          tags$li(strong("Transitivity:"), " Effect modifiers must be similarly distributed across the two trials."),
          tags$li(strong("Independence:"), " The two sets of evidence must not share patients (unless the multi-arm variance correction is applied).")
        ),

        h4("Limitations"),
        tags$ul(
          tags$li("This app does not test for inconsistency between direct and indirect evidence."),
          tags$li("Only a single indirect estimate is produced; network meta-analysis handles multiple comparisons."),
          tags$li("RD-based ITC is particularly sensitive to baseline risk heterogeneity across trials.")
        ),

        h4("References"),
        tags$ol(
          tags$li("Bucher HC, Guyatt GH, Griffith LE, Walter SD. The results of direct and indirect treatment comparisons in meta-analysis of randomized controlled trials.", em(" J Clin Epidemiol."), " 1997;50(6):683–691."),
          tags$li("Dias S, Welton NJ, Sutton AJ, Ades AE. NICE DSU Technical Support Document 2: A Generalised Linear Modelling Framework for Pair-wise and Network Meta-Analysis of Randomised Controlled Trials. 2011."),
          tags$li("Senn S. Cross-over Trials in Clinical Research. Chichester: Wiley; 1993. (Rho for multi-arm shared-control variance correction.)")
        )
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {

  # ---- Example loader ----
  observeEvent(input$example_select, {
    req(input$example_select)
    ex <- examples[[input$example_select]]

    updateSelectInput(session,  "measure",  selected = ex$measure)
    updateTextInput(session,    "name_A",   value    = ex$trt_names$A)
    updateTextInput(session,    "name_B",   value    = ex$trt_names$B)
    updateTextInput(session,    "name_C",   value    = ex$trt_names$C)
    updateSelectInput(session,  "dir_AC",   selected = ex$AC$dir)
    updateNumericInput(session, "est_AC",   value    = ex$AC$est)
    updateNumericInput(session, "lcl_AC",   value    = ex$AC$lcl)
    updateNumericInput(session, "ucl_AC",   value    = ex$AC$ucl)
    updateSelectInput(session,  "dir_BC",   selected = ex$BC$dir)
    updateNumericInput(session, "est_BC",   value    = ex$BC$est)
    updateNumericInput(session, "lcl_BC",   value    = ex$BC$lcl)
    updateNumericInput(session, "ucl_BC",   value    = ex$BC$ucl)

    if (grepl("Multi-arm", input$example_select)) {
      updateSliderInput(session,     "rho_multiarm",  value    = 0.5)
      updateRadioButtons(session,    "assump_multiarm", selected = "Yes")
      shinyalert("Multi-arm example loaded",
                 "The 'Multi-arm' assumption is set to Yes and rho = 0.5. Both estimates share the same control arm.",
                 type = "info")
    } else {
      updateRadioButtons(session, "assump_multiarm", selected = "No")
    }
  })

  output$example_desc <- renderText({ examples[[input$example_select]]$desc })

  # ---- Dynamic section labels using treatment names ----
  output$label_AC <- renderUI({
    paste0(
      if (nzchar(trimws(input$name_A))) trimws(input$name_A) else "A",
      " vs ",
      if (nzchar(trimws(input$name_C))) trimws(input$name_C) else "C"
    )
  })

  output$label_BC <- renderUI({
    paste0(
      if (nzchar(trimws(input$name_B))) trimws(input$name_B) else "B",
      " vs ",
      if (nzchar(trimws(input$name_C))) trimws(input$name_C) else "C"
    )
  })

  # ---- Assumptions ----
  assumptions_ok <- reactive({
    checks <- c(
      common_comp  = input$assump_common_comp == "Yes",
      same_outcome = input$assump_same_outcome == "Yes",
      transitivity = input$assump_transitivity == "Yes",
      independent  = input$assump_independent  == "Yes"
    )
    multiarm <- input$assump_multiarm == "Yes"
    list(
      ok       = all(checks),
      checks   = checks,
      multiarm = multiarm,
      rho      = if (multiarm) input$rho_multiarm else 0
    )
  })

  assump_ui_widget <- function(st) {
    failed <- names(st$checks)[!st$checks]
    if (length(failed) == 0) {
      tagList(
        tags$span(style = "color: green;", "✓ All assumptions satisfied."),
        if (st$multiarm) tags$p(tags$small("Multi-arm: rho adjustment active.")) else NULL
      )
    } else {
      tags$span(
        style = "color: darkorange;",
        paste0("⚠ ", length(failed), " assumption(s) not met: ", paste(failed, collapse = ", "))
      )
    }
  }

  output$assump_status        <- renderUI({ assump_ui_widget(assumptions_ok()) })
  output$assump_status_inline <- renderUI({ assump_ui_widget(assumptions_ok()) })

  observeEvent(input$btn_check, {
    st     <- assumptions_ok()
    failed <- names(st$checks)[!st$checks]

    if (length(failed) == 0 && !st$multiarm) {
      shinyalert(title = "Assumptions OK",
                 text  = "All key assumptions marked as Yes.",
                 type  = "success")
      return()
    }

    msg <- character(0)
    if (length(failed) > 0)
      msg <- c(msg, paste("Assumptions marked 'No':", paste(failed, collapse = ", ")))
    if (st$multiarm)
      msg <- c(msg, "Multi-arm flagged: variance will be adjusted using rho. Results are sensitive to the rho value.")

    shinyalert(
      title = "Warnings",
      text  = paste(msg, collapse = "\n\n"),
      type  = if (length(failed) > 0) "warning" else "info"
    )
  })

  # ---- Input parsing ----
  parsed_inputs <- reactive({
    measure  <- input$measure
    is_ratio <- measure %in% c("HR", "OR", "RR")

    # Silently wait for all fields to be populated
    req(input$est_AC, input$lcl_AC, input$ucl_AC,
        input$est_BC, input$lcl_BC, input$ucl_BC)

    if (is_ratio) {
      validate(
        need(isTRUE(input$est_AC > 0),           "A vs C estimate must be > 0."),
        need(isTRUE(input$lcl_AC > 0),           "A vs C lower CI must be > 0."),
        need(isTRUE(input$ucl_AC > 0),           "A vs C upper CI must be > 0."),
        need(isTRUE(input$est_BC > 0),           "B vs C estimate must be > 0."),
        need(isTRUE(input$lcl_BC > 0),           "B vs C lower CI must be > 0."),
        need(isTRUE(input$ucl_BC > 0),           "B vs C upper CI must be > 0."),
        need(isTRUE(input$est_AC < 1e6),         "A vs C estimate is implausibly large."),
        need(isTRUE(input$est_BC < 1e6),         "B vs C estimate is implausibly large."),
        need(isTRUE(input$lcl_AC < input$ucl_AC), "A vs C: lower CI must be < upper CI."),
        need(isTRUE(input$lcl_BC < input$ucl_BC), "B vs C: lower CI must be < upper CI.")
      )
    } else {
      validate(
        need(isTRUE(abs(input$est_AC) < 1e6),    "A vs C estimate is implausibly large."),
        need(isTRUE(abs(input$est_BC) < 1e6),    "B vs C estimate is implausibly large."),
        need(isTRUE(input$lcl_AC < input$ucl_AC), "A vs C: lower CI must be < upper CI."),
        need(isTRUE(input$lcl_BC < input$ucl_BC), "B vs C: lower CI must be < upper CI.")
      )
    }

    est_AC <- input$est_AC; lcl_AC <- input$lcl_AC; ucl_AC <- input$ucl_AC
    est_BC <- input$est_BC; lcl_BC <- input$lcl_BC; ucl_BC <- input$ucl_BC

    # Fix direction so all estimates are expressed as A vs C and B vs C
    if (input$dir_AC == "C_A") {
      inv    <- if (is_ratio) invert_ratio(est_AC, lcl_AC, ucl_AC) else invert_linear(est_AC, lcl_AC, ucl_AC)
      est_AC <- inv$est; lcl_AC <- inv$lower; ucl_AC <- inv$upper
    }
    if (input$dir_BC == "C_B") {
      inv    <- if (is_ratio) invert_ratio(est_BC, lcl_BC, ucl_BC) else invert_linear(est_BC, lcl_BC, ucl_BC)
      est_BC <- inv$est; lcl_BC <- inv$lower; ucl_BC <- inv$upper
    }

    # Post-inversion CI ordering check (catches e.g. swapped bounds before inversion)
    validate(
      need(isTRUE(lcl_AC < ucl_AC),
           "After direction correction, A vs C lower CI ≥ upper CI. Check your inputs."),
      need(isTRUE(lcl_BC < ucl_BC),
           "After direction correction, B vs C lower CI ≥ upper CI. Check your inputs.")
    )

    # Convert to calculation scale: log for ratios, raw for MD/RD
    if (is_ratio) {
      eff_AC <- log(est_AC); se_AC <- se_from_ci_logratio(lcl_AC, ucl_AC)
      eff_BC <- log(est_BC); se_BC <- se_from_ci_logratio(lcl_BC, ucl_BC)
    } else {
      eff_AC <- est_AC; se_AC <- se_from_ci_linear(lcl_AC, ucl_AC)
      eff_BC <- est_BC; se_BC <- se_from_ci_linear(lcl_BC, ucl_BC)
    }

    list(
      measure  = measure,
      is_ratio = is_ratio,
      AC = list(est = est_AC, lcl = lcl_AC, ucl = ucl_AC, eff = eff_AC, se = se_AC),
      BC = list(est = est_BC, lcl = lcl_BC, ucl = ucl_BC, eff = eff_BC, se = se_BC)
    )
  })

  # ---- Compute results on button click ----
  results <- eventReactive(input$btn_run, {
    st  <- assumptions_ok()
    pin <- parsed_inputs()

    # Warn (but do not block) on assumption violations
    failed <- names(st$checks)[!st$checks]
    if (length(failed) > 0) {
      shinyalert(
        title = "Assumption warning",
        text  = paste0(
          "The following assumptions are marked 'No':\n",
          paste(failed, collapse = ", "),
          "\n\nThe app will still compute, but results may be biased or invalid."
        ),
        type = "warning"
      )
    }

    if (pin$measure == "RD") {
      shinyalert(
        title = "RD caution",
        text  = "Risk Difference is sensitive to baseline event rates. Ensure control-arm rates are similar across both trials before trusting this result.",
        type  = "info"
      )
    }

    ind <- bucher_indirect(pin$AC$eff, pin$AC$se, pin$BC$eff, pin$BC$se, rho = st$rho)
    ci  <- ci95(ind$effect, ind$se)

    # Back-transform to original (ratio or linear) scale
    if (pin$is_ratio) {
      est_AB <- exp(ind$effect)
      lcl_AB <- exp(ci["lower"])
      ucl_AB <- exp(ci["upper"])
    } else {
      est_AB <- ind$effect
      lcl_AB <- ci["lower"]
      ucl_AB <- ci["upper"]
    }

    # Guard: extreme inputs can produce non-finite results
    validate(
      need(is.finite(est_AB) && is.finite(lcl_AB) && is.finite(ucl_AB),
           "Computed result is not finite. Please check that your input values are plausible.")
    )

    p_AB <- pval_z(ind$effect, ind$se)

    # Resolve treatment labels (fall back to A/B/C if blank)
    nm <- list(
      A = if (nzchar(trimws(input$name_A))) trimws(input$name_A) else "A",
      B = if (nzchar(trimws(input$name_B))) trimws(input$name_B) else "B",
      C = if (nzchar(trimws(input$name_C))) trimws(input$name_C) else "C"
    )

    df <- tibble(
      comparison = c(
        paste0(nm$A, " vs ", nm$C),
        paste0(nm$B, " vs ", nm$C),
        paste0("Indirect: ", nm$A, " vs ", nm$B)
      ),
      estimate = c(pin$AC$est, pin$BC$est, est_AB),
      lower    = c(pin$AC$lcl, pin$BC$lcl, unname(lcl_AB)),
      upper    = c(pin$AC$ucl, pin$BC$ucl, unname(ucl_AB)),
      p_value  = c(NA_real_, NA_real_, p_AB),
      type     = c("Direct", "Direct", "Indirect")
    )

    list(df = df, is_ratio = pin$is_ratio, measure = pin$measure, nm = nm)
  })

  # ---- Results table ----
  output$tbl_results <- renderTable({
    res <- results()
    req(res)

    res$df %>%
      mutate(
        estimate = round(estimate, 4),
        lower    = round(lower,    4),
        upper    = round(upper,    4),
        p_value  = vapply(p_value, function(p) if (is.na(p)) "—" else format_pval(p), character(1))
      ) %>%
      select(-type) %>%
      rename(
        Comparison   = comparison,
        Estimate     = estimate,
        `Lower 95%`  = lower,
        `Upper 95%`  = upper,
        `P-value`    = p_value
      )
  }, na = "—")

  # ---- Forest plot (reactive object reused for download) ----
  forest_plot <- reactive({
    res <- results()
    req(res)

    df <- res$df %>%
      mutate(comparison = factor(comparison, levels = rev(comparison)))

    # Blue for direct, red for indirect estimate
    make_plot <- function(xintercept) {
      ggplot(df, aes(y = comparison, x = estimate, xmin = lower, xmax = upper,
                     colour = type, shape = type)) +
        geom_vline(xintercept = xintercept, linetype = "dashed", colour = "grey60") +
        geom_errorbarh(height = 0.25) +
        geom_point(size = 3.5) +
        scale_colour_manual(values = c("Direct" = "#2C7BB6", "Indirect" = "#D7191C"),
                            name = NULL) +
        scale_shape_manual(values  = c("Direct" = 16,        "Indirect" = 18),
                           name = NULL) +
        labs(y = NULL) +
        theme_minimal(base_size = 13) +
        theme(legend.position  = "bottom",
              panel.grid.minor = element_blank())
    }

    if (res$is_ratio) {
      make_plot(1) + scale_x_log10() + labs(x = paste0(res$measure, " (log scale)"))
    } else {
      make_plot(0) + labs(x = res$measure)
    }
  })

  output$plot_forest <- renderPlot({ forest_plot() })

  # ---- Network diagram ----
  output$plot_network <- renderPlot({
    res <- results()
    req(res)
    nm <- res$nm

    par(mar = c(1, 1, 1, 1))
    plot(NA, xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")

    pts <- data.frame(
      node = c("A",   "B",   "C"),
      x    = c(2,     8,     5),
      y    = c(8,     8,     2),
      lbl  = c(nm$A, nm$B, nm$C),
      stringsAsFactors = FALSE
    )

    xy <- function(n) unlist(pts[pts$node == n, c("x", "y")])

    # Direct evidence edges (solid, blue)
    segments(xy("A")[1], xy("A")[2], xy("C")[1], xy("C")[2], lwd = 2.5, col = "#2C7BB6")
    segments(xy("B")[1], xy("B")[2], xy("C")[1], xy("C")[2], lwd = 2.5, col = "#2C7BB6")
    # Indirect path (dashed, red)
    segments(xy("A")[1], xy("A")[2], xy("B")[1], xy("B")[2], lwd = 1.5, lty = 2, col = "#D7191C")

    symbols(pts$x, pts$y, circles = rep(0.75, 3),
            inches = FALSE, add = TRUE, bg = "white", fg = "black")
    text(pts$x, pts$y, labels = pts$lbl, cex = 1.1, font = 2)

    legend("bottom",
           legend  = c("Direct evidence", "Indirect comparison"),
           lwd     = c(2.5, 1.5),
           lty     = c(1,   2),
           col     = c("#2C7BB6", "#D7191C"),
           bty     = "n",
           horiz   = TRUE)
  })

  # ---- Download handlers ----
  output$download_csv <- downloadHandler(
    filename = function() paste0("bucher_itc_", Sys.Date(), ".csv"),
    content  = function(file) {
      res <- results()
      req(res)
      out <- res$df %>%
        mutate(p_value = vapply(p_value,
                                function(p) if (is.na(p)) NA_character_ else format_pval(p),
                                character(1))) %>%
        select(-type) %>%
        rename(Comparison   = comparison,
               Estimate     = estimate,
               Lower_95pct  = lower,
               Upper_95pct  = upper,
               P_value      = p_value)
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() paste0("bucher_forest_", Sys.Date(), ".png"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = forest_plot(), width = 9, height = 5, dpi = 150)
    }
  )
}

shinyApp(ui, server)
