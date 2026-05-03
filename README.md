# Bucher's Method ITC — Shiny App

A web application for **adjusted indirect treatment comparison (ITC)** using Bucher's method. Use this when no head-to-head RCT exists between two treatments A and B, but both have been compared to a common comparator C.

---

## What it does

- Accepts effect estimates (HR, OR, RR, MD, or RD) with 95% CIs for two trials: A vs C and B vs C
- Computes the indirect estimate for A vs B with 95% CI and two-sided p-value
- Adjusts variance for shared-control multi-arm trials using a user-supplied correlation (rho)
- Produces a colour-coded forest plot and an evidence network diagram
- Warns when key methodological assumptions are violated (without blocking computation)

---

## Getting started

### Requirements

- R ≥ 4.1
- The following R packages:

```r
install.packages(c("shiny", "shinyalert", "ggplot2", "dplyr", "tibble"))
```

### Run locally

```r
shiny::runApp("app.R")
```

Or open `Butcher App V3.R` directly in RStudio and click **Run App**.

---

## How to use

1. **Assumptions tab** — Review and confirm the four key assumptions. Flag multi-arm design and set rho if applicable (typically 0.5 for a balanced 3-arm trial).
2. **Inputs tab** — Enter treatment names, choose an effect measure, and enter the estimate + 95% CI for each comparison. Use the example loader to explore pre-loaded published scenarios.
3. **Results tab** — View the indirect estimate table (with p-value), forest plot, and network diagram. Download results as CSV or PNG.
4. **About tab** — Statistical methodology, assumptions, limitations, and references.

---

## Statistical method

The indirect estimate on the log scale (for ratio measures):

```
log(HR_AB) = log(HR_AC) − log(HR_BC)
Var(AB)    = Var(AC) + Var(BC) − 2 · ρ · SE(AC) · SE(BC)
```

For linear measures (MD, RD) the same formula applies on the raw scale.

Standard errors are derived from 95% CIs assuming a symmetric normal distribution:

```
SE = (upper − lower) / 3.92
```

The 95% CI uses `estimate ± 1.96 × SE`. The p-value is a two-sided Z-test.

---

## Built-in examples

| Example | Measure |
|---|---|
| HIV: AZT vs Didanosine | HR |
| Stroke: Clopidogrel vs Aspirin | OR |
| Depression: Drug A vs Drug B | MD |
| Multi-arm trial (requires rho) | HR |

---

## Key assumptions (Bucher, 1997)

1. Both comparisons share the **same common comparator** C
2. **Outcome definitions** and follow-up windows are comparable
3. **Transitivity** — effect modifiers are similarly distributed across trials
4. **Independence** — no shared patients between the two evidence sets (unless multi-arm correction is applied)

---

## Files

| File | Purpose |
|---|---|
| `app.R` | Main Shiny app (standard deployment entry point) |
| `Butcher App V3.R` | Same app — working copy |
| `ButcherApp.bat` | Windows launcher (local use only) |

---

## References

1. Bucher HC, Guyatt GH, Griffith LE, Walter SD. The results of direct and indirect treatment comparisons in meta-analysis of randomized controlled trials. *J Clin Epidemiol.* 1997;50(6):683–691.
2. Dias S, Welton NJ, Sutton AJ, Ades AE. NICE DSU Technical Support Document 2: A Generalised Linear Modelling Framework for Pair-wise and Network Meta-Analysis of Randomised Controlled Trials. 2011.
3. Senn S. Cross-over Trials in Clinical Research. Chichester: Wiley; 1993.
