## ---------------------------------------------------------------------------
## Card (1995) - This script downloads the ORIGINAL data set from David Card's website
## (https://davidcard.berkeley.edu/data_sets.html)
## ---------------------------------------------------------------------------

## Opening Data --------------------------------------------------------------
library(readr)
library(car)

card <- read_csv(
  "your folder address/card1995.csv"
)

head(card)
summary(card)

## OLS wage equation ---------------------------------------------------------
ols_formula <- lwage ~ educ + exper + black + south + married + smsa
ols_fit <- lm(ols_formula, data = card)

print(summary(ols_fit))

## Main IV (2SLS) specification ----------------------------------------------
## Same wage equation, but educ is instrumented with nearc4 (grew up near a
## 4-year college)
iv_formula <- lwage ~ educ + exper + black + south + married + smsa |
  nearc4 + exper + black + south + married + smsa
iv_fit <- ivreg(iv_formula, data = card)

print(summary(iv_fit))

pe_vi <- lm(educ ~ nearc4 + exper + black + south + married + smsa,
            data = card)
linearHypothesis(pe_vi, c("nearc4 = 0"))
##------------------------------------------------------------------------------