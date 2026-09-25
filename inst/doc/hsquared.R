## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----f1-animal-path, echo = FALSE, fig.alt = "Flowchart of the default univariate animal model: pedigree and phenotype enter the hsquared formula, Julia REML returns variance components, heritability, and breeding values. Genomic and multivariate paths are not shown.", out.width = "100%"----
knitr::include_graphics("../man/figures/animal-model-path.svg")

## ----golden-data--------------------------------------------------------------
library(hsquared)

ped <- data.frame(
  id   = c("sire", "dam", "off1", "off2"),
  sire = c(NA, NA, "sire", "sire"),
  dam  = c(NA, NA, "dam", "dam")
)
dat <- data.frame(
  id     = c("sire", "dam", "off1", "off2"),
  sex    = c("m", "f", "m", "f"),
  weight = c(42, 38, 40, 37)
)

## ----golden-validate----------------------------------------------------------
hsquared(
  weight ~ sex + animal(1 | id, pedigree = ped),
  data = dat,
  control = hs_control(engine = "validate")
)

model_spec(
  weight ~ sex + animal(1 | id, pedigree = ped),
  data = dat
)

## ----golden-fit-path, eval = FALSE--------------------------------------------
# fit <- hsquared(
#   weight ~ sex + animal(1 | id, pedigree = ped),
#   data = dat,
#   control = hs_control(
#     engine_control = list(julia_project = "/path/to/HSquared.jl")
#   )
# )

## ----golden-fit, eval = FALSE-------------------------------------------------
# fit <- hsquared(
#   weight ~ sex + animal(1 | id, pedigree = ped),
#   data = dat
# )
# 
# fit_diagnostics(fit)
# heritability(fit)          # warns if the fit did not converge
# variance_components(fit)
# breeding_values(fit)

