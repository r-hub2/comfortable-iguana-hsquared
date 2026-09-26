# hsquared

<!-- ONE hex only: pkgdown uses man/figures/logo.* in the page header (drmTMB pattern).
     Do not also put a hex in this title or in a hero block — that made three. -->
<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/itchyshin/hsquared/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/itchyshin/hsquared/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<p class="hs-kicker">R interface · Julia engine · a twin, not a port</p>

<p class="hs-question">How much is genetic?</p>

<p class="hs-pitch">Formula first. The engine returns heritability, variance
components, and breeding values on the routes that are
<strong>covered</strong> — experimental 0.9.0 release, not production. CRAN
availability is tracked separately.</p>

<p class="hs-cta">
<a class="hs-btn hs-btn-primary" href="https://itchyshin.github.io/hsquared/articles/hsquared.html">Get started</a>
<a class="hs-btn" href="https://itchyshin.github.io/hsquared/articles/current-limits.html">Choose a model</a>
<a class="hs-btn" href="https://itchyshin.github.io/hsquared/articles/fitting-models.html">Fit a model</a>
</p>

<!-- Plain blockquote, never GitHub alert-callout syntax. pkgdown renders
     README.md through pandoc, which passes those bracketed alert markers
     through as literal text and leaks them onto the landing page. This wording
     carries the same weight on GitHub and on the site. -->
> <span class="hs-note-eyebrow">Before you report anything</span>
>
> **Warning — experimental 0.9.0 release.** Not production. CRAN availability
> is tracked separately. There are **7** R-public covered routes; fitting needs local Julia
> and an `HSquared.jl` checkout. Julia engine-covered ≠ R-public covered. Report
> point estimates only within a covered route. No interval is nominally
> coverage-calibrated; named univariate pedigree intervals are only classified
> as directional-conservative. See [Can I fit and report this?](https://itchyshin.github.io/hsquared/articles/current-limits.html).

## Your first analysis

<ol class="hs-workflow" aria-label="Five stages of an hsquared analysis">
<li><strong>1. Get started</strong><br><a href="https://itchyshin.github.io/hsquared/articles/hsquared.html">Install or validate a first animal-model formula</a>.</li>
<li><strong>2. Choose a model</strong><br><a href="https://itchyshin.github.io/hsquared/articles/current-limits.html">Check the route scope before fitting</a>.</li>
<li><strong>3. Fit</strong><br><a href="https://itchyshin.github.io/hsquared/articles/fitting-models.html">Run the smallest honest workflow</a>.</li>
<li><strong>4. Diagnose</strong><br><a href="https://itchyshin.github.io/hsquared/articles/visualizing-models.html">Read fit diagnostics before extracting results</a>.</li>
<li><strong>5. Report</strong><br><a href="https://itchyshin.github.io/hsquared/articles/current-limits.html">Use route-scoped point estimates only</a>.</li>
</ol>

## Start here

| If you want to… | Read this |
| --- | --- |
| check a formula without installing Julia | [Getting started](https://itchyshin.github.io/hsquared/articles/hsquared.html) (`engine = "validate"`) |
| decide whether a number may go in a paper | [Can I fit and report this?](https://itchyshin.github.io/hsquared/articles/current-limits.html) |
| see the seven R-public covered routes | [Model status](https://itchyshin.github.io/hsquared/articles/model-status.html) |
| set up the Julia engine | [Installation](#installation) |
| look up a function | [Reference](https://itchyshin.github.io/hsquared/reference/index.html) |
| understand evidence history without mistaking it for a release | [Progression & evidence](https://itchyshin.github.io/hsquared/articles/progression-evidence.html) |

## Quick start — no Julia required

This four-animal pedigree is a syntax demo, not a number for a paper.

```r
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

hsquared(
  weight ~ sex + animal(1 | id, pedigree = ped),
  data = dat,
  engine = "validate"
)
```

## Fitting — requires the Julia engine

With a local Julia, `JuliaCall`, and an `HSquared.jl` checkout (see
[Installation](#installation)):

```r
fit <- hsquared(
  weight ~ sex + animal(1 | id, pedigree = ped),
  data = dat
)
summary(fit)
```

Report only what [Can I fit and report this?](https://itchyshin.github.io/hsquared/articles/current-limits.html)
allows for the route you ran.

## Installation

```r
# install.packages("remotes")
remotes::install_github("itchyshin/hsquared")
```

Fitting also needs Julia, `JuliaCall`, and a local `HSquared.jl` tree.
`HSquared` is not in the Julia General registry — do not use
`Pkg.add("HSquared")` by name.

```r
install.packages("JuliaCall")
```

Clone or point at your `HSquared.jl` checkout per the getting-started guide.

## Authors

Shinichi Nakagawa ([ORCID 0000-0002-7765-5182](https://orcid.org/0000-0002-7765-5182),
itchyshin@gmail.com), Yefeng Yang
([ORCID 0000-0002-8610-4016](https://orcid.org/0000-0002-8610-4016)), Szymon Drobniak
([ORCID 0000-0001-8101-6247](https://orcid.org/0000-0001-8101-6247)).

## Twin

`hsquared` is the applied R interface.
[HSquared.jl](https://itchyshin.github.io/HSquared.jl/) is the sparse engine.
A twin, not a port.
