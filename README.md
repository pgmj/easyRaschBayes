# easyRaschBayes

<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/easyRaschBayes)](https://cran.r-project.org/package=easyRaschBayes)
[![Downloads](https://cranlogs.r-pkg.org/badges/easyRaschBayes?color=brightgreen)](https://CRAN.R-project.org/package=easyRaschBayes)
![Downloads Status](https://cranlogs.r-pkg.org/badges/grand-total/easyRaschBayes)
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<!-- badges: end -->

`easyRaschBayes` provides Rasch measurement theory analysis workflows for
Bayesian item response theory models fitted with
[brms](https://paulbuerkner.com/brms/). It supports dichotomous Rasch
models, partial credit and rating scale models, and a hurdle partial credit
model for data with many zero responses. All output exposes the full
posterior distribution.

The frequentist counterpart is [`easyRasch2`](https://pgmj.github.io/easyRasch2/),
which uses conditional maximum likelihood and simulation-based critical
values. Both packages include the same `phq9` example data, and their
vignettes analyse it side by side. For broader Rasch-analysis tutorials, see
the [vignette](https://pgmj.github.io/raschrvignette/RaschRvign.html) for the
archived package [`easyRasch`](https://pgmj.github.io/easyRasch/).

## Statement of need

`brms` fits Rasch-family models flexibly, with full posterior uncertainty
and with extensions such as person covariates, extra facets and custom
families (Bürkner, 2021). It does not provide the diagnostics a Rasch
analysis needs: item fit, local dependence, dimensionality, differential
item functioning, response category functioning, targeting and
reliability. `easyRaschBayes` provides these for `brms` models, with one
consistent interface, and reports each diagnostic with its posterior
uncertainty.

## Installation

Install the stable version from CRAN:

```r
install.packages("easyRaschBayes")
```

Or install the development version from GitHub:

```r
# install.packages("remotes") # if needed
remotes::install_github("pgmj/easyRaschBayes")
```

The development version changes several defaults and outputs, among them
conditional item fit given the total score and a new targeting plot. See
the [changelog](https://pgmj.github.io/easyRaschBayes/news/index.html).

The hurdle partial credit model needs a patched version of `brms`, since
CRAN `brms` 2.23 cannot fit custom ordinal families with item-specific
thresholds:

```r
remotes::install_github("rpsychologist/brms@fix-custom-ordinal-grouped-thres")
```

## Key design principles

- **Estimation**: models are fitted with `brms` and Stan. Item thresholds
  are read from the posterior draws in the same way by every function, for
  partial credit models (`thres(gr = item)`), rating scale models and
  dichotomous models with fixed or random item effects. Person locations
  are posterior means (EAP), with Warm's weighted likelihood estimates (WLE)
  as an alternative.
- **Inference**: posterior predictive checks. For Rasch-family models the
  total score is sufficient for the person location, and item fit and
  local dependence are computed conditional on it: item-restscore associations are compared with
  their expected values given the total scores (Kreiner, 2011), and infit
  follows the conditional infit of Christensen et al. (2013). Replicated
  data are drawn conditional on each person's total score. In null
  simulations the restscore and Q3 posterior predictive p-values were close
  to uniform, while infit stays conservative, so infit and restscore are
  meant to be used together. Posterior predictive p-values are not measures of the size of
  misfit.
- **Output**: tibbles for tables and `ggplot2` figures. Infit,
  item-restscore and Q3 are computed per draw by `*_statistic()` functions
  and summarised and plotted by `*_post()` functions, so the draws stay
  available for further analysis.
- **Naming**: snake_case, with the suffix `_hpcm` for functions for the
  hurdle partial credit model.

## Functions by domain

### Item fit

- `infit_statistic()` + `infit_post()`: conditional infit (and outfit) mean
  squares with posterior predictive p-values; `method = "theta"` for
  infit given each draw's person locations
- `item_restscore_statistic()` + `item_restscore_post()`: item-restscore
  association (Goodman-Kruskal's gamma) against its expected value given the
  total scores, with a max-T adjustment for the number of items
- `plot_icc()`: conditional item characteristic curves in class intervals
  of the total score
- `fit_statistic_pcm()`, `fit_statistic_rm()`: general posterior predictive
  item and person fit with user-supplied criteria

### Local dependence

- `q3_statistic()` + `q3_post()`: Yen's Q3 residual correlations given
  the total scores, against replicated data, with a max-T adjustment for
  the number of item pairs; `method = "theta"` for residuals given each
  draw's person locations

### Dimensionality

- `plot_residual_pca()`: first residual contrast, observed against
  replicated eigenvalue, with item loadings and locations

### Differential item functioning

- `dif_statistic()`: refits the model with item-specific group differences,
  estimates the difference in latent means (impact) separately, and
  summarises the posterior of each item's DIF
- `plot_icc()` with `dif_var`: conditional item characteristic curves by
  group, with partial gamma

### Response categories

- `plot_ipf()`: item category probability curves with credible bands
- `item_parameters()`: thresholds, with the posterior probability that
  each item's thresholds are disordered

### Targeting and reliability

- `plot_targeting()`: person-item targeting plot with the ranges where each
  response category is the most likely, threshold credible intervals, and
  the posterior probability of collapsed categories and disordered thresholds
- `RMUreliability()`: reliability from relative measurement uncertainty
  (Bignardi et al., 2026)

### Item and person parameters

- `item_parameters()`: item locations and thresholds with credible
  intervals
- `person_parameters()`: EAP and WLE person locations, and a sum score to
  WLE table
- `posterior_to_prior()`: informative priors from a fitted model, for use
  in a new model

### Hurdle partial credit model

- `hurdle_acat()` + `hurdle_acat_stanvars()`: custom `brms` family: a
  hurdle for any response above 0 and a partial credit model for the
  positive categories, each with its own person location
- `infit_statistic_hpcm()`, `q3_statistic_hpcm()`, `item_parameters_hpcm()`,
  `person_parameters_hpcm()`, `RMUreliability_hpcm()`, `plot_icc_hpcm()`,
  `plot_targeting_hpcm()`: the diagnostics above for each part of the model

### Data visualization

- `plot_tile()`, `plot_bars()`, `plot_stackedbars()`: response
  distributions from wide data

## Data

- `phq9`: PHQ-9 depression screener, 600 NHANES respondents selected to
  have few zero sum scores, with gender, age and education. The same data as
  in `easyRasch2`
- `phq9_nhanes`: PHQ-9, a random sample of 1000 NHANES respondents with the
  many zero responses of a general population sample, for the hurdle
  partial credit model

## Vignettes

- [Rasch partial credit model](https://pgmj.github.io/easyRaschBayes/articles/pcm-rasch-analysis.html),
  a full analysis of `phq9`
- [Hurdle partial credit model](https://pgmj.github.io/easyRaschBayes/articles/hpcm-rasch-analysis.html),
  `phq9_nhanes` analysed with an ordinary and a hurdle partial credit model

There is also a [comparison with frequentist Rasch tools in R](https://pgmj.github.io/rasch_bayes_comp.html).

## Example

```r
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(easyRaschBayes)

# Long format, one row per person and item. The acat family expects
# response categories starting at 1.
df <- phq9[, 1:9] %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response") %>%
  mutate(response = response + 1)

# Partial credit model
fit <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data   = df,
  family = acat,
  chains = 4, cores = 4
)

# Conditional item fit: item-restscore and infit
item_restscore_post(item_restscore_statistic(fit))
infit_post(infit_statistic(fit))

# Local dependence with Yen's Q3
q3_post(q3_statistic(fit))

# Response categories and targeting
plot_ipf(fit)
plot_targeting(fit)
```

## References

- Bignardi, G., Kievit, R., & Bürkner, P.-C. (2026). A general method for
  estimating reliability using Bayesian measurement uncertainty. *OSF*.
  <https://osf.io/preprints/psyarxiv/h54k8>
- Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
  with Bayesian item response models. *Journal of Intelligence, 8*(1).
  <https://doi.org/10.3390/jintelligence8010005>
- Bürkner, P.-C. (2021). Bayesian item response modeling in R with brms and
  Stan. *Journal of Statistical Software, 100*, 1–54.
  <https://doi.org/10.18637/jss.v100.i05>
- Christensen, K. B., Kreiner, S., & Mesbah, M. (Eds.) (2013). *Rasch
  models in health*. ISTE and Wiley.
- Kreiner, S. (2011). A note on item-restscore association in Rasch models.
  *Applied Psychological Measurement, 35*(7), 557–561.
  <https://doi.org/10.1177/0146621611410227>
- Magnus, B. E., & Garnier-Villarreal, M. (2022). A multidimensional
  zero-inflated graded response model for ordinal symptom data.
  *Psychological Methods, 27*(2), 261–279.
  <https://doi.org/10.1037/met0000395>
- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
  them? *Journal of Statistical Distributions and Applications, 7*(1), 5.
  <https://doi.org/10.1186/s40488-020-00108-7>

## Credits

This started as a small side project to update the code to assess item fit included in Bürkner (2020). With the help of Claude Opus things expanded to try to 
reproduce core Rasch analysis aspects with Bayesian methods. 
Most of the code in this package is produced by the LLM.

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed psychologist with a PhD in behavior analysis. He works as a research specialist at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 

## License

GPL (>= 3)
