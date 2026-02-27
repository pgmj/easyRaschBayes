# easyRaschBayes

<!-- badges: start -->
<!-- badges: end -->

**easyRaschBayes** provides functions to reproduce classic Rasch analysis
features using Bayesian item response theory (IRT) models fitted with
[brms](https://paulbuerkner.com/brms/). It supports both dichotomous
Rasch models and polytomous partial credit / rating scale models, and exposes
the full posterior distribution for all output.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("pgmj/easyRaschBayes")
```

## Main functions

| Function | Description |
|---|---|
| `dif_statistic()` | Differential Item Functioning (DIF) analysis |
| `fit_statistic_pcm()` | Posterior predictive item fit for polytomous models |
| `fit_statistic_rm()` | Posterior predictive item fit for dichotomous models |
| `infit_statistic()` | Conditional infit / outfit statistics |
| `item_restscore_statistic()` | Item–rest score associations with Goodman & Kruskal's gamma |
| `plot_residual_pca()` | Residual PCA contrast plot for dimensionality assessment |
| `q3_statistic()` | Yen's Q3 residual correlations for local dependence evaluation |
| `plot_ipf()` | Item category probability function curves |
| `plot_targeting()` | Person-item map (Wright map) |
| `RMUreliability()` | Reliability via Relative Measurement Uncertainty (RMU) |

## Usage

```r
library(brms)
library(easyRaschBayes)

# Fit a Bayesian Rasch model (dichotomous)
fit <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = my_data,
  family = bernoulli(),
  chains = 4, cores = 4
)

# Item fit
infit_statistic(fit)

# Person-item map
plot_targeting(fit)

# Local dependence
q3_statistic(fit)
```

## References

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS) with 
Bayesian Item Response Models. *Journal of Intelligence*, *8*(1). 
DOI: 10.3390/jintelligence8010005

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms and
Stan. *Journal of Statistical Software*, *100*, 1–54.
<https://doi.org/10.18637/jss.v100.i05>

## License

GPL (>= 3)
