# easyRaschBayes

**easyRaschBayes** provides functions to reproduce classic Rasch
analysis features using Bayesian item response theory (IRT) models
fitted with [brms](https://paulbuerkner.com/brms/). It supports both
dichotomous Rasch models and polytomous partial credit models, and
exposes the full posterior distribution for all output.

## Installation

Install the stable version from CRAN:

``` r
install.packages("easyRaschBayes")
```

Or install the development version from GitHub:

``` r
install.packages("remotes") # if you don't have `remotes` installed
remotes::install_github("pgmj/easyRaschBayes")
```

## Main functions

| Function                                                                                                    | Description                                                    |
|-------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------|
| [`dif_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/dif_statistic.md)                       | Differential Item Functioning (DIF) analysis                   |
| [`fit_statistic_pcm()`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_pcm.md)               | Posterior predictive item fit for polytomous models            |
| [`fit_statistic_rm()`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_rm.md)                 | Posterior predictive item fit for dichotomous models           |
| [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)                   | Conditional infit / outfit statistics                          |
| [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md) | Item–rest score associations with Goodman & Kruskal’s gamma    |
| [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/reference/plot_residual_pca.md)               | Residual PCA contrast plot for dimensionality assessment       |
| [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)                         | Yen’s Q3 residual correlations for local dependence evaluation |
| [`plot_ipf()`](https://pgmj.github.io/easyRaschBayes/reference/plot_ipf.md)                                 | Item category probability function curves                      |
| [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md)                     | Person-item map (Wright map)                                   |
| [`RMUreliability()`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md)                     | Reliability via Relative Measurement Uncertainty (RMU)         |

## Usage

``` r
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

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
with Bayesian Item Response Models. *Journal of Intelligence*, *8*(1).
<doi:10.3390/jintelligence8010005>

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
<https://doi.org/10.18637/jss.v100.i05>

## Credits

This started as a little side project to update the code posted by
Bürkner (2020) to assess item fit. With the help of Claude Opus 4.6 the
things expanded to try to reproduce the core Rasch analysis aspects.
Most of the code in this package is produced by the LLM.

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a
licensed psychologist with a PhD in behavior analysis. He works as a
research specialist at [Karolinska
Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research),
Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky:
  [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social)
- Mastodon: [@pgmj@scicomm.xyz](https://scicomm.xyz/@pgmj)

## License

GPL (\>= 3)
