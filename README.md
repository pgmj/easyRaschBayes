# easyRaschBayes

<!-- badges: start -->
[![CRAN Version](https://www.r-pkg.org/badges/version/easyRaschBayes)](https://cran.r-project.org/package=easyRaschBayes)
[![Downloads](https://cranlogs.r-pkg.org/badges/easyRaschBayes?color=brightgreen)](https://CRAN.R-project.org/package=easyRaschBayes)
![Downloads Status](https://cranlogs.r-pkg.org/badges/grand-total/easyRaschBayes)
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<!-- badges: end -->

**easyRaschBayes** provides functions to reproduce classic Rasch analysis
features using Bayesian item response theory (IRT) models fitted with
[brms](https://paulbuerkner.com/brms/). It supports both dichotomous
Rasch models and polytomous partial credit models, and exposes
the full posterior distribution for all output.

For more materials on Rasch analysis, see the 
[vignette](https://pgmj.github.io/raschrvignette/RaschRvign.html) for my 
(frequentist) package [`easyRasch`](https://pgmj.github.io/easyRasch/).

## Installation

Install the stable version from CRAN:

```r
install.packages("easyRaschBayes")
```

Or install the development version from GitHub:

```r
install.packages("remotes") # if you don't have `remotes` installed
remotes::install_github("pgmj/easyRaschBayes")
```

## Main functions

| Function | Description |
|---|---|
| `infit_statistic()` | Conditional infit (and outfit) statistics |
| `infit_post()` | Post-processing and plotting of infit statistics |
| `item_restscore_statistic()` | Item–rest score associations with Goodman & Kruskal's gamma |
| `item_restscore_post()` | Post-processing and plotting of item-restscore statistics |
| `plot_icc()` | Conditional item fit curves & DIF analysis |
| `plot_residual_pca()` | Residual PCA contrast plot for dimensionality assessment |
| `q3_statistic()` | Yen's Q3 residual correlations for local dependence evaluation |
| `q3_post()` | Post-processing and plotting of Q3 statistics |
| `plot_ipf()` | Item probability function curves |
| `plot_icc()` | Item category characteristics curves |
| `plot_targeting()` | Person-item map (Wright map) |
| `dif_statistic()` | Differential Item Functioning (DIF) analysis |
| `fit_statistic_pcm()` | Posterior predictive item/person fit for polytomous models |
| `fit_statistic_rm()` | Posterior predictive item/person fit for dichotomous models |
| `RMUreliability()` | Reliability via Relative Measurement Uncertainty (RMU) |
| `person_parameters()` | EAP and WLE latent scores and ordinal sum score transformation table |
| `item_parameters()` | Summarize item threshold locations |
| `plot_tile()` | Plot response data, see also `plot_bars()` and `plot_stackedbars()` |

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

# Conditional item infit
infit <- infit_statistic(fit)
infit_post(infit)

# Local dependence with Yen's Q3
q3_results <- q3_statistic(fit)
q3_post(q3_results)

# Response category functioning
plot_ipf(fit)

# Person-item map
plot_targeting(fit)
```

## References

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS) with 
Bayesian Item Response Models. *Journal of Intelligence*, *8*(1). 
<doi:10.3390/jintelligence8010005>

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms and
Stan. *Journal of Statistical Software*, *100*, 1–54.
<https://doi.org/10.18637/jss.v100.i05>

## Credits

This started as a side project to update the code to assess item fit included in Bürkner (2020).
With the help of Claude Opus 4.6 things expanded to try to 
reproduce core Rasch analysis steps. Most of the code in this package is
produced by the LLM.

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed psychologist with a PhD in behavior analysis. He works as a research specialist at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 

## License

GPL (>= 3)
