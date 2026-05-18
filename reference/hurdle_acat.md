# Hurdle Partial Credit Model Custom brms Family

A custom brms family implementing a hurdle Rasch partial credit model
(hPCM). A Bernoulli logit gate models \\P(Y = 0)\\; conditional on \\Y
\> 0\\, an acat-logit partial credit process governs transitions among
the positive categories \\1, \ldots, K - 1\\.

## Usage

``` r
hurdle_acat()
```

## Value

A `customfamily` object suitable for the `family` argument of
[`brm`](https://paulbuerkner.com/brms/reference/brm.html). The companion
stanvars (the Stan code for the custom lpmf) must be passed via the
`stanvars` argument; see
[`hurdle_acat_stanvars`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat_stanvars.md).

## Details

The two submodels can be interpreted as a *hurdle* (presence / absence
of a symptom or behaviour, sometimes called "susceptibility") and a
*partial credit* severity model (frequency / intensity given presence),
in the spirit of Magnus and Garnier-Villarreal (2022). Person random
effects on the two submodels can be modelled as correlated (via brms's
`(1 |g| id)` syntax), allowing the data to inform whether susceptibility
and severity are distinct latent constructs or essentially the same
trait.

## Usage


    library(brms)
    fit <- brm(
      bf(
        response | thres(gr = item) ~ 1 + (1 |g| id),
        hu ~ 0 + factor(item) + (1 |g| id)
      ),
      data     = dat,
      family   = hurdle_acat(),
      stanvars = hurdle_acat_stanvars(),
      ...
    )

## brms compatibility note

The custom family relies on brms's native ordinal infrastructure
(`specials = "ordinal"` + `thres(gr = item)`). On CRAN brms 2.23.x this
combination emits invalid Stan code for custom ordinal families with
grouped thresholds. A patched branch is available:


    devtools::install_github(
      "rpsychologist/brms@fix-custom-ordinal-grouped-thres"
    )

## Sign convention

The brms formula `hu ~ ... + (1 | id)` adds the person random effect to
the logit of `hu = P(Y = 0)`, so higher values of the person random
effect on `hu` mean *more* zeros (lower susceptibility). This is the
opposite sign of a Stan-style parameterisation in which higher
`theta_gate` means fewer zeros (higher susceptibility). The two
parameterisations imply identical likelihoods; only the sign of the
recovered correlation between the two person random effects is flipped
relative to a "susceptibility x severity" labelling. No inferences about
person ranking, gate probabilities, or severity probabilities are
affected.

## References

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`hurdle_acat_stanvars`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat_stanvars.md)
for the Stan function block,
[`infit_statistic_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic_hpcm.md)
for submodel-specific item infit,
[`q3_statistic_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic_hpcm.md)
for submodel-specific Q3 residual correlations.

## Author

Kristoffer Magnusson
