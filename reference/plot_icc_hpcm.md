# Item Characteristic Curves with Class Intervals for Hurdle PCM

Plots Item Characteristic Curves (ICCs) with class-interval overlays
separately for each submodel of a hurdle partial credit model fitted
with the
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
custom family. Returns two ggplot2 plots — one for the Bernoulli hurdle
on \\P(Y \> 0)\\, one for the conditional partial credit \\E(Y \mid Y \>
0)\\ — that can be inspected separately or combined with patchwork.
Conceptually the Bayesian / hurdle-PCM analogue of
[`plot_icc`](https://pgmj.github.io/easyRaschBayes/reference/plot_icc.md),
which only handles single-submodel ordinal / dichotomous IRT.

## Usage

``` r
plot_icc_hpcm(
  model,
  item_var = item,
  person_var = id,
  items = NULL,
  n_intervals = 5,
  theta_range = c(-4, 4),
  n_points = 200,
  center = TRUE,
  prob = 0.95,
  ncol = NULL,
  line_size = 0.8,
  ribbon_alpha = 0.3,
  point_size = 2.5,
  min_n = 5
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object using the
  [`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
  custom family.

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- items:

  An optional character vector of item names to plot. If `NULL` (the
  default), all items are plotted.

- n_intervals:

  Integer. The number of class intervals into which persons are binned
  within each submodel. Default is 5.

- theta_range:

  A numeric vector of length 2 specifying the range of theta for the
  expected curves. Default is `c(-4, 4)`.

- n_points:

  Integer. Number of evenly spaced theta values for computing the
  expected curves. Default is 200.

- center:

  Logical. If `TRUE` (the default), each submodel is recentered so that
  the mean item parameter is zero, matching
  [`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md)
  and
  [`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md).

- prob:

  Numeric in \\(0, 1)\\. Width of the credible interval ribbon around
  each expected curve. Default is 0.95.

- ncol:

  Integer. Number of columns in the faceted layout. If `NULL`, chosen
  automatically.

- line_size:

  Numeric. Line width for the expected curves. Default is 0.8.

- ribbon_alpha:

  Numeric in \\\[0, 1\]\\. Transparency of the credible interval
  ribbons. Default is 0.3.

- point_size:

  Numeric. Size of observed score points. Default is 2.5.

- min_n:

  Integer. Minimum number of observations required in a class interval
  (per item, per submodel) for the observed mean to be plotted. Default
  is 5.

## Value

A list with two elements:

- `hurdle`:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing per-item Bernoulli ICCs on the \\P(Y \> 0)\\ scale,
  with class-interval observed empirical hurdle-crossing rates overlaid.

- `pcm`:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing per-item conditional partial credit ICCs on the \\E(Y
  \mid Y \> 0)\\ scale, with class-interval observed mean severities
  overlaid (restricted to persons with \\Y \> 0\\ on that item).

## Details

**Hurdle submodel.** For each item the expected curve is \\P(Y\_{vi} \>
0 \mid \theta\_{hurdle}) = \mathrm{plogis} (\theta\_{hurdle} -
\delta\_{hurdle, i})\\, computed across posterior draws of the hurdle
item difficulty. Persons are binned into class intervals by their hurdle
sum score (count of items with \\Y \> 0\\); each interval is positioned
on the \\\theta\_{hurdle}\\ axis by inverting the posterior-mean
expected total \\\sum_i P(Y_i \> 0 \mid \theta\_{hurdle})\\, exactly
analogous to the procedure in
[`plot_icc`](https://pgmj.github.io/easyRaschBayes/reference/plot_icc.md).
The y-axis ranges from 0 to 1.

**Partial credit submodel.** For each item the expected curve is the
conditional expected severity \\E(Y\_{vi} \mid Y\_{vi} \> 0,
\theta\_{pcm}) = \sum\_{c = 1}^{K - 1} c \cdot P(Y\_{vi} = c \mid
Y\_{vi} \> 0)\\, with conditional PCM category probabilities computed
from the posterior draws of \\\tau\_{ik}\\. The y-axis ranges from 1 to
\\K - 1\\.

Persons are binned into class intervals by their **EAP**
\\\theta\_{pcm}\\ (taken from
[`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html)), not by a
sum score. This is intentional: under the hurdle PCM the sum of \\Y\\
over items with \\Y\_{vi} \> 0\\ is not a sufficient statistic for
\\\theta\_{pcm}\\ (because the number of contributing items varies
across persons; see
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)).
EAP-based binning sidesteps the inverse-CDF construction and avoids
privileging any particular summary score. Within each interval, the
observed mean response for an item is computed *only* over persons with
\\Y \> 0\\ on that item — i.e., the cells the PCM submodel actually
applies to.

**Uncertainty.** Each expected curve carries a credible interval ribbon
at width `prob`. Observed points have \\\pm 1.96\\SE\\ error bars where
SE is the within-bin sampling standard error of the mean (binomial for
the hurdle, ordinal for the partial credit submodel).

**Why two separate plots and not patchwork-combined.** The two submodels
have qualitatively different y-axes (probability vs. expected severity)
and different x-axis interpretations (presence trait vs. severity
trait). Returning them as a list preserves the option to inspect each in
isolation; combine them yourself with patchwork when a side-by-side or
stacked layout is wanted:


    plots <- plot_icc_hpcm(fit)
    patchwork::wrap_plots(plots$hurdle, plots$pcm, ncol = 1)

## References

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`plot_icc`](https://pgmj.github.io/easyRaschBayes/reference/plot_icc.md)
for the single-submodel ICC plot,
[`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md),
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md).
