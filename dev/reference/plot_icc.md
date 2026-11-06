# Item Characteristic Curves with Class Intervals

Plots Item Characteristic Curves (ICCs) for Bayesian Rasch-family models
fitted with brms. Each item panel shows the model-expected item score
curve (with a credible interval ribbon) overlaid with observed average
item scores computed within class intervals. Optionally, observed scores
can be split by a grouping variable to visually assess differential item
functioning (DIF).

## Usage

``` r
plot_icc(
  model,
  item_var = item,
  person_var = id,
  items = NULL,
  class_intervals = 5,
  theta_range = c(-4, 4),
  n_points = 200,
  center = TRUE,
  prob = 0.95,
  ncol = NULL,
  line_size = 0.8,
  ribbon_alpha = 0.3,
  point_size = 2.5,
  dif_var = NULL,
  dif_data = NULL,
  dif_labels = NULL,
  dif_stats = TRUE,
  min_n = 5,
  palette = NULL,
  method = c("quantile", "width", "score", "manual"),
  score_breaks = NULL,
  n_intervals = NULL
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object from an ordinal IRT model (e.g., `family = acat`) or a
  dichotomous model (`family = bernoulli()`).

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- items:

  An optional character vector of item names to plot. If `NULL` (the
  default), all items are plotted.

- class_intervals:

  Integer \>= 2. The number of class intervals for `method = "quantile"`
  and `method = "width"`. Default is 5. Ignored by `method = "score"`
  and `method = "manual"`.

- theta_range:

  A numeric vector of length 2 specifying the range of theta for the
  expected curve. Default is `c(-4, 4)`.

- n_points:

  Integer. Number of evenly spaced theta values for computing the
  expected curve. Default is 200.

- center:

  Logical. If `TRUE` (the default), the scale is recentered so the grand
  mean of item thresholds = 0, consistent with
  [`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).

- prob:

  Numeric in \\(0, 1)\\. Width of the credible interval ribbon around
  the expected curve. Default is 0.95.

- ncol:

  Integer. Number of columns in the faceted layout. If `NULL`, chosen
  automatically.

- line_size:

  Numeric. Line width for the expected curve. Default is 0.8.

- ribbon_alpha:

  Numeric in \\\[0, 1\]\\. Transparency of the credible interval ribbon.
  Default is 0.3.

- point_size:

  Numeric. Size of observed score points. Default is 2.5.

- dif_var:

  An optional unquoted variable name for a grouping variable to assess
  DIF visually. If supplied, observed scores are computed separately per
  group and coloured accordingly.

- dif_data:

  An optional data frame containing the DIF variable. Required when
  `dif_var` is specified and the variable was not part of the model
  formula (since brms drops unused variables from `model$data`). Must
  have the same rows and row order as the original model data.

- dif_labels:

  An optional character vector of labels for the DIF groups. If `NULL`,
  factor levels are used.

- dif_stats:

  Logical. If `TRUE` (the default when `dif_var` is specified), the
  partial gamma coefficient for each item is annotated in the plot
  panel. The partial gamma measures the strength of DIF, stratified by
  total score.

- min_n:

  Integer. Minimum number of persons required in a class interval for
  the observed mean to be plotted. Intervals with fewer persons are
  dropped to avoid unstable estimates. Default is 5.

- palette:

  An optional character vector of colors. If `NULL`, viridis colors are
  used.

- method:

  How persons are grouped along the sum score (see *Class intervals* in
  Details). `"quantile"` (the default) forms `class_intervals` groups of
  about equal numbers of persons. `"width"` forms `class_intervals`
  intervals of equal width over the observed sum-score range. `"score"`
  makes every observed sum score its own group. `"manual"` uses the
  groups given by `score_breaks`.

- score_breaks:

  Integer vector for `method = "manual"`: the sum scores at which a new
  group starts, in increasing order. Extreme scores are left out, so the
  first group starts at 1 and the last ends at the maximum score
  minus 1. For example, `score_breaks = c(3, 6)` with a maximum score of
  9 gives the groups 1-2, 3-5 and 6-8. Default `NULL`.

- n_intervals:

  Deprecated. Use `class_intervals`.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

This is the Bayesian analogue of `ICCplot()` from the iarm package,
using the class interval method.

**Expected curve:** For each item, the expected item score
\\E_i(\theta)\\ is computed for a dense grid of theta values using the
posterior draws of item thresholds. For the adjacent category (PCM/acat)
family: \$\$E_i(\theta) = \sum\_{c=0}^{K} c \cdot P\_{ic}(\theta)\$\$
where \\P\_{ic}\\ are the category probabilities. For dichotomous items,
\\E_i(\theta) = P_i(\theta)\\.

The posterior mean of \\E_i(\theta)\\ is plotted as a line, with a
credible interval ribbon showing the `prob` interval across posterior
draws.

**Class intervals:** Following the approach in
[`iarm::ICCplot()`](https://rdrr.io/pkg/iarm/man/ICCplot.html), persons
are binned into class intervals by their ordinal *sum score* (the
sufficient statistic in Rasch models), not by their EAP theta estimate.
The grouping rules are the same as in `RMitemICCPlot()` in the
easyRasch2 package, except that extreme sum scores are left out here.
With `method = "quantile"`, tied sum scores at the quantile boundaries
can give fewer groups than requested. Equal-width and manual grouping
can define groups that no person falls into. If fewer than two groups
can be formed, the quantile and width methods fall back to one group per
sum score. The plot caption states the grouping that was used. Within
each interval, the mean observed item response is computed and plotted
at the theta value corresponding to the mean sum score in that interval,
obtained by inverting the posterior-mean expected total score function.
Persons with extreme sum scores (0 and maximum) are excluded from the
class intervals, as their theta positions cannot be estimated. Sum
scores are only comparable across persons who answered the same items,
so, as in iarm, persons with missing responses are left out of the class
intervals and the partial gamma, with a message giving their number.

**Uncertainty around observed points:** For each class interval, the
standard error of the mean observed response is computed and displayed
as error bars showing \\\pm 1.96\\ SE. These reflect sampling
variability of the observed mean within each bin.

**DIF overlay:** When `dif_var` is specified, observed scores are
computed separately for each level of the DIF variable, producing
group-specific points connected by lines. Deviations between groups that
track alongside each other (but offset from the expected curve) suggest
uniform DIF. Crossing group lines suggest non-uniform DIF.

**Partial gamma:** When `dif_stats = TRUE` (default when DIF is
requested), the Goodman-Kruskal partial gamma coefficient is displayed
in each item panel. This measures the ordinal association between item
response and DIF group, stratified by total score. Values near 0
indicate no DIF; values near \\\pm 1\\ indicate strong DIF. Positive
values mean higher item scores, at the same total score, in the group
whose label sorts second (as
[`factor()`](https://rdrr.io/r/base/factor.html) orders the levels). The
computation reuses the concordant/discordant pair counting algorithm
from
[`gk_gamma`](https://pgmj.github.io/easyRaschBayes/dev/reference/gk_gamma.md).

## See also

[`plot_ipf`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md)
for category probability curves,
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md)
for person-item maps,
[`dif_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/dif_statistic.md)
for formal Bayesian DIF testing,
[`gk_gamma`](https://pgmj.github.io/easyRaschBayes/dev/reference/gk_gamma.md)
for the underlying gamma algorithm.

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# --- Partial Credit Model ---

df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data   = df_pcm,
  family = acat,
  chains = 4, cores = 2, iter = 1000 # use more iter and cores
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Basic ICC plot with 5 class intervals
plot_icc(fit_pcm)
#> Error: object 'fit_pcm' not found

# Select specific items, one group per sum score
plot_icc(fit_pcm, items = c("I1", "I2"), method = "score")
#> Error: object 'fit_pcm' not found

# Groups defined by hand: sum scores 1-3, 4 and 5-7
plot_icc(fit_pcm, method = "manual", score_breaks = c(4, 5))
#> Error: object 'fit_pcm' not found

# With DIF overlay and partial gamma annotation

# same dataset, adding a `gender` DIF variable
df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  mutate(gender = sample(c("M", "F"), nrow(.), TRUE)) %>% 
  pivot_longer(!c(id,gender), names_to = "item", values_to = "response")
  
plot_icc(fit_pcm, dif_var = gender, dif_data = df_pcm,
         dif_labels = c("Female", "Male"), class_intervals = 5)
#> Error: object 'fit_pcm' not found
# }
```
