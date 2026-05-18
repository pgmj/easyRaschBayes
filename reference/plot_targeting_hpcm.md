# Person-Item Targeting Plot for the Hurdle Partial Credit Model

Builds two person-item targeting plots — one for each submodel of a
hurdle partial credit model fitted with the
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
custom family. Each plot is a three-panel patchwork stack with the same
anatomy as
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md):
a person histogram, an inverted item-location histogram, and a
dot-and-whisker by item. Returning two plots — one for the presence
trait \\(\theta\_{hurdle}, \delta\_{hurdle})\\ and one for the severity
trait \\(\theta\_{pcm}, \tau\_{ik})\\ — reflects the fact that the two
submodels live on distinct latent scales and should be inspected on
their own axes.

## Usage

``` r
plot_targeting_hpcm(
  model,
  item_var = item,
  person_var = id,
  robust = FALSE,
  center = TRUE,
  sort_items = c("data", "location"),
  bins = 30,
  prob = 0.95,
  palette = NULL,
  person_fill = "#0072B2",
  threshold_fill = "#D55E00",
  height_ratios = c(3, 2, 5)
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

- robust:

  Logical. If `TRUE`, the central tendency / spread markers in the
  histogram panels use median \\\pm\\ MAD instead of mean \\\pm\\ SD.
  Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), each submodel is recentered so the
  mean item parameter is zero, with the corresponding person trait
  shifted by the same constant — matching
  [`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md)
  and
  [`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md).

- sort_items:

  One of `"data"` (default) or `"location"`, controlling the ordering of
  items in the dot-and-whisker panel.

- bins:

  Integer. Number of histogram bins. Default is 30.

- prob:

  Numeric in \\(0, 1)\\. Width of the credible interval shown as
  horizontal whiskers in the dot-and-whisker panel. Default is 0.95.

- palette:

  An optional character vector of colours for the threshold-category
  dot-and-whisker scale. If `NULL`, viridis is used. Applied to both
  submodels.

- person_fill:

  Fill colour for the person histograms. Default `"#0072B2"`.

- threshold_fill:

  Fill colour for the threshold histograms. Default `"#D55E00"`.

- height_ratios:

  A numeric vector of length 3 giving the relative heights of the
  (person, threshold, dot-and-whisker) panels. Default `c(3, 2, 5)`.

## Value

A list with two elements:

- `hurdle`:

  A patchwork object stacking the \\\theta\_{hurdle}\\ histogram, the
  \\\delta\_{hurdle}\\ histogram (inverted), and the per-item hurdle
  difficulty dot-and-whisker with the credible interval given by `prob`.

- `pcm`:

  A patchwork object with the same anatomy for the partial credit
  submodel: \\\theta\_{pcm}\\ histogram, PCM threshold histogram
  (inverted), and a per-item dot-and-whisker with one row per item and
  one coloured marker per threshold within item.

## Details

**Hurdle scale.** The presence person trait is taken as
\\\theta\_{hurdle, v} = -r\_{id}(v, \texttt{hu\\Intercept})\\ (the brms
random effect on `hu` with its sign flipped, so higher values mean
greater presence). Hurdle item difficulties are \\\delta\_{hurdle, i} =
b\_{hu\\factoritem, i}\\ directly (higher values = more zeros = harder
hurdle to cross). Under this convention, \\P(Y\_{vi} \> 0) =
\mathrm{plogis}(\theta\_{hurdle, v} - \delta\_{hurdle, i})\\, and the
histograms are directly comparable on a single x-axis.

**Partial credit scale.** The severity person trait is \\\theta\_{pcm,
v} = r\_{id}(v, \texttt{Intercept})\\, and per-item thresholds are
\\\tau\_{ik} = b\_{\texttt{Intercept}\[i, k\]}\\. The middle histogram
aggregates thresholds across items and threshold indices, exactly as in
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md).

**Independent centering per submodel.** When `center = TRUE`, the hurdle
is shifted by \\\overline{\delta\_{hurdle}}\\ and the PCM by
\\\overline{\tau}\\ (a different constant per submodel). The two
resulting x-axes therefore have a shared origin interpretation (mean
item parameter = 0) but cannot be combined on a single axis — these are
different latent traits.

**Combining the two plots.** Each list element is a valid patchwork
object, so they can be combined with the usual operators:


    plots <- plot_targeting_hpcm(fit)
    patchwork::wrap_plots(plots$hurdle, plots$pcm, ncol = 2)

For a tall, single-column layout, prefer `wrap_plots(..., ncol = 1)` so
each submodel keeps its own three-panel column.

## References

Wright, B. D. & Stone, M. H. (1979). *Best Test Design*. MESA Press.

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md)
for the single-submodel version,
[`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md),
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md).
