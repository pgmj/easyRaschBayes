# Next release (after 0.4.0): small fixes

## Count labels in plot_bars() and plot_stackedbars() (noted 2026-11-09)

`plot_bars()` (R/data_plots.R, geom_text after geom_col):

1. Labels overlap when panels are narrow. "n = 217" has a fixed size at a
   fixed height (`text_y = 6`). With `ncol >= 2`, the left item-label strip
   takes about half of each panel, the bars get narrower than a three-digit
   label, and neighbouring labels overlap. Seen with phq9 (600 persons) at
   7 x 7 inches, both 9 items in 2-3 columns and 4 items in 2 columns.
2. Poor contrast. The text uses reversed magma (`viridis_option = "A"`,
   `direction = -1`) on viridis "D" bars, so pink/orange text sits on dark
   purple/blue. Bars below 6% leave the light text on the white background
   (q4 "n = 29" nearly invisible).
3. `viridis_option` and `viridis_end` change only the text colour. The fill
   uses `scale_fill_viridis_d()` with defaults.

`plot_stackedbars()`: orange labels are hard to read on the yellow segments.

Possible fixes (behaviour changes, so not in 0.4.0):

- labels just above each bar in dark grey, without "n = ", or
- labels inside the bars in black or white chosen by fill lightness, with
  `check_overlap = TRUE`
- one palette for fill and text, so `viridis_option` affects the whole plot
