# Create an MD plot from expression data using base R

This function differs from
[`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
as it produces a plot in base R. Points are colored by adjusted p-value
(and logFC) significance, and the percentage and count of up/down
features are annotated on the plot area.

## Usage

``` r
plot_md2(
  df,
  x = "logCPM",
  y = "logFC",
  sig_col = "FDR",
  lab = NULL,
  fdr = 0.1,
  lfc = 0,
  label_sig = FALSE,
  n_labels = 10,
  annotate_counts = TRUE,
  add_grid = FALSE,
  xmax_label_offset = 0.85,
  ymax_label_offset = 0.8,
  ymin_label_offset = 0.8,
  lab_size = 1,
  lab_digits = 1,
  up_color = "red2",
  down_color = "royalblue2",
  nonde_color = "grey40",
  up_alpha = 1,
  down_alpha = 1,
  nonde_alpha = 1,
  up_size = 1,
  down_size = 1,
  nonde_size = 1,
  up_shape = 16,
  down_shape = 16,
  nonde_shape = ".",
  x_axis_limits = NULL,
  y_axis_limits = NULL,
  ...
)
```

## Arguments

- df:

  data.frame containing columns for average abundance, log-fold-change
  and adjusted p-values

- x:

  column in dataframe containing the logCPM data. Default "logCPM"

- y:

  column in dataframe containing the log-fold-change values. Default
  "logFC"

- sig_col:

  column in dataframe containing the adjusted p-values used to define
  significance. Default "FDR"

- lab:

  column in dataframe containing the labels to plot if label_sig = TRUE.
  Default NULL

- fdr:

  Significance threshold for features. Based on sig_col. Default 0.1

- lfc:

  logFC threshold for considering significant features. Features must
  have abs(logFC) greater than lfc. Default 0.

- label_sig:

  TRUE/FALSE. Label significant features with the values in `lab` using
  [`text()`](https://rdrr.io/r/graphics/text.html).

- n_labels:

  integer. Maximum number of features to label when label_sig = TRUE.
  Since base R cannot repel labels, only the `n_labels` significant
  features with the most extreme combined rank of adjusted p-value and
  absolute logFC are labeled. Default 10

- annotate_counts:

  TRUE/FALSE. Annotate the plot area with the percentage and count of
  up- and down-regulated features. Default TRUE

- add_grid:

  TRUE/FALSE. Draw a background grid with
  [`grid()`](https://rdrr.io/r/graphics/grid.html) beneath the points.
  Default FALSE

- xmax_label_offset:

  numeric. Value between 0 and 1 inclusive. Controls the x-position of
  the count labels as a fraction of the right x-axis limit

- ymax_label_offset:

  numeric. Value between 0 and 1 inclusive. Controls the y-position of
  the 'up' count label as a fraction of the upper y-axis limit

- ymin_label_offset:

  numeric. Value between 0 and 1 inclusive. Controls the y-position of
  the 'down' count label as a fraction of the lower y-axis limit

- lab_size:

  numeric. Size (cex) of the count labels when annotate_counts = TRUE.
  Default 1

- lab_digits:

  numeric. The number of digits used when rounding percentage values
  when annotate_counts = TRUE. Default 1

- up_color:

  Point color of the up-regulated features

- down_color:

  Point color of the down-regulated features

- nonde_color:

  Point color of the unperturbed features

- up_alpha:

  Point alpha of the up-regulated features

- down_alpha:

  Point alpha of the down-regulated features

- nonde_alpha:

  Point alpha of the unperturbed features

- up_size:

  Point size (cex) of the up-regulated features

- down_size:

  Point size (cex) of the down-regulated features

- nonde_size:

  Point size (cex) of the unperturbed features

- up_shape:

  Point shape (pch) of the up-regulated features

- down_shape:

  Point shape (pch) of the down-regulated features

- nonde_shape:

  Point shape (pch) of the unperturbed features

- x_axis_limits:

  numeric vector of length 2 giving the x-axis limits. Default (NULL)

- y_axis_limits:

  numeric vector of length 2 giving the y-axis limits. Default (NULL)

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). These
  override the defaults (e.g. `xlab`, `ylab`, `las`, `bty`, `main`)

## Value

Invisibly, a list with the plotted data and group counts

## Examples

``` r
plot_md2(GSE161650_de, fdr = 0.01, lfc = log2(2))
```
