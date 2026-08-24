# Default bin labels for an internal OSA call

`bins` and `bin_label` label the residual's bin axis, and
[`plot_resids`](https://chengmatt.github.io/SPoRC/dev/reference/plot_resids.md)
draws its second panel against them. Left empty the label columns never
reach the residual frame and that panel errors when it is printed, so
they are filled from the data here: the observed age bins for an age or
conditional age-at-length source and the length bins for a length
source. The observed bins are read off the observation array, since
ageing error can leave fewer of them than the model carries.

## Usage

``` r
osa_default_bins(data, comp_source, pop = FALSE, discard = FALSE)
```

## Arguments

- data:

  Data list from the fitted model.

- comp_source:

  Composition source name, as
  [`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
  takes it.

- pop, discard:

  Whether the source is population-specific or the discard stream, as
  [`get_osa`](https://chengmatt.github.io/SPoRC/dev/reference/get_osa.md)
  takes them.

## Value

List with `bins` and `bin_label`.
