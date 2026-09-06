# Number of observed bins a composition data source is recorded on

The `*_bins` arguments index into observed bins, so they need the bin
count of the array they will be applied to. That is normally read
straight off the supplied observation array, but a model with no data
for a data source can hand in an array with no dimensions at all, so the
model's own observed bin count stands in: the ageing error's
observed-age dimension for ages, and
[`obs_len_bins`](https://chengmatt.github.io/SPoRC/dev/reference/obs_len_bins.md)
for lengths.

## Usage

``` r
obs_bin_count(input_list, obs, dim_i, axis)
```

## Arguments

- input_list:

  Input list, used for the fallback.

- obs:

  The observation array for the data source, possibly dimensionless.

- dim_i:

  Integer. Which dimension of `obs` holds the bins.

- axis:

  Either `"age"` or `"len"`.

## Value

A single positive integer.
