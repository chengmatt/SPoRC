# Extend an array along a year dimension

This function takes an array and extends it along the specified year
dimension. The extension can either be filled with zeros or by repeating
the last year slice.

## Usage

``` r
extend_years(arr, n_years, yr_dim, fill = "zeros")
```

## Arguments

- arr:

  Array to extend. Can have any number of dimensions.

- n_years:

  Integer. The total number of years to extend the array to.

- yr_dim:

  Integer. The dimension of \`arr\` that corresponds to years.

- fill:

  Character or Numeric (scalar or array). How to fill the extended
  years: - \`"zeros"\`: fill with zeros - \`"last"\`: repeat the last
  year slice that is not a NaN or NA value. If all values are NA, then
  the array gets populated with NAs. - \`"mean"\`: take mean of time
  series - \`"F_pattern"\`: Used for fishery input sample sizes,
  dynamically fills in sample sizes based on fishing mortality values
  specified in a closed loop simulation - Numeric: Constant scalar or
  array

## Value

An array extended along the \`yr_dim\` dimension.
