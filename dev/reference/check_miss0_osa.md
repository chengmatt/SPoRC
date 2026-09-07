# Refuse the internal composition residuals for a zeros dropped logistic normal

Dropping the zeros leaves a different number of observations in each
cell, while the one step ahead residuals are packed as a fixed length
vector per cell, so the two cannot be combined.

## Usage

``` r
check_miss0_osa(input_list, like_vals, arg_name)
```

## Arguments

- input_list:

  Named list with `$data`.

- like_vals:

  Integer likelihood codes for this data source.

- arg_name:

  Name of the argument being checked, used in the message.

## Value

`NULL`, invisibly. Called for the error it raises.
