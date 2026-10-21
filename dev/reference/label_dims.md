# Name the dimensions of a report, parameter list or data list

Returns a copy of `x` with every array's dims named (`Pop`, `Region`,
`Year`, `Seas`, `Age`, `Len`, `Sex`, `Fleet`, ...) and labeled with
calendar years, model ages and length bins. `x` can be the report
(`obj$rep`), the parameter list (`input_list$par` or
`obj$env$parList()`) or the data list (`input_list$data`), and
[`reshape2::melt()`](https://rdrr.io/pkg/reshape2/man/melt.html) on a
labeled array returns those columns by name. Package functions read
these arrays by position, so pass them the unlabeled lists.

## Usage

``` r
label_dims(x, data = x, labels = list())
```

## Arguments

- x:

  Report, parameter or data list from a SPoRC model.

- data:

  The model's data list, `input_list$data`, which gives the years, ages
  and length bins. Defaults to `x`, for labeling the data list itself.

- labels:

  Optional named list of labels for `Pop`, `Region`, `Seas`, `Sex`,
  `Fish_Fleet` or `Srv_Fleet`, one label per population, region and so
  on (e.g. `list(Sex = c("Female", "Male"))`). Unlabeled dims are
  numbered from one.

## Value

`x`, with `dimnames` on its arrays and `names` on its vectors. Values
are unchanged. Scalars, lists, vectors of indices and arrays of a shape
this function does not describe are returned as they were. Dims that
count blocks, parameters or selectivity bins (`Block`, `Sel_Par`, `Bin`,
...) keep any names they had and are otherwise numbered from one, as is
any dim whose length differs from its label set, such as observed ages
pooled differently from the model ages. Year dims past the last model
year continue the calendar years, so `NAA` ends at the terminal year
plus one.

## See also

Other Utility:
[`fit_model()`](https://chengmatt.github.io/SPoRC/dev/reference/fit_model.md),
[`get_logistN_Sigma()`](https://chengmatt.github.io/SPoRC/dev/reference/get_logistN_Sigma.md),
[`get_par_est_info()`](https://chengmatt.github.io/SPoRC/dev/reference/get_par_est_info.md),
[`post_optim_sanity_checks()`](https://chengmatt.github.io/SPoRC/dev/reference/post_optim_sanity_checks.md),
[`rho_trans()`](https://chengmatt.github.io/SPoRC/dev/reference/rho_trans.md),
[`set_data_indicator_unused()`](https://chengmatt.github.io/SPoRC/dev/reference/set_data_indicator_unused.md)

## Examples

``` r
if (FALSE) { # \dontrun{
labels <- list(Sex = c("Female", "Male"))
rep <- label_dims(obj$rep, input_list$data, labels)
par <- label_dims(input_list$par, input_list$data, labels)
dat <- label_dims(input_list$data, labels = labels)
rep$SSB[1, , "2020"]
ssb_df <- reshape2::melt(rep$SSB) # columns Pop, Region, Year, value
} # }
```
