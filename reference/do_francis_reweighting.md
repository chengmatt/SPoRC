# Get Francis Weights

Get Francis Weights

## Usage

``` r
do_francis_reweighting(data, rep, age_labels, len_labels, year_labels)
```

## Arguments

- data:

  List of data inputs

- rep:

  Report file list

- age_labels:

  Age labels

- len_labels:

  Length labels

- year_labels:

  Year labels

## Value

A list object of francis weights (note that it will be NAs for some, if
using jnt composition approaches - i.e., only uses one dimension), as
well as a dataframe of francis mean fits

## Details

Function to get francis weights. Used inside the wrapper function
run_francis(), or can be defined by the user as a loop to extract
Francis weights (see example).

## See also

Other Francis Reweighting:
[`run_francis()`](https://chengmatt.github.io/SPoRC/reference/run_francis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
for(j in 1:5) {

  if(j == 1) { # reset weights at 1
    data$Wt_FishAgeComps[] <- 1
    data$Wt_FishLenComps[] <- 1
    data$Wt_SrvAgeComps[] <- 1
    data$Wt_SrvLenComps[] <- 1
  } else {
    data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
    data$Wt_FishLenComps[] <- wts$new_fish_len_wts
    data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
    data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
  }

  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )

  rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report
  wts <- do_francis_reweighting(data = data, rep = rep, age_labels = 2:31,
                                len_labels = seq(41, 99, 2), year_labels = 1960:2024)
}
} # }
```
