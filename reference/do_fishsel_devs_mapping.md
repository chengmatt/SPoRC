# Helper function to set up fishery selectivity deviations mapping

Helper function to set up fishery selectivity deviations mapping

## Usage

``` r
do_fishsel_devs_mapping(
  input_list,
  fish_sel_devs_spec,
  fishsel_devs_shared_ages
)
```

## Arguments

- input_list:

  Input list

- fish_sel_devs_spec:

  Character vector specifying fishery selectivity deviations
  parameterization

- fishsel_devs_shared_ages:

  List object for specifying which ages are shared when selectivity
  deviations are semi-parametric (e.g., list(1:5, 6:10, 11:30) specifies
  that ages 1-5, 6-10, and 11-30 have the same deviations.)
