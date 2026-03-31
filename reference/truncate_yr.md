# Helper function to truncate data years, parameters, and mapping to conduct retrospective diagnostics. Called within do_retrospective function.

Helper function to truncate data years, parameters, and mapping to
conduct retrospective diagnostics. Called within do_retrospective
function.

## Usage

``` r
truncate_yr(j, data, parameters, mapping)
```

## Arguments

- j:

  The years to truncate from the terminal year

- data:

  Data list used for the RTMB model

- parameters:

  Parameter list used for the RTMB model

- mapping:

  Mapping list used for the RTMB model

## Value

List of data, parameters, and mapping that have truncated dimensions
from the original data, parameters, and mapping list

## Examples

``` r
if (FALSE) { # \dontrun{
retro_list <- retro_truncate_year(j = 0, data, parameters, mapping) # does not remove any data
retro_list <- retro_truncate_year(j = 1, data, parameters, mapping) # removes last year of data
} # }
```
