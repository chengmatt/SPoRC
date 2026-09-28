# Shorten a dsem to fewer grid years

A retrospective peel drops years from every array, and the dsem has to
follow: the arrows do not change, but the cells they land on do. A
linked series' cells move whenever its array has a dim after the year
dim, so they are rebuilt from where that series sits rather than
trimmed.

## Usage

``` r
peel_dsem_years(data, parameters, mapping, n_grid_yrs)
```

## Arguments

- data:

  Data list holding the dsem fields.

- parameters:

  Parameter list, with its deviation arrays already peeled.

- mapping:

  Map list.

- n_grid_yrs:

  Rows the grid should keep.

## Value

List with the peeled `data`, `parameters` and `mapping`.
