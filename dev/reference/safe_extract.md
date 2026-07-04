# Safely extract a named element from a TMB report object

Returns the named element if it exists and is non-`NULL`; returns `0`
otherwise. Used to guard against missing or `NULL` report fields when
accumulating likelihood components.

## Usage

``` r
safe_extract(obj, name)
```

## Arguments

- obj:

  Named list, typically `obj$rep` from a fitted RTMB model.

- name:

  Character string. Name of the element to extract.

## Value

The value of `obj[[name]]` if present and non-`NULL`; `0` otherwise.
