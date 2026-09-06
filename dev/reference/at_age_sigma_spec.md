# Should an at-age data source's observation error parameter be estimated?

A data source nobody fits has nothing to inform its standard deviation,
and a data source taking its error from reported standard errors alone
has no parameter to read, so both are kept fixed whatever the spec says.

## Usage

``` r
at_age_sigma_spec(spec, form, any_used)
```

## Arguments

- spec:

  `"est"` or `"fix"` as supplied by the caller.

- form:

  The data source's `sigma_form`.

- any_used:

  `TRUE` when any fleet fits this data source.

## Value

`"est"` or `"fix"`.
