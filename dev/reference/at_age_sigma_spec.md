# Should an at-age stream's observation error parameter be estimated?

A stream nobody fits has nothing to inform its standard deviation, and a
stream taking its error from reported standard errors alone has no
parameter to read, so both are held fixed whatever the spec says.

## Usage

``` r
at_age_sigma_spec(spec, form, any_used)
```

## Arguments

- spec:

  `"est"` or `"fix"` as supplied by the caller.

- form:

  The stream's `sigma_form`.

- any_used:

  `TRUE` when any fleet fits this stream.

## Value

`"est"` or `"fix"`.
