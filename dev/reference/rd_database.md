# The package's Rd database, wherever it is being read from

Three places this gets called from and each needs a different route. An
installed package answers to its name. A source tree answers to its
root, which is not the working directory when the caller is a vignette
or a test, so the root is walked up to. And
[`pkgload::load_all`](https://pkgload.r-lib.org/reference/load_all.html)
shadows the installed help without building its index, which is why the
name route is tried and allowed to fail rather than relied on.

## Usage

``` r
rd_database()
```

## Value

A named list of parsed Rd, empty when none can be found.
