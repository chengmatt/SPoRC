# Argument descriptions from a function's own help page

Reads the `\arguments` section of the installed (or source) Rd, so the
reference and
[`?Setup_Mod_Rec`](https://chengmatt.github.io/SPoRC/dev/reference/Setup_Mod_Rec.md)
cannot disagree. One Rd item may document several arguments at once
(`"a,b,c"`), and each gets the shared text.

## Usage

``` r
rd_argument_text(topic, db)
```

## Arguments

- topic:

  Function name.

- db:

  Rd database, from [`Rd_db`](https://rdrr.io/r/tools/Rdutils.html).

## Value

Named character vector of descriptions, empty when the topic is absent.
