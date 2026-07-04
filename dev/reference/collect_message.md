# Append a message to the global messages list

Concatenates its arguments into a single string and appends the result
to `messages_list` in the calling environment via `<<-`. Used internally
to accumulate validation and setup messages for deferred display.

## Usage

``` r
collect_message(...)
```

## Arguments

- ...:

  Character strings passed to `paste(..., sep = "")`.

## Value

`NULL` invisibly. Side effect: `messages_list` is updated in the parent
environment.
