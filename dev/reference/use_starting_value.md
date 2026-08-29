# Substitute a user-supplied starting value, checking its shape first

Every `Setup_Mod_*` stage builds a default for each parameter and then
lets `starting_values` replace it. The replacement used to be taken as
given. A value of the wrong shape is not rejected by that: it is carried
into the objective, read position by position, and indexes past its own
end somewhere far from the argument that caused it. What comes back is
RTMB's `'*this' is not a valid 'advector'`, which names nothing useful.

## Usage

``` r
use_starting_value(default, starting_values, nm)
```

## Arguments

- default:

  The parameter as the stage built it.

- starting_values:

  The user's list of starting values.

- nm:

  Name of the parameter.

## Value

`starting_values[[nm]]` when it is supplied and correctly shaped,
otherwise `default`.

## Details

The default already carries the shape the model expects, so it is the
reference the supplied value is measured against.
