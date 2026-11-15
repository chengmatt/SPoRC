# Move a set of tag vectors through one matrix

Every vector is a cohort's tags by region, moved one cohort at a time,
`t(tag_vec) %*% move_mat` for each.

## Usage

``` r
move_tag_vectors(tag_vecs, move_mat)
```

## Arguments

- tag_vecs:

  List of numeric vectors `[region]`.

- move_mat:

  Matrix `[region_from, region_to]`.

## Value

List of the moved vectors `[region_to]`, in the same order.
