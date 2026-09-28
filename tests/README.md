# SPoRC tests

Read this before you change a test to make it pass.

Tests in SPoRC serve two rather different purposes, and the distinction matters because it
determines what to do when one of them fails. Most tests check a relationship that must hold
whatever the particular numbers happen to be, such as whether the model recovers the
parameters used to simulate its own data, so that a failure points fairly directly at what
changed. The regression tests, in contrast, compare a fit against values stored from an
earlier run, and a failure there indicates only that the fit has moved, without saying
whether that movement was intended. In the following, we describe each kind in turn, and
then give the steps to follow when a stored comparison fails.

## Tests that check a relationship

These tests assert a property that any correct implementation must satisfy, rather than a
particular value, and so a failure identifies the property that no longer holds.

- `test-sim_selftest_*.R` simulate data from known parameters, refit the model, and check
  that the estimates recover those parameters. A failure indicates that the model no longer
  recovers the parameters it generated.
- `test-integration_*.R` compare two routes to the same quantity that must agree, for
  example an operating model against the estimation model. A failure indicates that the two
  have diverged, although not which of them is at fault.
- `expect_jnLL_decomposes()` (see `helper-jnll_decomposition.R`) checks that the likelihood
  components sum to the reported total, so that a failure indicates a component missing from
  the sum, or counted twice.
- `test-setup_*.R` and `test-utils_*.R` check input validation, dimension handling, and
  mapping. A failure indicates that a setup function accepts something it should reject, or
  rejects something it should accept.
- `test-model_population_dynamics.R` writes its expectations as arithmetic, for example
  `40 * exp(-0.3)`, such that each assertion has its own derivation and can be checked by
  hand.

## Regression tests

The following tests instead compare against stored numeric vectors:

- `test-regression_dusky.R`
- `test-regression_ebs_pollock_sgl.R`
- `test-regression_sabie_sgl.R`
- `test-regression_sabie_three_rg.R`
- `test-refpts_sgl_rg_spr.R`

These stored values are output from earlier SPoRC fits of these assessments, taken from runs
we had checked and trusted at the time. They are neither hand-derived nor imported from an
external model, and so they hold no authority of their own; their purpose is to detect a
change that moves a fitted result without anyone intending it. A failure is therefore
ambiguous, and should be resolved in the following order.

1. Assume first that you have introduced a bug, which is by far the most common explanation.
   Something in the likelihood, the population dynamics, or the setup path has changed the
   fit without your intending it.
2. If the numerical change really was intended, regenerate the vectors deliberately, confirm
   that the new fit is one you would defend, and record the change in `NEWS.md` together with
   the reason for it.
3. Do not regenerate because the difference appeared small. A change of a few percent in SSB
   is precisely what a real bug tends to look like.

A failure that appears on one platform only, and that sits in the last digit or two, is
generally optimizer or BLAS sensitivity rather than a change in the code. Such comparisons
should be given enough tolerance to absorb that variation, and no stored comparison should be
tightened to `tolerance = 0`.

## Adding a new model or assessment

We recommend a self-test, or a comparison against another route through the package that must
agree, in preference to a new stored vector, because both of those remain informative when
they fail. A fit should be stored only when the point of the test is that this particular
configuration continues to produce this particular answer. If you do store one, add a header
comment giving the run the values came from, and add the file to the list above.
