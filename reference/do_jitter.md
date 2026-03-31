# Run Jitter Analysis

Run Jitter Analysis

## Usage

``` r
do_jitter(
  data,
  parameters,
  mapping,
  random = NULL,
  sd,
  n_jitter,
  n_newton_loops,
  do_par,
  n_cores,
  par_vec = NULL
)
```

## Arguments

- data:

  Data list to make obj

- parameters:

  Parameter list to make obj

- mapping:

  Mapping list to make obj

- random:

  Character of random effects

- sd:

  sd for jitter (additive)

- n_jitter:

  Number of jitters to do

- n_newton_loops:

  Number of newton loops to do

- do_par:

  Whether to do paralleizaiton or not (boolean)

- n_cores:

  Number of cores to use

- par_vec:

  Vector of parameter starting values to use for jitter analysis. The
  default of this is NULL (jitters the starting value of the model). If
  a vector is provided, the jitter is initialized at the MLE parameters

## Value

Dataframe of jitter values

## See also

Other Model Diagnostics:
[`do_likelihood_profile()`](https://chengmatt.github.io/SPoRC/reference/do_likelihood_profile.md),
[`do_retrospective()`](https://chengmatt.github.io/SPoRC/reference/do_retrospective.md),
[`do_runs_test()`](https://chengmatt.github.io/SPoRC/reference/do_runs_test.md),
[`get_catch_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_catch_fits_plot.md),
[`get_comp_prop()`](https://chengmatt.github.io/SPoRC/reference/get_comp_prop.md),
[`get_idx_fits()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits.md),
[`get_idx_fits_plot()`](https://chengmatt.github.io/SPoRC/reference/get_idx_fits_plot.md),
[`get_model_rep_from_mcmc()`](https://chengmatt.github.io/SPoRC/reference/get_model_rep_from_mcmc.md),
[`get_nLL_plot()`](https://chengmatt.github.io/SPoRC/reference/get_nLL_plot.md),
[`get_osa()`](https://chengmatt.github.io/SPoRC/reference/get_osa.md),
[`get_retrospective_plot()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_plot.md),
[`get_retrospective_relative_difference()`](https://chengmatt.github.io/SPoRC/reference/get_retrospective_relative_difference.md),
[`plot_resids()`](https://chengmatt.github.io/SPoRC/reference/plot_resids.md)

## Examples

``` r
if (FALSE) { # \dontrun{
   library(ggplot2)
   # get jitter values
   jit <- do_jitter(data = data,
                 parameters = parameters,
                 mapping = mapping,
                 random = NULL,
                 sd = 0.1,
                 n_jitter = 100,
                 n_newton_loops = 3,
                 do_par = TRUE,
                 n_cores = 8)

   # get proportion converged
   prop_converged <- jit %>%
   filter(Year == 1, Type == 'Recruitment') %>%
     summarize(prop_conv = sum(Hessian) / length(Hessian))

   # get final model results
   final_mod <- reshape2::melt(sabie_rtmb_model$rep$SSB) %>%
   rename(Region = Var1, Year = Var2) %>%
   mutate(Type = 'SSB') %>%
     bind_rows(reshape2::melt(sabie_rtmb_model$rep$Rec) %>%
     rename(Region = Var1, Year = Var2) %>% mutate(Type = 'Recruitment'))

   # comparison of SSB and recruitment
  ggplot() +
    geom_line(jit, mapping = aes(x = Year + 1959, y = value,
     group = jitter, color = Hessian), lwd = 1) +
    geom_line(final_mod, mapping = aes(x = Year + 1959, y = value),
    color = "black", lwd = 1.3 , lty = 2) +
    facet_grid(Type~Region, scales = 'free',
               labeller = labeller(Region = function(x) paste0("Region ", x),
    Type = c("Recruitment" = "Age 2 Recruitment (millions)",
     "SSB" = 'SSB (kt)'))) +
    labs(x = "Year", y = "Value") +
    theme_bw(base_size = 20) +
    scale_color_manual(values = c("red", 'grey')) +
    geom_text(data = jit %>% filter(Type == 'SSB', Year == 1, jitter == 1),
              aes(x = Inf, y = Inf, label =
              paste("Proportion Converged: ",
              round(prop_converged$prop_conv, 3))),
              hjust = 1.1, vjust = 1.9, size = 6, color = "black")

   # compare jitter of max gradient and hessian PD
   ggplot(jit, aes(x = jitter, y = jnLL,
   color = Max_Gradient, shape = Hessian)) +
     geom_point(size = 5, alpha = 0.3) +
     geom_hline(yintercept = min(sabie_rtmb_model$rep$jnLL),
     lty = 2, size = 2, color = "blue") +
     facet_wrap(~Hessian, labeller = labeller(
       Hessian = c("FALSE" = "non-PD Hessian", "TRUE" = 'PD Hessian')
     )) +
     scale_color_viridis_c() +
     theme_bw(base_size = 20) +
     theme(legend.position = "bottom") +
     guides(color = guide_colorbar(barwidth = 15, barheight = 0.5)) +
     labs(x = 'Jitter') +
     geom_text(data = jit %>% filter(Hessian == TRUE,
     Year == 1, jitter == 1),
               aes(x = Inf, y = Inf, label =
               paste("Proportion Converged: ",
               round(prop_converged$prop_conv, 3))),
               hjust = 1.1, vjust = 1.9, size = 6, color = "black")
} # }
```
