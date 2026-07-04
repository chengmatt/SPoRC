library(devtools)
library(pkgdown)
library(roxygen2)
library(attachment)
library(usethis)
library(here)

# Authentication -----------------------------------------------
# usethis::create_github_token()
# gitcreds::gitcreds_set()
Sys.unsetenv("GITHUB_PAT") # unset before building vignettes

# Code Coverage -----------------------------------------------------------
# covr::package_coverage(quiet = FALSE, type = 'all')

# Documentation -----------------------------------------------------------
desc::desc_set_version("1.2.0.9000")

# Build Package -----------------------------------------------------------
document() # document functions
roxygenise() # make sure functions have roxygen documentation
att_amend_desc(update = TRUE)

# Dependency adjustments
desc::desc_del_dep("compResidual", "Imports")
desc::desc_del_dep("remotes", "Imports")
desc::desc_del_dep("compResidual", "Suggests")

# Add suggested packages
usethis::use_package("igraph", type = "Suggests")
usethis::use_package("splines2", type = "Suggests")

# pkgdown Site ------------------------------------------------------------
pkgdown::clean_cache()
pkgdown::clean_site(force = TRUE)
pkgdown::check_pkgdown()
pkgdown::build_news()
pkgdown::build_site()
servr::httd("docs")

# Build Ignore ------------------------------------------------------------
use_build_ignore("dev")
usethis::use_build_ignore("_pkgdown.yml")
usethis::use_build_ignore(".claude")

# Check, Build & Install --------------------------------------------------
test()
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning")
devtools::check(args = c("--no-tests"))
build()
install()
unloadNamespace("SPoRC")


# Vignettes ---------------------------------------------------------------
# build_vignettes()

# Build Vignettes ---------------------------------------------------------

# usethis::use_vignette("a_model_dimensions")
# usethis::use_vignette("b_model_parameters")
# usethis::use_vignette("c_model_equations")
# usethis::use_vignette("d_model_report")
# usethis::use_vignette("e_single_region_sablefish_case_study")
# usethis::use_vignette("f_single_region_ebs_pollock_case_study")
# usethis::use_vignette("g_spatial_sablefish_case_study")
# usethis::use_vignette("h_closed_loop_simulations")
# usethis::use_vignette("i_reference_points")
# usethis::use_vignette("j_starting_mapping")
# usethis::use_vignette("k_defining_priors")
# usethis::use_vignette("l_simulation_testing")
# usethis::use_vignette("m_simulation_dimensions")
# usethis::use_vignette("n_single_region_ebs_pollock_randomeff_case_study")
# usethis::use_vignette("o_get_started")
# usethis::use_vignette("p_single_region_dusky_alt_mp_testing")
# usethis::use_vignette("q_movement_param")
# usethis::use_vignette("r_natal-homing-pop-lrgr-rg")
# usethis::use_vignette("s_discard_retention")
# usethis::use_vignette('t_model_options')


# Integration Tests --------------------------------------------------------------
# devtools::test()

# Misc One-Time Setup -----------------------------------------------------
# usethis::use_news_md()
# usethis::edit_git_ignore()
# usethis::use_github_action("test-coverage")
# usethis::use_github_action("check-standard")
# usethis::use_github_action("pkgdown")
# usethis::use_citation()

# Undo dev/ build ignore:
# rbuildignore <- readLines(".Rbuildignore")
# rbuildignore <- rbuildignore[!grepl("^\\^dev\\$", rbuildignore)]
# writeLines(rbuildignore, ".Rbuildignore")
