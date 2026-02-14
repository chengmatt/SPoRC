library(devtools)
library(pkgdown)
library(roxygen2)
library(attachment)
library(usethis)
library(here)

# Build Package -----------------------------------------------------------
att_amend_desc(update = TRUE)
document() # document functions
roxygenise() # make sure functions have roxygen documentation
desc::desc_del_dep("compResidual", "Imports") # move to suggests
desc::desc_set_dep("compResidual", "Suggests") # move to suggests
desc::desc_del_dep("remotes", "Imports") # remove remotes
Sys.unsetenv("GITHUB_PAT") # may need to unset to build vignettes
build_site(examples = FALSE)
use_build_ignore("dev") # ignore dev folder
usethis::use_build_ignore("_pkgdown.yml") # ignore pkgdown.yml
check() # check package stuff
Sys.unsetenv("GITHUB_PAT") # may need to unset to build vignettes
build() # build package
install() # install locally
unloadNamespace('SPoRC')

# Do tests
# devtools::test()

# Gitignore
# usethis::edit_git_ignore()
# unignore dev folder
# rbuildignore <- readLines(".Rbuildignore")
# rbuildignore <- rbuildignore[!grepl("^\\^dev\\$", rbuildignore)]
# writeLines(rbuildignore, ".Rbuildignore")

# Add news
# usethis::use_news_md() # add news md

# Unit Tests --------------------------------------------------------------
# usethis::use_testthat()
# usethis::use_test("dusky_rtmb")
# usethis::use_test("sabie_sgl_rtmb")
# usethis::use_test("ebs_pol_sgl_rtmb")
# usethis::use_test("sabie_three_rg_rtmb")
# usethis::use_test("sgl_rg_simple_sim_test")

# Build Vignettes ---------------------------------------------------------
# build_vignettes() # build vignettes
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
