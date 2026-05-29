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
covr::package_coverage(quiet = FALSE, type = 'all')

# Documentation -----------------------------------------------------------
att_amend_desc(update = TRUE)
document()
roxygenise()

# Dependency adjustments
desc::desc_del_dep("compResidual", "Imports")
desc::desc_del_dep("remotes", "Imports")
desc::desc_del_dep("compResidual", "Suggests")

# pkgdown Site ------------------------------------------------------------
pkgdown::clean_cache()
pkgdown::clean_site(force = TRUE)
pkgdown::check_pkgdown()
pkgdown::build_news()
pkgdown::build_site()
servr::httd("docs")

# One-time setup
# usethis::use_pkgdown_github_pages()
# pkgdown::build_site(preview = TRUE)

# Build Ignore ------------------------------------------------------------
use_build_ignore("dev")
usethis::use_build_ignore("_pkgdown.yml")
usethis::use_build_ignore(".claude")

# Check, Build & Install --------------------------------------------------
test()
devtools::check(args = c("--no-tests"))
build()
install()
unloadNamespace("SPoRC")

# Vignettes ---------------------------------------------------------------
# build_vignettes()

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

# usethis::use_testthat()
# usethis::use_test("dusky_rtmb")
# usethis::use_test("sabie_sgl_rtmb")
# usethis::use_test("ebs_pol_sgl_rtmb")
# usethis::use_test("sabie_three_rg_rtmb")
# usethis::use_test("sgl_rg_simple_sim_test")
# usethis::use_test("sgl_rg_spr_sabie_test")
# usethis::use_test("sgl_rg_bh_msy_sabie_test")
# usethis::use_test("mlt_rg_global_spr_sabie_test")
# usethis::use_test("mlt_rg_global_bh_msy_sabie_test")
# usethis::use_test("mlt_rg_local_bh_msy_sabie_test")
# usethis::use_test("mlt_rg_pop_seas_global_spr_test")
# usethis::use_test("mlt_rg_pop_seas_local_bh_msy_test")
# usethis::use_test("custom_distributions")
# usethis::use_test("custom_simulate_distributions")
# usethis::use_test("diag_util_func_dusky")
# usethis::use_test('3d_precision_matrix')
# usethis::use_test('logistN_utils')
# usethis::use_test('release_tag_attr')


# Misc One-Time Setup -----------------------------------------------------
# usethis::use_news_md()
# usethis::edit_git_ignore()
# usethis::use_github_action("test-coverage")

# Undo dev/ build ignore:
# rbuildignore <- readLines(".Rbuildignore")
# rbuildignore <- rbuildignore[!grepl("^\\^dev\\$", rbuildignore)]
# writeLines(rbuildignore, ".Rbuildignore")
