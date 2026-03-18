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

# Documentation -----------------------------------------------------------
att_amend_desc(update = TRUE)
document()
roxygenise()

# Dependency adjustments
desc::desc_del_dep("compResidual", "Imports")
desc::desc_set_dep("compResidual", "Suggests")
desc::desc_del_dep("remotes", "Imports")

# pkgdown Site ------------------------------------------------------------
pkgdown::clean_cache()
pkgdown::check_pkgdown()
pkgdown::build_site()
pkgdown::build_news()
servr::httd("docs")

# One-time setup
# usethis::use_pkgdown_github_pages()
# pkgdown::build_site(preview = TRUE)

# Build Ignore ------------------------------------------------------------
use_build_ignore("dev")
usethis::use_build_ignore("_pkgdown.yml")

# Check, Build & Install --------------------------------------------------
check()
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


# Unit Tests --------------------------------------------------------------
# devtools::test()

# usethis::use_testthat()
# usethis::use_test("dusky_rtmb")
# usethis::use_test("sabie_sgl_rtmb")
# usethis::use_test("ebs_pol_sgl_rtmb")
# usethis::use_test("sabie_three_rg_rtmb")
# usethis::use_test("sgl_rg_simple_sim_test")


# Misc One-Time Setup -----------------------------------------------------
# usethis::use_news_md()
# usethis::edit_git_ignore()

# Undo dev/ build ignore:
# rbuildignore <- readLines(".Rbuildignore")
# rbuildignore <- rbuildignore[!grepl("^\\^dev\\$", rbuildignore)]
# writeLines(rbuildignore, ".Rbuildignore")
