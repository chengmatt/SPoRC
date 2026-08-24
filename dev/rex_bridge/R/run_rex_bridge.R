# Evaluate the SPoRC rex sole model at the SS3 maximum likelihood estimate and compare
# every bridged quantity with the assessment's report. Run from the package root.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("tests/testthat/helper-bridge_goa_rex.R")
dat <- readRDS("dev/rex_bridge/output/rex_bridge_data.rds")

# quantities the builder left as raw tables
dat$catch_se_value <- unique(na.omit(as.vector(dat$catch_se)))[1]
dat$growth_A1 <- dat$ctl$Growth_Age_for_L1; dat$growth_A2 <- dat$ctl$Growth_Age_for_L2
dat$mat$first_mature_age <- dat$ctl$First_Mature_Age
va <- dat$var_adj; dat$var_adj_len <- rep(1, 3); dat$var_adj_age <- rep(1, 3)
names(va) <- tolower(names(va))
for(i in seq_len(nrow(va))) { if(va$factor[i] == 4) dat$var_adj_len[va$fleet[i]] <- va$value[i]; if(va$factor[i] == 5) dat$var_adj_age[va$fleet[i]] <- va$value[i] }
qp <- dat$q$parms; dat$q$prior_mean <- qp[grep("LnQ_base_NonEastern", rownames(qp)), "PRIOR"]; dat$q$prior_sd <- qp[grep("LnQ_base_NonEastern", rownames(qp)), "PR_SD"]
cat("A1/A2:", dat$growth_A1, dat$growth_A2, " first mature age:", dat$mat$first_mature_age, " catch se:", dat$catch_se_value, "\n")
cat("var adj len:", dat$var_adj_len, " age:", dat$var_adj_age, " q prior:", dat$q$prior_mean, dat$q$prior_sd, "\n")

input <- suppressMessages(suppressWarnings(build_goa_rex_input(dat)))
input <- seed_goa_rex_mle(input, dat)
obj <- fit_model(input$data, input$par, input$map, do_optim = FALSE, silent = TRUE)
r <- obj$rep
yrs <- dat$years; n_yrs <- length(yrs); n_ages <- length(dat$ages)
pct <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)

cat("\n=== growth (area 1 female): SPoRC vs SS3 ===\n")
g <- dat$ss3$growth[[1]][[1]]
cmp <- data.frame(age = dat$ages, Len_Beg_ss3 = g$Len_Beg, Len_Beg = r$mean_LAA_spawn[1, 1, 1, 1, , 1], Len_Mid_ss3 = g$Len_Mid, Len_Mid = r$mean_LAA_srv[1, 1, 1, 1, , 1, 1],
                  SD_Mid_ss3 = g$SD_Mid, SD_Mid = r$sd_LAA_srv[1, 1, 1, 1, , 1, 1], Wt_Beg_ss3 = g$Wt_Beg, Wt_Beg = r$WAA[1, 1, 1, 1, , 1], Wt_Mid_ss3 = g$Wt_Mid, Wt_Mid = r$WAA_fish[1, 1, 1, 1, , 1, 1])
print(round(cmp[c(1:4, 19:21), ], 5))
cat("max % gap Len_Beg:", signif(pct(cmp$Len_Beg, cmp$Len_Beg_ss3), 3), " Len_Mid:", signif(pct(cmp$Len_Mid, cmp$Len_Mid_ss3), 3), " SD_Mid:", signif(pct(cmp$SD_Mid, cmp$SD_Mid_ss3), 3),
    " Wt_Beg:", signif(pct(cmp$Wt_Beg, cmp$Wt_Beg_ss3), 3), " Wt_Mid:", signif(pct(cmp$Wt_Mid, cmp$Wt_Mid_ss3), 3), "\n")
for(a in 1:2) for(s in 1:2) { g <- dat$ss3$growth[[a]][[s]]; cat(sprintf("  area %d sex %d: Len_Mid %.2g%%  SD_Mid %.2g%%  Wt_Mid %.2g%%\n", a, s, pct(r$mean_LAA_srv[1, a, 1, 1, , s, 1], g$Len_Mid), pct(r$sd_LAA_srv[1, a, 1, 1, , s, 1], g$SD_Mid), pct(r$WAA_fish[1, a, 1, 1, , s, 1], g$Wt_Mid))) }

cat("\n=== ALK (mid season, area 1 female) ===\n")
alk_ss3 <- dat$ss3$ALK[, , "Seas: 1 Sub_Seas: 2 Morph: 1"][rev(seq_len(n_lens_alk <- dim(dat$ss3$ALK)[1])), ] # SS3 prints lengths largest first
alk_sporc <- r$SizeAgeTrans[1, 1, 1, 1, , , 1]
cat("max abs gap:", signif(max(abs(alk_sporc - alk_ss3)), 3), "\n")

cat("\n=== selectivity (fishery, ages 0..20) ===\n")
for(s in 1:2) cat(sprintf("  sex %d max abs gap: %.3g\n", s, max(abs(r$fish_sel[1, 1, 1, 1, , s, 1] - dat$ss3$sel[[1]][[s]]))))
for(sf in 1:2) for(s in 1:2) cat(sprintf("  survey %d sex %d max abs gap: %.3g\n", sf, s, max(abs(r$srv_sel[1, 1, 1, 1, , s, sf] - dat$ss3$sel[[sf + 1]][[s]]))))

cat("\n=== numbers at age (start of year) ===\n")
for(a in 1:2) for(s in 1:2) cat(sprintf("  area %d sex %d: max %% gap %.3g (1982 age 0: %.1f vs %.1f; age 20: %.1f vs %.1f)\n", a, s,
  pct(r$NAA[1, a, 1:n_yrs, 1, , s], dat$ss3$NAA[a, , , s]), r$NAA[1, a, 1, 1, 1, s], dat$ss3$NAA[a, 1, 1, s], r$NAA[1, a, 1, 1, n_ages, s], dat$ss3$NAA[a, 1, n_ages, s]))

cat("\n=== SSB / recruitment / total biomass ===\n")
for(a in 1:2) cat(sprintf("  area %d: SSB %% gap %.3g  Rec %% gap %.3g  Bio_all %% gap %.3g\n", a,
  pct(r$SSB[1, a, 1:n_yrs], dat$ss3$SSB[, a]), pct(r$Rec[1, a, 1:n_yrs], dat$ss3$Rec[, a]),
  pct(sapply(1:n_yrs, function(y) sum(r$NAA[1, a, y, 1, , ] * r$WAA[1, a, y, 1, , ])), dat$ss3$Bio_all[, a])))
cat("  SSB 1982:", r$SSB[1, , 1], " vs", dat$ss3$SSB[1, ], "\n")

cat("\n=== catch and indices ===\n")
cat("  catch % gap:", signif(pct(r$PredCatch[1, 1, , 1, 1], dat$ss3$dead_B), 3), "\n")
cp <- dat$ss3$cpue
for(sf in 1:2) { ci <- cp[cp$Fleet == sf + 1, ]; cat(sprintf("  survey %d index %% gap: %.3g\n", sf, pct(r$PredSrvIdx[1, sf, match(ci$Yr, yrs), 1, sf], ci$Exp))) }

cat("\n=== likelihood components ===\n")
cat("  SPoRC: catch", round(sum(r$Catch_nLL), 4), " srv idx", round(sum(r$SrvIdx_nLL), 4), " fish len", round(sum(r$FishLenComps_nLL), 3), " srv len", round(sum(r$SrvLenComps_nLL), 3),
    " fish age", round(sum(r$FishAgeComps_nLL), 3), " srv caal", round(sum(r$Srv_caal_nLL), 3), " rec", round(sum(r$Rec_nLL), 3), " init", round(sum(r$Init_Rec_nLL), 3), " q prior", round(sum(r$srv_q_nLL), 4), "\n")
print(dat$ss3$likelihoods)
cat("\nmax |gradient| at the SS3 MLE:", signif(max(abs(obj$gr(obj$par))), 3), "\n")
g <- obj$gr(obj$par); nm <- names(obj$par); big <- order(-abs(g))[1:8]; print(data.frame(par = nm[big], grad = signif(g[big], 3)))
saveRDS(list(input = input, rep = r), "dev/rex_bridge/output/rex_sporc_at_mle.rds")
