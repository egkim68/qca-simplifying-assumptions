# analyze.R -- every number, table and figure input of the paper, from the raw simulation output.
# Usage: Rscript analyze.R [cores]   (after simulate.R main and simulate.R extension; writes results/)
suppressMessages(source("simulate.R"))
cores <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else max(1, parallel::detectCores() - 1)
OUT <- "results"; dir.create(OUT, showWarnings = FALSE)
RES <- list(); save <- function(x, f) { RES[[sub("\\.csv$", "", f)]] <<- x; invisible(x) }   # collected, written at the end
readall <- function(dir, pat) do.call(rbind, lapply(list.files(dir, pat, full.names = TRUE), read.csv))
auc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y)); n0 <- as.numeric(sum(!y)); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
keys <- c("condition_id", "shape", "n", "N", "dist", "noise", "incl.cut", "family")
ty <- c("SA", "csa", "noncsa", "unt", "ten", "d1", "d2", "SAc")

d <- readall("raw/main", "^cell_"); tm <- readall("raw/main", "^terms_")
d$family <- ifelse(d$shape == "N1", "null", "structured")
d0_total <- nrow(d); d0_excluded <- d$status != "ok"
cat("studies:", nrow(d), " excluded:", sum(d$status != "ok"), "\n")
d <- d[d$status == "ok", ]; main <- d$incl.cut == 0.8
tm <- merge(tm, unique(d[, keys]), by = "condition_id")
tp <- tm[tm$solution == "pars" & tm$family == "structured" & tm$incl.cut == 0.8, ]
tp$ds <- paste(tp$condition_id, tp$rep)

# ---- 1. precision of assumption types (Table 2, Figure 2, H1-H4) ----
okc <- c(rbind(c(ty, "unt2", "ten2"), paste0(c(ty, "unt2", "ten2"), "_ok")))
cell <- aggregate(d[, okc], d[, keys], sum)
pr <- function(k, m) ifelse(m > 0, k / m, NA)
cell <- within(cell, {
  p_SA <- pr(SA_ok, SA); p_csa <- pr(csa_ok, csa); p_noncsa <- pr(noncsa_ok, noncsa)
  p_unt <- pr(unt_ok, unt); p_ten <- pr(ten_ok, ten); p_d1 <- pr(d1_ok, d1); p_d2 <- pr(d2_ok, d2); p_SAc <- pr(SAc_ok, SAc)
  h1 <- ifelse(csa >= 30 & noncsa >= 30, p_noncsa - p_csa, NA)
  h2 <- ifelse(unt >= 30 & ten >= 30, p_ten - p_unt, NA)
  h3 <- p_SAc - p_SA; retained <- pr(SAc, SA)
  h4 <- ifelse(d1 >= 30 & d2 >= 30, p_d1 - p_d2, NA)
})
save(cell, "cells.csv")
pool <- function(z, by) {
  g <- aggregate(z[, c(rbind(ty, paste0(ty, "_ok")))], z[by], sum)
  with(g, data.frame(g[by], prec_all = SA_ok / SA, csa_share = csa / SA, prec_csa = csa_ok / csa, prec_noncsa = noncsa_ok / noncsa,
                     unt_share = unt / SA, prec_unt = unt_ok / unt, prec_ten = ten_ok / ten, d2_share = d2 / SA,
                     prec_d1 = d1_ok / d1, prec_d2 = d2_ok / d2, prec_clean = SAc_ok / SAc, retained = SAc / SA))
}
S <- d[main & d$family == "structured", ]
save(pool(S, "noise"), "pooled_by_noise.csv")
save(pool(d[main, ], c("family", "noise")), "pooled_by_family_noise.csv")
save(pool(S, c("N", "noise")), "pooled_by_N_noise.csv")
nec <- aggregate(cbind(unt2, unt2_ok, ten2, ten2_ok, unt) ~ noise, S, sum)
save(within(nec, { prec_unt2 <- unt2_ok / unt2; prec_ten2 <- ten2_ok / ten2 }), "necessity_consistency_only.csv")

# ---- 2. terms: sufficiency, termwise B&T, support (Figures 3, 4a; H5) ----
save(aggregate(cbind(pars_all_suff, pars_all_bt, clean_all_suff, clean_all_bt) ~ noise, S, mean), "solution_level.csv")
tp$sgrp <- cut(tp$support_rows, c(0, 2, 4, Inf), labels = c("1-2", "3-4", "5+"))
ts <- aggregate(cbind(sufficient, bt_correct) ~ sgrp + noise, tp, mean)
ts$terms <- aggregate(sufficient ~ sgrp + noise, tp, length)$sufficient; save(ts, "terms_by_support.csv")
tcl <- tm[tm$family == "structured" & tm$incl.cut == 0.8, ]
save(aggregate(cbind(sufficient, bt_correct, support_rows, rem_share) ~ solution + noise, tcl, mean), "terms_pars_vs_clean.csv")
h5 <- do.call(rbind, lapply(split(tp[tp$noise == 0, ], tp$n[tp$noise == 0]), function(z)
  data.frame(n = z$n[1], low_terms = sum(z$sgrp == "1-2"), low_suff = mean(z$sufficient[z$sgrp == "1-2"]),
             high_terms = sum(z$sgrp == "5+"), high_suff = if (any(z$sgrp == "5+")) mean(z$sufficient[z$sgrp == "5+"]) else NA)))
h5$ok <- with(h5, ifelse(low_terms >= 100 & high_terms >= 100, low_suff <= 0.60 & high_suff >= 0.85, NA))
save(h5, "h5_terms.csv")
nt <- tm[tm$family == "null", ]   # random control
save(do.call(rbind, lapply(split(nt, nt$solution), function(z) data.frame(solution = z$solution[1],
  sufficient = mean(z$sufficient), bt_correct = mean(z$bt_correct), terms = nrow(z)))), "random_control_terms.csv")

# ---- 3. predicting term sufficiency (H6, Figure 4b): half split by study, bootstrap of test studies ----
fits <- list(raw_cov = sufficient ~ raw_cov, support = sufficient ~ log(support_rows), rem_share = sufficient ~ rem_share,
             full = sufficient ~ log(support_rows) + rem_share + log(size))
set.seed(20261112)
h6 <- do.call(rbind, lapply(sort(unique(tp$noise)), function(nz) {
  z <- tp[tp$noise == nz & !is.na(tp$raw_cov), ]
  ids <- unique(z$ds); trn <- sample(ids, length(ids) / 2)
  A <- z[z$ds %in% trn, ]; B <- z[!z$ds %in% trn, ]
  sc <- lapply(fits, function(f) predict(glm(f, binomial, A), B))
  idx <- split(seq_len(nrow(B)), B$ds); u <- names(idx)
  bs <- replicate(2000, { i <- unlist(idx[sample(u, replace = TRUE)], use.names = FALSE)
    auc(sc$full[i], B$sufficient[i]) - auc(sc$raw_cov[i], B$sufficient[i]) })
  data.frame(noise = nz, test_terms = nrow(B), t(sapply(sc, auc, y = B$sufficient)),
             diff = auc(sc$full, B$sufficient) - auc(sc$raw_cov, B$sufficient), lo = quantile(bs, .025), hi = quantile(bs, .975))
}))
save(h6, "h6_auc.csv")
q <- tp[!is.na(tp$raw_cov), ]
save(do.call(rbind, lapply(sort(unique(q$noise)), function(nz) { z <- q[q$noise == nz, ]; p <- z[z$rem_share > 0, ]
  data.frame(noise = nz, share_rho0 = mean(z$rem_share == 0), auc_rem_pos = auc(-p$rem_share, p$sufficient), auc_raw_pos = auc(p$raw_cov, p$sufficient))
})), "remainder_share_restricted.csv")

# ---- 4. hypothesis decisions (Table 3) ----
SC <- cell[cell$family == "structured" & cell$incl.cut == 0.8, ]
share <- function(v, thr) c(mean(v >= thr, na.rm = TRUE), sum(!is.na(v)))
H <- rbind(H1 = share(SC$h1, 0.15), H2 = share(SC$h2, 0.10), H3 = share(SC$h3, 0.05), H4 = share(SC$h4, 0.15))
hyp <- data.frame(hypothesis = paste0("H", 1:6),
  estimate = c(H[, 1], mean(h5$ok, na.rm = TRUE), min(h6$diff)), usable = c(H[, 2], sum(!is.na(h5$ok)), nrow(h6)),
  second = c(NA, NA, median(SC$retained, na.rm = TRUE), NA, NA, min(h6$lo)),
  decision = c(H["H1", 1] >= 0.80, H["H2", 1] >= 0.75, H["H3", 1] >= 0.75 & median(SC$retained, na.rm = TRUE) >= 0.60,
               H["H4", 1] >= 0.80, all(h5$ok, na.rm = TRUE), min(h6$diff) >= 0.05 && min(h6$lo) > 0),
  share_gain_positive = c(NA, NA, mean(SC$h3 > 0, na.rm = TRUE), NA, NA, NA))
save(hyp, "hypotheses.csv")

# ---- 5. sensitivity checks (Section 3.6) ----
qq <- tp[!is.na(tp$raw_cov), ]
within_ <- do.call(rbind, lapply(c(0, 0.1), function(nz) do.call(rbind, lapply(c("shape", "n", "N", "dist"), function(v) {
  z <- qq[qq$noise == nz, ]
  do.call(rbind, lapply(split(z, z[[v]]), function(g) data.frame(noise = nz, factor = v, level = as.character(g[[v]][1]),
    terms = nrow(g), auc_rem = auc(-g$rem_share, g$sufficient))))
}))))
save(within_[within_$terms >= 500, ], "within_levels.csv")
held <- do.call(rbind, lapply(c(0, 0.1), function(nz) do.call(rbind, lapply(c("shape", "n", "N", "dist"), function(v) {
  z <- qq[qq$noise == nz, ]
  do.call(rbind, lapply(sort(unique(z[[v]])), function(lv) {
    A <- z[z[[v]] != lv, ]; B <- z[z[[v]] == lv, ]
    data.frame(noise = nz, held_out = v, level = as.character(lv), test_terms = nrow(B),
               auc_raw = auc(predict(glm(fits$raw_cov, binomial, A), B), B$sufficient),
               auc_support_model = auc(predict(glm(fits$full, binomial, A), B), B$sufficient))
  }))
}))))
save(held, "held_out.csv")

# Regenerated studies (identical data from the saved random-number streams)
g <- conditions_main(); gD <- g[g$shape != "N1" & g$incl.cut == 0.8, ]
pp <- function(z, k) sum(z[[paste0(k, "_ok")]], na.rm = TRUE) / sum(z[[k]], na.rm = TRUE)
# (a) consistency cut-offs on the same studies (Supplementary Table S3)
cut_cells <- gD$condition_id[gD$shape %in% c("D2", "D3") & gD$n %in% 5:6 & gD$N == 40 & gD$dist == "uniform"]
cuts <- do.call(rbind, parallel::mclapply(cut_cells, function(id) do.call(rbind, lapply(c(0.75, 0.8, 1), function(ct) {
  cond <- gD[gD$condition_id == id, ]; cond$incl.cut <- ct
  x <- run_cell(cond, 60261110 + id)$ds; x$cut <- ct; x[x$status == "ok", ] })), mc.cores = cores))
cs <- do.call(rbind, lapply(split(cuts, list(cuts$noise, cuts$cut)), function(z) data.frame(noise = z$noise[1], cut = z$cut[1],
  studies = nrow(z), t(sapply(ty, pp, z = z)))))
save(cs[order(cs$noise, cs$cut), ], "cutoff_same_studies.csv")
# (b) alternative criterion: terms scored against the generating conjunctions, without noise
alt <- do.call(rbind, parallel::mclapply(gD$condition_id[gD$noise == 0], function(id)
  run_cell(gD[gD$condition_id == id, ], 60261110 + id, ref = "gen")$terms, mc.cores = cores))
w <- aggregate(bt_correct ~ condition_id + rep + solution, alt, all)
save(data.frame(aggregate(bt_correct ~ solution, alt, mean), solution_all = aggregate(bt_correct ~ solution, w, mean)$bt_correct),
     "generating_conjunctions.csv")
# (c) outcome prevalence among remainders (base rate of 0.275)
base <- do.call(rbind, parallel::mclapply(gD$condition_id, function(id) {
  cond <- gD[gD$condition_id == id, ]; XNOISE <<- 0
  RNGkind("L'Ecuyer-CMRG"); set.seed(60261110 + id); s <- .Random.seed; R <- R1 <- 0
  for (r in 1:500) {
    assign(".Random.seed", s, envir = .GlobalEnv); s <- parallel::nextRNGStream(s)
    tr <- gen_departure(cond$shape, cond$n); tt <- truth_from_cases(gen_cases(tr$f, cond$n, cond$N, cond$dist, cond$noise), cond$n)
    if (length(tt$T) && length(tt$R)) { R <- R + length(tt$R); R1 <- R1 + sum(tr$f[tt$R + 1L]) }
  }
  data.frame(noise = cond$noise, R = R, R1 = R1) }, mc.cores = cores))
b <- aggregate(cbind(R, R1) ~ noise, base, sum); b$prev_remainders <- b$R1 / b$R; save(b, "base_rate.csv")

# ---- 6. errors in the conditions (Supplementary Table S4, Section 3.7) ----
x <- readall("raw/extension", "^cell_"); xt <- readall("raw/extension", "^terms_")
m0 <- d[d$family == "structured" & d$noise == 0 & main, ]; m0$xnoise <- 0
xk <- c("condition_id", "xnoise", c(rbind(ty, paste0(ty, "_ok"))))
dx <- rbind(m0[, xk], x[x$status == "ok", xk])
save(do.call(rbind, lapply(split(dx, dx$xnoise), function(z) data.frame(xnoise = z$xnoise[1], studies = nrow(z),
  t(sapply(ty, pp, z = z)), unt_share = sum(z$unt) / sum(z$SA)))), "condition_errors.csv")
xt <- merge(xt, unique(x[, c("condition_id", "xnoise")]))
tt0 <- rbind(transform(tm[tm$family == "structured" & tm$noise == 0 & tm$incl.cut == 0.8, c("solution", "sufficient", "raw_cov", "rem_share")], xnoise = 0),
             xt[, c("solution", "sufficient", "raw_cov", "rem_share", "xnoise")])
save(do.call(rbind, lapply(split(tt0, tt0$xnoise), function(z) { p <- z[z$solution == "pars" & !is.na(z$raw_cov), ]
  data.frame(xnoise = z$xnoise[1], auc_rem = auc(-p$rem_share, p$sufficient), auc_raw = auc(p$raw_cov, p$sufficient)) })),
  "condition_errors_terms.csv")

# ---- 7. six published data sets (Table 4, Supplementary Table S5) ----
sets <- list(
  Krook = list(file = "Krook.tab", outcome = "WNP", conds = c("ES", "QU", "WS", "WM", "LP")),
  Hino = list(file = "HC.tab", outcome = "VOTE", conds = c("FOREIGN", "UNEMP", "CONV", "PRES80")),
  Lipset = list(file = "LR.tab", outcome = "SURV", conds = c("DEV", "URB", "LIT", "IND", "STB"),
                prep = function(d) data.frame(   # crisp thresholds of Rihoux and De Meur (2009)
                  DEV = as.integer(d$DEV >= 550), URB = as.integer(d$URB >= 50), LIT = as.integer(d$LIT >= 75),
                  IND = as.integer(d$IND >= 30), STB = as.integer(d$STB < 10), SURV = as.integer(d$SURV > 0))),
  SchneiderSadowski = list(file = "d.education.tab", outcome = "O", conds = paste0("NPM", 1:6)),
  Cragun = list(file = "d.tumorscreen.tab", outcome = "HPF", conds = c("CA", "AR", "RR", "DR", "DC", "PR")),
  Borradaile = list(file = "borradaile2021.tab", outcome = "ENCR", conds = c("EASE", "MOTI", "MOTG", "COMMS", "RISK")))
illus <- function(nm, s, ic) {
  dd <- read.table(file.path("data", s$file), header = TRUE); if (!is.null(s$prep)) dd <- s$prep(dd)
  n <- length(s$conds); cases <- data.frame(row = as.integer(as.matrix(dd[, s$conds]) %*% (2L^((n - 1):0))), y = dd[[s$outcome]])
  tt <- truth_from_cases(cases, n, ic)
  sol <- minimize_sop(tt$T, n, dc = tt$R, names = s$conds); SA <- intersect(covered_rows(sol$terms, n), tt$R)
  csa <- SA %in% (if (length(tt$Tneg)) parsimonious_assumptions(tt$Tneg, tt$R, n) else integer(0))
  nt <- necessity_table(cases, n); nec <- nt[!is.na(nt$inclN) & nt$inclN >= 0.9 & !is.na(nt$RoN) & nt$RoN >= 0.6, ]
  unt_R <- untenable_rows(tt$R, nec, n)
  solc <- minimize_sop(tt$T, n, dc = setdiff(tt$R, c(tt$R[unt_R], SA[csa])), names = s$conds)
  pos <- cases$row[cases$y == 1]
  terms <- function(tr, lab) do.call(rbind, lapply(seq_len(nrow(tr)), function(t) { rows <- cube_rows(tr$value[t], tr$mask[t], n)
    data.frame(dataset = nm, incl.cut = ic, solution = lab, term = format_sop(tr[t, , drop = FALSE], n, s$conds),
               support_rows = sum(rows %in% tt$T), raw_cov = round(mean(pos %in% rows), 3), rem_share = round(mean(rows %in% tt$R), 3)) }))
  list(summary = data.frame(dataset = nm, incl.cut = ic, cases = nrow(dd), parsimonious = sol$expression, SA = length(SA),
         CSA = sum(csa), untenable = sum(SA %in% tt$R[unt_R]), distance2plus = sum(vapply(SA, min_dist, 0, T = tt$T) >= 2),
         cleaned = solc$expression),
       terms = rbind(terms(sol$terms, "parsimonious"), terms(solc$terms, "cleaned")))
}
il <- unlist(lapply(names(sets), function(nm) lapply(c(0.8, 1), function(ic) illus(nm, sets[[nm]], ic))), recursive = FALSE)
save(do.call(rbind, lapply(il, `[[`, "summary")), "published_summary.csv")
save(do.call(rbind, lapply(il, `[[`, "terms")), "published_terms.csv")

# ---- 8. write one file per table and figure of the paper, plus the numbers quoted in the text ----
r3 <- function(x) round(x, 3)
w <- function(x, f) write.csv(x, file.path(OUT, f), row.names = FALSE)
pn <- RES$pooled_by_noise; rc <- RES$pooled_by_family_noise; rc <- rc[rc$family == "null", ]
w(data.frame(noise = c(pn$noise, "random control"), all = r3(c(pn$prec_all, rc$prec_all)), csa = r3(c(pn$prec_csa, rc$prec_csa)),
             untenable = r3(c(pn$prec_unt, rc$prec_unt)), distance2 = r3(c(pn$prec_d2, rc$prec_d2))), "table2.csv")
w(transform(RES$hypotheses, estimate = r3(estimate), second = r3(second), share_gain_positive = r3(share_gain_positive)), "table3.csv")
ps <- RES$published_summary; pt <- RES$published_terms
ps <- ps[ps$incl.cut == 0.8, ]; p8 <- pt[pt$incl.cut == 0.8 & pt$solution == "parsimonious", ]
weak <- do.call(rbind, lapply(split(p8, p8$dataset), function(z) z[order(z$support_rows, -z$rem_share)[1], ]))
w(data.frame(dataset = ps$dataset, cases = ps$cases, SA = ps$SA, CSA = ps$CSA, untenable = ps$untenable,
             cleaning_changed = ps$parsimonious != ps$cleaned, weakest_term = weak$term[match(ps$dataset, weak$dataset)],
             weakest_rows = weak$support_rows[match(ps$dataset, weak$dataset)], weakest_raw_cov = weak$raw_cov[match(ps$dataset, weak$dataset)]), "table4.csv")
w(data.frame(noise = pn$noise, other_than_csa = r3(pn$prec_noncsa), tenable = r3(pn$prec_ten), distance1 = r3(pn$prec_d1), cleaned = r3(pn$prec_clean)), "tableS1.csv")
w(transform(RES$held_out, auc_raw = r3(auc_raw), auc_support_model = r3(auc_support_model)), "tableS2.csv")
w(data.frame(RES$cutoff_same_studies[, c("noise", "cut", "studies")], r3(RES$cutoff_same_studies[, ty])), "tableS3.csv")
ce <- RES$condition_errors; w(data.frame(condition_error = ce$xnoise, all = r3(ce$SA), csa = r3(ce$csa), untenable = r3(ce$unt), distance2 = r3(ce$d2)), "tableS4.csv")
w(pt[pt$incl.cut == 0.8, c("dataset", "solution", "term", "support_rows", "raw_cov", "rem_share")], "tableS5.csv")
cs8 <- RES$cells; cs8 <- cs8[cs8$family == "structured" & cs8$incl.cut == 0.8, ]
w(cs8[, c("condition_id", "shape", "n", "N", "dist", "noise", "h1", "h2")], "figure2.csv")
w(RES$terms_pars_vs_clean[, c("solution", "noise", "sufficient", "bt_correct")], "figure3.csv")
w(RES$terms_by_support, "figure4a.csv")
w(RES$h6_auc, "figure4b.csv")
# numbers quoted in the text (long format: section, quantity, group, value)
tx <- function(sec, q, g, v) data.frame(section = sec, quantity = q, group = as.character(g), value = r3(v))
br <- RES$base_rate; nn <- RES$necessity_consistency_only; sl <- RES$solution_level; rct <- RES$random_control_terms
rs <- RES$remainder_share_restricted; wl <- RES$within_levels; gc <- RES$generating_conjunctions; h5 <- RES$h5_terms
bn <- RES$pooled_by_N_noise; cet <- RES$condition_errors_terms
w(rbind(
  tx("3.1", "excluded studies (share)", "all", sum(d0_excluded) / d0_total),
  tx("3.1", "outcome among remainders (base rate)", br$noise, br$prev_remainders),
  tx("3.1", "CSA share of SAs", pn$noise, pn$csa_share),
  tx("3.1", "CSA precision by cases", paste0("N=", bn$N, ", noise=", bn$noise), bn$prec_csa),
  tx("3.1", "other SA precision by cases", paste0("N=", bn$N, ", noise=", bn$noise), bn$prec_noncsa),
  tx("3.2", "untenable share of SAs", pn$noise, pn$unt_share),
  tx("3.2", "untenable precision, consistency-only necessity", nn$noise, nn$prec_unt2),
  tx("3.2", "tenable precision, consistency-only necessity", nn$noise, nn$prec_ten2),
  tx("3.2", "untenable SAs, consistency-only vs main criterion (ratio)", nn$noise, nn$unt2 / nn$unt),
  tx("3.3", "share of whole solutions with all terms sufficient", c(paste("pars", sl$noise), paste("clean", sl$noise)), c(sl$pars_all_suff, sl$clean_all_suff)),
  tx("3.3", "share of whole solutions with all terms termwise B&T", c(paste("pars", sl$noise), paste("clean", sl$noise)), c(sl$pars_all_bt, sl$clean_all_bt)),
  tx("3.3", "random control: sufficient terms", rct$solution, rct$sufficient),
  tx("3.3", "random control: termwise B&T terms", rct$solution, rct$bt_correct),
  tx("3.3", "generating conjunctions: correct terms", gc$solution, gc$bt_correct),
  tx("3.3", "generating conjunctions: whole solutions", gc$solution, gc$solution_all),
  tx("3.4", "SAs at distance 2 or more (share)", pn$noise, pn$d2_share),
  tx("3.5", "H5: sufficient terms with 1-2 rows", h5$n, h5$low_suff),
  tx("3.5", "H5: sufficient terms with 5+ rows", h5$n, h5$high_suff),
  tx("3.5", "terms with no unobserved rows (share)", rs$noise, rs$share_rho0),
  tx("3.5", "AUC remainder share, terms with unobserved rows", rs$noise, rs$auc_rem_pos),
  tx("3.5", "AUC raw coverage, terms with unobserved rows", rs$noise, rs$auc_raw_pos),
  tx("3.6", "AUC remainder share within design levels (min, max)", rep(c("noise 0", "noise 0.10"), each = 2),
     c(range(wl$auc_rem[wl$noise == 0]), range(wl$auc_rem[wl$noise == 0.1]))),
  tx("3.7", "untenable share of SAs under condition errors", ce$xnoise, ce$unt_share),
  tx("3.7", "AUC remainder share under condition errors", cet$xnoise, cet$auc_rem),
  tx("3.7", "AUC raw coverage under condition errors", cet$xnoise, cet$auc_raw)), "text_numbers.csv")
print(hyp)
