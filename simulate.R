# simulate.R -- data generation, QCA analysis and scoring of every simulated study.
# Usage: Rscript simulate.R main [cores]        (440 cells x 500 studies -> raw/main/)
#        Rscript simulate.R extension [cores]   (192 cells x 500 studies -> raw/extension/)
# Truth-table rows are integers 0..2^n - 1 (first condition = most significant bit).

# ---- exact Boolean minimization (Quine-McCluskey prime implicants + exact cover) ----
# An implicant is a pair (value, mask): bits set in mask are eliminated; value is zero there.

.bit <- function(i, n) bitwShiftL(1L, n - i)   # bit for condition i (1-based)

# All prime implicants of the function whose on-set plus don't-care set is `pts`.
prime_implicants <- function(pts, n) {
  pts <- sort(unique(as.integer(pts)))
  if (length(pts) == 0) return(data.frame(value = integer(0), mask = integer(0)))
  size <- bitwShiftL(1L, n)
  cur <- list("0" = pts)          # mask -> values present at that mask
  primes_v <- integer(0); primes_m <- integer(0)
  while (length(cur) > 0) {
    nxt <- list()
    for (mk in names(cur)) {
      mask <- as.integer(mk)
      vals <- cur[[mk]]
      present <- logical(size); present[vals + 1L] <- TRUE
      used <- logical(size)
      for (b in 0:(n - 1)) {
        bit <- bitwShiftL(1L, b)
        if (bitwAnd(mask, bit) != 0L) next
        lo <- vals[bitwAnd(vals, bit) == 0L]
        hi <- lo + bit
        ok <- present[hi + 1L]
        if (any(ok)) {
          lo <- lo[ok]; hi <- hi[ok]
          used[lo + 1L] <- TRUE; used[hi + 1L] <- TRUE
          key <- as.character(bitwOr(mask, bit))
          nxt[[key]] <- c(nxt[[key]], lo)
        }
      }
      keep <- vals[!used[vals + 1L]]
      primes_v <- c(primes_v, keep); primes_m <- c(primes_m, rep(mask, length(keep)))
    }
    cur <- lapply(nxt, function(v) sort(unique(v)))
  }
  data.frame(value = primes_v, mask = primes_m)
}

.literals <- function(mask, n) {
  n - vapply(mask, function(m) sum(bitwAnd(m, bitwShiftL(1L, 0:(n - 1))) != 0L), 0)
}

# Exact minimum cover: fewest terms first, then fewest literals.
.exact_cover <- function(cov, cost) {
  # cov: logical matrix, rows = prime implicants, cols = minterms to cover
  best <- list(k = Inf, lit = Inf, sel = integer(0))
  rec <- function(rows_ok, cols_left, sel, k, lit) {
    # Reductions: essential rows, dominated columns, dominated rows.
    repeat {
      changed <- FALSE
      cl <- which(cols_left)
      if (length(cl) == 0) break
      sub <- cov[, cl, drop = FALSE] & rows_ok
      cnt <- colSums(sub)
      if (any(cnt == 0)) return(invisible())
      ess <- unique(apply(sub[, cnt == 1, drop = FALSE], 2, which))
      if (length(ess)) {
        for (r in ess) {
          sel <- c(sel, r); k <- k + 1L; lit <- lit + cost[r]
          cols_left <- cols_left & !cov[r, ]
        }
        rows_ok[ess] <- FALSE
        if (k > best$k) return(invisible())
        next
      }
      # column dominance: if rows(i) contain rows(j), covering j covers i
      if (length(cl) > 1) {
        inter <- crossprod(sub)            # |rows(i) & rows(j)|
        dom <- inter == matrix(cnt, length(cl), length(cl), byrow = TRUE) # rows(j) subset of rows(i)?
        # dom[i, j]: rows(j) within rows(i); drop i if some other j qualifies
        diag(dom) <- FALSE
        tie <- dom & t(dom)
        dom[tie & upper.tri(dom)] <- FALSE   # identical columns: keep one
        drop_col <- rowSums(dom) > 0
        if (any(drop_col)) { cols_left[cl[drop_col]] <- FALSE; changed <- TRUE }
      }
      # row dominance: row a dropped if row b covers a superset at no more cost
      ro <- which(rows_ok)
      if (length(ro) > 1) {
        cl <- which(cols_left)
        rs <- cov[ro, cl, drop = FALSE]
        rsize <- rowSums(rs)
        if (any(rsize == 0)) { rows_ok[ro[rsize == 0]] <- FALSE; changed <- TRUE; next }
        inter <- tcrossprod(rs)             # |cols(a) & cols(b)|
        sub_ab <- inter == matrix(rsize, length(ro), length(ro))  # cols(a) within cols(b)
        cheaper <- outer(cost[ro], cost[ro], ">=")                # cost(a) >= cost(b)
        d <- sub_ab & cheaper; diag(d) <- FALSE
        eq <- d & t(d); d[eq & lower.tri(d)] <- FALSE              # identical rows: keep one
        drop_row <- rowSums(d) > 0
        if (any(drop_row)) { rows_ok[ro[drop_row]] <- FALSE; changed <- TRUE }
      }
      if (!changed) break
    }
    if (!any(cols_left)) {
      if (k < best$k || (k == best$k && lit < best$lit))
        best <<- list(k = k, lit = lit, sel = sel)
      return(invisible())
    }
    # lower bound: greedy set of columns that share no covering row
    lb <- 0L; blocked <- logical(nrow(cov))
    for (j in which(cols_left)) {
      r <- cov[, j] & rows_ok
      if (!any(r & blocked)) { lb <- lb + 1L; blocked <- blocked | r }
    }
    if (k + lb > best$k) return(invisible())
    # branch on the column with the fewest candidate rows
    cl <- which(cols_left)
    counts <- colSums(cov[rows_ok, cl, drop = FALSE])
    if (any(counts == 0)) return(invisible())
    j <- cl[which.min(counts)]
    cand <- which(cov[, j] & rows_ok)
    cand <- cand[order(-colSums(t(cov[cand, , drop = FALSE]) & cols_left), cost[cand])]
    for (r in cand) {
      rows_ok[r] <- FALSE   # covers using r are all explored in this branch
      rec(rows_ok, cols_left & !cov[r, ], c(sel, r), k + 1L, lit + cost[r])
    }
    invisible()
  }
  rec(rep(TRUE, nrow(cov)), rep(TRUE, ncol(cov)), integer(0), 0L, 0)
  best
}

#' Minimal sum-of-products for a Boolean function
#'
#' Finds a minimal sum-of-products expression (fewest terms, then fewest
#' literals) covering every on-set minterm, optionally using don't-care
#' minterms.
#'
#' @param on Integer vector of on-set minterms (rows with outcome 1).
#' @param n Number of conditions (variables).
#' @param dc Integer vector of don't-care minterms (may be covered or not).
#' @param names Optional condition names, used to print the expression.
#' @return A list with `terms` (data frame of value/mask), `n_terms`,
#'   `n_literals` and `expression`.
#' @examples
#' minimize_sop(c(3, 6, 7, 15), n = 4)
#' minimize_sop(c(3, 6, 7, 15), n = 4, dc = 2)
#' @export
minimize_sop <- function(on, n, dc = integer(0), names = NULL) {
  on <- sort(unique(as.integer(on))); dc <- setdiff(as.integer(dc), on)
  if (is.null(names)) names <- LETTERS[seq_len(n)]
  if (length(on) == 0)
    return(structure(list(terms = data.frame(value = integer(0), mask = integer(0)),
                          n_terms = 0L, n_literals = 0L, expression = "0"),
                     class = "sop"))
  pis <- prime_implicants(c(on, dc), n)
  cov <- outer(seq_len(nrow(pis)), on, function(i, m)
    bitwAnd(m, bitwNot(pis$mask[i])) == pis$value[i])
  cost <- .literals(pis$mask, n)
  # drop primes that cover no on-set minterm
  keep <- rowSums(cov) > 0
  pis <- pis[keep, , drop = FALSE]; cov <- cov[keep, , drop = FALSE]; cost <- cost[keep]
  best <- .exact_cover(cov, cost)
  terms <- pis[best$sel, , drop = FALSE]
  terms <- terms[order(-cost[best$sel], terms$value), , drop = FALSE]
  rownames(terms) <- NULL
  structure(list(terms = terms, n_terms = nrow(terms),
                 n_literals = as.integer(sum(.literals(terms$mask, n))),
                 expression = format_sop(terms, n, names)),
            class = "sop")
}

#' Write implicants as a readable expression
#'
#' Present conditions are written as the name, absent ones with a tilde,
#' joined by `*`, and terms joined by `+` (the usual QCA notation).
#' @param terms Data frame with columns `value` and `mask`.
#' @param n Number of conditions.
#' @param names Condition names.
#' @export
format_sop <- function(terms, n, names = LETTERS[seq_len(n)]) {
  if (nrow(terms) == 0) return("0")
  out <- vapply(seq_len(nrow(terms)), function(t) {
    parts <- character(0)
    for (i in seq_len(n)) {
      b <- .bit(i, n)
      if (bitwAnd(terms$mask[t], b) != 0L) next
      parts <- c(parts, if (bitwAnd(terms$value[t], b) != 0L) names[i] else paste0("~", names[i]))
    }
    if (length(parts) == 0) "1" else paste(parts, collapse = "*")
  }, "")
  paste(out, collapse = " + ")
}

#' @export

covered_rows <- function(terms, n) {
  m <- 0:(2^n - 1); out <- logical(2^n)
  for (t in seq_len(nrow(terms))) out <- out | (bitwAnd(m, bitwNot(terms$mask[t])) == terms$value[t])
  which(out) - 1L
}
parsimonious_assumptions <- function(T, R, n) {
  if (length(T) == 0 || length(R) == 0) return(integer(0))
  intersect(covered_rows(minimize_sop(T, n, dc = R)$terms, n), R)
}

# ---- true functions ----
.bits_of <- function(n) { m <- 0:(2^n - 1)
  sapply(seq_len(n), function(i) as.integer(bitwAnd(m, bitwShiftL(1L, n - i)) != 0L)) }
conj_rows <- function(vars, vals, X) {
  hit <- rep(TRUE, nrow(X))
  for (j in seq_along(vars)) hit <- hit & (X[, vars[j]] == vals[j])
  hit
}
# D1-D3: disjunctions of random conjunctions, redrawn until prevalence is 0.25-0.35. D4: at least 3 of 4.
# Returns f (logical, 2^n) and gen, the generating conjunctions as (value, mask) implicants.
gen_departure <- function(shape, n, max_tries = 10000) {
  X <- .bits_of(n); full <- 2^n - 1
  cube <- function(vars, vals) { b <- bitwShiftL(1L, n - vars); data.frame(value = sum(b * vals), mask = bitwAnd(full, bitwNot(sum(b)))) }
  for (try in seq_len(max_tries)) {
    if (shape == "D4") {
      vars <- sample(n, 4); f <- rowSums(X[, vars, drop = FALSE]) >= 3
      return(list(f = f, gen = do.call(rbind, lapply(combn(vars, 3, simplify = FALSE), cube, vals = rep(1L, 3)))))
    }
    lens <- switch(shape, D1 = 2L, D2 = sample(2:4, 2, TRUE), D3 = sample(3:min(5, n), 4, TRUE))
    f <- logical(2^n); gen <- NULL
    for (len in lens) {
      vars <- sample(n, len); vals <- sample(0:1, len, TRUE)
      f <- f | conj_rows(vars, vals, X); gen <- rbind(gen, cube(vars, vals))
    }
    if (mean(f) >= 0.25 && mean(f) <= 0.35) return(list(f = f, gen = gen))
  }
  stop("Could not reach the prevalence band for ", shape, " with n = ", n, ".")
}
gen_null <- function(n) list(f = runif(2^n) < 0.5, gen = NULL)   # N1: random control

# ---- cases, truth table, checks ----
XNOISE <- 0   # probability that each condition is recorded wrongly (extension only)
gen_cases <- function(f, n, N, dist = c("uniform", "clustered"), noise = 0, beta = 0.7) {
  dist <- match.arg(dist)
  w <- if (dist == "uniform") rep(1, 2^n) else exp(beta * (rowSums(.bits_of(n)) - n / 2))
  r <- sample(0:(2^n - 1), N, replace = TRUE, prob = w)
  y <- xor(f[r + 1L], runif(N) < noise)
  if (XNOISE > 0) r <- bitwXor(r, as.integer(matrix(runif(N * n) < XNOISE, N, n) %*% (2L^((n - 1):0))))
  data.frame(row = r, y = as.integer(y))
}
truth_from_cases <- function(cases, n, incl.cut = 0.8) {
  cnt <- tabulate(cases$row + 1L, 2^n); pos <- tabulate(cases$row[cases$y == 1] + 1L, 2^n)
  obs <- which(cnt > 0) - 1L; incl <- pos[obs + 1L] / cnt[obs + 1L]
  list(T = obs[incl >= incl.cut - 1e-12], F = obs[incl < incl.cut - 1e-12],
       Tneg = obs[(1 - incl) >= incl.cut - 1e-12], R = setdiff(0:(2^n - 1), obs))
}
# inclN = sum(min(X,Y)) / sum(Y); RoN = sum(1 - X) / sum(1 - min(X,Y)), for every literal (condition k at value v).
necessity_table <- function(cases, n) {
  bits <- bitwShiftL(1L, (n - 1):0); Y <- cases$y
  out <- expand.grid(k = seq_len(n), v = 0:1); out$inclN <- out$RoN <- NA_real_
  for (i in seq_len(nrow(out))) {
    X <- as.integer((bitwAnd(cases$row, bits[out$k[i]]) != 0L) == (out$v[i] == 1))
    out$inclN[i] <- if (sum(Y)) sum(pmin(X, Y)) / sum(Y) else NA
    d <- sum(1 - pmin(X, Y)); out$RoN[i] <- if (d > 0) sum(1 - X) / d else NA
  }
  out
}
untenable_rows <- function(R, nec, n) {
  if (!length(R) || !nrow(nec)) return(rep(FALSE, length(R)))
  bits <- bitwShiftL(1L, (n - 1):0)
  vapply(R, function(r) any(vapply(seq_len(nrow(nec)), function(i) as.integer(bitwAnd(r, bits[nec$k[i]]) != 0L) != nec$v[i], TRUE)), TRUE)
}
cube_rows <- function(v, m, n) { r <- 0:(2^n - 1); r[bitwAnd(r, bitwNot(m)) == v] }
# termwise B&T criterion: the term's literals are a subset of those of some reference implicant.
bt_correct <- function(tv, tm, ref) any(bitwAnd(ref$mask, bitwNot(tm)) == 0L & bitwAnd(ref$value, bitwNot(tm)) == tv)
min_dist <- function(r, T) min(vapply(T, function(t) sum(as.integer(intToBits(bitwXor(r, t)))[1:16]), 0))

# ---- design ----
conditions_main <- function() {
  s <- expand.grid(shape = paste0("D", 1:4), n = 4:7, N = c(20, 40, 80), dist = c("uniform", "clustered"),
                   noise = c(0, 0.05, 0.10, 0.20), incl.cut = 0.8, stringsAsFactors = FALSE)
  c0 <- expand.grid(shape = "N1", n = 4:7, N = c(20, 40, 80), dist = c("uniform", "clustered"),
                    noise = 0, incl.cut = 0.8, stringsAsFactors = FALSE)
  sens <- expand.grid(shape = c("D2", "D3"), n = 5:6, N = 40, dist = "uniform",
                      noise = c(0, 0.05, 0.10, 0.20), incl.cut = c(1, 0.75), stringsAsFactors = FALSE)
  g <- rbind(s, c0, sens); g$condition_id <- seq_len(nrow(g)); g
}
conditions_extension <- function() {
  g <- expand.grid(shape = paste0("D", 1:4), n = 4:7, N = c(20, 40, 80), dist = c("uniform", "clustered"),
                   noise = 0, xnoise = c(0.05, 0.10), incl.cut = 0.8, stringsAsFactors = FALSE)
  g$condition_id <- seq_len(nrow(g)); g
}

# ---- one simulated study ----
# ref = "pi" scores terms against the prime implicants of the true function; "gen" against its generating conjunctions.
run_rep <- function(cond, rep, ref = "pi") {
  n <- cond$n
  tr <- if (cond$shape == "N1") gen_null(n) else gen_departure(cond$shape, n)
  cases <- gen_cases(tr$f, n, cond$N, cond$dist, cond$noise)
  tt <- truth_from_cases(cases, n, cond$incl.cut)
  ds <- data.frame(rep = rep, n_T = length(tt$T), n_F = length(tt$F), n_R = length(tt$R), prevalence = mean(tr$f), status = "ok")
  if (!length(tt$T) || !length(tt$R)) { ds$status <- "degenerate"; return(list(ds = ds, terms = NULL)) }
  sol <- minimize_sop(tt$T, n, dc = tt$R)
  SA <- intersect(covered_rows(sol$terms, n), tt$R)
  csa <- SA %in% (if (length(tt$Tneg)) parsimonious_assumptions(tt$Tneg, tt$R, n) else integer(0))
  nt <- necessity_table(cases, n)
  nec_real <- nt[!is.na(nt$inclN) & nt$inclN >= 0.9 & !is.na(nt$RoN) & nt$RoN >= 0.6, ]
  nec_cons <- nt[!is.na(nt$inclN) & nt$inclN >= 0.9, ]          # consistency only (Section 3.2)
  unt_R <- untenable_rows(tt$R, nec_real, n); unt <- SA %in% tt$R[unt_R]; unt2 <- SA %in% tt$R[untenable_rows(tt$R, nec_cons, n)]
  dist <- vapply(SA, min_dist, 0, T = tt$T)
  truly <- tr$f[SA + 1L]
  solc <- minimize_sop(tt$T, n, dc = setdiff(tt$R, c(tt$R[unt_R], SA[csa])))   # cleaned solution
  SAc <- intersect(covered_rows(solc$terms, n), tt$R)
  refs <- if (ref == "gen") tr$gen else if (any(tr$f)) prime_implicants(which(tr$f) - 1L, n) else data.frame(value = integer(0), mask = integer(0))
  pos_cases <- cases$row[cases$y == 1]
  term_rec <- function(terms, label) {
    if (!nrow(terms)) return(NULL)
    do.call(rbind, lapply(seq_len(nrow(terms)), function(t) {
      rows <- cube_rows(terms$value[t], terms$mask[t], n)
      data.frame(rep = rep, solution = label, term = t, size = length(rows),
                 support_rows = sum(rows %in% tt$T), raw_cov = if (length(pos_cases)) mean(pos_cases %in% rows) else NA,
                 rem_share = mean(rows %in% tt$R), sufficient = all(tr$f[rows + 1L]),
                 bt_correct = bt_correct(terms$value[t], terms$mask[t], refs))
    }))
  }
  tp <- term_rec(sol$terms, "pars"); tc <- term_rec(solc$terms, "clean")
  s <- function(sel) c(sum(sel), sum(truly[sel]))
  ds <- cbind(ds, data.frame(
    SA = length(SA), SA_ok = sum(truly), csa = s(csa)[1], csa_ok = s(csa)[2], noncsa = s(!csa)[1], noncsa_ok = s(!csa)[2],
    unt = s(unt)[1], unt_ok = s(unt)[2], ten = s(!unt)[1], ten_ok = s(!unt)[2],
    unt2 = s(unt2)[1], unt2_ok = s(unt2)[2], ten2 = s(!unt2)[1], ten2_ok = s(!unt2)[2],
    d1 = s(dist == 1)[1], d1_ok = s(dist == 1)[2], d2 = s(dist >= 2)[1], d2_ok = s(dist >= 2)[2],
    SAc = length(SAc), SAc_ok = sum(tr$f[SAc + 1L]),
    pars_all_suff = all(tp$sufficient), pars_all_bt = all(tp$bt_correct),
    clean_all_suff = if (is.null(tc)) NA else all(tc$sufficient), clean_all_bt = if (is.null(tc)) NA else all(tc$bt_correct)))
  list(ds = ds, terms = rbind(tp, tc))
}

# ---- one design cell: 500 studies, each with its own L'Ecuyer-CMRG stream ----
run_cell <- function(cond, seed, reps = 500, ...) {
  XNOISE <<- if (is.null(cond$xnoise)) 0 else cond$xnoise
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed); s <- .Random.seed
  res <- vector("list", reps)
  for (r in seq_len(reps)) {
    assign(".Random.seed", s, envir = .GlobalEnv); s <- parallel::nextRNGStream(s)
    res[[r]] <- run_rep(cond, r, ...)
  }
  cols <- unique(unlist(lapply(res, function(x) names(x$ds))))
  ds <- do.call(rbind, lapply(res, function(x) { for (cc in setdiff(cols, names(x$ds))) x$ds[[cc]] <- NA; x$ds[cols] }))
  tm <- do.call(rbind, lapply(res, `[[`, "terms"))
  list(ds = cbind(cond[rep(1, nrow(ds)), ], ds), terms = if (is.null(tm)) NULL else cbind(condition_id = cond$condition_id, tm))
}

if (sys.nframe() == 0) {
  a <- commandArgs(TRUE); which <- a[1]
  cores <- if (length(a) > 1) as.integer(a[2]) else max(1, parallel::detectCores() - 1)
  g <- if (which == "main") conditions_main() else conditions_extension()
  base <- if (which == "main") 60261110 else 70261110
  out <- file.path("raw", which); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  todo <- g$condition_id[!file.exists(file.path(out, sprintf("cell_%03d.csv", g$condition_id)))]
  invisible(parallel::mclapply(todo, function(id) {
    x <- run_cell(g[g$condition_id == id, ], base + id)
    write.csv(x$ds, file.path(out, sprintf("cell_%03d.csv", id)), row.names = FALSE)
    if (!is.null(x$terms)) write.csv(x$terms, file.path(out, sprintf("terms_%03d.csv", id)), row.names = FALSE)
  }, mc.cores = cores, mc.preschedule = FALSE))
}
