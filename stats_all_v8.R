## =====================================================================
## Kirs et al. — Hawaiian streams: ALL statistical analyses
## Reproduces the tests behind Sections 3.2–3.6, Table 4, and Tables S6–S10.
##
## INPUT FILES (place in the working directory):
##   data-final-5d81e002.xlsx        (4 header rows, then 80 data rows)
##   landuse-clean-c615a1e2.xlsx      (NOAA C-CAP land cover, % per watershed)
##
## REQUIREMENTS:  install.packages("readxl")     # base R is used otherwise
##
## DESIGN NOTES
##   * 80 samples = 20 streams x 4 temporal replicates. Streams are the
##     clustering/independent unit (G = 20). All model-based inference uses
##     cluster-robust (sandwich) SEs clustered on stream with a t(G-1)=t(19)
##     reference.
##   * Concentrations are log10-transformed. Left-censored (below-detection)
##     values are handled by Tobit MLE for analytes with non-detects; E. coli
##     and enterococci (no non-detects) use OLS.
##   * Watershed variables (land cover, cesspool density) are constant within a
##     stream and are therefore tested at the STREAM level (n = 20) by Spearman
##     rank correlation, to avoid pseudoreplication.
##   * Land cover is from the C-CAP file (NOT the older columns embedded in
##     data-final).
##
## The script prints every result to the console; it does not write files.
## =====================================================================

library(readxl)
set.seed(2025)

DATA <- "data-final-5d81e002.xlsx"
LU   <- "landuse-clean-c615a1e2.xlsx"

## ------------------------------------------------------------------ ##
## 1. LOAD DATA                                                        ##
## ------------------------------------------------------------------ ##
raw <- as.data.frame(read_excel(DATA, col_names = FALSE, .name_repair = "minimal"))
raw <- raw[-(1:4), ]                       # drop 4 grouped header rows -> 80 rows
gv  <- function(j) raw[[j]]                # get raw column by 1-based index

## 1-based column positions in data-final
C_SITE <- 1; C_ISLAND <- 2
C_TEMP <- 8; C_SAL <- 9; C_COND <- 10; C_DO <- 11; C_PH <- 13; C_TURB <- 14
C_RAIN3 <- 15; C_RAIN6 <- 16; C_RAIN12 <- 17; C_RAIN24 <- 18
C_ECOLI <- 20; C_ENT <- 21; C_CPERF <- 22; C_FPLUS <- 23; C_SOMATIC <- 24
C_SAL_P <- 25; C_CAMP <- 26; C_CJE <- 27
C_CRASS_PA <- 35; C_HF183_PA <- 36
C_CRASS_C <- 39; C_HF183_C <- 40
C_AREA <- 46; C_NCESS <- 47

site   <- as.integer(factor(as.character(gv(C_SITE))))            # 20 stream clusters
stream_id <- as.character(gv(C_SITE))
island <- factor(trimws(as.character(gv(C_ISLAND))),
                 levels = c("Oahu", "Maui", "Kauai", "Hawaii"))
num  <- function(j) suppressWarnings(as.numeric(as.character(gv(j))))
zsc  <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)

## ---- parsers ----
parse_conc <- function(x) {                # ">x"->x ; "<x"->x/2 ; blank/ND->NA
  if (is.na(x)) return(NA_real_)
  s <- trimws(as.character(x))
  if (s == "" || tolower(s) %in% c("na","nan","nd","none")) return(NA_real_)
  if (grepl("^>", s)) return(as.numeric(gsub("[^0-9.]", "", s)))
  if (grepl("^<", s)) return(as.numeric(gsub("[^0-9.]", "", s)) / 2)
  suppressWarnings(as.numeric(gsub("[^0-9.]", "", s)))
}
log10_sub <- function(x, floor = 1) {      # log10 with DL/2 for '<' and 0 (predictor use)
  v <- parse_conc(x); if (is.na(v)) return(NA_real_)
  if (v > 0) log10(v) else log10(floor / 2)
}
detected_any_volume <- function(x) {       # pathogen positive at any filtered volume
  if (is.na(x)) return(0L)
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("na","nan","nd","none","","0")) return(0L)
  as.integer(grepl("[1-9]", s))
}
is_pos <- function(x) if (is.na(x)) 0L else as.integer(tolower(trimws(as.character(x))) == "pos")

## ---- C-CAP land cover, mapped to stream id ----
lu_raw <- as.data.frame(read_excel(LU, .name_repair = "minimal"))
lu_classes <- names(lu_raw)[-1]                              # 9 C-CAP classes (data order)
lu_name2id <- c(Heeia="1",Waiahole="2",Punaluu="3",Kaluanui="4",Anahulu="5",
                Paukauila="6",Kiikii="7",Moanalua="8",Keehi="9",Honokohau="M1",
                Honolua="M2",Waiehu="M3",Iao="M4",Maliko="M5",Waikomo="K1",
                Niumalu="K2",Nawiliwili="K3",Hanalei="K4",Wailuku="H1",Honolii="H2")
lu_by_id <- setNames(lapply(seq_len(nrow(lu_raw)),
                            function(i) as.numeric(lu_raw[i, -1])),
                     lu_name2id[as.character(lu_raw[[1]])])
lc_col <- function(cls) sapply(stream_id, function(s) lu_by_id[[s]][match(cls, lu_classes)])

## ------------------------------------------------------------------ ##
## 2. ESTIMATORS (cluster-robust)                                      ##
## ------------------------------------------------------------------ ##
tp <- function(b, se, G) 2 * pt(-abs(b / se), df = G - 1)     # two-sided t(G-1) p-value

## ---- OLS with cluster-robust SEs ----
ols_cr <- function(X, y, cl) {
  ok <- complete.cases(X, y); X <- X[ok, , drop = FALSE]; y <- y[ok]; cl <- cl[ok]
  n <- nrow(X); k <- ncol(X); XtXi <- solve(crossprod(X))
  b <- XtXi %*% crossprod(X, y); u <- as.vector(y - X %*% b)
  G <- length(unique(cl)); meat <- matrix(0, k, k)
  for (c in unique(cl)) { s <- crossprod(X[cl == c, , drop = FALSE], u[cl == c]); meat <- meat + tcrossprod(s) }
  V <- (G/(G-1))*((n-1)/(n-k)) * (XtXi %*% meat %*% XtXi)
  se <- sqrt(diag(V)); list(beta = as.vector(b), se = se, p = tp(b, se, G), G = G, ncens = 0)
}
## ---- Tobit (left-censored Gaussian) MLE + cluster-robust SEs ----
tobit_nll <- function(par, X, y, L, cens) {
  k <- ncol(X); b <- par[1:k]; sig <- exp(par[k+1]); xb <- as.vector(X %*% b)
  ll <- numeric(length(y))
  u <- cens == 0; ll[u]  <- dnorm(y[u], xb[u], sig, log = TRUE)
  c2 <- cens == 1; ll[c2] <- pnorm((L[c2] - xb[c2]) / sig, log.p = TRUE)
  -sum(ll)
}
tobit_scores <- function(par, X, y, L, cens) {
  k <- ncol(X); b <- par[1:k]; sig <- exp(par[k+1]); xb <- as.vector(X %*% b)
  S <- matrix(0, nrow(X), k + 1)
  u <- cens == 0; r <- (y[u] - xb[u]) / sig
  S[u, 1:k] <- (r / sig) * X[u, , drop = FALSE]; S[u, k+1] <- r^2 - 1
  cc <- cens == 1; a <- (L[cc] - xb[cc]) / sig; lam <- dnorm(a) / pmax(pnorm(a), 1e-12)
  S[cc, 1:k] <- -(lam / sig) * X[cc, , drop = FALSE]; S[cc, k+1] <- -a * lam
  S
}
tobit_cr <- function(X, y, L, cens, cl) {
  ok <- !is.na(y); X <- X[ok,,drop=FALSE]; y <- y[ok]; L <- L[ok]; cens <- cens[ok]; cl <- cl[ok]
  n <- nrow(X); k <- ncol(X)
  start <- c(as.vector(solve(crossprod(X), crossprod(X, y))), log(sd(y)))
  fit <- optim(start, tobit_nll, X = X, y = y, L = L, cens = cens,
               method = "BFGS", hessian = TRUE, control = list(maxit = 500))
  Ainv <- solve(fit$hessian); S <- tobit_scores(fit$par, X, y, L, cens)
  G <- length(unique(cl)); meat <- matrix(0, k+1, k+1)
  for (c in unique(cl)) { s <- colSums(S[cl == c,,drop=FALSE]); meat <- meat + tcrossprod(s) }
  V <- (G/(G-1))*((n-1)/(n-(k+1))) * (Ainv %*% meat %*% Ainv)
  se <- sqrt(diag(V))[1:k]; b <- fit$par[1:k]
  list(beta = b, se = se, p = tp(b, se, G), G = G, ncens = sum(cens == 1))
}
## ---- ridge-stabilised logistic IRLS + cluster-robust SEs (+ optional score WCB) ----
irls_ridge <- function(X, y, l2 = 1e-4, it = 200) {
  k <- ncol(X); b <- rep(0, k)
  for (i in seq_len(it)) {
    eta <- pmin(pmax(as.vector(X %*% b), -30), 30); p <- 1/(1+exp(-eta)); W <- p*(1-p)
    H <- -(t(X) %*% (X * W)) - l2 * diag(k); g <- t(X) %*% (y - p) - l2 * b
    bn <- b - solve(H, g); if (max(abs(bn - b)) < 1e-9) { b <- bn; break }; b <- bn
  }
  as.vector(b)
}
logit_cr_wcb <- function(x, y, cl, B = 0) {          # univariate: y ~ z(x); B>0 adds score WCB
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; cl <- cl[ok]
  z <- (x - mean(x)) / sd(x); X <- cbind(1, z); n <- nrow(X); k <- 2
  b <- irls_ridge(X, y); p <- 1/(1+exp(-pmin(pmax(as.vector(X %*% b), -30), 30))); W <- p*(1-p)
  Ainv <- solve(t(X) %*% (X * W) + 1e-4 * diag(k)); g <- (y - p) * X
  cls <- sort(unique(cl)); G <- length(cls); idx <- lapply(cls, function(c) which(cl == c))
  corr <- (G/(G-1))*((n-1)/(n-k)); meat <- matrix(0, k, k)
  for (ii in idx) { s <- colSums(g[ii,,drop=FALSE]); meat <- meat + tcrossprod(s) }
  V <- corr * (Ainv %*% meat %*% Ainv); se <- sqrt(V[2,2]); t_obs <- b[2] / se
  out <- c(beta = b[2], p_t = tp(b[2], se, G))
  if (B > 0) {                                        # score-based wild cluster bootstrap-t
    tstar <- numeric(B)
    for (bi in seq_len(B)) {
      w <- sample(c(-1,1), G, replace = TRUE); gw <- g
      for (j in seq_along(idx)) gw[idx[[j]], ] <- g[idx[[j]],,drop=FALSE] * w[j]
      bd <- Ainv %*% colSums(gw); ms <- matrix(0,k,k)
      for (ii in idx) { s <- colSums(gw[ii,,drop=FALSE]); ms <- ms + tcrossprod(s) }
      Vs <- corr * (Ainv %*% ms %*% Ainv); ses <- sqrt(Vs[2,2]); tstar[bi] <- if (ses>0) bd[2]/ses else 0
    }
    out <- c(out, p_wcb = (sum(abs(tstar) >= abs(t_obs)) + 1) / (B + 1))
  }
  out
}
## ---- rank helpers, permutation Spearman, permutation Kruskal-Wallis, BH-FDR ----
rank_avg <- function(a) { r <- rank(a, ties.method = "average"); r }
spearman_perm <- function(x, y, nperm = 20000) {
  ok <- !is.na(x) & !is.na(y); rx <- rank_avg(x[ok]); ry <- rank_avg(y[ok])
  rho <- cor(rx, ry); cnt <- sum(replicate(nperm, abs(cor(rx, sample(ry))) >= abs(rho) - 1e-12))
  c(rho = rho, p = (cnt + 1) / (nperm + 1))
}
kw_perm <- function(y, g, nperm = 10000) {
  ok <- !is.na(y); y <- y[ok]; g <- droplevels(g[ok])
  H <- unname(kruskal.test(y, g)$statistic); cnt <- 1L
  for (i in seq_len(nperm)) if (unname(kruskal.test(y, sample(g))$statistic) >= H - 1e-9) cnt <- cnt + 1L
  c(H = H, p_perm = cnt / (nperm + 1))
}
bh <- function(p) p.adjust(p, "BH")

## ------------------------------------------------------------------ ##
## 3. BUILD ANALYSIS VARIABLES                                         ##
## ------------------------------------------------------------------ ##
conc <- list(                              # median-usable parsed concentrations
  "E. coli"=vapply(gv(C_ECOLI),parse_conc,numeric(1)),
  "Enterococci"=vapply(gv(C_ENT),parse_conc,numeric(1)),
  "C. perfringens"=vapply(gv(C_CPERF),parse_conc,numeric(1)),
  "F+ coliphage"=vapply(gv(C_FPLUS),parse_conc,numeric(1)),
  "Somatic coliphage"=vapply(gv(C_SOMATIC),parse_conc,numeric(1)))
logind <- list(                            # log10 (DL/2) predictor form of each indicator/marker
  "E. coli"=vapply(gv(C_ECOLI),log10_sub,numeric(1)),
  "Enterococci"=vapply(gv(C_ENT),log10_sub,numeric(1)),
  "C. perfringens"=vapply(gv(C_CPERF),log10_sub,numeric(1)),
  "F+ coliphage"=vapply(gv(C_FPLUS),log10_sub,numeric(1)),
  "Somatic coliphage"=vapply(gv(C_SOMATIC),log10_sub,numeric(1)),
  "crAssphage"=vapply(gv(C_CRASS_C),log10_sub,numeric(1)),
  "HF183"=vapply(gv(C_HF183_C),log10_sub,numeric(1)))
patho <- list(
  "Salmonella spp."=vapply(gv(C_SAL_P),detected_any_volume,integer(1)),
  "Campylobacter spp."=vapply(gv(C_CAMP),detected_any_volume,integer(1)),
  "C. jejuni"=vapply(gv(C_CJE),detected_any_volume,integer(1)))
marker_pa <- list(
  "crAssphage"=vapply(gv(C_CRASS_PA),is_pos,integer(1)),
  "HF183"=vapply(gv(C_HF183_PA),is_pos,integer(1)))
any_marker <- as.integer(marker_pa[["crAssphage"]] | marker_pa[["HF183"]])
cess_density <- num(C_NCESS) / num(C_AREA)

cat("Loaded", nrow(raw), "samples;", length(unique(site)), "streams;",
    "islands:", paste(table(island), collapse="/"), "\n")

## ------------------------------------------------------------------ ##
## 4. Section 3.2 — ISLAND COMPARISONS                                 ##
##    concentrations -> Kruskal-Wallis; detection -> Fisher's exact    ##
## ------------------------------------------------------------------ ##
cat("\n================ 3.2  ISLAND COMPARISONS ================\n")
cat("-- Kruskal-Wallis across islands (FIB/coliphage concentration) --\n")
for (nm in names(conc)) {
  r <- kw_perm(conc[[nm]], island)
  cat(sprintf("  %-18s H=%5.2f  perm p=%.4f\n", nm, r["H"], r["p_perm"]))
}
cat("-- Fisher's exact across islands (pathogen & marker detection) --\n")
for (nm in names(patho)) {
  d <- patho[[nm]]; p <- fisher.test(table(factor(d,c(1,0)), island))$p.value
  cat(sprintf("  %-18s Fisher p=%.4f\n", nm, p))
}
for (nm in names(marker_pa)) {
  d <- marker_pa[[nm]]; p <- fisher.test(table(factor(d,c(1,0)), island))$p.value
  cat(sprintf("  %-18s Fisher p=%.4f\n", nm, p))
}

## ------------------------------------------------------------------ ##
## 5. Section 3.3 — INDICATOR -> PATHOGEN (Table 4 / Table S7)         ##
##    univariate cluster-robust logistic + score WCB                   ##
## ------------------------------------------------------------------ ##
cat("\n================ 3.3  INDICATOR/MARKER -> PATHOGEN (Table 4/S7) ================\n")
cat(sprintf("%-18s%-20s%8s%9s%10s\n","indicator/marker","pathogen","beta/SD","p(t19)","p(wboot)"))
for (pn in names(patho)) for (inm in names(logind)) {
  r <- logit_cr_wcb(logind[[inm]], as.numeric(patho[[pn]]), site, B = 9999)
  cat(sprintf("%-18s%-20s%8.2f%9.3f%10.3f\n", inm, pn, r["beta"], r["p_t"], r["p_wcb"]))
}
cat("E. coli–enterococci Spearman rho =",
    round(cor(logind[["E. coli"]], logind[["Enterococci"]], method="spearman", use="complete.obs"),2),"\n")

## ------------------------------------------------------------------ ##
## 6. Section 3.4 — FIO CONCENTRATION vs MARKER STATUS (Figure 4)      ##
##    cluster-robust OLS of log10(indicator) on any-marker (present/absent)
##    below-detection coliphage/C.perfringens excluded (per Methods)   ##
## ------------------------------------------------------------------ ##
cat("\n================ 3.4  FIO CONCENTRATION vs ANY-MARKER (Figure 4) ================\n")
for (inm in c("E. coli","Enterococci","C. perfringens","F+ coliphage","Somatic coliphage")) {
  y <- log10(conc[[inm]]); x <- any_marker
  ok <- is.finite(y)                                   # drop non-detects (NA) & non-positive
  fit <- ols_cr(cbind(1, x[ok]), y[ok], site[ok])
  cat(sprintf("  %-18s higher-in-marker+ beta=%+.2f  p=%.3f  (n=%d)\n",
              inm, fit$beta[2], fit$p[2], sum(ok)))
}

## ------------------------------------------------------------------ ##
## 7. Section 3.5 — MARKER -> PATHOGEN DECOUPLING (Table S8)           ##
##    markers as predictors of pathogen detection: concentration & P/A ##
## ------------------------------------------------------------------ ##
cat("\n================ 3.5  MARKER -> PATHOGEN DECOUPLING (Table S8) ================\n")
cat("(a) marker log10 concentration as predictor (univariate logistic):\n")
for (mk in c("crAssphage","HF183")) for (pn in names(patho)) {
  r <- logit_cr_wcb(logind[[mk]], as.numeric(patho[[pn]]), site, B = 0)
  cat(sprintf("  %-11s -> %-20s p=%.3f\n", mk, pn, r["p_t"]))
}
cat("(b) marker presence/absence as predictor:\n")
for (mk in c("crAssphage","HF183")) for (pn in names(patho)) {
  r <- logit_cr_wcb(as.numeric(marker_pa[[mk]]), as.numeric(patho[[pn]]), site, B = 0)
  cat(sprintf("  %-11s -> %-20s p=%.3f\n", mk, pn, r["p_t"]))
}

## ------------------------------------------------------------------ ##
## 8. Section 3.6 — ENVIRONMENTAL / LAND-USE DRIVERS                   ##
## ------------------------------------------------------------------ ##
cat("\n================ 3.6  ENVIRONMENTAL / LAND-USE DRIVERS ================\n")

## 8a. Table S9 — MULTIVARIABLE concentration models (per +1 SD)
##     predictors: COND, TURB, RAIN, OSDS(cesspool), DEV, AG, GRASS (C-CAP)
Z <- list(COND = zsc(num(C_COND)),
          TURB = zsc(log10(num(C_TURB) + 0.1)),
          RAIN = zsc(log1p(num(C_RAIN24))),
          OSDS = zsc(log10(cess_density + 0.1)),
          DEV  = zsc(lc_col("Developed Land")),
          AG   = zsc(lc_col("Agricultural Land")),
          GRASS= zsc(lc_col("Grassland")))
PRED <- c("COND","TURB","RAIN","OSDS","DEV","AG","GRASS")
PRED_HF <- c("COND","OSDS","DEV")                        # HF183: reduced (91% censored)
Xmat <- function(p) cbind(1, do.call(cbind, Z[p]))
tobit_parse <- function(j, floor = 1) {                  # -> (ylog, Llog, cens)
  y<-rep(NA,80);L<-rep(NA,80);cens<-rep(FALSE,80)
  for (i in seq_len(80)) {
    s <- trimws(as.character(gv(j)[i]))
    if (is.na(gv(j)[i]) || s=="" || tolower(s)%in%c("na","nan","nd","none")) next
    if (grepl("^<",s)) { dl<-as.numeric(gsub("[^0-9.]","",s)); cens[i]<-TRUE; L[i]<-log10(dl); y[i]<-L[i] }
    else { v<-as.numeric(gsub("[^0-9.]","",s)); if (v>0){y[i]<-log10(v);L[i]<-y[i]} else {cens[i]<-TRUE;L[i]<-log10(floor);y[i]<-L[i]} }
  }
  list(y=y,L=L,cens=cens)
}
cat("-- Table S9: multivariable concentration ~ environment (beta per +1 SD) --\n")
cat(sprintf("%-18s%-6s%8s%8s%8s %s\n","response","pred","beta","se","p","ncens"))
runS9 <- function(label, kind, j, preds) {
  X <- Xmat(preds)
  if (kind == "ols") {
    y <- vapply(gv(j), function(x){v<-parse_conc(x); if(is.na(v)||v<=0) NA_real_ else log10(v)}, numeric(1))
    fit <- ols_cr(X, y, site)
  } else { m <- tobit_parse(j); fit <- tobit_cr(X, m$y, m$L, m$cens, site) }
  for (pi in seq_along(preds))
    cat(sprintf("%-18s%-6s%8.3f%8.3f%8.3f %s\n", label, preds[pi],
                fit$beta[pi+1], fit$se[pi+1], fit$p[pi+1], fit$ncens))
}
runS9("E. coli","ols",C_ECOLI,PRED); runS9("Enterococci","ols",C_ENT,PRED)
runS9("C. perfringens","tobit",C_CPERF,PRED); runS9("F+ coliphage","tobit",C_FPLUS,PRED)
runS9("Somatic coliphage","tobit",C_SOMATIC,PRED); runS9("crAssphage","tobit",C_CRASS_C,PRED)
runS9("HF183 (91% cens)","tobit",C_HF183_C,PRED_HF)

## 8b. Table S10 — pathogen detection ~ physico-chemical (univariate logistic, per +1 SD)
cat("\n-- Table S10a: pathogen detection ~ physico-chemical (logistic, per +1 SD) --\n")
physchem <- list(Conductivity=num(C_COND), Turbidity=num(C_TURB), pH=num(C_PH),
                 "Dissolved oxygen"=num(C_DO), Temperature=num(C_TEMP), "24-h rainfall"=num(C_RAIN24))
for (en in names(physchem)) {
  line <- sprintf("  %-16s", en)
  for (pn in names(patho)) {
    r <- logit_cr_wcb(physchem[[en]], as.numeric(patho[[pn]]), site, B = 0)
    line <- paste0(line, sprintf("  %s: b=%+.2f p=%.3f", substr(pn,1,4), r["beta"], r["p_t"]))
  }
  cat(line, "\n")
}

## 8c. Table S10 — pathogen detection ~ watershed land cover (STREAM-level Spearman, n=20)
cat("\n-- Table S10b: pathogen detection-rate ~ land cover (Spearman, n=20) + BH-FDR --\n")
streams <- sort(unique(stream_id))
stream_rate <- function(v) sapply(streams, function(s) mean(v[stream_id == s]))
stream_val  <- function(v) sapply(streams, function(s) mean(v[stream_id == s], na.rm = TRUE))
lc_stream <- function(cls) sapply(streams, function(s) lu_by_id[[s]][match(cls, lu_classes)])
land_pred <- c("Developed Land","Agricultural Land","Grassland","Forest Land")
prate <- lapply(patho, stream_rate)
allp <- c(); tags <- c()
for (cls in c(land_pred)) {
  xs <- lc_stream(cls); line <- sprintf("  %-18s", cls)
  for (pn in names(patho)) {
    sp <- spearman_perm(xs, prate[[pn]])
    line <- paste0(line, sprintf("  %s rho=%+.2f p=%.3f", substr(pn,1,4), sp["rho"], sp["p"]))
    allp <- c(allp, sp["p"]); tags <- c(tags, paste(cls, pn))
  }
  cat(line, "\n")
}
## cesspool density (stream level)
xs <- stream_val(cess_density); line <- sprintf("  %-18s", "Cesspool density")
for (pn in names(patho)) { sp <- spearman_perm(xs, prate[[pn]]); line <- paste0(line, sprintf("  %s rho=%+.2f p=%.3f", substr(pn,1,4), sp["rho"], sp["p"])) }
cat(line, "\n")
cat(sprintf("\n  Benjamini-Hochberg FDR across the %d land-cover tests: min q = %.3f (none < 0.05 => land-cover associations not robust)\n",
            length(allp), min(bh(allp))))

cat("\n================ DONE ================\n")
