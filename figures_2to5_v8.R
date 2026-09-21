## =====================================================================
## Kirs et al. — Hawaiian streams: code for Figures 2, 3, 4, and 5
##   Fig 2  Indicator bacteria & human markers by island (boxplots, log10)
##   Fig 3  Pathogen & human-marker detection rates by island (grouped bars)
##   Fig 4  FIB/coliphage concentrations by human-marker status (boxplots)
##   Fig 5  Spearman heatmap: microbial indicators x physico-chemical variables
##
## INPUT : data-final-5d81e002.xlsx  (4 header rows, then 80 data rows)
## NEEDS : install.packages(c("readxl","ggplot2","patchwork"))
## Each section writes a 600-dpi PNG + TIFF.
## =====================================================================

library(readxl); library(ggplot2); library(patchwork)
set.seed(42)
DATA <- "data-final-5d81e002.xlsx"

island_levels <- c("Oahu","Maui","Kauai","Hawaii")
island_labels <- c("O‘ahu","Maui","Kaua‘i","Hawai‘i")
wong <- c(Oahu="#0072B2", Maui="#E69F00", Kauai="#009E73", Hawaii="#CC79A7")

## 1-based columns
COL <- c(island=2, Ecoli=20, ENT=21, CPERF=22, Fplus=23, Somatic=24,
         Salmonella=25, Campylobacter=26, Cjejuni=27,
         crass_pa=35, hf183_pa=36, crass_c=39, hf183_c=40,
         temp=8, cond=10, DO=11, pH=13, turb=14, rain24=18)

parse_conc <- function(x) {                      # ">x"->x ; "<x"->x/2 ; blank/ND->NA
  if (is.na(x)) return(NA_real_)
  s <- trimws(as.character(x))
  if (s=="" || tolower(s) %in% c("na","nan","nd","none")) return(NA_real_)
  if (grepl("^>",s)) return(as.numeric(gsub("[^0-9.]","",s)))
  if (grepl("^<",s)) return(as.numeric(gsub("[^0-9.]","",s))/2)
  suppressWarnings(as.numeric(gsub("[^0-9.]","",s)))
}
is_pos <- function(x) if (is.na(x)) 0L else as.integer(tolower(trimws(as.character(x)))=="pos")
det_any <- function(x){ if(is.na(x)) return(0L); s<-tolower(trimws(as.character(x)))
  if(s%in%c("na","nan","nd","none","","0")) 0L else as.integer(grepl("[1-9]",s)) }
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
floor_pos <- function(v) ifelse(!is.na(v) & v<=0, 0.5, v)     # keep zeros visible on log axis

raw <- as.data.frame(read_excel(DATA, col_names=FALSE, .name_repair="minimal"))[-(1:4),]
col <- function(k) raw[[COL[[k]]]]
island <- factor(trimws(as.character(col("island"))), levels=island_levels)
site   <- as.character(raw[[1]])

## =====================================================================
## FIGURE 2 — concentrations by island (indicators top, markers bottom)
## =====================================================================
long <- function(cols, labs, parse=TRUE)
  do.call(rbind, Map(function(cc,lb) data.frame(
    analyte=lb, island=island,
    value=if(parse) floor_pos(vapply(col(cc),parse_conc,numeric(1))) else num(col(cc))),
    cols, labs))
fib <- long(c("Ecoli","ENT","CPERF"), c("E. coli","Enterococci","C. perfringens"))
fib$analyte <- factor(fib$analyte, levels=c("E. coli","Enterococci","C. perfringens"))
mk  <- long(c("crass_c","hf183_c"), c("crAssphage","HF183"))
mk$analyte  <- factor(mk$analyte, levels=c("crAssphage","HF183"))
stv <- data.frame(analyte=factor(c("E. coli","Enterococci","C. perfringens"),
                                 levels=levels(fib$analyte)), stv=c(410,130,50))
panel <- function(df, ylab, stv_df=NULL) {
  p <- ggplot(df, aes(island, value)) +
    geom_jitter(width=0.15, height=0, size=0.7, colour="grey45", alpha=0.5)
  if (!is.null(stv_df)) p <- p + geom_hline(data=stv_df, aes(yintercept=stv),
                                            colour="#B2182B", linetype="dotted", linewidth=0.6)
  p + geom_boxplot(aes(fill=island), outlier.shape=NA, alpha=0.72, colour="black", linewidth=0.4) +
    facet_wrap(~analyte, nrow=1) + scale_y_log10() +
    scale_x_discrete(limits=island_levels, labels=island_labels) +
    scale_fill_manual(values=wong, guide="none") + labs(x=NULL, y=ylab) +
    theme_classic(base_size=11) +
    theme(strip.background=element_blank(), strip.text=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=8))
}
fig2 <- panel(fib,"Concentration (MPN or CFU / 100 mL)", stv) /
        panel(mk, "Concentration (gc / 100 mL)")
ggsave("figure2_island_concentrations.png", fig2, width=9, height=7, dpi=600)
ggsave("figure2_island_concentrations.tif", fig2, width=9, height=7, dpi=600)

## =====================================================================
## FIGURE 3 — pathogen & human-marker detection rate by island (grouped bars)
## =====================================================================
d <- data.frame(island=island)
for (k in c("Salmonella","Campylobacter","Cjejuni")) d[[k]] <- vapply(col(k), det_any, integer(1))
d$crass <- vapply(col("crass_pa"), is_pos, integer(1))
d$hf183 <- vapply(col("hf183_pa"), is_pos, integer(1))
orgs <- list("Salmonella spp."="Salmonella","Campylobacter spp."="Campylobacter",
             "C. jejuni"="Cjejuni","crAssphage"="crass","HF183"="hf183")
det <- do.call(rbind, lapply(names(orgs), function(lb) data.frame(
  organism=lb, island=island_levels,
  pct=sapply(island_levels, function(i) 100*mean(d[[orgs[[lb]]]][d$island==i])))))
det$organism <- factor(det$organism, levels=names(orgs))
det$island   <- factor(det$island, levels=island_levels)
## italicise the taxa (parsed labels); markers upright
xlab_expr <- c("Salmonella spp."="italic('Salmonella')~'spp.'",
               "Campylobacter spp."="italic('Campylobacter')~'spp.'",
               "C. jejuni"="italic('C. jejuni')","crAssphage"="'crAssphage'","HF183"="'HF183'")
fig3 <- ggplot(det, aes(organism, pct, fill=island)) +
  geom_col(position=position_dodge(width=0.8), width=0.75, colour="black", linewidth=0.2) +
  scale_fill_manual(values=wong, labels=island_labels, name="Island") +
  scale_x_discrete(labels=function(v) parse(text=xlab_expr[v])) +
  scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(x=NULL, y="Samples positive (%)") + theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave("figure3_detection_by_island.png", fig3, width=8, height=4.5, dpi=600)
ggsave("figure3_detection_by_island.tif", fig3, width=8, height=4.5, dpi=600)

## =====================================================================
## FIGURE 4 — FIB/coliphage concentration by human-marker status
##   any-marker = crAssphage OR HF183 positive. Non-detects dropped for the
##   boxplots. p-values are from cluster-robust OLS of log10(conc) on marker
##   status (see stats_all_v8.R, Section 3.4).
## =====================================================================
anymk <- as.integer(vapply(col("crass_pa"),is_pos,integer(1)) |
                    vapply(col("hf183_pa"),is_pos,integer(1)))
ols_cr_p <- function(y, x, cl) {                 # cluster-robust p for slope of y~x
  ok <- is.finite(y); y<-y[ok]; x<-x[ok]; cl<-cl[ok]
  X<-cbind(1,x); n<-nrow(X); k<-2; XtXi<-solve(crossprod(X))
  b<-XtXi%*%crossprod(X,y); u<-as.vector(y-X%*%b); G<-length(unique(cl)); meat<-matrix(0,k,k)
  for(c in unique(cl)){s<-crossprod(X[cl==c,,drop=FALSE],u[cl==c]);meat<-meat+tcrossprod(s)}
  V<-(G/(G-1))*((n-1)/(n-k))*(XtXi%*%meat%*%XtXi); se<-sqrt(V[2,2])
  2*pt(-abs(b[2]/se), df=G-1)
}
inds <- list("E. coli"="Ecoli","Enterococci"="ENT","C. perfringens"="CPERF",
             "F+ coliphage"="Fplus","Somatic coliphage"="Somatic")
d4 <- do.call(rbind, lapply(names(inds), function(lb){
  v <- vapply(col(inds[[lb]]), parse_conc, numeric(1))
  data.frame(indicator=lb, value=v, marker=factor(ifelse(anymk==1,"Marker +","Marker −"),
             levels=c("Marker −","Marker +")), site=site)
}))
d4$indicator <- factor(d4$indicator, levels=names(inds))
pann <- do.call(rbind, lapply(names(inds), function(lb){
  v <- vapply(col(inds[[lb]]), parse_conc, numeric(1))
  data.frame(indicator=lb, p=ols_cr_p(log10(v), anymk, as.integer(factor(site))))
}))
pann$indicator <- factor(pann$indicator, levels=names(inds))
pann$lab <- paste0("p = ", formatC(pann$p, format="f", digits=3))
fig4 <- ggplot(subset(d4, is.finite(value) & value>0), aes(marker, value)) +
  geom_boxplot(outlier.shape=NA, width=0.6, fill="grey85", colour="grey30") +
  geom_jitter(width=0.15, height=0, size=1, alpha=0.5, colour="#2166AC") +
  facet_wrap(~indicator, nrow=1, scales="free_y") + scale_y_log10() +
  geom_text(data=pann, aes(x=1.5, y=Inf, label=lab), vjust=1.4, size=3, inherit.aes=FALSE) +
  labs(x=NULL, y="Concentration (MPN, CFU or PFU / 100 mL)") +
  theme_classic(base_size=11) +
  theme(strip.background=element_blank(), strip.text=element_text(size=9),
        axis.text.x=element_text(size=8))
ggsave("figure4_fib_by_marker.png", fig4, width=11, height=4, dpi=600)
ggsave("figure4_fib_by_marker.tif", fig4, width=11, height=4, dpi=600)

## =====================================================================
## FIGURE 5 — Spearman heatmap: microbial indicators x physico-chemical
##   permutation Spearman (10,000 relabelings); asterisks * <.05 ** <.01 *** <.001
## =====================================================================
micro <- list("E. coli"="Ecoli","Enterococci"="ENT","C. perfringens"="CPERF",
              "F+ coliphage"="Fplus","Somatic coliphage"="Somatic","crAssphage"="crass_c")
env   <- list("Temperature"="temp","Conductivity"="cond","Dissolved oxygen"="DO",
              "pH"="pH","Turbidity"="turb","Rainfall (24 h)"="rain24")
Mvals <- lapply(micro, function(cc) vapply(col(cc), parse_conc, numeric(1)))
Evals <- lapply(env,   function(cc) num(col(cc)))
spearman_perm <- function(x, y, nperm=10000) {
  ok <- !is.na(x)&!is.na(y); rx<-rank(x[ok]); ry<-rank(y[ok]); rho<-cor(rx,ry)
  perm <- replicate(nperm, cor(rx, sample(ry)))
  c(rho=rho, p=(sum(abs(perm) >= abs(rho)-1e-12)+1)/(nperm+1))
}
res <- expand.grid(micro=names(micro), env=names(env), stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE)
res$rho <- NA; res$p <- NA
for (i in seq_len(nrow(res))) { r<-spearman_perm(Mvals[[res$micro[i]]], Evals[[res$env[i]]]); res$rho[i]<-r["rho"]; res$p[i]<-r["p"] }
star <- ifelse(res$p<0.001,"***",ifelse(res$p<0.01,"**",ifelse(res$p<0.05,"*","")))
res$label <- ifelse(star=="", sprintf("%.2f",res$rho), sprintf("%.2f\n%s",res$rho,star))
res$micro <- factor(res$micro, levels=names(micro)); res$env <- factor(res$env, levels=names(env))
res$txt <- ifelse(abs(res$rho)>0.5, "white", "black")
to_expr <- function(v) vapply(as.character(v), function(s)
  if (s %in% c("E. coli","C. perfringens")) sprintf("italic('%s')",s) else sprintf("'%s'",s), character(1))
fig5 <- ggplot(res, aes(env, micro, fill=rho)) +
  geom_tile(colour="white", linewidth=1) +
  geom_text(aes(label=label, colour=txt), size=3, lineheight=0.85) +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F7", high="#8B0A2A", midpoint=0,
                       limits=c(-0.7,0.7), oob=scales::squish, name="Spearman's ρ") +
  scale_colour_identity() +
  scale_y_discrete(limits=rev(names(micro)), labels=function(v) parse(text=to_expr(v))) +
  scale_x_discrete(labels=function(v) parse(text=to_expr(v))) +
  labs(x=NULL, y=NULL) + coord_fixed() + theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1), panel.grid=element_blank())
ggsave("figure5_micro_env_heatmap.png", fig5, width=7.2, height=5.6, dpi=600)
ggsave("figure5_micro_env_heatmap.tif", fig5, width=7.2, height=5.6, dpi=600)

cat("Saved figures 2, 3, 4, 5 (png + tif).\n")
