options(digits=17)
o <- normalizePath(commandArgs(trailingOnly=TRUE)[1], mustWork=TRUE)
r44 <- "closure/inputs"
r63 <- "models/flower"
library(edgeR); library(limma)
stopifnot(as.character(getRversion())=="4.6.1",as.character(packageVersion("edgeR"))=="4.10.7",
          as.character(packageVersion("limma"))=="3.68.5")
source(file.path(r63,"scripts/model_functions.R"))
gene_import <- readRDS(file.path(r44,"flower_gene_import.rds"))
counts <- gene_import$counts; tpm <- gene_import$abundance
A <- readRDS(file.path(r63,"inputs/filtered_model_inputs.rds"))
meta <- A$meta; runs <- meta$run; group <- meta$group
ann <- read.delim(file.path(r63,"context_inputs/flower/gene_filter.tsv"))
stopifnot(identical(dimnames(counts),dimnames(tpm)),ncol(counts)==12,
          all(ann$gene_id %in% rownames(counts)))
len <- gene_import$length
scaled <- tximport:::makeCountsFromAbundance(counts,tpm,NULL,countsFromAbundance="scaledTPM")
keep <- function(m) ann$gene_id[rowSums(m[ann$gene_id,runs,drop=FALSE]>=10)>=4]
domains <- list(A=rownames(A$counts), B1=rownames(A$counts), B2=keep(counts), C=keep(scaled))
base <- read.delim(file.path(o,"flower/gene_results_all.tsv"))
res <- list(A=base)
res$C <- fit_counts(scaled[domains$C,runs],group)$results
# B1/B2 are fitted only when the original exported effective lengths are readable.
if(!is.null(len)){
  stopifnot(identical(dimnames(counts), dimnames(len)))
  for(sc in c("B1","B2")){
    ids <- domains[[sc]]
    raw <- counts[ids,runs]; norm <- len[ids,runs]
    stopifnot(all(is.finite(norm)),all(norm>0))
    norm <- norm/exp(rowMeans(log(norm)))
    normCounts <- raw/norm
    eff <- calcNormFactors(normCounts)*colSums(normCounts)
    y <- scaleOffset(DGEList(counts=raw,group=group),log(sweep(norm,2,eff,"*")))
    fit <- glmQLFit(y,model.matrix(~group),robust=TRUE,legacy=FALSE,keep.unit.mat=TRUE)
    qt <- topTags(glmQLFTest(fit,coef=2),n=Inf,sort.by="none")$table
    tr <- topTags(glmTreat(fit,coef=2,lfc=log2(1.5)),n=Inf,sort.by="none")$table
    res[[sc]] <- data.frame(gene_id=ids,edgeR_logFC=qt$logFC,QL_BH=qt$FDR,TREAT_BH=tr$FDR,voom_BH=NA_real_)
  }
}
stopifnot(!is.null(len),length(domains$B1)==12141,length(domains$B2)==12167)
for(sc in c("B1","B2")) write.csv(res[[sc]],file.path(o,"round71b/task1",paste0(sc,"_all_genes.csv")),row.names=FALSE)
sym <- setNames(ann$symbol,ann$gene_id)
sel <- base$gene_id[base$QL_BH < .05]
stopifnot(length(sel)==13,all(c("Nat8l","Cdkn1a") %in% sym[sel]))
write.csv(res$C,file.path(o,"round71b/task1/C_all_genes.csv"),row.names=FALSE)
rows <- list(); sums <- list(); fields <- c("QL_BH","TREAT_BH","voom_BH")
for(sc in names(domains)){
  z <- res[[sc]]
  s <- data.frame(scheme=sc,retained=length(domains[[sc]]),n_domain_intersect_A=length(intersect(domains$A,domains[[sc]])))
  for(f in fields){
    available <- !is.null(z) && any(is.finite(z[[f]]))
    hits <- if(available) z$gene_id[which(z[[f]] < .05)] else character()
    overlap <- intersect(hits,base$gene_id[base[[f]]<.05])
    s[[paste0(f,"_discoveries")]] <- if(available) length(hits) else NA_integer_
    s[[paste0(f,"_A_overlap_n")]] <- if(available) length(overlap) else NA_integer_
    s[[paste0(f,"_A_overlap_symbols")]] <- if(available) paste(sym[overlap],collapse="|") else NA_character_
  }
  sums[[sc]] <- s
  tab <- data.frame(scheme=sc,gene_id=sel,symbol=unname(sym[sel]),passed_filter=sel %in% domains[[sc]])
  for(gr in levels(group)){
    for(m in c("unscaled_fragments","lengthScaledTPM")){
      x <- if(m=="unscaled_fragments") counts[sel,runs[group==gr]] else A$counts[sel,runs[group==gr]]
      tab[[paste(gr,m,"min",sep="_")]] <- apply(x,1,min)
      tab[[paste(gr,m,"max",sep="_")]] <- apply(x,1,max)
    }
  }
  for(f in c("edgeR_logFC",fields)) tab[[f]] <- if(is.null(z)) NA_real_ else z[[f]][match(sel,z$gene_id)]
  tab$status <- ifelse(!tab$passed_filter,"filtered_out",if(is.null(z)) "missing_original_length_matrix" else "computed")
  rows[[sc]] <- tab
}
write.csv(do.call(rbind,rows),file.path(o,"tximport_scale_sensitivity_13genes.csv"),row.names=FALSE)
write.csv(do.call(rbind,sums),file.path(o,"tximport_scale_global_summary.csv"),row.names=FALSE)
g <- data.frame(gene_id=ann$gene_id,symbol=ann$symbol)
for(sc in names(domains)) g[[sc]] <- g$gene_id %in% domains[[sc]]
write.csv(g,file.path(o,"round71b/task1/filter_domains_all.csv"),row.names=FALSE)
print(do.call(rbind,sums))
