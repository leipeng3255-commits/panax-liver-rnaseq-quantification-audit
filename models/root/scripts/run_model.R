# Adapter for a filtered model-input tier reproduction. No original project path.
# Rscript may escape spaces in its internal --file argument. The wrapper passes
# its resolved package directory explicitly, rather than parsing that argument.
bundle_env <- Sys.getenv('D03_REPRO_BUNDLE','')
stopifnot(nzchar(bundle_env))
bundle <- normalizePath(bundle_env,mustWork=TRUE)
lib_env <- Sys.getenv('D03_R_LIBS','')
if(nzchar(lib_env)) .libPaths(c(strsplit(lib_env,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite);library(digest)})
required <- c('edgeR','limma','statmod','locfit','jsonlite','digest')
versions <- as.list(vapply(required,function(p)as.character(packageVersion(p)),character(1)))
args <- commandArgs(trailingOnly=TRUE)
if(identical(args,'--check')){
  cat(toJSON(list(status='DEPENDENCIES_AVAILABLE',R_version=R.version.string,
    package_versions=versions,library_paths=.libPaths(),platform=R.version$platform),
    pretty=TRUE,auto_unbox=TRUE),'\n')
  quit(save='no',status=0)
}
stopifnot(length(args)==2L,args[1]=='--output')
P <- normalizePath(args[2],mustWork=TRUE)
allowed_start <- c('RUN_CLAIM.json','execution.log')
stopifnot(file.exists(file.path(P,'RUN_CLAIM.json')),
          all(list.files(P,all.files=TRUE,no..=TRUE) %in% allowed_start))
claim <- fromJSON(file.path(P,'RUN_CLAIM.json'))
stopifnot(identical(claim$bundle, bundle),identical(claim$tier,'filtered model-input tier'))
options(digits=17);set.seed(20261030)
write_tsv<-function(x,name)write.table(x,file.path(P,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
input_path <- file.path(bundle,'inputs','filtered_model_inputs.rds')
input <- readRDS(input_path)
stopifnot(identical(names(input),c('counts','annotation','meta')))
cts <- input$counts;anno <- input$annotation;meta <- input$meta;runs <- meta$run
stopifnot(is.matrix(cts),is.numeric(cts),all(is.finite(cts)),all(cts>=0),
          identical(dim(cts),c(13318L,6L)),!anyDuplicated(rownames(cts)),
          identical(rownames(cts),anno$gene_id),identical(colnames(cts),runs),
          identical(as.character(meta$group),rep(c('Model','PNS'),each=3)),
          identical(levels(meta$group),c('Model','PNS')),!anyDuplicated(meta$biosample))
y<-calcNormFactors(DGEList(counts=cts,group=meta$group),method='TMM')
design<-model.matrix(~group,meta)
stopifnot(qr(design)$rank==2,nrow(design)-qr(design)$rank==4)
GMT <- file.path(bundle,'resources','mh.all.v2025.1.Mm.symbols.gmt')
write_tsv(data.frame(run=runs,design,check.names=FALSE),'design_matrix.tsv')
write_tsv(data.frame(run=runs,y$samples,check.names=FALSE),'reconstructed_y_samples.tsv')
writeLines(capture.output(sessionInfo()),file.path(P,'sessionInfo.txt'))
write_json(list(tier='filtered model-input tier',R_version=R.version.string,
  package_versions=versions,library_paths=.libPaths(),input_sha256=digest(file=input_path,algo='sha256'),
  gmt_sha256=digest(file=GMT,algo='sha256'),seed=20261030L,
  revalidated_filtering_or_tximport_or_FASTQ=FALSE),file.path(P,'RUN_ENVIRONMENT.json'),
  pretty=TRUE,auto_unbox=TRUE,digits=16)
# This file contains the original round40 core, with only its GMT path adapted.
source(file.path(bundle,'scripts','statistical_core.R'),local=TRUE,chdir=FALSE,keep.source=FALSE)
writeLines(capture.output(sessionInfo()),file.path(P,'sessionInfo_complete.txt'))
