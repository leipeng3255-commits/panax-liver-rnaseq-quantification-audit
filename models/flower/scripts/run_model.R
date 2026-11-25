# Filtered model-input adapter; no original study path is used at execution.
be<-Sys.getenv('D03_REPRO_BUNDLE','');stopifnot(nzchar(be));bundle<-normalizePath(be,mustWork=TRUE)
le<-Sys.getenv('D03_R_LIBS','');if(nzchar(le)) .libPaths(c(strsplit(le,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite);library(digest)})
required<-c('edgeR','limma','statmod','locfit','jsonlite','digest')
versions<-as.list(vapply(required,function(p)as.character(packageVersion(p)),character(1)))
args<-commandArgs(trailingOnly=TRUE)
if(identical(args,'--check')){
 cat(toJSON(list(status='DEPENDENCIES_AVAILABLE',R_version=R.version.string,package_versions=versions,library_paths=.libPaths(),platform=R.version$platform),pretty=TRUE,auto_unbox=TRUE),'\n');quit(save='no',status=0)
}
stopifnot(length(args)==2L,args[1]=='--output');P<-normalizePath(args[2],mustWork=TRUE)
stopifnot(file.exists(file.path(P,'RUN_CLAIM.json')),all(list.files(P,all.files=TRUE,no..=TRUE)%in%c('RUN_CLAIM.json','execution.log')))
claim<-fromJSON(file.path(P,'RUN_CLAIM.json'));stopifnot(identical(claim$bundle,bundle),identical(claim$tier,'filtered model-input tier'))
options(digits=17);set.seed(20261030)
write_tsv<-function(x,n)write.table(x,file.path(P,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
ip<-file.path(bundle,'inputs','filtered_model_inputs.rds');s<-readRDS(ip)
stopifnot(identical(names(s),c('counts','annotation','meta')),identical(dim(s$counts),c(12141L,8L)),
 is.matrix(s$counts),is.numeric(s$counts),all(is.finite(s$counts)),all(s$counts>=0),
 identical(rownames(s$counts),s$annotation$gene_id),identical(colnames(s$counts),s$meta$run),
 identical(as.character(s$meta$group),rep(c('DM','PNF'),each=4)),identical(levels(s$meta$group),c('DM','PNF')),
 !anyDuplicated(rownames(s$counts)),!anyDuplicated(s$meta$run),!anyDuplicated(s$meta$biosample))
R40<-file.path(bundle,'context_inputs','root')
for(n in c('gene_filter.tsv','CONTEXT_CONTRACT.json'))stopifnot(file.copy(file.path(bundle,'context_inputs','flower',n),file.path(P,n),overwrite=FALSE))
GMT<-file.path(bundle,'resources','mh.all.v2025.1.Mm.symbols.gmt')
write_tsv(data.frame(gene_id=rownames(s$counts),s$counts,check.names=FALSE),'model_input_counts.tsv')
writeLines(capture.output(sessionInfo()),file.path(P,'sessionInfo.txt'))
write_json(list(tier='filtered model-input tier',R_version=R.version.string,package_versions=versions,
 library_paths=.libPaths(),input_sha256=digest(file=ip,algo='sha256'),gmt_sha256=digest(file=GMT,algo='sha256'),
 original_filter_tximport_or_reads_revalidated=FALSE,root_context_inputs_are_fixed_diagnostics_not_new_root_results=TRUE,
 measurement_trace_scope='Low-count short-transcript estimates remain a biological interpretation limit; computational reproduction does not resolve fragment assignment.'),file.path(P,'RUN_ENVIRONMENT.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
source(file.path(bundle,'scripts','model_functions.R'),local=TRUE,chdir=FALSE)
source(file.path(bundle,'scripts','statistical_core.R'),local=TRUE,chdir=FALSE)
# Adapter-only normalization/design exports, after the unchanged original core.
write_tsv(data.frame(run=s$meta$run,fit$y$samples,check.names=FALSE),'reconstructed_y_samples.tsv')
write_tsv(data.frame(run=s$meta$run,fit$design,check.names=FALSE),'reconstructed_design.tsv')
writeLines(capture.output(sessionInfo()),file.path(P,'sessionInfo_complete.txt'))
