fit<-fit_counts(s$counts,s$meta$group)
# Numerical tail identity; this is not a null-calibration check.
stopifnot(max(abs(2*pt(-abs(fit$limma$t[,2]),df=fit$limma$df.total)-fit$results$voom_P))<1e-12)
out<-merge(s$annotation,fit$results,by='gene_id',sort=FALSE);out<-out[match(rownames(s$counts),out$gene_id),];write_tsv(out,'gene_results_all.tsv')
write_tsv(data.frame(gene_id=rownames(fit$voom$E),fit$voom$E,check.names=FALSE),'voom_expression.tsv');write_tsv(data.frame(gene_id=rownames(s$counts),fit$voom$weights,check.names=FALSE),'voom_weights.tsv')
gmt<-strsplit(readLines(GMT),'\t');sets<-setNames(lapply(gmt,function(z)unique(z[-(1:2)])),sapply(gmt,'[',1));stopifnot(length(sets)==50)
u<-which(s$annotation$unique_symbol_in_domain);symbols<-s$annotation$symbol[u];indexes<-lapply(sets,function(x)which(symbols%in%x));ok<-lengths(indexes)>=15&lengths(indexes)<=500
write_tsv(data.frame(gene_set=names(indexes),source_size=lengths(sets),matched=lengths(indexes),eligible=ok),'pathway_coverage.tsv');indexes<-indexes[ok];stopifnot(length(indexes)>0)
write_tsv(do.call(rbind,lapply(names(indexes),function(n)data.frame(gene_set=n,gene_id=s$annotation$gene_id[u[indexes[[n]]]],symbol=symbols[indexes[[n]]]))),'pathway_members.tsv')
cam<-run_camera(fit$voom,fit$design,u,indexes,'full');write_tsv(cam,'pathway_CAMERA_all.tsv')
loo<-list();cams<-list()
for(j in seq_len(8)){
 inc<-setdiff(seq_len(8),j);f<-fit_counts(s$counts[,inc,drop=FALSE],s$meta$group[inc],FALSE)
 loo[[j]]<-data.frame(omitted=s$meta$run[j],f$results);cams[[j]]<-run_camera(f$voom,f$design,u,indexes,paste0('omit_',s$meta$run[j]));cat('finished omission',s$meta$run[j],'\n')
}
write_tsv(do.call(rbind,loo),'gene_leave_one_all.tsv');write_tsv(do.call(rbind,cams),'pathway_CAMERA_leave_one.tsv')
root<-read.delim(file.path(R40,'gene_results_all.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
domain<-read.delim(file.path(R40,'gene_annotation.tsv'),stringsAsFactors=FALSE)
root_filter<-read.delim(file.path(R40,'gene_filter.tsv'),stringsAsFactors=FALSE)
flower_filter<-read.delim(file.path(P,'gene_filter.tsv'),stringsAsFactors=FALSE)
context<-fromJSON(file.path(P,'CONTEXT_CONTRACT.json'),simplifyVector=FALSE)
shared<-assemble_shared(root,out,domain,root_filter,flower_filter,context)
# Root-only pointwise bound is fixed before flower values: every BY value must be 1.
stopifnot(all(shared$conditional_BY==1))
write_tsv(shared,'cross_context_shared_direction.tsv')
write_tsv(as.data.frame(table(shared$coverage_status),stringsAsFactors=FALSE),'cross_context_coverage_summary.tsv')

pdf(file.path(P,'model_diagnostics.pdf'),width=10,height=8,useDingbats=FALSE);par(mfrow=c(2,2));plotQLDisp(fit$fit);plotMD(fit$ql);hist(out$QL_P,breaks=30,main='All QL P values',xlab='Conditional P');plot(out$edgeR_logFC,out$voom_logFC,pch=16,cex=.2,col='#24758c55',xlab='edgeR log2FC',ylab='voom log2FC');abline(0,1);dev.off()
saveRDS(fit,file.path(P,'model_objects.rds'))
write_json(list(stage='CONDITIONAL_EXPLORATORY_MODELS_COMPLETE_PENDING_INDEPENDENT_QA',genes=nrow(out),samples=8,edgeR_QL_BH_lt05=sum(out$QL_BH<.05),voom_BH_lt05=sum(out$voom_BH<.05),TREAT_BH_lt05=sum(out$TREAT_BH<.05),camera_BH_lt05=as.list(tapply(cam$BH<.05,cam$correlation_setting,sum)),cross_context_fixed_universe=nrow(shared),cross_context_both_modeled=sum(shared$coverage_status=='BOTH_MODELED'),shared_diagnostic_BY_lt05=sum(shared$conditional_BY<.05),shared_diagnostic_BH_lt05=sum(shared$conditional_BH<.05),cross_context_discovery_status='DISABLED_ROOT_ONLY_BOUND_PRECLUDES_DISCOVERY',calibrated_replicability_FDR_established=FALSE,biological_independence_externally_verified=FALSE,causal_or_material_equivalence=FALSE),file.path(P,'RESULTS_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(P,'sessionInfo_fit.txt'))
