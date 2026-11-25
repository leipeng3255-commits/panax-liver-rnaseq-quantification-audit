fit<-glmQLFit(y,design,robust=TRUE,legacy=FALSE,keep.unit.mat=TRUE)
qt<-glmQLFTest(fit,coef=2);tr<-glmTreat(fit,coef=2,lfc=log2(1.5))
qtb<-topTags(qt,n=Inf,sort.by='none')$table;trb<-topTags(tr,n=Inf,sort.by='none')$table
v<-voom(y,design,plot=FALSE);lf<-eBayes(lmFit(v,design),robust=TRUE)
lt<-topTable(lf,coef=2,number=Inf,sort.by='none',confint=TRUE)
stopifnot(identical(rownames(qtb),anno$gene_id),identical(rownames(lt),anno$gene_id))
out<-data.frame(anno,edgeR_logFC=qtb$logFC,logCPM=qtb$logCPM,QL_F=qtb$F,QL_P=qtb$PValue,QL_BH=qtb$FDR,TREAT_P=trb$PValue,TREAT_BH=trb$FDR,voom_logFC=lt$logFC,voom_CI_low=lt$CI.L,voom_CI_high=lt$CI.R,voom_t=lt$t,voom_P=lt$P.Value,voom_BH=lt$adj.P.Val,QL_var_post=fit$s2.post)
write_tsv(out,'gene_results_all.tsv')
write_tsv(data.frame(gene_id=rownames(v$E),v$E,check.names=FALSE),'voom_expression.tsv')
write_tsv(data.frame(gene_id=anno$gene_id,v$weights,check.names=FALSE),'voom_weights.tsv')
# Some limma versions do not preserve rownames on the weights matrix.
if(!file.exists(file.path(P,'voom_weights.tsv'))||nrow(read.delim(file.path(P,'voom_weights.tsv')))!=nrow(anno))stop('weights export failure')
# Interpret the 50 sets against the actual unique-mapped tested gene background.
gmt<-strsplit(readLines(GMT),'\t')
sets<-setNames(lapply(gmt,function(z)unique(z[-(1:2)])),sapply(gmt,'[',1));stopifnot(length(sets)==50,!anyDuplicated(names(sets)))
u<-which(anno$unique_symbol_in_domain);symbols<-anno$symbol[u]
indexes<-lapply(sets,function(s)which(symbols%in%s)); eligible<-lengths(indexes)>=15&lengths(indexes)<=500
coverage<-data.frame(gene_set=names(sets),source_size=lengths(sets),matched=lengths(indexes),eligible=eligible)
write_tsv(coverage,'pathway_coverage.tsv');indexes<-indexes[eligible]
write_tsv(do.call(rbind,lapply(names(indexes),function(n)data.frame(gene_set=n,gene_id=anno$gene_id[u[indexes[[n]]]],symbol=symbols[indexes[[n]]]))),'pathway_members.tsv')
get_camera<-function(vo,dd,scenario){
 all<-list(); vu<-vo[u,]
 for(tag in c('estimated','zero','0.01','0.05')){
  rho<-switch(tag,estimated=NA_real_,zero=0,`0.01`=.01,`0.05`=.05)
  cam<-camera(vu,indexes,dd,contrast=2,inter.gene.cor=rho,allow.neg.cor=FALSE,sort=FALSE)
  rr<-if('Correlation'%in%names(cam))cam$Correlation else rep(rho,nrow(cam))
  all[[tag]]<-data.frame(scenario=scenario,correlation_setting=tag,gene_set=rownames(cam),NGenes=cam$NGenes,correlation=rr,VIF=1+(cam$NGenes-1)*pmax(rr,0),Direction=cam$Direction,P=cam$PValue,BH=cam$FDR)
 }
 do.call(rbind,all)
}
cam<-get_camera(v,design,'full');write_tsv(cam,'pathway_CAMERA_all.tsv')
# Refit every prespecified deletion, with the same gene domain.
loo<-list();cams<-list()
for(j in seq_along(runs)){
 inc<-setdiff(seq_along(runs),j);m<-meta[inc,,drop=FALSE];de<-model.matrix(~group,m)
 yy<-calcNormFactors(DGEList(counts=cts[,inc,drop=FALSE],group=m$group),method='TMM')
 ff<-glmQLFit(yy,de,robust=TRUE,legacy=FALSE);tt<-topTags(glmQLFTest(ff,coef=2),n=Inf,sort.by='none')$table
 vv<-voom(yy,de,plot=FALSE)
 loo[[j]]<-data.frame(omitted=runs[j],gene_id=rownames(tt),logFC=tt$logFC,P=tt$PValue,BH=tt$FDR)
 cams[[j]]<-get_camera(vv,de,paste0('omit_',runs[j]));cat('finished omission',runs[j],'\n')
}
lo<-do.call(rbind,loo);write_tsv(lo,'gene_leave_one_all.tsv');write_tsv(do.call(rbind,cams),'pathway_CAMERA_leave_one.tsv')
lo_fc<-sapply(loo,function(x)x$logFC);stopifnot(all(vapply(loo,function(x)identical(x$gene_id,anno$gene_id),logical(1))))
rob<-data.frame(gene_id=anno$gene_id,min_logFC=apply(lo_fc,1,min),max_logFC=apply(lo_fc,1,max),same_sign_n=rowSums(sign(lo_fc)==sign(qtb$logFC)),BH_lt05_n=rowSums(sapply(loo,function(x)x$BH)<.05))
write_tsv(rob,'gene_leave_one_summary.tsv')
pdf(file.path(P,'model_diagnostics.pdf'),width=10,height=8,useDingbats=FALSE);par(mfrow=c(2,2))
plotQLDisp(fit);plotMD(qt,status=ifelse(qtb$FDR<.05,sign(qtb$logFC),0),main='edgeR QL: PNS - Model')
hist(qtb$PValue,breaks=30,main='All QL P values',xlab='Model-based P')
plot(out$edgeR_logFC,out$voom_logFC,pch=16,cex=.25,col='#24758c55',xlab='edgeR log2 FC',ylab='voom log2 FC',main='Same data, two models');abline(0,1,col='grey50');dev.off()
summary<-list(stage='CONDITIONAL_EXPLORATORY_MODEL_COMPLETE',genes=nrow(out),samples=6,edgeR_QL_BH_lt05=sum(out$QL_BH<.05),edgeR_TREAT_1_5fold_BH_lt05=sum(out$TREAT_BH<.05),voom_BH_lt05=sum(out$voom_BH<.05),both_models_BH_lt05=sum(out$QL_BH<.05&out$voom_BH<.05),effect_spearman=cor(out$edgeR_logFC,out$voom_logFC,method='spearman'),same_sign_n=sum(sign(out$edgeR_logFC)==sign(out$voom_logFC)),QL_significant_stable_sign_all_6=sum(out$QL_BH<.05&rob$same_sign_n==6),QL_significant_retained_BH_all_6=sum(out$QL_BH<.05&rob$BH_lt05_n==6),pathway_background_genes=length(u),pathways_tested=length(indexes),camera_BH_lt05=as.list(tapply(cam$BH<.05,cam$correlation_setting,sum)),source_reported_independence_not_externally_verified=TRUE,causal_or_restoration_claim=FALSE)
write_json(summary,file.path(P,'RESULTS_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
saveRDS(list(y=y,fit=fit,ql=qt,treat=tr,voom=v,limma=lf,design=design,annotation=anno),file.path(P,'model_objects.rds'))
cat(toJSON(summary,auto_unbox=TRUE,digits=12),'\n')
