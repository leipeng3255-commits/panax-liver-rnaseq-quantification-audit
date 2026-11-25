fit_counts <- function(cts,group,with_treat=TRUE) {
 group<-factor(group,levels=c('DM','PNF'));stopifnot(!anyNA(group),all(table(group)>=3))
 design<-model.matrix(~group);stopifnot(qr(design)$rank==2)
 y<-edgeR::calcNormFactors(edgeR::DGEList(counts=cts,group=group),method='TMM')
 fit<-edgeR::glmQLFit(y,design,robust=TRUE,legacy=FALSE,keep.unit.mat=TRUE)
 qt<-edgeR::glmQLFTest(fit,coef=2);q<-edgeR::topTags(qt,n=Inf,sort.by='none')$table
 v<-limma::voom(y,design,plot=FALSE);lm<-limma::eBayes(limma::lmFit(v,design),robust=TRUE)
 l<-limma::topTable(lm,coef=2,number=Inf,sort.by='none',confint=TRUE)
 stopifnot(identical(rownames(q),rownames(cts)),identical(rownames(l),rownames(cts)))
 out<-data.frame(gene_id=rownames(cts),edgeR_logFC=q$logFC,logCPM=q$logCPM,QL_F=q$F,QL_P=q$PValue,QL_BH=q$FDR,
  voom_logFC=l$logFC,voom_CI_low=l$CI.L,voom_CI_high=l$CI.R,voom_t=l$t,voom_P=l$P.Value,voom_BH=l$adj.P.Val)
 if(with_treat){tr<-edgeR::topTags(edgeR::glmTreat(fit,coef=2,lfc=log2(1.5)),n=Inf,sort.by='none')$table;out$TREAT_P<-tr$PValue;out$TREAT_BH<-tr$FDR}
 list(results=out,y=y,fit=fit,ql=qt,voom=v,limma=lm,design=design)
}
run_camera <- function(voom,design,u,indexes,scenario) {
 do.call(rbind,lapply(c('estimated','zero','0.01','0.05'),function(tag){
  rho<-switch(tag,estimated=NA_real_,zero=0,`0.01`=.01,`0.05`=.05)
  a<-limma::camera(voom[u,],indexes,design,contrast=2,inter.gene.cor=rho,allow.neg.cor=FALSE,sort=FALSE)
  cr<-if('Correlation'%in%names(a))a$Correlation else rep(rho,nrow(a))
  data.frame(scenario=scenario,correlation_setting=tag,gene_set=rownames(a),NGenes=a$NGenes,correlation=cr,
   VIF=1+(a$NGenes-1)*pmax(cr,0),Direction=a$Direction,P=a$PValue,BH=a$FDR)
 }))
}

# Symmetric two-sided P and its signed t must describe the same coefficient.
shared_direction <- function(p1,t1,p2,t2) {
 stopifnot(length(p1)>0,length(p1)==length(p2),length(t1)==length(p1),length(t2)==length(p2),
           all(is.finite(c(p1,p2,t1,t2))),all(p1>=0&p1<=1),all(p2>=0&p2<=1),
           all(p1[t1==0]==1),all(p2[t2==0]==1))
 # Direct lower/upper tails avoid cancellation for large negative t.
 up1<-ifelse(t1>=0,p1/2,1-p1/2);up2<-ifelse(t2>=0,p2/2,1-p2/2)
 dn1<-ifelse(t1<=0,p1/2,1-p1/2);dn2<-ifelse(t2<=0,p2/2,1-p2/2)
 data.frame(p_shared=pmin(1,2*pmin(pmax(up1,up2),pmax(dn1,dn2))),
  sign_relation=ifelse(sign(t1)==sign(t2)&t1!=0,ifelse(t1>0,'both_up','both_down'),'opposite_or_zero'))
}
assemble_shared <- function(root,flower,domain,root_filter,flower_filter,context) {
 stopifnot(context$root$dataset=='PRJNA1462623',context$flower$dataset=='PRJNA1545393',
  context$root$taxon_id==10090,context$flower$taxon_id==10090,
  context$root$annotation=='GENCODE_M37',context$flower$annotation=='GENCODE_M37',
  context$root$contrast=='PNS_minus_Model',context$flower$contrast=='PNF_minus_DM')
 valid_ids<-function(x){stopifnot(is.character(x),length(x)>0,!anyNA(x),!anyDuplicated(x),all(grepl('^ENSMUSG[0-9]{11}\\.[0-9]+$',x)))}
 valid_ids(domain$gene_id);stopifnot(all(domain$gene_type=='protein_coding'))
 for(flt in list(root_filter,flower_filter)){
  valid_ids(flt$gene_id);stopifnot(setequal(flt$gene_id,domain$gene_id),is.logical(flt$included),!anyNA(flt$included))
 }
 needed<-c('gene_id','voom_logFC','voom_CI_low','voom_CI_high','voom_t','voom_P')
 validate_result<-function(r,flt){
  stopifnot(all(needed%in%names(r)));valid_ids(r$gene_id)
  stopifnot(setequal(r$gene_id,flt$gene_id[flt$included]))
  for(n in needed[-1])stopifnot(is.numeric(r[[n]]),all(is.finite(r[[n]])))
  stopifnot(all(r$voom_P>=0&r$voom_P<=1),all(sign(r$voom_logFC)==sign(r$voom_t)),
   all(r$voom_CI_low<=r$voom_logFC),all(r$voom_CI_high>=r$voom_logFC),all(r$voom_P[r$voom_t==0]==1))
 }
 validate_result(root,root_filter);validate_result(flower,flower_filter)
 ids<-domain$gene_id;ir<-match(ids,root$gene_id);ifl<-match(ids,flower$gene_id)
 hasr<-!is.na(ir);hasf<-!is.na(ifl);both<-hasr&hasf
 status<-ifelse(both,'BOTH_MODELED',ifelse(hasr,'FLOWER_FILTERED',ifelse(hasf,'ROOT_FILTERED','BOTH_FILTERED')))
 out<-data.frame(gene_id=ids,coverage_status=status,p_shared=rep(1,length(ids)),sign_relation='NOT_ESTIMATED_IN_BOTH',
  p_assignment=ifelse(both,'TWO_MODEL_DIRECTIONAL_CONJUNCTION','ONE_FOR_UNTESTED_NOT_BIOLOGICAL_NULL'))
 for(prefix in c('root','flower')){
  r<-if(prefix=='root')root else flower;i<-if(prefix=='root')ir else ifl
  for(n in needed[-1])out[[paste0(prefix,'_',sub('voom_','',n))]]<-r[[n]][i]
 }
 if(any(both)){
  z<-shared_direction(root$voom_P[ir[both]],root$voom_t[ir[both]],flower$voom_P[ifl[both]],flower$voom_t[ifl[both]])
  out$p_shared[both]<-z$p_shared;out$sign_relation[both]<-z$sign_relation
 }
 out$conditional_BY<-p.adjust(out$p_shared,'BY');out$conditional_BH<-p.adjust(out$p_shared,'BH')
 out$interpretation<-'CONDITIONAL_MODEL_DIAGNOSTIC_NOT_CALIBRATED_BIOLOGICAL_REPLICATION'
 out
}
