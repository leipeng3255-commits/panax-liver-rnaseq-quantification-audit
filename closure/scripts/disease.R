o <- normalizePath(commandArgs(trailingOnly=TRUE)[1], mustWork=TRUE)
r <- file.path(o,'root')
f <- 'closure/inputs'
p <- 'models/flower'
library(edgeR);library(limma)
stopifnot(as.character(getRversion())=='4.6.1',as.character(packageVersion('edgeR'))=='4.10.7',as.character(packageVersion('limma'))=='3.68.5')
sessionInfo()
tab <- function(x) read.delim(x,check.names=FALSE)
mat <- function(file){a<-tab(file);stopifnot(!anyDuplicated(a$gene_id));z<-as.matrix(a[,-1]);rownames(z)<-a$gene_id;z}
tp <- mat(file.path(f,'import_abundance.tsv.gz'))
raw <- mat(file.path(f,'import_counts.tsv.gz'))
stopifnot(identical(dimnames(raw),dimnames(tp)))
s <- readRDS(file.path(p,'inputs/filtered_model_inputs.rds'))
j <- read.csv(file.path(f,'flower_samples.csv.gz'))
groups <- split(j$run,j$condition)
groups <- lapply(groups,sort)
a <- tab(file.path(p,'context_inputs/root/gene_annotation.tsv'));rownames(a)<-a$gene_id
g <- read.csv(file.path(f,'flower_descriptive.csv.gz'),check.names=FALSE)
z <- read.csv(file.path(f,'root_descriptive.csv.gz'));rownames(z)<-z$gene_id
sdir<-file.path(o,'round71b/task2')
rr<-readRDS(file.path(f,'root_disease.rds'))
len<-mat(file.path(f,'import_length.tsv.gz'))
stopifnot(identical(dimnames(len),dimnames(raw)),all(is.finite(len)),all(len>0))
stopifnot(as.character(packageVersion('tximport'))=='1.40.0')
A<-tximport:::makeCountsFromAbundance(raw,tp,len,'lengthScaledTPM')
stopifnot(isTRUE(all.equal(A[rownames(s$counts),colnames(s$counts)],s$counts,tolerance=1e-12,check.attributes=FALSE)))
gmt<-strsplit(readLines(file.path(p,'resources/mh.all.v2025.1.Mm.symbols.gmt')),'\t')
fit<-function(counts,lo,hi,tag,len=NULL){
 runs<-c(sort(lo),sort(hi));grp<-factor(rep(c(0,1),c(length(lo),length(hi))));X<-model.matrix(~grp)
 m<-counts[a$gene_id,runs,drop=FALSE];m<-m[rowSums(m>=10)>=length(lo),,drop=FALSE]
 if(is.null(len)){
   y<-calcNormFactors(DGEList(m,group=grp));v<-voom(y,X,plot=FALSE)
 }else{
   L<-len[rownames(m),runs];N<-L/exp(rowMeans(log(L)));adj<-m/N;eff<-colSums(adj)*calcNormFactors(adj)
   full<-sweep(N,2,eff,'*');y<-scaleOffset(DGEList(m,group=grp),log(full))
   v<-voom(m,X,lib.size=eff,offset.prior=log(N),plot=FALSE)
   stopifnot(max(abs(v$E-log2((m+.5)/(full+1)*1e6)))<1e-11)
 }
 q<-glmQLFit(y,X,robust=TRUE,legacy=FALSE,keep.unit.mat=TRUE)
 z<-topTags(glmQLFTest(q,coef=2),n=Inf,sort.by='none')$table
 t<-eBayes(lmFit(v,X),robust=TRUE)$t[,2]
 stats<-data.frame(gene_id=rownames(m),edgeR_logFC=z$logFC,voom_t=t)
 u<-which(a[rownames(m),'unique_symbol_in_domain']);symbols<-a[rownames(m)[u],'symbol']
 ix<-setNames(lapply(gmt,function(gs) which(symbols%in%gs[-(1:2)])),vapply(gmt,'[','',1))
 eligible<-lengths(ix)>=15 & lengths(ix)<=500
 cc<-do.call(rbind,lapply(c('estimated','zero'),function(k){
  x<-camera(v[u,],ix[eligible],X,contrast=2,inter.gene.cor=if(k=='zero') 0 else NA_real_,allow.neg.cor=FALSE,sort=FALSE)
  x<-x[match(names(ix),rownames(x)),,drop=FALSE]
  data.frame(study=if(tag=='root_A') 'root' else 'flower',scale=tag,correlation=k,gene_set=names(ix),NGenes=lengths(ix),eligible=eligible,direction=x$Direction,p=x$PValue,BH=x$FDR,significant=x$FDR<.05,overlap_other_study=NA,status=ifelse(eligible,'computed','outside_15_500'))
 }))
 write.csv(stats,file.path(sdir,paste0(tag,'_allgenes.csv')),row.names=FALSE)
 cat(tag,nrow(m),paste(runs,collapse='|'),'\n')
 list(stats=stats,cam=cc)
}
rlo<-rr$samples$run[rr$samples$group=='Control'];rhi<-rr$samples$run[rr$samples$group=='Model']
stopifnot(length(rlo)==3,length(rhi)==3,length(groups$NC)==4,length(groups$DM)==4)
root<-fit(rr$length_scaled_estimated_counts,rlo,rhi,'root_A')
models<-list(flower_A=fit(A,groups$NC,groups$DM,'flower_A'),flower_B2=fit(raw,groups$NC,groups$DM,'flower_B2',len))
orig<-read.csv(file.path(sdir,'pre71b_cross.csv'));rows<-list(orig[is.finite(orig$rho),])
for (tag in names(models)){
 rt<-root$stats;ft<-models[[tag]]$stats;ids<-intersect(rt$gene_id,ft$gene_id);ir<-match(ids,rt$gene_id);ig<-match(ids,ft$gene_id)
 expr<-rowMeans(log2(1+tp[ids,]));low<-expr<=median(expr)
 for(nm in c('edgeR_logFC','voom_t')) for(t in c('all','lower','higher')){
 k<-switch(t,all=rep(TRUE,length(ids)),lower=low,higher=!low)
 rows[[length(rows)+1]]<-data.frame(contrast='disease',statistic=nm,stratum=t,n_genes=sum(k),rho=cor(rt[[nm]][ir][k],ft[[nm]][ig][k],method='spearman'),domain='intersection of respective coding disease model domains',count_scale=paste('root_A vs',tag))
 }
}
cross<-do.call(rbind,rows)
old<-read.csv(file.path(sdir,'pre71b_camera.csv'));old<-old[old$scale=='scaledTPM_secondary',]
cc<-do.call(rbind,c(list(root$cam),lapply(models,'[[','cam'),list(old)))
summ<-list()
for(tag in c(names(models),'scaledTPM_secondary')) for(k in c('estimated','zero')){
 ai<-which(cc$scale=='root_A' & cc$correlation==k);bi<-which(cc$scale==tag & cc$correlation==k);bi<-bi[match(cc$gene_set[ai],cc$gene_set[bi])]
 both<-cc$significant[ai] & cc$significant[bi];same<-both & cc$direction[ai]==cc$direction[bi]
 cc$overlap_other_study[bi]<-both
 field<-paste0('overlap_',tag);if(!field%in%names(cc)) cc[[field]]<-NA;cc[[field]][ai]<-both
 summ[[length(summ)+1]]<-data.frame(flower_scale=tag,correlation=k,root_eligible=sum(cc$eligible[ai]),flower_eligible=sum(cc$eligible[bi]),root_hits=sum(cc$significant[ai],na.rm=TRUE),flower_hits=sum(cc$significant[bi],na.rm=TRUE),overlap=sum(both,na.rm=TRUE),same_direction=sum(same,na.rm=TRUE),overlap_names=paste(cc$gene_set[ai][which(both)],collapse='|'),same_direction_names=paste(cc$gene_set[ai][which(same)],collapse='|'))
}
stopifnot(nrow(cc)==400,!anyDuplicated(cc[,c('scale','correlation','gene_set')]),all(is.finite(cross$rho)))
cs<-do.call(rbind,summ);write.csv(cs,file.path(sdir,'camera_overlap.csv'),row.names=FALSE)
write.csv(cross,file.path(o,'cross_context_snr_calibration.csv'),row.names=FALSE)
write.csv(cc,file.path(o,'disease_positive_control_camera.csv'),row.names=FALSE)
print(cross[cross$contrast=='disease' & cross$stratum=='all',]);print(cs[,1:8]);sessionInfo()
