#!/usr/bin/Rscript --vanilla



combine_cov<-function(sample.infile=NULL, cov.file=NULL, cov.comb, last.autosome=22) {
  names=NULL;  

  if(!is.null(sample.infile)) {
    if(!is.null(cov.file)) {
      print("Error:  Both sample batch file and and individual sample have been input.");
      quit(10);
    }
    sample.info=read.table(sample.infile, stringsAsFactors=FALSE);
    nsamps=nrow(sample.info);
    for(i in 1:nsamps) {
      dd=read.table(sample.info[i,3]);
      # note!, expects coverage files to be 4 columns.  4th column (target.name) not used--could be removed
      colnames(dd)=c("chrstr", "start", "stop", "target.name", "cov", "zz", "target.len", "zz1");
      dd=dd[!duplicated(dd[, c("chrstr", "start", "stop")]),];
      if(i==1) {
        covs=dd[, c("chrstr", "start", "stop", "target.len", "cov")];
      } else {
        covs=merge(covs, dd[, c("chrstr", "start", "stop", "target.len", "cov")], by=c("chrstr", "start", "stop", "target.len"), sort=FALSE);
      }
    }
    names=make.names(sample.info[,1]);
    colnames(covs)[5:ncol(covs)]=names;
  } else if (!is.null(cov.file)) {
    covs=read.table(cov.file);
    colnames(covs)=c("chrstr", "start", "stop", "target.name", "cov", "zz", "target.len", "zz1"); 
    covs=covs[, c("chrstr", "start", "stop", "target.len", "cov")];
    names=make.names(basename(cov.file));
    colnames(covs)[5]=names;
    names=c(names);
  } else {
    print("Error: missing input file");
    quit(10);
  }
  dd=NULL;
  num=covs[, "chrstr"];
  if(grepl("chr", covs[1, "chrstr"])) {	
    num=substring(covs[, "chrstr"], 4);
  }
  num=as.character(num);
  num[num=="X"]<-last.autosome+1;
  num[num=="Y"]<-last.autosome+2;
  chr=as.numeric(num);
  covs=cbind( chr, covs);
  covs=covs[order(covs$chr, covs$start),];                          #re-sort after merge
  write.table(covs, file=cov.comb, quote=F, row.names=F, col.names=T, sep=",");
  return(names);
}


mask_low_coverage<-function(ctlmedfile, combfile, min_cov, last.autosome=22) {
  ctlmed=read.csv(ctlmedfile);
  comb=read.csv(combfile);
  merged=merge(ctlmed, comb, by=infocols, sort=FALSE);
  merged=merged[merged$med*100/merged$target.len>min_cov & merged$chr<=last.autosome,];					#assumes 100bp reads
  merged=merged[order(merged$chr, merged$start),];
  write.table(merged[, colnames(comb)], file=paste(combfile, ".low_cov_flt.csv", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",");
}


get.gc.content<-function(fasta=NULL, bedfile=NULL, gcfile) {
  if(is.null(fasta) || is.null(bedfile)) {
    print("Error:  missing reference.fa or target bedfile");
    quit(10);
  }
  cmd=paste(BEDTools,"nucBed -fi ", fasta, " -bed ", bedfile, " > ", gcfile, sep="");
  print(cmd);
  system(cmd); 
}

select.normals.get.med<-function(combfile, ctlcombfile, afdata=NULL, min.corr, min.normals, min.cov, last.autosome, medfile) {
	
  
  casecov=read.csv(combfile, stringsAsFactors=FALSE);
  casecov=casecov[casecov$chr<=last.autosome,];
  if(ncol(casecov)>6) {
    print("select normals only implemented for single-sample mode.");
	quit(10);
  } 
  case=colnames(casecov)[6];
  ctlcov=read.csv(ctlcombfile, stringsAsFactors=FALSE);
  ctlnames=colnames(ctlcov)[6:ncol(ctlcov)];
	
  if(!is.null(afdata)) {									#calculate correlation for targets with mean vaf in (0.40, 0.6)
    load(afdata);
	afd=af.list[[1]][af.list[[1]]$GT=="0/1",];
	ag=setNames(aggregate(maxaf~chr+start+stop, data=afd, FUN="mean"), c("chr", "start", "stop", "meanaf"));
	ag1=setNames(aggregate(maxaf~chr+start+stop, data=afd, FUN="length"), c("chr", "start", "stop", "ct"));
	ag=merge(ag, ag1);
	ag=merge( ag, casecov);
	ag=ag[ ag$meanaf<0.6 ,]; 
  } else {													#correlation for all targets
	ag=casecov;
  }
  ag=merge(ag, ctlcov);
  corrmat=cor(ag[, case], ag[, ctlnames]);
  keepctl=colnames(corrmat[1,corrmat[1,]>min.corr, drop=FALSE]);
  if(length(keepctl)<min.normals) {
    print("Error:  Low correlation between case and control coverage profiles.  Try reducing min.correlation or min.normals.");
	print(corrmat);
	quit(10);
  }
  ag=casecov=NULL;
  ctlcov=ctlcov[, c(infocols, keepctl)]; 
  write.table(ctlcov, file=ctlcombfile, row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE);
	
  counts=ctlcov[, keepctl, drop=FALSE];
  dd=ctlcov[, infocols];
  total=apply(counts, FUN=sum, MARGIN=2);
  geomean=exp(mean(log(total)));
  for(i in 1:ncol(counts)) {
	counts[,i]=counts[,i]/total[i]*geomean;													#normalize each sample by geometric mean
  }
  med=apply(counts, FUN=median, MARGIN=1);
  min=apply(counts, 1, FUN=min);
  max=apply(counts,1, FUN=max);
  min=log2((min+big.epsilon)/(med+big.epsilon));
  max=log2((max+big.epsilon)/(med+big.epsilon));
  dd=cbind(dd, med, min, max);
  write.table(dd, file=medfile, sep=",", row.names=FALSE, col.names=TRUE, quote=F);
  dd=dd[dd$med*100/dd$target.len>min.cov & dd$chr<=last.autosome,];							#masks low coverage here
  write.table(dd, file=paste(medfile, ".low_cov_flt.csv", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",");

}
	

		

gc_vaf_correct<-function(combfile, medfile, normfile, gcfile, names, frac.normal, last.autosome=22, afdata=NULL) {


  gc=read.table(file=gcfile, header=TRUE, comment.char="", stringsAsFactors=FALSE);
  gc=gc[, c(1:3,grep("gc", colnames(gc)))];
  colnames(gc)=c("chr", "start", "stop", "pct_gc");
  gc=gc[!duplicated(gc[, c("chr", "start", "stop")]),];
  dd=read.csv(combfile, stringsAsFactors=FALSE);
  auto=(as.data.frame(dd)$chr<=last.autosome);
  med=read.csv(medfile, stringsAsFactors=FALSE);
  med=med[, c(infocols, "med")];
  colnames(med)[6]="ctl.med";
  data=as.data.frame(merge(dd, med, by=infocols));
  data=data[order(data$chr, data$start),];
  data=data[, c(names, "ctl.med"), drop=F];
  Mtest=log2(data+epsilon);		#add small constant to avoid taking log(0), etc

  Mnorm=dd[, infocols];

  gc=as.data.frame(merge(Mnorm, gc), by=c("chr", "start", "stop"));
  gc=gc[order(gc$chr, gc$start),];

  load(afdata);

  for(i in 1:length(names)) {

    case=names[i];
    log.ratio=Mtest[, case]-Mtest[, "ctl.med"];
    ls=loess(log.ratio[auto]~gc[auto,]$pct_gc);
    pred=predict(ls, gc$pct_gc);
    log.ratio=log.ratio-pred;
    afd=af.list[[i]][af.list[[i]]$GT=="0/1",];
    ag=setNames(aggregate(maxaf~chr+start+stop, data=afd, FUN="mean"), c("chr", "start", "stop", "meanaf"));
    ag1=setNames(aggregate(maxaf~chr+start+stop, data=afd, FUN="length"), c("chr", "start", "stop", "ct"));
    ag=merge(ag, ag1);
    ag=merge(cbind(Mnorm, log.ratio), ag);
    ag=ag[ag$ct>1 & ag$meanaf<0.6 & ag$chr<=last.autosome,];						#recenter based on median ration in regions without LOH
    log.ratio=log.ratio-median(ag$log.ratio);
    Mnorm=cbind(Mnorm, log.ratio);
  }
  colnames(Mnorm)[6:ncol(Mnorm)]=names;
  write.table(Mnorm, file=normfile, quote=F, row.names=F, col.names=T, sep=",");
}



inv_set_norm_gc_loess<-function(combfile, medfile, normfile, gcfile, names, frac.normal, last.autosome=22) {

  gc=read.table(file=gcfile, header=TRUE, comment.char="", stringsAsFactors=FALSE);
  gc=gc[, c(1:3,grep("gc", colnames(gc)))];
  colnames(gc)=c("chr", "start", "stop", "pct_gc");
  gc=gc[!duplicated(gc[, c("chr", "start", "stop")]),];
  dd=read.csv(combfile, stringsAsFactors=FALSE);
  auto=(as.data.frame(dd)$chr<=last.autosome);
  med=read.csv(medfile, stringsAsFactors=FALSE);
  med=med[, c(infocols, "med")];
  colnames(med)[6]="ctl.med";
  data=as.data.frame(merge(dd, med, by=infocols));
  data=data[order(data$chr, data$start),];
  data=data[, c(names, "ctl.med"), drop=F];
  Mtest=log2(data+epsilon);		#add small constant to avoid taking log(0), etc

  count=floor(frac.normal*sum(auto));
  Mnorm=dd[, infocols];
  gc=as.data.frame(merge(Mnorm, gc), by=c("chr", "start", "stop"));
  gc=gc[order(gc$chr, gc$start),];

  for(case in names) {

    log.ratio=Mtest[, case]-Mtest[, "ctl.med"];
    ls=loess(log.ratio[auto]~gc[auto,]$pct_gc);
    pred=predict(ls, gc$pct_gc);
    log.ratio=log.ratio-pred;
    M1=cbind(case=Mtest[, "ctl.med"]+log.ratio, Mtest[,"ctl.med"]);
    M1.auto=M1[auto,];		                #normalize each individual case against control median

    total=nrow(M1.auto);																							#adapted from GRSN to find rank-invariant set
    idx<-1:total;   subSet<-1:total;
    discardNumber <- (total - count) / 10
    while (TRUE)	{
      total <- floor(max(total - discardNumber, count))
      M2 <- apply(M1.auto[idx, ], 2, rank)
      V2 <- apply(M2, 1, var)
      subSet <- order(V2, decreasing=FALSE)[1:total]     
      idx <- idx[subSet]
      if (total == count) break
    }

    ls=loess(M1.auto[idx, "case"]~M1.auto[idx, 2]);
    pred=predict(ls, M1[, 2]);
    pred[is.na(pred)]<-M1[is.na(pred),2];
    log.ratio=M1[,1]-pred;
    Mnorm=cbind(Mnorm, log.ratio);
  }
  colnames(Mnorm)[6:ncol(Mnorm)]=names;
  write.table(Mnorm, file=normfile, quote=F, row.names=F, col.names=T, sep=",");
}


  
annotate.with.genes<-function(norm.cov.file, target.bedfile, annot.bedfile) {
  ncov=read.csv(norm.cov.file, stringsAsFactors=FALSE);
  if(!is.null(annot.bedfile)){
    key=read.table(pipe(paste(BEDTools, "intersectBed -wa -u -a ", annot.bedfile, " -b ", target.bedfile, " | ", BEDTools, "closestBed -a ", target.bedfile, " -b stdin -t first ", sep="")));  
    key=key[, c(1:3,8)];
    colnames(key)=c("chrstr", "start", "stop", "gene");
    ncov=merge(ncov, key, by=c("chrstr", "start", "stop"), all.x=TRUE, sort=FALSE);
  } else {
    ncov=cbind(ncov, gene=ncov$chrstr);
  }
  ncov=ncov[order(ncov$chr, ncov$start),];
  write.table(ncov, file=norm.cov.file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",");
}

merge.score.call<-function(depthlist, aflist, final.list, j, p.cov.cutoff, p.af.cutoff, af.threshold=0.58, coverage.min.ratio=0.75, coverage.max.ratio=1.25, maxhet=0.98) {
  load(depthlist);
  afdata=NULL;
  ncontrols=5;
  depth=depth.list[[j]];
  if(!is.null(aflist)) {
    load(aflist);
    afdata=af.list[[j]];
  }
  load(final.list);

  if(!is.null(afdata)) {

    aggregate(maxaf~chr+start+stop, data=afdata, FUN=function(x) c(mn =sum(x[x<maxhet]), n=length(x[x<maxhet]), maxaf=max(x), minaf=min(x))) ->ag
    mat=matrix(lapply(ag[,4], "[", 1), ncol=4, byrow=FALSE)
    ag=as.data.frame(cbind(ag[,-4], mat))
    colnames(ag)=c("chr", "start", "stop", "sum", "ct", "maxmaf", "minmaf")
    for(col in c("sum", "ct", "maxmaf", "minmaf")) {
      ag[, col]=unlist(ag[,col]);
    }
    mm=merge(depth, ag, by=c("chr", "start", "stop"), all.x=TRUE)
    ag2=setNames(aggregate(ct~segno, data=mm, FUN=sum, na.rm=TRUE), c("segno", "seghetct"));
    ag3=setNames(aggregate(sum~segno, data=mm, FUN=sum, na.rm=TRUE), c("segno", "mu"));
    ag3=merge(ag2, ag3, by="segno")
    ag3[, "mu"]=ag3[, "mu"]/ag3[, "seghetct"]
    mm1=merge(mm, ag3, by="segno", all.x=TRUE);
    mm1=mm1[, !(colnames(mm1) %in% c("sum", "ct"))]
    mm1=cbind(mm1, p.af.t=rep(NA, nrow(mm1)));
    mm1[is.na(mm1$seghetct), "seghetct"]<-0;
    mm1=mm1[, c(colnames(depth),  "seghetct", "mu", "p.af.t", "minmaf", "maxmaf")];
  } else {
    mm1=cbind(depth, maxmaf=rep(NA, nrow(depth)), minmaf=rep(NA, nrow(depth)), seghetct=rep(NA, nrow(depth)), mu=rep(NA, nrow(depth)), p.af.t=rep(NA, nrow(depth)));
  }		   
  cname=colnames(mm1)[6];
  mm1=mm1[, c("chr", "chrstr", "start", "stop",  "target.len", cname, "med", "min", "max", "gene", "segno", "seglogratio", "seghetct", "mu", "p.af.t", "minmaf", "maxmaf")]
  mm1=mm1[order(mm1$chr, mm1$start),];

  final=mm1;  mm1=NULL;
  final=cbind(final, p.cov.np=rep(0, nrow(final)), p.cov.t=rep(0, nrow(final)), p.combined=rep(0, nrow(final)), cn=rep(0, nrow(final)), finalcn=rep(0, nrow(final)));

  for(seg in unique(final$segno)) {
        
    dd1=final[final$segno==seg,];
	log2ratio=dd1[,6];
	if(dd1$seglogratio[1]>log2(coverage.max.ratio)) {
	  ct1=sum(log2ratio>dd1$max);
	  p.cov.np=pbinom(nrow(dd1)-ct1, size=nrow(dd1), prob=1/(1+ncontrols));
      p.cov.t=safe.t.test(log2ratio, alternative="greater");          
	} else if (dd1$seglogratio[1]<log2(coverage.min.ratio)){
	  ct1=sum(log2ratio<dd1$min);
	  p.cov.np=pbinom(nrow(dd1)-ct1, size=nrow(dd1), prob=1/(1+ncontrols));
      p.cov.t=safe.t.test(log2ratio, alternative="less");
	} else {
	  p.cov.np=1;
	  p.cov.t=1;
	}
	p.af.t=dd1$p.af.t[1];
    if(!is.na(p.cov.np) && !is.na(p.cov.t)) {
	  fishsum=-2*(log(p.cov.t)+log(p.cov.np));
	  p.combined=pchisq(fishsum, df=4, lower.tail=FALSE);
	} else {
	  p.combined=min(p.cov.np, p.cov.t, na.rm=TRUE);
	}
	final[final$segno==seg, ]$p.cov.np=p.cov.np;
	final[final$segno==seg, ]$p.cov.t=p.cov.t;
	final[final$segno==seg, ]$p.combined=p.combined;

	paf=final[final$segno==seg, "p.af.t"][1];
	pcov=final[final$segno==seg, "p.cov.t"][1];
	mu=final[final$segno==seg, "mu"][1];
    pcomb=final[final$segno==seg, "p.combined"][1];
    seglogratio=final[final$segno==seg,  "seglogratio"][1]
	cn=cn1=2;

    if(!is.na(pcov) && pcov<p.cov.cutoff) {
	  cn=2*2^seglogratio;
	  if(is.na(mu) ||  mu>af.threshold) {
		cn1=2*2^seglogratio;
	  }
    }
    final[final$segno==seg, "cn"]<-cn;
    final[final$segno==seg, "finalcn"]<-cn1;

  }

  depth.list.final[[j]]<-final;
  save(depth.list.final, file=final.list);

}


plot.cn1106<-function(dd, restrict=NULL, thin=NULL, title) {

  layout(c(1,2), widths=c(1), height=c(1,1));
  dd=dd[order(dd$chr, dd$start),];
  if(!is.null(thin)) {
    ss=seq(from=1, to=nrow(dd), by=thin);
    dd=dd[ss,];
  }
  gg=dd$genenum;
  gg1=c(0, gg[-length(gg)]);
  cc=dd$chr;
  cc1=c(0, cc[-length(cc)])
  dd=cbind(dd, diff=gg1-gg, index=c(1:nrow(dd)), chrdiff=cc1-cc);
  vlines=dd[dd$diff!=0, "index"]
  darklines=dd[dd$chrdiff<0, "index"];
  chrlabs=dd[dd$chrdiff!=0, "chr"];
  labs=dd[dd$diff!=0, "gene"];
  dd=cbind(dd, log2ratio=dd[,6]);
  dd=cbind(dd, col=ifelse(dd$chr %% 2 ==0, "black", "grey"));
 
  dd$col=as.character(dd$col);
  par(mar=c(0.5,3,1,1),oma=c(0.5,1,0.5,1));
  if(is.null(restrict)) {
    plot(dd$log2ratio, axes=F, main=title, ylim=c((min(min(dd$log2ratio), log2(0.25)))-0.5, (max(max(dd$log2ratio), log2(4)))+0.5), cex=0.5, col="black"); #col=dd$col);
  } else {
	dd$log2ratio=pmin(dd$log2ratio, log2(restrict[2]));       #clip log ratio
    dd$log2ratio=pmax(dd$log2ratio, log2(restrict[1]));
    plot(dd$log2ratio, axes=FALSE, main=title, ylim=c( log2(restrict[1]),  log2(restrict[2])), cex=0.5, col="black");		
  }
  axis(2, labels=NA, tck=-0.02, cex.axis=0.75);
  axis(2, lwd=0, line=-0.4, las=1, cex.axis=0.75);
  box();

  mtext(side=2, "log2ratio", line=2);
  dd$cn=dd$cn/2;
  dd$finalcn=dd$finalcn/2;
  if(!is.null(restrict)) {
    dd$cn=pmin(dd$cn, restrict[2]);
    dd$cn=pmax(dd$cn, restrict[1]);
    dd$finalcn=pmin(dd$finalcn, restrict[2]);
    dd$finalcn=pmax(dd$finalcn, restrict[1]);
  }
  points(log2(dd$cn), col="orange", cex=0.5);
  points(log2(dd$finalcn), col="red", cex=0.5);
  

  for(j in darklines) {abline(v=j, col="darkgrey")}
  ytext=c(max(max(dd$log2ratio), log2(4) ));
   if(!is.null(restrict)) {
     ytext=c( log2(restrict[2]-1) );
  }
  xtext=0.5*(darklines+c(darklines[-1], nrow(dd)));
  text(xtext, ytext, labels=paste("chr",chrlabs), srt=90);
  abline(h=log2(0.75), col="green");
  abline(h=log2(1.25), col="green");
  abline(h=log2(0.5), col="orange");
  abline(h=log2(1.5), col="orange");

  hetct=pmin(dd$seghetct, 5);
	par(mar=c(3,3,0,1));
  #plot(dd$mu, axes=F, col=(hetct+1), ylim=c(0,1), cex=0.5);
  plot(dd$maxmaf, axes=F, col="grey", ylim=c(0,1), cex=0.3);
  #points(dd$maxmaf, col="grey", cex=0.5);
  points(dd$minmaf, col="grey", cex=0.3);
  points(dd$mu, col=(hetct+1), cex=0.5);
 
  abline(h=0.5, col="green");
  abline(h=0.58, col="green");
  pos=0.5*(vlines+c(vlines[-1], nrow(dd)));
  alt=c(TRUE, FALSE);
  alt1=c(FALSE, TRUE);
  axis(1, at=pos[alt], labels=labs[alt], las=2, cex.axis=0.5, tck=-0.005, mgp=c(3,0.3, 0));
  if(length(labs)>1) {
    axis(1, at=pos[alt1], labels=labs[alt1], las=2, cex.axis=0.5, tck=-0.005, mgp=c(3,2.2,0));
  }
	
  # axis(1, at=0.5*(vlines+c(vlines[-1], nrow(dd))),labels=labs,  las=2, cex.axis=0.5, tck=-0.005)
  axis(2, labels=NA, tck=-0.02, cex.axis=0.75);
  axis(2, lwd=0, line=-0.4, las=1, cex.axis=0.75);
  box();
  mtext(side=2, "maxaf", line=2);
  for(j in darklines) {abline(v=j, col="darkgrey")}
  legend("bottomleft", legend=c("0", "1", "2","3", "4", "5+"), col=c(1:6), pch=15);
}



safe.t.test<-function(...) {
    obj=try(t.test(...), silent=TRUE);
    if (is(obj, "try-error")) return(NA) else return(obj$p.value);
}

segment_norm<-function( norm.cov.file, range.file, outfile, names, minwid=3, segalpha=0.01, sdundo=2) {
  
  nsamps=length(names);
  depth.list=vector("list", nsamps);
  ncov=read.csv(norm.cov.file,  stringsAsFactors=FALSE);
  rg=read.csv(range.file,  stringsAsFactors=FALSE);
  ncov=merge(ncov, rg, by=infocols, sort=FALSE);
  ncov=ncov[order(ncov$chr, ncov$start),];
  ncov=cbind(ncov, gene=rep(0, nrow(ncov)));
  
  CNAobj=CNA(ncov[, names, drop=FALSE], ncov$chr, 0.5*(ncov$start+ncov$stop));
  smCNAobj=smooth.CNA(CNAobj);
  seg=segment(smCNAobj, undo.splits="sdundo", p.method="perm", alpha=segalpha, min.width=minwid, undo.SD=sdundo);
  
  for(j in 1:nsamps) {
    out=seg$out;
    out=out[out$ID==paste("Sample.",j, sep=""),];
    
    test=ncov[, c(infocols, names[j], "med", "min", "max", "gene")];
    test=cbind(test, segno=rep(0,nrow(test)), seglogratio=rep(0,nrow(test)));
    for(i in 1:nrow(out)) {
      sub=(test$chr==out$chrom[i])& ((test$start>=out$loc.start[i] & test$start<=out$loc.end[i])|(test$stop>=out$loc.start[i] & test$stop<=out$loc.end[i])| (test$start<out$loc.start[i] & test$stop>=out$loc.start[i]));
      nr=length(test[sub, "chr"]);
      test[sub, c("segno", "seglogratio")]<-cbind(rep(i, nr), rep(out$seg.mean[i], nr));
    }
    depth.list[[j]]<-test;
  }
  save(depth.list, file=outfile);
}


af.to.exons<-function(exon.target, snp.wingspan, last.autosome=22, sample.infile=NULL, vcf.file=NULL, min.cov=20, outfile, var.bedfile, vcftype, source_dir) {
  
  af.list=NULL;
  if(snp.wingspan>0) {
    bed=read.table(exon.target, stringsAsFactors=FALSE);
    bed[,2]=pmax(1, bed[,2]-snp.wingspan);
    bed[,3]=bed[,3]+snp.wingspan;
    write.table(bed, file=var.bedfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t");
  } else {
    var.bedfile=exon.target;
  }
  vcfs=c();
  if(!is.null(sample.infile)) {
    sample.info=read.table(sample.infile, stringsAsFactors=FALSE);
    vcfs=sample.info[,2];
  } else if (!is.null(vcf.file)) {
    vcfs=c(vcf.file);
  } else {
    print("Error: no allele frequency data.");
    quit(10);
  }
  af.list<-vector("list", length(vcfs));
  for(j in 1:length(vcfs)) {
    if(vcftype=="SAMTOOLS") {
  	  print(paste(BEDTools, "intersectBed -wa -u -a ", vcfs[j], " -b ", var.bedfile, " | ", BEDTools, "closestBed -a stdin -b ", exon.target, " -t first | perl ", source_dir, "/parse_vcf_st.pl ", min.cov, sep=""));
        vafs=read.csv(pipe(paste(BEDTools, "intersectBed -wa -u -a ", vcfs[j], " -b ", var.bedfile, " | ", BEDTools, "closestBed -a stdin -b ", exon.target, " -t first | perl ", source_dir, "/parse_vcf_st.pl ", min.cov, sep="")), stringsAsFactors=FALSE);
      } else if(vcftype=="VARSCAN") {
        print(paste( "cut -f 1-10 ", vcfs[j], " | intersectBed -wa -u -a stdin -b ", var.bedfile, " | ", BEDTools, "closestBed -a stdin -b ", exon.target, " -t first | perl ", source_dir, "/parse_vcf_gi1vs.pl ", min.cov, sep=""));
        vafs=read.csv(pipe(paste("cut -f 1-10 ", vcfs[j],  " | intersectBed -wa -u -a stdin  -b ", var.bedfile, " | ", BEDTools, "closestBed -a stdin -b ", exon.target, " -t first | perl ", source_dir, "/parse_vcf_gi1vs.pl ", min.cov, sep="")), stringsAsFactors=FALSE);
      } else {
		print(paste(BEDTools, "intersectBed -wa -u -a ", vcfs[j], " -b ", var.bedfile, " | ", BEDTools, "closestBed -a stdin -b ", exon.target, " -t first | perl ", source_dir, "/parse_vcf_gi1.pl ", min.cov, sep=""));
		vafs=read.csv(pipe(paste(BEDTools, "intersectBed -wa -u -a ", vcfs[j], " -b ", var.bedfile, " | ", BEDTools, "closestBed -a stdin -b ", exon.target, " -t first | perl ", source_dir, "/parse_vcf_gi1.pl ", min.cov, sep="")), stringsAsFactors=FALSE);
      }
          
      chrnum=vafs[ , "chr"];
      if(grepl("chr", vafs[1, "chr"])) {
        chrnum=substring(vafs[, "chr"], 4);
    }
    chrnum[chrnum=="X"]<-last.autosome+1;
    chrnum[chrnum=="Y"]<-last.autosome+2;
    chrnum=as.numeric(chrnum);
    vafs=cbind(vafs, chrnum, maxaf=pmax(vafs$maf, 1-vafs$maf));
    vafs=vafs[vafs$depth*vafs$maxaf<(vafs$depth-3),];
    vafs=vafs[vafs$GT=="0/1" | vafs$GT=="1/1",];
    af.list[[j]]<-vafs;
  } 
  save(af.list, file=outfile)
  
}

merge_depth_af<-function( depthlist, aflist, outlist, p.cov.cutoff, p.af.cutoff, useaf, af.threshold=0.58, coverage.min.ratio=0.75, coverage.max.ratio=1.25) { 

  load(depthlist);
  #if(useaf) { load(aflist);}
  nsamps=length(depth.list);
  depth.list.final<- vector("list", nsamps);
  save(depth.list.final, file=outlist);

  for(j in 1:nsamps) {
    if(useaf) {
      merge.score.call(depthlist, aflist, outlist, j, p.cov.cutoff, p.af.cutoff, af.threshold=0.58, coverage.min.ratio=0.75, coverage.max.ratio=1.25, maxhet=0.98);
    } else {
      merge.score.call(depthlist, NULL, outlist, j, p.cov.cutoff, p.af.cutoff, af.threshold=0.58, coverage.min.ratio=0.75, coverage.max.ratio=1.25, maxhet=0.98);
    }
  }
}


plot_depth_af.1213<-function(names, target.bedfile, tempfile, plotfile, regstarts, regstops, restrict, thin=NULL, restrictbed=NULL, png=FALSE) {
  
  if(!png) {
    pdf(file=paste(plotfile, ".pdf", sep="") ,width=10.5, height=7.5, paper="USr");
    par(ask=F, oma=c(0,0,0,0));
  
    if(!is.null(restrictbed)) {
      subspace=read.table(pipe(paste(BEDTools, "intersectBed -wa -u -a ",restrictbed, " -b ", target.bedfile, " | ", BEDTools, "closestBed -a ", target.bedfile, " -b stdin -t first ", sep="")));  
      subspace=subspace[!duplicated(subspace[, 1:3])];
    }
    nsamps=length(names);
  
    load(tempfile);	#regstarts=c(1,4,7,10,14, 18);	#regstops=c(3,6,9,13, 17, 24);
  
    for(i in 1:nsamps) {
      print(i);
    
      for( reg in 1:length(regstarts)) {
        dd=depth.list.final[[i]];
        if(!is.null(restrictbed)) {
          dd=dd[paste(dd[,1], dd[,2], dd[,3]) %in% paste(subspace[,1], subspace[,2], subspace[,3]),];
        }
        dd=dd[dd$chr>=regstarts[reg] & dd$chr<=regstops[reg],];
        dd=cbind(dd, genenum=rep(0, nrow(dd)));
        dd$genenum=as.numeric(as.factor(dd$gene));
        plot.cn1106(dd, restrict, thin, names[i]);
      }
    }
    dev.off()
  }
}



write_outfile<-function(names, outlist, outfile, sigoutfile, p.af.cutoff=0.05, p.cov.cutoff=0.05) {
  
  nsamps=length(names);
  allresults=NULL;
  sigresults=NULL;
  load(outlist);
  for(j in 1:nsamps) {
    name=names[j];
    dd=depth.list.final[[j]];
    results=cbind(sample=rep(name, nrow(dd)), dd[, c("chr", "chrstr", "start", "stop", "gene", "segno", "seghetct", "mu", "p.af.t", "p.cov.np", "p.cov.t", "p.combined", "cn", "finalcn")]);
    minag=aggregate(start~chr+chrstr+gene+segno, data=results, min)
    maxag=aggregate(stop~chr+chrstr+gene+segno, data=results, max)
    ag=merge(maxag, minag)
    results=results[, c("sample", "chr", "chrstr", "gene", "segno", "seghetct", "mu", "p.af.t", "p.cov.np", "p.cov.t", "p.combined", "cn", "finalcn")]
    results=results[!duplicated(results),]
    results=merge(results, ag);
    results=results[, c("sample", "chr", "chrstr", "start", "stop", "gene", "segno", "seghetct", "mu", "p.af.t", "p.cov.np", "p.cov.t", "p.combined", "cn", "finalcn")];
    results=results[order(results$chr, results$start),] 
    allresults=rbind(allresults, results);
    sigresults=rbind(sigresults, results[results[, "cn"]!=2,]);
  }
  write.table(allresults, outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",");
  write.table(sigresults, sigoutfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",");
}





    







	
	
