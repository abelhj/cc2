#!/usr/bin/Rscript --vanilla


source('/gscuser/habel/programs/cc/cnv_fns0725.R'); 
library(DNAcopy);
library(ff);
library(ffbase);
library("optparse");

option_list<-list(
  #make_option("--test_info", default=NULL, help="info file for test samples [default %default]"),
  make_option("--test_cov", default=NULL, help="single test coverage file [default %default]"),
  make_option("--test_vcf", default=NULL, help="single test vcf_file [default %default]"),
  make_option("--ctl_info", default=NULL, help="info file for ctl samples [default %default]"),
  make_option("--calc_case_coverage", default=FALSE, type="logical", help="calculate cov for cases [default %default]"),
  make_option("--calc_ctl_coverage", default=FALSE, type="logical", help="calculate cov for cases [default %default]"),
  make_option("--use_af", default=TRUE, type="logical", help="use allele frequency information [default %default]"),
  make_option("--min_coverage", default=50, type="integer", help="minimum coverage depth [default %default]"),
  make_option("--plot.only", default=FALSE, type="logical", help="[default %default]"),
  make_option("--outdir", default=NULL, help="output directory[default %default]"),
  make_option("--test_cov_dir", default=NULL, help="directory for test sample coverage (required only if calc_case_coverage=TRUE) [default %default]"),
  make_option("--ctl_cov_dir", default=NULL, help="directory for test sample coverage (required only if calc_ctl_coverage=TRUE) [default %default]"),
  make_option("--target.bedfile", default=NULL, help="bedfile used for coverage calculation [default %default]"),
  make_option("--annot.bedfile", default=NULL, help="bedfile used for annotation with gene names (optional) [default %default]"),
  make_option("--snp.wingspan", default=0, type="integer", help="distance to extend bedfile to look for intronic snvs [default %default]"),
  make_option("--plot.bedfile", default=NULL, help="bedfile specifying plot region [default %default]"),
  make_option("--outprefix", default=NULL, help="prefix for output files [default %default]"),
  make_option("--fa", default=NULL, help="fasta file [default %default]"),
  make_option("--plot.starts", default="1", type="character", help="comma separated list of start chr for plots (e.g., to plot all chr on same page, starts=1, stops=22.  To plot one chr/page, starts=1,2,3...22, stops=1,2,3...22), [default %default]"),
  make_option("--plot.stops", default="22",type="character", help="[default %default]"),
  make_option("--p.af.cutoff", default=1e-4, type="double", help="[default %default]"),
  make_option("--p.cov.cutoff", default=1e-4, type="double", help="[default %default]"),
  make_option("--af.threshold", default=0.58, type="double", help="[default %default]"),
  make_option("--coverage.min.ratio", default=0.75, type="double", help="[default %default]"),
  make_option("--coverage.max.ratio", default=1.25, type="double", help="[default %default]"),
  make_option("--minwid", default=5, type="integer", help="min number of regions in segment [default %default]"),
  make_option("--segalpha", default=0.01, type="double", help="alpha for CBS [default %default]"),
  make_option("--sdundo", default=2, type="double", help="sdundo for CBS [default %default]"),
  make_option("--restrict", default=NULL, type="character", help="restict y-axis in plot, [default %default]"), 
  make_option("--last.autosome", default=22, type="integer", help="number of last autosome [default %default]"),
  make_option("--frac.normal", default=0.5, type="double", help="fraction of targeted space assumed normal [default %default]"),
  make_option("--thin.rate", default=NULL, type="integer", help="thinning rate for plot"),
  make_option("--plot.id", default="", type="character", help="new plot label"),
  make_option("--png", default=FALSE, type="logical", help="plot is .png file [default %default]"),
  make_option("--vcftype", default=NULL, type="character", help="TUMOR, NORMAL, VARSCAN, SAMTOOLS  or NULL (single sample VCF)"),
  make_option("--min.var.cov", default=50, type="integer", help="min coverage for including AF [default %default]"),
  make_option("--vafs_normalize", default=FALSE, type="logical", help="use vafs to determine copy neutral regions for normalization"),
  make_option("--min.normal.corr", default=0.70, type="double", help="min correlation between case and control coverage [default %default]"),
  make_option("--min.num.normals", default=3, type="integer", help="min number of controls to pool")
);
 
opt=parse_args(OptionParser(option_list=option_list));
vafs_normalize=opt$vafs_normalize;
if(vafs_normalize && !opt$use_af) {
  print("Can't do VAF-based normalization with use_af=FALSE.");
  quit(10);
}

if(is.null(opt$test_info) && is.null(opt$test_cov)) {
  print("test_info and test cov file both missing");
  quit(10);
}
if(is.null(opt$ctl_info)) {
  print("ctl info file missing.")
  quit(10);
}
if(is.null(opt$outdir)) {
  print("must specify outdir.");
  quit(10);
}
if(is.null(opt$target.bedfile)) {
  print("target bedfile missing.");
  quit(10);
}
if(is.null(opt$fa)) {
  print("fasta file missing.")
  quit(10);
}
plot.starts=as.numeric(unlist(strsplit(opt$plot.starts, split=",")))
plot.stops=as.numeric(unlist(strsplit(opt$plot.stops, split=",")))
restrict=opt$restrict;
if(!is.null(restrict)) {
  restrict=as.numeric(unlist(strsplit(opt$restrict, split=",")));
  print(restrict);
}

tempdir="";
sample_names=c();
control_names=c();
epsilon=0.25;
big.epsilon=2;
infocols=c("chr", "chrstr", "start", "stop", "target.len");
BEDTools="";

ctlcomb=ctlmedfile=combfile=combfile.flt=ctlmedfile.flt=ctlcombfile.flt=mrnfile=temp1=temp2=temp3=NULL;



dir.create(opt$outdir);
tempdir=paste(opt$outdir, "/temp_", opt$outprefix, sep="");
dir.create(tempdir);

ctlcombfile=paste(tempdir, '/ctl_comb.csv', sep="");
ctlmedfile=paste(tempdir, '/ctl.med.csv', sep="");
combfile=paste(tempdir, '/comb.csv', sep="");
combfile.flt=paste(combfile, ".low_cov_flt.csv", sep="");
ctlmedfile.flt=paste(ctlmedfile, ".low_cov_flt.csv", sep="");
ctlcombfile.flt=paste(ctlcombfile, ".low_cov_flt.csv", sep="");
mrnfile=paste(tempdir, '/test.mrn.csv', sep="");
gcfile=paste(tempdir, '/gc.txt', sep="");

temp1=paste(tempdir, "/temp1.RData", sep="");
temp2=paste(tempdir, "/temp2.RData", sep="");
temp3=paste(tempdir, "/temp3.RData", sep="");
 
  
print("combining case coverage");
case.names=combine_cov(sample.infile=opt$test_info, cov.file=opt$test_cov, cov.comb=combfile, last.autosome=opt$last.autosome);

if(opt$plot.only) {
  plot_depth_af.1213( names=case.names, opt$target.bedfile, tempfile=temp3, plotfile=paste(opt$outdir, "/cc.", opt$outprefix, opt$plot.id, sep=""), plot.starts, plot.stops, restrict=restrict, thin=opt$thin.rate, restrictbed=opt$plot.bedfile, opt$png);
} else {
  print("combining control coverage");
  combine_cov(sample.infile=opt$ctl_info, cov.comb=ctlcombfile, cov.file=NULL, last.autosome=opt$last.autosome);
  if(opt$use_af) {
    print("calculating vafs");
    af.to.exons(exon.target=opt$target.bedfile, snp.wingspan=opt$snp.wingspan, last.autosome=opt$last.autosome, sample.infile=opt$test_info, vcf.file=opt$test_vcf, min.cov=opt$min.var.cov, outfile=temp2, var.bedfile=paste(tempdir, "/temp.bed", sep=""), vcftype=opt$vcftype);
  } 
  select.normals.get.med(combfile, ctlcombfile, temp2, opt$min.normal.corr, opt$min.num.normals, opt$min_coverage, last.autosome=opt$last.autosome, ctlmedfile);
    
  print("masking low coverage regions");
  mask_low_coverage(ctlmedfile, combfile, opt$min_coverage, opt$last.autosome);
  
  print("calculating gc content and normalizing coverage data");
  get.gc.content(opt$fa, opt$target.bedfile, gcfile);
 

  if(!vafs_normalize) {
    print("gc correction and rank invariant set normalizaton");
    inv_set_norm_gc_loess(combfile.flt, ctlmedfile.flt, mrnfile, gcfile, case.names, opt$frac.normal, opt$last.autosome);
  } else {
    print("gc correction and VAF-based normalization");
    gc_vaf_correct(combfile.flt, ctlmedfile.flt, mrnfile, gcfile, case.names, opt$last.autosome, afdata=temp2);
  }
 
  annotate.with.genes(mrnfile, opt$target.bedfile, opt$annot.bedfile);
  print("segmentation")
  segment_norm(mrnfile, ctlmedfile.flt, temp1, case.names, minwid=opt$minwid, segalpha=opt$segalpha, sdundo=opt$sdundo);
  
  print("merging coverage data with vafs");
  merge_depth_af( temp1, temp2, temp3, opt$p.cov.cutoff, opt$p.af.cutoff, opt$use_af, opt$af.threshold, opt$coverage.min.ratio, opt$coverage.max.ratio); 
  print("plotting");
  plot_depth_af.1213( names=case.names, opt$target.bedfile, tempfile=temp3, plotfile=paste(opt$outdir, "/cc.", opt$outprefix, opt$plot.id,  sep=""), plot.starts, plot.stops, restrict=restrict, thin=opt$thin.rate, restrictbed=opt$plot.bedfile, opt$png);
  print("writing output file");
  write_outfile(names=case.names, temp3, outfile=paste(opt$outdir, "/cc.results.", opt$outprefix, ".csv", sep=""), sigoutfile=paste(opt$outdir, "/cc.sig.results.", opt$outprefix, ".csv", sep=""), opt$p.af.cutoff, opt$p.cov.cutoff);
}




	

