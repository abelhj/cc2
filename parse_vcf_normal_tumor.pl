#!/usr/bin/perl                                                                                                                                                       

use strict;
my $tn=$ARGV[0];
my $mindepth=$ARGV[1];
my %nt=('A'=>1,
        'C'=>1,
        'T'=>1,
        'G'=>1);
print STDOUT "chr,pos,depth,GT,maf,locus\n";
while(my $line=<STDIN>) {
    chomp($line);
    if(!(substr($line, 0,1) eq "#")) {
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $Nsample, $Tsample, $chr0, $start0, $stop0, $locus)=split(/\s+/, $line);
        if($nt{$alt}==1 && $nt{$ref}==1) {
	    my $depth=0;
	    my $gt=0;
	    my $maf=0;
	    #my $nref=0;
	    #my $nalt=0;
	    my ($nreff, $nrefr, $naltf, $naltr)=();
	    my @keys=split(":", $format);
	    my @vals=();
	    if($tn eq "TUMOR") {
		@vals=split(":", $Tsample);
	    } elsif ($tn eq "NORMAL") {
		@vals=split(":", $Nsample);
	    }
	    for(my $i=0; $i<@keys; $i++) {
		if($keys[$i] eq "DP") {
		    $depth=$vals[$i];
		} elsif($keys[$i] eq "GT") {
		    $gt=$vals[$i];
		} elsif($keys[$i] eq "FA") {
		    #($nreff, $nrefr, $naltf, $naltr)=split(",", $vals[$i]);
		    #$maf=($naltr+$naltf+0.5)/($naltr+$naltf+$nreff+$nrefr+0.5);
            $maf=$vals[$i];
		}
	    }
	    if(($depth>$mindepth) && !($gt eq "0/0")) {
		print STDOUT "$chr,$pos,$depth,$gt,$maf,$locus\n";
	    }
	}
    }
}
