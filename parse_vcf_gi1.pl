#!/usr/bin/perl                                                                                                                                                       

use strict;
my $mindepth=$ARGV[0];
my %nt=('A'=>1,
        'C'=>1,
        'T'=>1,
        'G'=>1);
print STDOUT "chr,pos,depth,GT,maf,start,stop\n";
while(my $line=<STDIN>) {
    chomp($line);
    if(!(substr($line, 0,1) eq "#")) {
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample, $chr0, $start0, $stop0, $locus)=split(/\s+/, $line);
        if($nt{$alt}==1 && $nt{$ref}==1) {
	    my $depth=0;
	    my $gt=0;
	    my $maf=0;
	    my $nref=0;
	    my $nalt=0;
	    my @keys=split(":", $format);
	    my @vals=split(":", $sample);
	    for(my $i=0; $i<@keys; $i++) {
		if($keys[$i] eq "DP") {
		    $depth=$vals[$i];
		} elsif($keys[$i] eq "GT") {
		    $gt=$vals[$i];
		} elsif($keys[$i] eq "AD") {
		    ($nref, $nalt)=split(",", $vals[$i]);
		    $maf=($nalt+0.5)/($nalt+$nref+0.5);
		}
	    }
	    if($depth>$mindepth) {
		print STDOUT "$chr,$pos,$depth,$gt,$maf,$start0,$stop0\n";
	    }
	}
    }
}
