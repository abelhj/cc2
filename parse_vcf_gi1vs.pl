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
        #my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample, $chr0, $start0, $stop0, $locus)=split(/\s+/, $line);
	my ($chr, $pos, $ref, $alt, $j0, $j1, $j2, $j3, $nref, $nalt, $chr0, $start0, $stop0, $locus)=split(/\s+/, $line);
        if($nt{$alt}==1 && $nt{$ref}==1) {
	    my $depth=$nref+$nalt;
	    my $maf=($nalt+0.5)/($nalt+$nref+0.5);
	    my $gt="0/1";
	    if($maf<0.03 || $maf>0.97) {
		$gt="1/1";
	    }
	    if($depth>$mindepth) {
		print STDOUT "$chr,$pos,$depth,$gt,$maf,$start0,$stop0\n";
	    }
	}
    }
}
