use strict;
use warnings;

my %comp = ("A" => "T",
	    "C" => "G",
	    "G" => "C",
	    "T" => "A",
	    "N" => "N");

my ($inputf,$outputf,$path_to_samtools) = @ARGV;
die "perl $0 [input bam file] [output bam file] <path to samtools>\n" 
    if @ARGV < 2;


if(!defined($path_to_samtools))
{
    $path_to_samtools = "";
}
else
{
    $path_to_samtools .= "/"; 
}

open(DATA,"$path_to_samtools"."samtools view -h $inputf|") or 
    die "cannot open $inputf\n";
open(OUT,"|$path_to_samtools"."samtools view -S -b - > $outputf") or
    die "cannot create $outputf\n";

while(<DATA>)
{
    chomp;
    next if length($_) == 0;
    if(substr($_,0,1) eq '@')
    {
	print OUT $_,"\n";
	next;
    }

    #Split
    my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$RNEXT,$PNEXT,$TLEN,$SEQ,$QUAL,@other) = split(/\t/,$_);
    
    #If not R2, print it out and skip.
    if( ($FLAG & 128) == 0)
    {
	print OUT $_,"\n";
	next;
    }
    #R2
    if( ($FLAG & 16) == 0)#Forward strand
    {       
	$FLAG += 16;#Flip
    }
    else#Revsersed strand
    {
	$FLAG -= 16;#Flip
    }
    #Reverse complement sequence and revserse quality
    #$SEQ =~ tr/[a-z]/[A-Z]/;
    #$SEQ =~ tr/[ACGT]/[TGCA]/;
    
    print OUT join("\t",($QNAME,$FLAG,$RNAME,$POS,$MAPQ,
			 $CIGAR,$RNEXT,$PNEXT,$TLEN,
			 reverse($SEQ),reverse($QUAL),@other)),"\n";   
}
