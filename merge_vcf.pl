#!usr/perl -w

use strict;
use warnings;

$ENV{_JAVA_OPTIONS}="-Xmx8g";

my @list;

open (LIST, "$ARGV[0]");

while (my $l = <LIST>)
{
	chomp ($l);
	my @input = split ("\t", $l);
	push @list, "$input[0]/$input[1].recalibrated_snps_indels.vcf";
}
close LIST;

my $vcf_file = "--variant $list[0]";

for (my $i=1; $i<=$#list; $i++)
{
	$vcf_file = "$vcf_file --variant $list[$i]";	
}

print "Merging gVCF files\n";

my $cmd = "java -jar /scratch/ratan/humans/code/GATK/GenomeAnalysisTK.jar -T CombineGVCFs -R /scratch/ratan/humans/resource/human_g1k_v37_decoy.fasta -o $ARGV[1].cohort.g.vcf $vcf_file";
print "$cmd\n";
system ($cmd);
