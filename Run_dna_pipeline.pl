#!/usr/bin/perl -w

use strict;
use warnings;
use FindBin qw($Bin);
use POSIX qw(strftime);

#Initialzing the variables

my ($Sample,$Path,$f1,$f2,$f3,$f4,$f5,$f6,$f7);

#Inputing the config file

my $Config_File=$ARGV[0];
open (CONFIG,$Config_File) or die "Can't open $Config_File\t$!\n";

while (my $line=<CONFIG>)
{
	chomp($line);
	my (@info)=split(/\t/,$line,999);
	if ($info[0] eq "SAMPLE_NAME") { $Sample = $info[1]; } 
	if ($info[0] eq "WORKING_DIR") { $Path = $info[1]; }
	if ($info[0] eq "SPLIT_FLAG") { $f1 = $info[1]; }
	if ($info[0] eq "ADAPTRIM_FLAG") { $f2 = $info[1]; }
	if ($info[0] eq "ALIGN_FLAG") { $f3 = $info[1]; }
	if ($info[0] eq "SAMBAMBA_FLAG") { $f4 = $info[1]; }
	if ($info[0] eq "SAMTOOLS_FLAG") { $f5 = $info[1]; }
	if ($info[0] eq "VCF_FLAG") { $f6 = $info[1]; }
	if ($info[0] eq "STATS_FLAG") { $f7 = $info[1]; }
}
close(CONFIG);

my $SCRIPT_PATH = $Bin;

if ($f1==1) 
{
	system("perl $SCRIPT_PATH/dna_1_split.pl $Config_File"); 
}

if ($f2==1) 
{
	system("perl $SCRIPT_PATH/dna_2_trimming.pl $Config_File");
}

if ($f3==1)
{
	system("perl $SCRIPT_PATH/dna_3_sam_generation.pl $Config_File");
	system("perl $SCRIPT_PATH/dna_3.5_merge_sam.pl $Config_File");
}

if ($f4==1)
{
	system("perl $SCRIPT_PATH/dna_4_sambamba.pl $Config_File");
}

if ($f5==1)
{
	system("perl $SCRIPT_PATH/dna_4_samtools.pl $Config_File");
}

if ($f6==1)
{
	system("perl $SCRIPT_PATH/dna_5_gvcf_generation.pl $Config_File");
}

if ($f7==1)
{
	system("perl $SCRIPT_PATH/dna_report_1_fastq_summary.pl $Config_File");
	system("perl $SCRIPT_PATH/dna_report_2_trim_summary.pl $Config_File");
	system("perl $SCRIPT_PATH/dna_report_3_nbase.pl $Config_File");
        system("perl $SCRIPT_PATH/dna_report_4_alignment_summary.pl $Config_File");
        system("perl $SCRIPT_PATH/dna_report_5_insert_mapping.pl $Config_File");
        system("perl $SCRIPT_PATH/dna_report_6_graph.pl $Sample $Path");
	system("perl $SCRIPT_PATH/dna_report_7_final_report.pl $Config_File");
}
