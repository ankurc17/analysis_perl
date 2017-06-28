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
	if ($info[0] eq "SPLIT_FLAG") { $f1 = $info[1]; }
	if ($info[0] eq "RNA_STAR") { $f2 = $info[1]; }
	if ($info[0] eq "MERGE_FLAG") { $f3 = $info[1]; }
	if ($info[0] eq "SAMBAMBA_FLAG") { $f4 = $info[1]; }
}
close(CONFIG);

my $SCRIPT_PATH = $Bin;

if ($f1==1) 
{
	system("perl $SCRIPT_PATH/rna_1_split.pl $Config_File"); 
}

if ($f2==1) 
{
	system("perl $SCRIPT_PATH/rna_2_aligner_star.pl $Config_File");
}

if ($f3==1)
{
	system("perl $SCRIPT_PATH/rna_3_merge_sam.pl $Config_File");
}

if ($f4==1)
{
	system("perl $SCRIPT_PATH/rna_4_sambamba.pl $Config_File");
}
