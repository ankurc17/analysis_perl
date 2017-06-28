#!usr/perl -w

use strict;
use warnings;
use FindBin qw($Bin);

my $file = $ARGV[0];
my $path = $ARGV[1];
my $fastq_path = $ARGV[2];

open (CFG, "$file") or die "Can't open that file!\n";
while (my $l = <CFG>)
{
	chomp ($l);
	my ($line, $type) = split ("\t", $l);

	my $working_dir = "$path/$line";
	
	system ("mkdir $working_dir");

	my ($r1, $r2, $inputpath);

	if(defined $path){ $inputpath = "$fastq_path/"; }
	else { print "Provide the path for fastq files. Exiting the process\n"; exit; }

	opendir(INDIR,$inputpath) or die ("couldn't open dir");
	my @files = readdir INDIR;

	foreach my $r(@files){
		chomp($r);
		if($r=~m/$line/){
			if($r=~m/_1.fastq.gz$/){ $r1 = "$inputpath/$r"; }
			elsif($r=~m/_2.fastq.gz$/){ $r2 = "$inputpath/$r"; }
		}
	}

	open (SHELL, ">>$working_dir/$line.sh");	
	
	print SHELL "#!/bin/bash\n#$ -cwd\n#$ -V\n\n## Join the standard error and the standard output into 1 file output\n#\$ -j y\n\n#\$ -q normal.q\n#\$ -pe smp 8\n";
	
	if ($type eq "DNA") { 
		print SHELL "perl /$Bin/Run_dna_pipeline.pl $working_dir/config.txt\n";	
		system ("cp /$Bin/dna_template.txt $working_dir/config.txt");
	}
	else
	{
		print SHELL "perl /$Bin/Run_rna_pipeline.pl $working_dir/config.txt\n";
		system ("cp /$Bin/rna_template.txt $working_dir/config.txt");
	}

	my ($count, $len_r1, $len_r2) = 0;
	
	my $value1 = qx (zcat $r1 | head -2 | sed -n '2p' | wc -L);
	my $value2 = qx (zcat $r2 | head -2 | sed -n '2p' | wc -L);

	my $len;
	if ($value1 > $value2) { $len  = $value1+1; } else { $len = $value2+1;}
	open (INDCFG, ">>$working_dir/config.txt") or die "Can't open $working_dir/config.txt\n";
	print INDCFG "WORKING_DIR\t$working_dir\nREAD1\t$r1\nREAD2\t$r2\nSAMPLEID\t$line\n";
	print INDCFG "R1LENGTH\t$len\nR2LENGTH\t$len";
		
	if (-e "$working_dir/config.txt") { print "File $working_dir/config.txt created\n"; }
}


close (CFG);
