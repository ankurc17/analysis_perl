#!usr/perl -w

use strict;
use warnings;
use Data::Dumper;

my @list = ("clean","splitters","discordant");

my ($sampleid,$working_dir);
my %insert = ();
my %mapping = ();

my ($size, $quality,$readper_clean, $readper_discordant, $readper_splitters);

open (INFP, "$ARGV[0]");

while (my $l = <INFP>)
{
        chomp ($l);
        my @input = split ("\t", $l);
        if ($input[0] eq "SAMPLEID") { $sampleid = $input[1]; }
        elsif ($input[0] eq "WORKING_DIR") { $working_dir = $input[1]; }
}

my $temp1 = 0;
my $temp = 0;

my $insert_file = '';
my $mapping_file = '';

for (my $i=0; $i<=$#list; $i++)
{
	$insert_file = "$working_dir/$sampleid.$list[$i].InsertSizeDistribution.txt";

	open (ISD, "$insert_file") or die "Can't open $insert_file\n";
	while (my $line1 = <ISD>)
	{
	        chomp($line1);
	        my @line_data1=split(/\t/,$line1);
	        $temp1 += $line_data1[1];
	}
	
	close ISD;
	
	open (ISD, "$insert_file") or die "Can't open $insert_file\n";
	while (my $line1 = <ISD>)
	{
        	chomp($line1);
		my @line_data1=split(/\t/,$line1);	
		my $b = ($line_data1[1]/$temp1)*100;
		$insert{$line_data1[0]}{$i} = "$line_data1[0]\t$b";
	}
	close ISD;
	$temp1 = 0;
	$mapping_file = "$working_dir/$sampleid.$list[$i].MappingQualityDistribution.txt";

	open(MQD, "$mapping_file") or die "Can't open $mapping_file\n";
	while (my $line1 = <MQD>)
	{
		chomp($line1);
		my @line_data1=split(/\t/,$line1);
		if ($line_data1[0]=~ m/[0-9]+/)
		{
			$temp+=$line_data1[1];
		}
	}
	close MQD;

	open(MQD, "$mapping_file") or die "Can't open $mapping_file\n";
	while (my $line1 = <MQD>)
	{
		chomp($line1);
		my @line_data1=split(/\t/,$line1);
		my $b2 = ($line_data1[1]/$temp)*100;
		$mapping{$line_data1[0]}{$i} = "$line_data1[0]\t$b2";
	}
	close MQD;
	$temp = 0;
}

open (OUTFP, ">$working_dir/$sampleid.insert.txt");
print OUTFP "Insert Size\tReadPercentage_Clean\tReadPercentage_Splitters\tReadPercentage_Discordant\n";

foreach my $key1 (sort {$a <=> $b} keys %insert) {
	print OUTFP "$key1\t";
	for (my $i=0; $i<=$#list; $i++) {
		if (exists($insert{$key1}{$i})) {
			my @a= split("\t",$insert{$key1}{$i});
			if ($i==2) { print OUTFP "$a[1]"; }
			else {print OUTFP "$a[1]\t"; }
		}
	}
	print OUTFP "\n";
}
close OUTFP;

open (OUTFP, ">$working_dir/$sampleid.mapping.txt");
print OUTFP "MQ\tReadPercentage_Clean\tReadPercentage_Splitters\tReadPercentage_Discordant\n";
foreach my $key1 ( sort {$a <=> $b} keys %mapping) {
	print OUTFP "$key1\t";
	for (my $i=0; $i<=$#list; $i++) {
		if (exists($mapping{$key1}{$i})) {
			my @a= split("\t",$mapping{$key1}{$i});
			if ($i==2) { print OUTFP "$a[1]"; }
			else {print OUTFP "$a[1]\t"; }
		}
	}
	print OUTFP "\n";
}
close OUTFP;
