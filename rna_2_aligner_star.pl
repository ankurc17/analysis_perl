#!usr/perl -w

use strict;
use warnings;
use Parallel::ForkManager;

print "Using Star Aligner for RNA data\n";
star($ARGV[0]);
print "Star Aligner finished\n";
#########################################################################################################################################################################################

sub star 
{
	my ($config) = @_;

	#Variables initiated
	my ($sampleid, $working_dir, $star, $ref, $genome_dir, $threads);
	open (CFG, $config) or die "Can't open config file\n";
	while (my $l = <CFG>)
	{
		chomp ($l);
		my @input = split ("\t", $l);
		if ($input[0] eq "SAMPLEID") { $sampleid = $input[1]; }
		if ($input[0] eq "WORKING_DIR") { $working_dir = $input[1]; }
		if ($input[0] eq "STAR") { $star = $input[1]; }
		if ($input[0] eq "REF") { $ref = $input[1]; }
		if ($input[0] eq "INDEX_DIR") { $genome_dir = $input[1]; }
		if ($input[0] eq "TOTAL_THREADS") { $threads = $input[1]; }
	}
	close CFG;

	my $total_sub_fastq;
	open (INFP, "$working_dir/$sampleid.totalsubfastq.txt") or die "Can't open total sub fastq file\n";
        while(my $r = <INFP>) {
		chomp($r);
		$total_sub_fastq = $r;
	}
	close(INFP);
	print "File with sub-file count number opened successfully and now starting alignment\n";
	close INFP;

	my $manager = Parallel::ForkManager->new($threads);
	for (my $fno = 1; $fno<=$total_sub_fastq; $fno++)
	{
		run_star($star,$fno,$working_dir."/".$sampleid."_1_".$fno.".fastq.gz",$working_dir."/".$sampleid."_2_".$fno.".fastq.gz",$sampleid,$ref,$genome_dir);
		$manager->finish;
	}
	$manager->wait_all_children;
}
###############################################################################################################################################################################################
sub run_star
{
	my ($tool,$i,$r1,$r2,$sample,$reference,$directory)= @_;
#	print "$tool --genomeDir $directory --runThreadN 5 --outSAMstrandField intronMotif --outFileNamePrefix $sample.$i. --readFilesIn $r1 $r2 --readFilesCommand zcat\n";
	system ("$tool --genomeDir $directory --runThreadN 5 --outSAMstrandField intronMotif --outFileNamePrefix $sample.$i. --readFilesIn $r1 $r2 --readFilesCommand zcat");
}
