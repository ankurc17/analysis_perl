#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;

print "Starting alignment operation\n";
run_alignment($ARGV[0]);
print "Alignment operation finished successfully\n";

################################################################################################################################################

sub run_alignment {
	my ($config) = @_;
	my @info = ();
	my($ALIGNER,$bwa_path,$samblaster,$ref,$working_dir,$threads,$sampleid) = '';
	my $line ='';
	open (INFP, $config) or die "Can't open $config";
	print "\nInput parameter for Alignment step from config file....\n";
	print "---------------------------------------------------------\n";
	while ($line = <INFP>) {
		chomp ($line);
		@info = split(/\t/,$line);
		if($info[0] eq "BWA") {$bwa_path = $info[1]; }
		elsif($info[0] eq "SAMBLASTER") { $samblaster = $info[1]; }
		elsif($info[0] eq "WORKING_DIR") { $working_dir = $info[1]; }
		elsif($info[0] eq "SAMPLEID") { $sampleid = $info[1]; }
		elsif($info[0] eq "TOTAL_THREADS") { $threads = $info[1]; }
		elsif($info[0] eq "REF") { $ref = $info[1]; }
	}

	print "Aligner: $ALIGNER\n";
	print "Aligner Path: $bwa_path\n";
	print "Working Dir: $working_dir\n";
	print "Sample Name: $sampleid\n";
	print "Reference: $ref\n";

	my $TOTAL_SUB_FASTQ;
	my $manager;
	my $total_subfastq_file = $working_dir."/"."$sampleid.totalsubfastq.txt";
	
	open(INFP,"<$total_subfastq_file") || die "Could not open $total_subfastq_file input file\n";
	while(my $r = <INFP>) {
		chomp($r);
		$TOTAL_SUB_FASTQ = $r;
	}
	close(INFP);

	print "File with sub-file count number opened successfully and now starting alignment\n";
	
	$manager = Parallel::ForkManager->new($threads);
	for(my $fno = 1; $fno <= $TOTAL_SUB_FASTQ; $fno++) {
		$manager->start and next;
		run_samblaster($bwa_path,$samblaster,$working_dir."/".$sampleid."_1_".$fno."_val_1.fq.gz",$working_dir."/".$sampleid."_2_".$fno."_val_2.fq.gz",$sampleid,"$working_dir/$sampleid.$fno.overall.sam",$ref,"$working_dir/$sampleid.$fno.discordant.sam","$working_dir/$sampleid.$fno.splitters.sam","$working_dir/$sampleid.$fno.unmapped.fq","$working_dir/$sampleid.$fno.clean.sam");
		$manager->finish;
	}
	$manager->wait_all_children;
}
#############################################################################################################################################################################################################

sub run_samblaster {
	my ($bwa,$samblaster,$fastq_1,$fastq_2,$sample,$sam,$index,$discordant_sam,$splitters_sam,$unmapped_fastq,$clean_sam) = @_;
	my $label = "\@RG\tID:$sample\tSM:$sample\tLB:lib\tPU:1";
	print "Running alignment for $fastq_1 $fastq_2\n";
	system("$bwa mem -t 5 -Y -R '$label' $index $fastq_1 $fastq_2 > $sam");
	system("$samblaster --addMateTags -d $discordant_sam -s $splitters_sam -u $unmapped_fastq -e -i $sam -o $clean_sam");
}

#############################################################################################################################################################################################################
