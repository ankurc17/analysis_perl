#!usr/perl -w

use strict;
use warnings;
use Parallel::ForkManager;
print "Starting adapter trimming operation";

remove_adapter($ARGV[0]);

print "Adapter trimming sucessfully completed\n";
################################################################################################################################################

sub remove_adapter {
	my ($config) = @_;
	my @info = ();
	my($path,$working_dir,$sampleid,$threads) = '';
	my $line ='';
	open (INFP, $config) or die "Can't open $config";
	print "\nInput parameter for adaptor step from config file....\n";
	print "---------------------------------------------------------\n";
	while ($line = <INFP>) {
		chomp ($line);
		@info = split(/\t/,$line);
		if($info[0] eq "TRIM_GALORE") { $path = $info[1]; }
		elsif($info[0] eq "WORKING_DIR") { $working_dir = $info[1]; }
		elsif($info[0] eq "SAMPLEID") { $sampleid = $info[1]; }
		elsif($info[0] eq "TOTAL_THREADS") { $threads = $info[1]; }
	}
	my (@a,$i);
	
	my $TOTAL_SUB_FASTQ;
        my $total_subfastq_file = $working_dir."/"."$sampleid.totalsubfastq.txt";

        open(INFP,"<$total_subfastq_file") || die "Could not open $total_subfastq_file input file\n";
        while(my $r = <INFP>) {
                chomp($r);
                $TOTAL_SUB_FASTQ = $r;
        }
        close(INFP);

       print "File with sub-file count number opened successfully\n";

	print "Running adapter trimming program\n";
	print "Program Path: $path\n";
	print "working Dir: $working_dir\n";
	print "Sample Name: $sampleid\n";
	print "Adapter Threads: $threads\n";
	print "Total Sub fastq: $TOTAL_SUB_FASTQ\n";
	print "\n";
#$OUTPUT_PATH."/".$sampleid."_1_".$split_fastq_number.".fastq.gz";
	my $manager = Parallel::ForkManager->new($threads);
	for(my $fno=1; $fno<=$TOTAL_SUB_FASTQ; $fno++) {
		$manager->start and next;
		fastq_mcf_trimming($path,$working_dir."/".$sampleid."_1_".$fno.".fastq.gz",$working_dir."/".$sampleid."_2_".$fno.".fastq.gz",$working_dir);	
	$manager->finish;
        }
        $manager->wait_all_children;
        
	my $to_remove = $working_dir."/"."*skip";
#	system("rm $to_remove");

	print "Completed trimming of fastq files...\n";
}

#########################################################################################################################################################
###                                      Run adapter trimming program
##########################################################################################################################################################
sub fastq_mcf_trimming {
	my($TRIM_GALORE,$fastq_1,$fastq_2, $working_dir) = @_;
	print "Running adapter trimming on $fastq_1 and $fastq_2\n";
	print "$TRIM_GALORE -q 20 --gzip -o $working_dir --paired $fastq_1 $fastq_2\n";
	system("$TRIM_GALORE -q 20 --gzip -o $working_dir --paired $fastq_1 $fastq_2");
}
