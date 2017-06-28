#!usr/perl -w

use strict;
use warnings;
use Parallel::ForkManager;

print "Working on Split operation\n";
fastq_split ($ARGV[0]);
print "Splitting Finished\n";

##############################################################################################################################################################################

sub fastq_split
{
	my ($input_file) = @_;
	print "Input is $input_file\n";
	my ($read1, $read2, $splitnumber, $maxlength, $working_dir, $sampleid, $r1len, $r2len, $threads);
	open (INFP, "$input_file") or die  "Can't open config at split step\n";
	while (my $l = <INFP>)
	{
		chomp ($l);
		my @records = split("\t", $l);
		if ($records[0] eq "READ1") { $read1 = $records[1]; }
		if ($records[0] eq "READ2") { $read2 = $records[1]; }
		if ($records[0] eq "READSPERFASTQ") { $splitnumber = $records[1]; }
		if ($records[0] eq "SAMPLEID") { $sampleid = $records[1]; }
		if ($records[0] eq "WORKING_DIR") { $working_dir = $records[1]; }
	}
	close (INFP);

	print "Read1:$read1\nRead2:$read2\nReadsPerFastQ:$splitnumber\nSampleID:$sampleid\nWorking_dir:$working_dir\n";
	my $total_sub_fastq = split_fastq_file($read1,$read2,$sampleid,$splitnumber,$working_dir);	
	system ("chmod 777 $working_dir/*.gz");
}	

#####################################################################################################################################################################################
sub split_fastq_file {
	my ($read1,$read2,$sampleid,$totalreadsfastq,$output_path) = @_;

	my $split_fastq_number = 1;
        my $total = 0;

        
	my $to_remove1 = $output_path."/".$sampleid."_1_*.fastq.gz";
        my $to_remove2 = $output_path."/".$sampleid."_2_*.fastq.gz";

        system("rm $to_remove1 $to_remove2");

        if ($read1=~ m/\.gz/) {
		open(INFP1,"gunzip -c $read1 |") || die "Could not open $read1 input fastq file\n";
	        open(INFP2,"gunzip -c $read2 |") || die "Could not open $read2 input fastq file\n"; 
	}
	else {
		open(INFP1,"<$read1") || die "Could not open $read1 input fastq file\n";
        	open(INFP2,"<$read2") || die "Could not open $read2 input fastq file\n";
	}

	my $outFile1 = $output_path."/".$sampleid."_1_".$split_fastq_number.".fastq.gz";
	my $outFile2 = $output_path."/".$sampleid."_2_".$split_fastq_number.".fastq.gz";

	open (OUTFP1, "| gzip -c > $outFile1");
	open (OUTFP2, "| gzip -c > $outFile2");

	print "Splitting files, files opened properly\n";
	my ($row1_1,$row2_1,$row3_1,$row4_1,$row1_2,$row2_2,$row3_2,$row4_2,$row4_2_1,$row4_1_1,$row2_1_1,$row2_2_1);
        while($row1_1=<INFP1>) {
                $row2_1=<INFP1>; $row3_1=<INFP1>; $row4_1=<INFP1>;
                $row1_2=<INFP2>; $row2_2=<INFP2>; $row3_2=<INFP2>; $row4_2=<INFP2>;

                chomp($row1_1); chomp($row2_1); chomp($row3_1); chomp($row4_1);
                chomp($row1_2); chomp($row2_2); chomp($row3_2); chomp($row4_2);

                if($total % $totalreadsfastq == 0 && $total > 0) {
                        close(OUTFP1); close(OUTFP2);
                        $split_fastq_number++; #Open the next file

                        $outFile1 = $output_path."/".$sampleid."_1_".$split_fastq_number.".fastq.gz";
                        $outFile2 = $output_path."/".$sampleid."_2_".$split_fastq_number.".fastq.gz";

                        open (OUTFP1, "| gzip -c > $outFile1");
                        open (OUTFP2, "| gzip -c > $outFile2");
		}

                print OUTFP1 "$row1_1\n$row2_1\n$row3_1\n$row4_1\n";
                print OUTFP2 "$row1_2\n$row2_2\n$row3_2\n$row4_2\n";

                $total++;
        }
        close(INFP1); close(INFP2); close(OUTFP1); close(OUTFP2); 

        print "Splitting of fastq file finished, total $split_fastq_number sub-fastq files generated\n";

	my $total_split_file = $output_path."/"."$sampleid.totalsubfastq.txt";

        open(OUTFP,">$total_split_file") || die "Could not open $total_split_file\n";
        print OUTFP $split_fastq_number,"\n";
        close(OUTFP);

        return($split_fastq_number);
}
########################################################################################################################################################################################
