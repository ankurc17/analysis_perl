#!usr/perl -w

use strict;
use warnings;
use Parallel::ForkManager;

print "Working on Adap Trim FastQC Stats Generation\n";
fastq_summary ($ARGV[0]);
print "Adap Trim FastQC Finished\n";

##############################################################################################################################################################################

sub fastq_summary
{
	my ($input_file) = @_;
	print "Input is $input_file\n";
	my ($read1, $read2, $splitnumber, $maxlength, $working_dir, $sampleid, $r1len, $r2len, $threads);
	open (INFP, "$input_file") or die  "Can't open config at split step\n";
	while (my $l = <INFP>)
	{
		chomp ($l);
		my @records = split("\t", $l);
		if ($records[0] eq "SAMPLEID") { $sampleid = $records[1]; }
		if ($records[0] eq "WORKING_DIR") { $working_dir = $records[1]; }
		if ($records[0] eq "R1LENGTH") { $r1len = $records[1]; }
		if ($records[0] eq "R2LENGTH") { $r2len = $records[1]; }
		if ($records[0] eq "TOTAL_THREADS") { $threads = $records[1]; }
	}
	close (INFP);

	my $total_sub_fastq;
	my $total_subfastq_file = $working_dir."/"."$sampleid.totalsubfastq.txt";
#Ankur_1_299_val_1.fq.gz
	open(INFP,"<$total_subfastq_file") || die "Could not open $total_subfastq_file input file\n";
        while(my $r = <INFP>) { chomp($r); $total_sub_fastq = $r; }
	close(INFP);

	print "File with sub-file count number opened successfully\n";
	
	print "SampleID:$sampleid\nWorking_dir:$working_dir\nRead1Length:$r1len\nRead2Length:$r2len\n";
	
	if ($r1len > $r2len) { $maxlength = $r1len; }
	else { $maxlength = $r2len; }
	my $manager = Parallel::ForkManager->new($threads);
	for(my $fno=1; $fno<=$total_sub_fastq; $fno++) {
		$manager->start and next;
		calculate_fastq_summary("$working_dir/$sampleid\_1_$fno\_val_1.fq.gz","$working_dir/$sampleid\_2_$fno\_val_2.fq.gz",$fno,$working_dir,$maxlength,$sampleid);
		$manager->finish;
	}
	$manager->wait_all_children;
	print "Summary Generation Finished.\n";
	combine_summary ( $total_sub_fastq,$working_dir,$sampleid,$maxlength);

	for(my $fno = 1; $fno <= $total_sub_fastq; $fno++) {
		if ((!-z "$working_dir/$sampleid.AT.FastqSummary.txt") || (!-z "$working_dir/$sampleid.FastqSummary.txt"))
                {
			system ("rm $working_dir/$sampleid\_1_$fno\_val_1.fq.gz $working_dir/$sampleid\_2_$fno\_val_2.fq.gz");
			system ("rm $working_dir/$sampleid\_1_$fno.fastq.gz $working_dir/$sampleid\_2_$fno.fastq.gz");
                }
                else
                {
                        print "Files for are empty. Can't delete the fastq files\n";
                }
        }
}	

#####################################################################################################################################################################################

sub calculate_fastq_summary {
	my ($r1_file,$r2_file,$fno,$working_dir,$readlength,$sampleid) = @_;
	my $QUALITY_ENCODING = 33;
	
	#---------------------------Initializing variables-----------------------------------------------------
	my (@sample_qual_R1, @sample_qual_R2, @sample_base_R1, @sample_base_R2, @length_freq_R1, @length_freq_R2);
	my ($r1_R1, $r2_R1, $r3_R1, $r4_R1, $r1_R2, $r2_R2, $r3_R2, $r4_R2, $total_gc_R1, $total_gc_R2, $i);
	my (@read_R1, @read_R2, @quality_R1, @quality_R2, $total_read_qual_R1, $total_read_qual_R2);
	
	for(my $i=0;$i<$readlength;$i++) {
		$sample_qual_R1[0][$i] = 0; $sample_qual_R1[1][$i] = 0; $sample_qual_R1[2][$i] = 0; $sample_qual_R1[3][$i] = 0;
		$sample_qual_R2[0][$i] = 0; $sample_qual_R2[1][$i] = 0; $sample_qual_R2[2][$i] = 0; $sample_qual_R2[3][$i] = 0;
		$sample_base_R1[0][$i] = 0; $sample_base_R1[1][$i] = 0; $sample_base_R1[2][$i] = 0; $sample_base_R1[3][$i] = 0; $sample_base_R1[4][$i] = 0; 
		$sample_base_R2[0][$i] = 0; $sample_base_R2[1][$i] = 0; $sample_base_R2[2][$i] = 0; $sample_base_R2[3][$i] = 0; $sample_base_R2[4][$i] = 0; 	
		$length_freq_R1[$i] = 0; $length_freq_R2[$i] = 0;
	}

	my $total_data_R1 = 0;
	my $total_data_R2 = 0;
	my $total_qual_R1 = 0;
	my $total_qual_R2 = 0;
	my @gc_fraction_R1 = (0,0,0,0,0,0,0,0,0,0);
	my @gc_fraction_R2 = (0,0,0,0,0,0,0,0,0,0);
	my $total_reads = 0;
	my @average_read_qual_R1 = (0,0,0,0);
	my @average_read_qual_R2 = (0,0,0,0);
	
	#-------------------------------------------------------------------------------------------------------
	open(INFP1,"gunzip -c $r1_file |") || die "Could not open input file $r1_file\n";
	open(INFP2,"gunzip -c $r2_file |") || die "Could not open input file $r2_file\n";
	
	while($r1_R1 = <INFP1>) {
		$r2_R1 = <INFP1>; $r3_R1 = <INFP1>; $r4_R1 = <INFP1>;
		$r1_R2 = <INFP2>; $r2_R2 = <INFP2>; $r3_R2 = <INFP2>; $r4_R2 = <INFP2>;
		chomp($r2_R1); chomp($r4_R1);
		chomp($r2_R2); chomp($r4_R2);
		
		$total_reads++;
		
		@read_R1 = split("",$r2_R1);	#Read1
		@read_R2 = split("",$r2_R2);	#Read2
		@quality_R1 = split("",$r4_R1);	#Quality1
		@quality_R2 = split("",$r4_R2);	#Quality2
		
		$total_data_R1 = $total_data_R1 + $#read_R1 + 1;
		$total_data_R2 = $total_data_R2 + $#read_R2 + 1;
		$total_read_qual_R1 = 0; $total_read_qual_R2 = 0;
		$total_gc_R1 = 0; $total_gc_R2 = 0;
		
		$length_freq_R1[$#read_R1]++;
		$length_freq_R2[$#read_R2]++;
		
		for($i=0;$i<=$#read_R1;$i++) {
			$total_read_qual_R1 = $total_read_qual_R1 + (ord($quality_R1[$i])-$QUALITY_ENCODING);
			if(ord($quality_R1[$i])-$QUALITY_ENCODING < 10) {
				$sample_qual_R1[0][$i]++;
			}
			elsif(ord($quality_R1[$i])-$QUALITY_ENCODING < 20) {
				$sample_qual_R1[1][$i]++;
			}
			elsif(ord($quality_R1[$i])-$QUALITY_ENCODING < 30) {
				$sample_qual_R1[2][$i]++;
			}
			else {
				$sample_qual_R1[3][$i]++;
			}
			
			if($read_R1[$i] eq "A") {
				$sample_base_R1[0][$i]++;
			}
			elsif($read_R1[$i] eq "C") {
				$sample_base_R1[1][$i]++;
				$total_gc_R1++;
			}
			elsif($read_R1[$i] eq "G") {
				$sample_base_R1[2][$i]++;
				$total_gc_R1++;
			}
			elsif($read_R1[$i] eq "T") {
				$sample_base_R1[3][$i]++;
			}
			else {
				$sample_base_R1[4][$i]++;
			}
		}
		
		if(int($total_read_qual_R1/($#read_R1+1)) < 10) {
			$average_read_qual_R1[0]++;
		}
		elsif(int($total_read_qual_R1/($#read_R1+1)) < 20) {
			$average_read_qual_R1[1]++;
		}
		elsif(int($total_read_qual_R1/($#read_R1+1)) < 30) {
			$average_read_qual_R1[2]++;
		}
		else {
			$average_read_qual_R1[3]++;
		}
		
		for($i=0;$i<=$#read_R2;$i++) {
			$total_read_qual_R2 = $total_read_qual_R2 + (ord($quality_R2[$i])-$QUALITY_ENCODING);
			if(ord($quality_R2[$i])-$QUALITY_ENCODING < 10) {
				$sample_qual_R2[0][$i]++;
			}
			elsif(ord($quality_R2[$i])-$QUALITY_ENCODING < 20) {
				$sample_qual_R2[1][$i]++;
			}
			elsif(ord($quality_R2[$i])-$QUALITY_ENCODING < 30) {
				$sample_qual_R2[2][$i]++;
			}
			else {
				$sample_qual_R2[3][$i]++;
			}
			
			if($read_R2[$i] eq "A") {
				$sample_base_R2[0][$i]++;
			}
			elsif($read_R2[$i] eq "C") {
				$sample_base_R2[1][$i]++;
				$total_gc_R2++;
			}
			elsif($read_R2[$i] eq "G") {
				$sample_base_R2[2][$i]++;
				$total_gc_R2++;
			}
			elsif($read_R2[$i] eq "T") {
				$sample_base_R2[3][$i]++;
			}
			else {
				$sample_base_R2[4][$i]++;
			}
		}
		
		if(int($total_read_qual_R2/($#read_R2+1)) < 10) {
			$average_read_qual_R2[0]++;
		}
		elsif(int($total_read_qual_R2/($#read_R2+1)) < 20) {
			$average_read_qual_R2[1]++;
		}
		elsif(int($total_read_qual_R2/($#read_R2+1)) < 30) {
			$average_read_qual_R2[2]++;
		}
		else {
			$average_read_qual_R2[3]++;
		}
		
		$gc_fraction_R1[int($total_gc_R1/($#read_R1+1)*10)]++;
		$gc_fraction_R2[int($total_gc_R2/($#read_R2+1)*10)]++;
		
		$total_qual_R1 = $total_qual_R1 + $total_read_qual_R1;
		$total_qual_R2 = $total_qual_R2 + $total_read_qual_R2;
	}
	close(INFP1); close(INFP2);
	
	my @total_bases_per_position = ();
	for($i=0;$i<$readlength;$i++) {
		$total_bases_per_position[0][$i] = 0;
		$total_bases_per_position[1][$i] = 0;
	}
		
	for($i=0;$i<$readlength;$i++) {
		if($length_freq_R1[$i] > 0) {
			for(my $k=0; $k<=$i; $k++) {
				$total_bases_per_position[0][$k]+=$length_freq_R1[$i];
			}
		}
		
		if($length_freq_R2[$i] > 0) {
			for(my $k=0; $k<=$i; $k++) {
				$total_bases_per_position[1][$k]+=$length_freq_R2[$i];
			}
		}
	}

	#-------------------------------------------------------------------------------------------------------
	#					Summary file
	#-------------------------------------------------------------------------------------------------------
	my $summary_file = "$working_dir/$sampleid.$fno".".AT.summary.txt";
	open(OUTFP1,">$summary_file") || die "Could not open $summary_file file\n";
	
	print OUTFP1 "$total_reads\t$total_data_R1\t$total_data_R2\t$total_qual_R1\t$total_qual_R2\n";
	print OUTFP1 $average_read_qual_R1[0],"\t",$average_read_qual_R2[0],"\n";
	print OUTFP1 $average_read_qual_R1[1],"\t",$average_read_qual_R2[1],"\n";
	print OUTFP1 $average_read_qual_R1[2],"\t",$average_read_qual_R2[2],"\n";
	print OUTFP1 $average_read_qual_R1[3],"\t",$average_read_qual_R2[3],"\n";
	
	for($i=0;$i<10;$i++) {
		print OUTFP1 $i*10,"\t",$gc_fraction_R1[$i],"\t",$gc_fraction_R2[$i],"\n";
	}
	
	for($i=0;$i<$readlength;$i++) {
		print OUTFP1 $i+1;
		print OUTFP1 "\t",$sample_qual_R1[0][$i],"\t",$sample_qual_R1[1][$i],"\t",$sample_qual_R1[2][$i],"\t",$sample_qual_R1[3][$i];
		print OUTFP1 "\t",$sample_qual_R2[0][$i],"\t",$sample_qual_R2[1][$i],"\t",$sample_qual_R2[2][$i],"\t",$sample_qual_R2[3][$i];
		print OUTFP1 "\t",$sample_base_R1[0][$i],"\t",$sample_base_R1[1][$i],"\t",$sample_base_R1[2][$i],"\t",$sample_base_R1[3][$i];
		print OUTFP1 "\t",$sample_base_R2[0][$i],"\t",$sample_base_R2[1][$i],"\t",$sample_base_R2[2][$i],"\t",$sample_base_R2[3][$i];
		print OUTFP1 "\t",$total_bases_per_position[0][$i],"\t",$total_bases_per_position[1][$i];
		print OUTFP1 "\t",$length_freq_R1[$i],"\t",$length_freq_R2[$i];
		print OUTFP1 "\n";
	}
	close(OUTFP1);
}

############################################################################################################################################################################################
sub combine_summary {
	my($total_files,$working_dir,$sampleid,$readlength) = @_;
	
	my @total_reads = (0,0,0); #Read1 & Read2
	my (@base_composition_R1, @base_composition_R2, @base_qual_composition_R1, @base_qual_composition_R2);
	my (@read_gc_content,@read_qual, $i,@records,$total_records,$line,$input_file,$k,@a);
	my @total_base_comp = ();
	my @total_base_qual = ();
	my @total_bases_per_position = ();
	my @total_quality = (0,0);
	my @read_length_R1 = ();
	my @read_length_R2 = ();
	my @length_freq_R1 = ();
	my @length_freq_R2 = ();

	for($i=0; $i<$readlength; $i++) {
		$base_composition_R1[0][$i] = 0; $base_composition_R1[1][$i] = 0; $base_composition_R1[2][$i] = 0; $base_composition_R1[3][$i] = 0;
		$base_qual_composition_R1[0][$i] = 0; $base_qual_composition_R1[1][$i] = 0; $base_qual_composition_R1[2][$i] = 0; $base_qual_composition_R1[3][$i] = 0;
		$base_composition_R2[0][$i] = 0; $base_composition_R2[1][$i] = 0; $base_composition_R2[2][$i] = 0; $base_composition_R2[3][$i] = 0;
		$base_qual_composition_R2[0][$i] = 0; $base_qual_composition_R2[1][$i] = 0; $base_qual_composition_R2[2][$i] = 0; $base_qual_composition_R2[3][$i] = 0;
		if($i < 10) {
			$read_gc_content[0][$i] = 0; $read_gc_content[1][$i] = 0;
		}
		if($i < 4) {
			$read_qual[0][$i] = 0; $read_qual[1][$i] = 0;
			$total_base_qual[0][$i] = 0; $total_base_qual[1][$i] = 0;
			$total_base_comp[0][$i] = 0; $total_base_comp[1][$i] = 0;
		}
		
		$total_bases_per_position[0][$i] = 0; $total_bases_per_position[1][$i] = 0;
		$read_length_R1[$i] = 0; $read_length_R2[$i] = 0;
	}
	
	for($i=0; $i<$readlength+20; $i++) {
		for($k=0; $k<25;$k++) {
			$records[$i][$k] = 0;
		}
	}
	
	my $total_read_length_R1 = 0;
	my $total_read_length_R2 = 0;

	open(OUTFP,">$working_dir/$sampleid.AT.Read_Qual_File.txt") || die "Could not open $working_dir/$sampleid.AT.Read_Qual_File.txt output file\n";
	print OUTFP "Percentage\tQuality\n";
	for(my $fno = 1; $fno <= $total_files; $fno++) {
		$input_file = "$working_dir/$sampleid".".$fno."."AT.summary.txt";
		@records = ();
		$total_records = 0;
		open(INFP,"<$input_file") || die "Could not open $input_file input file\n";
		while($line = <INFP>) {
			chomp($line);
			@a=split(/\t/,$line);
			for($i=0; $i<=$#a;$i++) {
				$records[$total_records][$i] = $a[$i];
			}
			$total_records++;
		}
		close(INFP);
		
		system("rm $input_file");	#delete the temp statistics file
		
		#Total reads, total data R1, total data R2
		$total_reads[0]+=$records[0][0]; $total_reads[1]+=$records[0][1]; $total_reads[2]+=$records[0][2];
		
		#Total base quality
		$total_quality[0]+=$records[0][3]; $total_quality[1]+=$records[0][4];
		
		#Average read quality R1, R2 (Q10, Q20, Q30, moreEqQ30)
		for($i = 0; $i < 4; $i++) {
			$read_qual[0][$i]+=$records[$i+1][0]; $read_qual[1][$i]+=$records[$i+1][1];
		}
		
		#For Read Quality Box-Plot
		print OUTFP $records[1][0]/$records[0][0]*100,"\tLessQ10_R1\n";
		print OUTFP $records[1][1]/$records[0][0]*100,"\tLessQ10_R2\n";
		print OUTFP $records[2][0]/$records[0][0]*100,"\tLessQ20_R1\n";
		print OUTFP $records[2][1]/$records[0][0]*100,"\tLessQ20_R2\n";
		print OUTFP $records[3][0]/$records[0][0]*100,"\tLessQ30_R1\n";
		print OUTFP $records[3][1]/$records[0][0]*100,"\tLessQ30_R2\n";
		print OUTFP $records[4][0]/$records[0][0]*100,"\tMoreEqQ30_R1\n";
		print OUTFP $records[4][1]/$records[0][0]*100,"\tMoreEqQ30_R2\n";
		
		#For GC-Content
		for($i = 0; $i < 10; $i++) {
			$read_gc_content[0][$i]+=$records[$i+5][1]; $read_gc_content[1][$i]+=$records[$i+5][2];
		}
		
		#For base composition and quality composition
		for($i = 0; $i < $readlength; $i++) {
			$base_qual_composition_R1[0][$i]+=$records[$i+15][1]; $base_qual_composition_R1[1][$i]+=$records[$i+15][2];
			$base_qual_composition_R1[2][$i]+=$records[$i+15][3]; $base_qual_composition_R1[3][$i]+=$records[$i+15][4];
			$base_qual_composition_R2[0][$i]+=$records[$i+15][5]; $base_qual_composition_R2[1][$i]+=$records[$i+15][6];
			$base_qual_composition_R2[2][$i]+=$records[$i+15][7]; $base_qual_composition_R2[3][$i]+=$records[$i+15][8];
			
			$base_composition_R1[0][$i]+=$records[$i+15][9];  $base_composition_R1[1][$i]+=$records[$i+15][10];
			$base_composition_R1[2][$i]+=$records[$i+15][11]; $base_composition_R1[3][$i]+=$records[$i+15][12];
			$base_composition_R2[0][$i]+=$records[$i+15][13]; $base_composition_R2[1][$i]+=$records[$i+15][14];
			$base_composition_R2[2][$i]+=$records[$i+15][15]; $base_composition_R2[3][$i]+=$records[$i+15][16];
			
			$total_base_comp[0][0]+=$records[$i+15][9];   $total_base_comp[1][0]+=$records[$i+15][13];
			$total_base_comp[0][1]+=$records[$i+15][10];  $total_base_comp[1][1]+=$records[$i+15][14];
			$total_base_comp[0][2]+=$records[$i+15][11];  $total_base_comp[1][2]+=$records[$i+15][15];
			$total_base_comp[0][3]+=$records[$i+15][12];  $total_base_comp[1][3]+=$records[$i+15][16];
			
			$total_base_qual[0][0]+=$records[$i+15][1]; $total_base_qual[1][0]+=$records[$i+15][5];
			$total_base_qual[0][1]+=$records[$i+15][2]; $total_base_qual[1][1]+=$records[$i+15][6];
			$total_base_qual[0][2]+=$records[$i+15][3]; $total_base_qual[1][2]+=$records[$i+15][7];
			$total_base_qual[0][3]+=$records[$i+15][4]; $total_base_qual[1][3]+=$records[$i+15][8];
			
			$total_bases_per_position[0][$i]+=$records[$i+15][17]; $total_bases_per_position[1][$i]+=$records[$i+15][18];
			$read_length_R1[$i]+=(($i+1)*$records[$i+15][19]); $read_length_R2[$i]+=(($i+1)*$records[$i+15][20]);
			$total_read_length_R1+=(($i+1)*$records[$i+15][19]); $total_read_length_R2+=(($i+1)*$records[$i+15][20]);
			$length_freq_R1[$i]+=$records[$i+15][19]; $length_freq_R2[$i]+=$records[$i+15][20];	
		}
	}
	close(OUTFP);
	
	open(STATFP,">$working_dir/$sampleid.AT.FastqSummary.txt") || die "Could not open statistics output file\n";

	print STATFP "SAMPLEID\t$sampleid\n";
	print STATFP "TOTAL_PAIRED_READS\t$total_reads[0]\n";
	print STATFP "TOTAL_READS\t",2*$total_reads[0],"\n";

	my $scale ;
	my $type;
	if(($total_reads[1]+$total_reads[2]) >= 1000000000) {	#Check for total data
		$scale = 1000000000;
		$type = "Gb";
	}
	else {
		$scale = 1000000;
		$type = "Mb";
	}

	print STATFP "TOTAL_DATA\t",($total_reads[1]+$total_reads[2])/$scale,"($type)\n";
	print STATFP "AVERAGE_READ_LENGTH\t",($total_read_length_R1+$total_read_length_R2)/(2*$total_reads[0]),"\n";
	print STATFP "TOTAL_A\t",($total_base_comp[0][0]+$total_base_comp[1][0])/$scale,"($type)\t",($total_base_comp[0][0]+$total_base_comp[1][0])/($total_reads[1]+$total_reads[2])*100,"(%)\n";
	print STATFP "TOTAL_C\t",($total_base_comp[0][1]+$total_base_comp[1][1])/$scale,"($type)\t",($total_base_comp[0][1]+$total_base_comp[1][1])/($total_reads[1]+$total_reads[2])*100,"(%)\n";
	print STATFP "TOTAL_G\t",($total_base_comp[0][2]+$total_base_comp[1][2])/$scale,"($type)\t",($total_base_comp[0][2]+$total_base_comp[1][2])/($total_reads[1]+$total_reads[2])*100,"(%)\n";
        print STATFP "TOTAL_T\t",($total_base_comp[0][3]+$total_base_comp[1][3])/$scale,"($type)\t",($total_base_comp[0][3]+$total_base_comp[1][3])/($total_reads[1]+$total_reads[2])*100,"(%)\n";
	my $total_n=(($total_reads[1]+$total_reads[2]) - ($total_base_comp[0][0]+$total_base_comp[1][0] + $total_base_comp[0][1]+$total_base_comp[1][1] + $total_base_comp[0][2]+$total_base_comp[1][2] + $total_base_comp[0][3]+$total_base_comp[1][3]))/$scale;
	my $total_n_per=100 - ($total_base_comp[0][0]+$total_base_comp[1][0] + $total_base_comp[0][1]+$total_base_comp[1][1] + $total_base_comp[0][2]+$total_base_comp[1][2] + $total_base_comp[0][3]+$total_base_comp[1][3])/($total_reads[1]+$total_reads[2])*100;
	print STATFP "TOTAL_N\t",$total_n,"($type)\t",$total_n_per,"(%)\n";
	print STATFP "GC_PERCENTAGE\t",($total_base_comp[0][1]+$total_base_comp[1][1]+$total_base_comp[0][2]+$total_base_comp[1][2])/($total_reads[1]+$total_reads[2])*100,"\n";
	
	print STATFP "AVERAGE_BASE_QUALITY\t",($total_quality[0]+$total_quality[1])/($total_reads[1]+$total_reads[2]),"\n";
	print STATFP "TOTAL_DATA_MORE_THAN_Q30\t",($total_base_qual[0][3]+$total_base_qual[1][3])/($total_reads[1]+$total_reads[2])*100,"\n";
	print STATFP "TOTAL_DATA_MORE_THAN_Q20\t",($total_base_qual[0][3]+$total_base_qual[1][3]+$total_base_qual[0][2]+$total_base_qual[1][2])/($total_reads[1]+$total_reads[2])*100,"\n";
	print STATFP "TOTAL_DATA_MORE_THAN_Q10\t",($total_base_qual[0][3]+$total_base_qual[1][3]+$total_base_qual[0][2]+$total_base_qual[1][2]+$total_base_qual[0][1]+$total_base_qual[1][1])/($total_reads[1]+$total_reads[2])*100,"\n";


	if($total_reads[1] >= 1000000000 || $total_reads[2] >= 1000000000) {
		$type = "Gb";
		$scale = 1000000000;
	}
	else {
		$type = "Mb";
		$scale = 1000000;
	}


	print STATFP "READ1_AVG_BASE_QUALITY\t",$total_quality[0]/$total_reads[1],"\n";
	print STATFP "READ2_AVG_BASE_QUALITY\t",$total_quality[1]/$total_reads[2],"\n";
	print STATFP "READ1_TOTAL_DATA\t",$total_reads[1]/$scale,"($type)\n";
	print STATFP "READ2_TOTAL_DATA\t",$total_reads[2]/$scale,"($type)\n";
	print STATFP "READ1_AVG_READ_LENGTH\t",$total_read_length_R1/$total_reads[0],"\n";
	print STATFP "READ2_AVG_READ_LENGTH\t",$total_read_length_R2/$total_reads[0],"\n";
	
	print STATFP "READ1_TOTAL_A\t",$total_base_comp[0][0]/$scale,"($type)\n";
	print STATFP "READ1_TOTAL_C\t",$total_base_comp[0][1]/$scale,"($type)\n";
	print STATFP "READ1_TOTAL_G\t",$total_base_comp[0][2]/$scale,"($type)\n";
	print STATFP "READ1_TOTAL_T\t",$total_base_comp[0][3]/$scale,"($type)\n";
	my $r1_total_n = ($total_reads[1] - ($total_base_comp[0][0]+$total_base_comp[0][1]+$total_base_comp[0][2]+$total_base_comp[0][3]))/$scale;
	print STATFP "READ1_TOTAL_N\t",$r1_total_n,"($type)\n";

	print STATFP "READ2_TOTAL_A\t",$total_base_comp[1][0]/$scale,"($type)\n";
	print STATFP "READ2_TOTAL_C\t",$total_base_comp[1][1]/$scale,"($type)\n";
        print STATFP "READ2_TOTAL_G\t",$total_base_comp[1][2]/$scale,"($type)\n";
        print STATFP "READ2_TOTAL_T\t",$total_base_comp[1][3]/$scale,"($type)\n";
	my $r2_total_n = ($total_reads[2] - ($total_base_comp[1][0]+$total_base_comp[1][1]+$total_base_comp[1][2]+$total_base_comp[1][3]))/$scale;
	print STATFP "READ2_TOTAL_N\t",$r2_total_n,"($type)\n";
	print STATFP "READ1_TOTAL_DATA_MORE_THAN_Q30\t",$total_base_qual[0][3]/$total_reads[1]*100,"(%)\n";
	print STATFP "READ1_TOTAL_DATA_MORE_THAN_Q20\t",($total_base_qual[0][3]+$total_base_qual[0][2])/$total_reads[1]*100,"(%)\n";
	print STATFP "READ1_TOTAL_DATA_MORE_THAN_Q10\t",($total_base_qual[0][3]+$total_base_qual[0][2]+$total_base_qual[0][1])/$total_reads[1]*100,"(%)\n";

	print STATFP "READ2_TOTAL_DATA_MORE_THAN_Q30\t",$total_base_qual[1][3]/$total_reads[2]*100,"(%)\n";
	print STATFP "READ2_TOTAL_DATA_MORE_THAN_Q20\t",($total_base_qual[1][3]+$total_base_qual[1][2])/$total_reads[2]*100,"(%)\n";
	print STATFP "READ2_TOTAL_DATA_MORE_THAN_Q10\t",($total_base_qual[1][3]+$total_base_qual[1][2]+$total_base_qual[1][1])/$total_reads[2]*100,"(%)\n";

	open(OUTFP1,">$working_dir/$sampleid.AT.BaseComposition_R1.txt"); open(OUTFP2,">$working_dir/$sampleid.AT.BaseComposition_R2.txt");
	open(OUTFP3,">$working_dir/$sampleid.AT.BaseQualComposition_R1.txt"); open(OUTFP4,">$working_dir/$sampleid.AT.BaseQualComposition_R2.txt");
	
	print OUTFP1 "BasePosition\tPerA\tPerC\tPerG\tPerT\tPerN\n"; print OUTFP2 "BasePosition\tPerA\tPerC\tPerG\tPerT\tPerN\n";
	print OUTFP3 "BasePosition\tLessQ10\tLessQ20\tLessQ30\tMoreEqQ30\n"; print OUTFP4 "BasePosition\tLessQ10\tLessQ20\tLessQ30\tMoreEqQ30\n";
	for($i = 0; $i< $readlength; $i++) {
		if($total_bases_per_position[0][$i] > 0) {
			my $d=100-($base_composition_R1[0][$i]+$base_composition_R1[1][$i]+$base_composition_R1[2][$i]+$base_composition_R1[3][$i])/$total_bases_per_position[0][$i]*100;
			print OUTFP1 $i+1,"\t",$base_composition_R1[0][$i]/$total_bases_per_position[0][$i]*100,"\t",$base_composition_R1[1][$i]/$total_bases_per_position[0][$i]*100,"\t",$base_composition_R1[2][$i]/$total_bases_per_position[0][$i]*100,"\t",$base_composition_R1[3][$i]/$total_bases_per_position[0][$i]*100,"\t$d\n";
		}
		else {
			print OUTFP1 $i+1,"\t0\t0\t0\t0\t0\n";
		}

		if($total_bases_per_position[1][$i] > 0) {
			my $e=100-($base_composition_R2[0][$i]+$base_composition_R2[1][$i]+$base_composition_R2[2][$i]+$base_composition_R2[3][$i])/$total_bases_per_position[1][$i]*100;
			print OUTFP2 $i+1,"\t",$base_composition_R2[0][$i]/$total_bases_per_position[1][$i]*100,"\t",$base_composition_R2[1][$i]/$total_bases_per_position[1][$i]*100,"\t",$base_composition_R2[2][$i]/$total_bases_per_position[1][$i]*100,"\t",$base_composition_R2[3][$i]/$total_bases_per_position[1][$i]*100,"\t$e\n";
		}
		else {
			print OUTFP2 $i+1,"\t0\t0\t0\t0\t0\n";
		}

		if($total_bases_per_position[0][$i] > 0) {
			print OUTFP3 $i+1,"\t",$base_qual_composition_R1[0][$i]/$total_bases_per_position[0][$i]*100,"\t",$base_qual_composition_R1[1][$i]/$total_bases_per_position[0][$i]*100,"\t",$base_qual_composition_R1[2][$i]/$total_bases_per_position[0][$i]*100,"\t",$base_qual_composition_R1[3][$i]/$total_bases_per_position[0][$i]*100,"\n";
		}
		else {
			print OUTFP3 $i+1,"\t0\t0\t0\t0\n";
		}
		
		if($total_bases_per_position[1][$i] > 0) {
			print OUTFP4 $i+1,"\t",$base_qual_composition_R2[0][$i]/$total_bases_per_position[1][$i]*100,"\t",$base_qual_composition_R2[1][$i]/$total_bases_per_position[1][$i]*100,"\t",$base_qual_composition_R2[2][$i]/$total_bases_per_position[1][$i]*100,"\t",$base_qual_composition_R2[3][$i]/$total_bases_per_position[1][$i]*100,"\n";
		}
		else {
			print OUTFP4 $i+1,"\t0\t0\t0\t0\n";
		}
	}
	close(OUTFP1); close(OUTFP2); close(OUTFP3); close(OUTFP4);
	

	open(OUTFP1,">$working_dir/$sampleid.AT.ReadGCContent.txt");
	for($i=0; $i < 10; $i++) {
		print OUTFP1 $i*10,"-",($i+1)*10,"\t",$read_gc_content[0][$i]/$total_reads[0]*100,"\t",$read_gc_content[1][$i]/$total_reads[0]*100,"\n";
	}
	close(OUTFP1);

	open(OUTFP1,">$working_dir/$sampleid.AT.ReadLength.txt");
	print OUTFP1 "Readlength\tRead1\tRead2\n";
	for($i=0; $i < $readlength; $i++) {
		print OUTFP1 $i+1,"\t",$length_freq_R1[$i]/$total_reads[0]*100,"\t",$length_freq_R2[$i]/$total_reads[0]*100,"\n";
	}
	close(OUTFP1);
}
