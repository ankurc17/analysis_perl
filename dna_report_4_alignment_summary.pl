#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use FindBin qw($Bin);


print "Starting SAM processing operation";
my ($tot_hardclipped_set, $total_SE_failed_alignment, $total_SE_unaligned, $total_SE_passed_alignment, $total_PE_unaligned, $total_PE_crossmapped, $total_PE_passed_alignment, $total_PE_failed_alignment, $total_SE_aligned, $total_PE_aligned) = (0,0,0,0,0,0,0,0,0,0);

my %crossmapped_hash = (); my %chrom_pass_alignment_freq = ();
my @overall_mapq = (); my @passed_mapq = (); my @failed_mapq = ();
my @overall_insz = (); my @passed_insz = (); my @failed_insz = ();
my %hardclipped_hash = ();

split_sam($ARGV[0]);
print "SAM processing completed\n";
################################################################################################################################################
sub split_sam {
	my ($config) = @_;
	my($working_dir,$sampleid,$threads,$samtools) = '';
	my @info = ();
	my $line = '';
	open (INFP, $config) or die "Can't open $config";
	print "\nInput parameter for filter alignment step from config file....\n";
	print "----------------------------------------------------------------\n";
	while ($line = <INFP>) {
		chomp ($line);
		@info = split(/\t/,$line);
		if ($info[0] eq "WORKING_DIR") {$working_dir = $info[1]; }
		elsif($info[0] eq "SAMPLEID") {$sampleid = $info[1]; }
		elsif($info[0] eq "TOTAL_THREADS") { $threads = $info[1]; }
		elsif($info[0] eq "SAMTOOLS_PATH") { $samtools = $info[1]; }
	}
	
	print "Working Dir: $working_dir\n";
	print "Sample Name: $sampleid\n";
	print "Threads: $threads\n";
	print "Samtools Path: $samtools\n";

	my $TOTAL_SUB_FASTQ = 0;
	my $total_subfastq_file = "$working_dir/$sampleid.totalsubfastq.txt";        
        open(INFP,"<$total_subfastq_file") || die "Could not open $total_subfastq_file input file\n";
        while(my $r = <INFP>) { chomp($r); $TOTAL_SUB_FASTQ = $r; } close(INFP);

	print "File with sub-file count number opened successfully and now starting alignment\n";
	$TOTAL_SUB_FASTQ = 1;
	print "$TOTAL_SUB_FASTQ\n";
        my $manager = Parallel::ForkManager->new($threads);
        for(my $fno = 1; $fno <= $TOTAL_SUB_FASTQ; $fno++) {
                $manager->start and next;
#		generate_stats("overall",$working_dir,$fno,$sampleid,$samtools); 
		generate_stats("discordant",$working_dir,$fno,$sampleid,$samtools);
		generate_stats("clean",$working_dir,$fno,$sampleid,$samtools);
		generate_stats("splitters",$working_dir,$fno,$sampleid,$samtools);
                $manager->finish;
        }
        $manager->wait_all_children;

#	prepare_stats($sampleid,$working_dir,$TOTAL_SUB_FASTQ,"overall");
	prepare_stats($sampleid,$working_dir,$TOTAL_SUB_FASTQ,"discordant");
	prepare_stats($sampleid,$working_dir,$TOTAL_SUB_FASTQ,"clean");
	prepare_stats($sampleid,$working_dir,$TOTAL_SUB_FASTQ,"splitters");
	print "Alignment stats calculation completed\n";
	do_cleanup ($sampleid, $working_dir,$TOTAL_SUB_FASTQ);

}
################################################################################################################################################
sub do_cleanup
{
	my($sampleid,$working_dir,$TOTAL_SUB_FASTQ) = @_;
	for (my $i = 1; $i <=$TOTAL_SUB_FASTQ; $i++)
	{
	        my $hardclipped_stat_file = "$working_dir/$sampleid.$i.*.hardclipped.txt";
	        my $crossmapped_stat_file = "$working_dir/$sampleid.$i.*.crossmapped.txt";
		my $overall = "$working_dir/$sampleid.$i.overall.sam";
		my $discordant = "$working_dir/$sampleid.$i.discordant.sam";
		my $clean = "$working_dir/$sampleid.$i.clean.sam";
		my $splitters = "$working_dir/$sampleid.$i.splitters.sam";
		system ("rm $hardclipped_stat_file $crossmapped_stat_file $overall $discordant $clean $splitters");
	}
}
################################################################################################################################################
sub generate_stats
{
	my($type,$working_dir,$fno,$sampleid,$samtools) = @_;
	print "Running stats generation for $type\n";
	
#	open(INFP,"<$working_dir/$sampleid.$type.sam") || die "Could not open $working_dir/$sampleid.$fno.$type.sam file";
	open(INFP,"<$working_dir/$sampleid.$fno.$type.sam") || die "Could not open $working_dir/$sampleid.$fno.$type.sam file";
	my $r;
	while($r = <INFP>) {
		if($r =~ /^\@/) { }
		else { last; }
	}

	my $line1 = $r; chomp($line1); my @l1 = split(/\t/,$line1);
	my @arr = (); my $tot = 0;
	my $id = $l1[0];
	my $total_paired_reads = 0;
	my $total_se_reads = 0;

	$arr[$tot] = $line1; $tot++;
	my $flag = 0;
	my @l;
	my @temp = ();

	for(my $i=0; $i<=1000; $i++) { $overall_insz[$i] = 0; $passed_insz[$i] = 0; $failed_insz[$i] = 0; }
	for(my $i=0; $i<=60; $i++)   { $overall_mapq[$i] = 0; $passed_mapq[$i] = 0; $failed_mapq[$i] = 0; }


	my $ERROR_FLAG = 0;

	while(1) {
		my $line = <INFP>;
	        if(!defined($line)) { $flag = 1; }
		else {	
			chomp($line);
			@l = split(/\t/,$line);
		}

		if($flag == 0 && $l[0] eq $id) { $arr[$tot] = $line; $tot++; }
		else {
			$total_paired_reads++;
			if($tot < 2) { 
				$ERROR_FLAG = 1;	#Only one end read found
				$total_se_reads++;
			}
			elsif($tot > 2) {
				$tot_hardclipped_set++;
				@temp = (); my $new_tot = 0;
				for(my $i=0; $i<$tot; $i++) {
					my @b = split(/\t/,$arr[$i]);
					if(!($b[5] =~ m/H/)) {	$temp[$new_tot] = $arr[$i]; $new_tot++; }
					else { 
						if(!exists($hardclipped_hash{$b[2]})) { $hardclipped_hash{$b[2]} = 1; }
						else { $hardclipped_hash{$b[2]}++; }
					}
				}
				write_to_files($new_tot,2,@temp);
			}
			else {
				write_to_files($tot,1,@arr);
			}
			
			if($flag == 1) { last; }
			else { @arr = (); $tot = 0; $arr[$tot] = $line; $id = $l[0]; $tot++; }
		}
	}

	close(INFP);

	my $alignment_stat_file1 = "$working_dir/$sampleid.$fno.$type.alignmentsummary.txt";
	my $all_mapping_quality_file = "$working_dir/$sampleid.$fno.$type.allmappingquality.txt";
	my $all_insertsize_dist_file = "$working_dir/$sampleid.$fno.$type.allinsertsize.distribution.txt";
	my $hardclipped_stat_file = "$working_dir/$sampleid.$fno.$type.hardclipped.txt";
	my $crossmapped_stat_file = "$working_dir/$sampleid.$fno.$type.crossmapped.txt"; 

	open(OUTFP,">$all_insertsize_dist_file");
	my ($overall_insz_mean,$passed_insz_mean,$failed_insz_mean, $overall_tot, $passed_tot, $failed_tot) = (0,0,0,0,0,0);
	for(my $i=0; $i<= 1000; $i++) {
		print OUTFP "$i\t$overall_insz[$i]\t$passed_insz[$i]\t$failed_insz[$i]\n";
		$overall_insz_mean += ($i*$overall_insz[$i]); $passed_insz_mean += ($i*$passed_insz[$i]); $failed_insz_mean += ($i*$failed_insz[$i]);
		$overall_tot+=$overall_insz[$i]; $passed_tot+=$passed_insz[$i]; $failed_tot+=$failed_insz[$i];
	}

	close(OUTFP);
	if ($overall_tot eq 0) { $overall_insz_mean = 0; }
	else { $overall_insz_mean = $overall_insz_mean/$overall_tot; }

	if ($passed_tot eq 0) { $passed_insz_mean = 0; }
	else { $passed_insz_mean=$passed_insz_mean/$passed_tot; }

	if ($failed_tot eq 0) { $failed_insz_mean = 0; }
	else { $failed_insz_mean=$failed_insz_mean/$failed_tot; }

	open(OUTFP,">$all_mapping_quality_file");
	my ($overall_mapq_mean,$passed_mapq_mean,$failed_mapq_mean, $overall_tot_mapq, $passed_tot_mapq, $failed_tot_mapq) = (0,0,0,0,0,0);
	for(my $i=0; $i<= 60; $i++) {
	        print OUTFP "$i\t$overall_mapq[$i]\t$passed_mapq[$i]\t$failed_mapq[$i]\n";
	        $overall_mapq_mean += ($i*$overall_mapq[$i]); $passed_mapq_mean += ($i*$passed_mapq[$i]); $failed_mapq_mean += ($i*$failed_mapq[$i]);
	        $overall_tot_mapq+=$overall_mapq[$i]; $passed_tot_mapq+=$passed_mapq[$i]; $failed_tot_mapq+=$failed_mapq[$i];
	}
	close(OUTFP);
	
	if ($overall_tot_mapq eq 0) { $overall_mapq_mean = 0; }
	else { $overall_mapq_mean = $overall_mapq_mean/$overall_tot_mapq; }

	if ($passed_tot_mapq eq 0) { $passed_mapq_mean = 0; }
	else { $passed_mapq_mean=$passed_mapq_mean/$passed_tot_mapq; }

	if ($failed_tot_mapq eq 0) { $failed_mapq_mean = 0; }
	else { $failed_mapq_mean=$failed_mapq_mean/$failed_tot_mapq; }

	open(OUTFP,">$hardclipped_stat_file");
	foreach my $id (keys %hardclipped_hash) {
		print OUTFP "$id\t$hardclipped_hash{$id}\n";
	}
	close(OUTFP);

	open(OUTFP,">$crossmapped_stat_file");
	foreach my $id (keys %crossmapped_hash) {
		print OUTFP "$id\t$crossmapped_hash{$id}\n";
	}
	close(OUTFP);

	my $total_reads = 2*$total_paired_reads+$total_se_reads;

	open(statFp,">$alignment_stat_file1") || die "Could not open $alignment_stat_file1 file\n";
	print statFp "TOTAL_READS\t",$total_reads,"\n";
	print statFp "TOTAL_PE_READS\t$total_paired_reads\t(",2*$total_paired_reads/$total_reads*100,"%)\n";
	print statFp "TOTAL_SE_READS\t$total_se_reads\t(",2*$total_se_reads/$total_reads*100,"%)\n";
	print statFp "TOTAL_ALIGNED\t",2*$total_PE_aligned+$total_SE_aligned,"\t(",(2*$total_PE_aligned+$total_SE_aligned)/$total_reads*100,"%)\n";
	print statFp "TOTAL_SE_ALIGNMENT\t$total_SE_aligned\t(";
	if($total_se_reads == 0) { print statFp "0%)\n"; } else { print statFp $total_SE_aligned/$total_se_reads*100,"%)\n"; }
	print statFp "TOTAL_PE_ALIGNMENT\t$total_PE_aligned\t(";
	if($total_paired_reads == 0) { print statFp "0%)\n"; } else { print statFp $total_PE_aligned/$total_paired_reads*100,"%)\n"; }
	print statFp "TOTAL_PASSED_ALIGNMENT\t",$total_PE_passed_alignment+$total_SE_passed_alignment,"\t(",($total_PE_passed_alignment+$total_SE_passed_alignment)/$total_reads*100,"%)\n";
	print statFp "TOTAL_FAILED_ALIGNMENT\t",$total_PE_failed_alignment+$total_SE_failed_alignment,"\t(",($total_PE_failed_alignment+$total_SE_failed_alignment)/$total_reads*100,"%)\n";
	print statFp "TOTAL_CROSSMAPPED\t$total_PE_crossmapped","\t(",$total_PE_crossmapped/$total_reads*100,"%)\n";
	print statFp "TOTAL_UNALIGNED\t",2*$total_PE_unaligned+$total_SE_unaligned,"\t(",(2*$total_PE_unaligned+$total_SE_unaligned)/$total_reads*100,"%)\n";
	print statFp "TOTAL_SE_UNALIGNED\t$total_SE_unaligned\n";
	print statFp "TOTAL_PE_UNALIGNED\t$total_PE_unaligned\n";
	print statFp "TOTAL_HARDCLIPPED_SET\t$tot_hardclipped_set\n";
	print statFp "ERROR_FLAG\t$ERROR_FLAG\n";
	print statFp "OVERALL_MEAN_INSERTSIZE\t$overall_insz_mean\t$overall_tot\n";
	print statFp "PASSED_MEAN_INSERTSIZE\t$passed_insz_mean\t$passed_tot\n";
	print statFp "FAILED_MEAN_INSERTSIZE\t$failed_insz_mean\t$failed_tot\n";
	print statFp "OVERALL_MEAN_MAPQ\t$overall_mapq_mean\t$overall_tot_mapq\n";
	print statFp "PASSED_MEAN_MAPQ\t$passed_mapq_mean\t$passed_tot_mapq\n";
	print statFp "FAILED_MEAN_MAPQ\t$failed_mapq_mean\t$failed_tot_mapq\n";
	
}

#####################################################################################################################################################################################
sub write_to_files {
	my($totrec,$fg,@records) = @_;

	if($totrec == 2) {
		my @read1 = split(/\t/,$records[0]); my @read2 = split(/\t/,$records[1]);
		if($read1[2] ne "*" && $read2[2] ne "*") { #both reads are aligned
			$total_PE_aligned++;
			$overall_mapq[$read1[4]]++; $overall_mapq[$read2[4]]++;

			if(abs($read1[8]) > 1000) { $overall_insz[1000]+=2; }
			else { $overall_insz[abs($read1[8])]+=2; }

			if($read1[2] eq $read2[2]) {
				if($read1[4] >= 29 && abs($read1[8]) >= 100 && abs($read1[8]) <= 1000) { 
					$total_PE_passed_alignment++; $passed_mapq[$read1[4]]++;
					if(abs($read1[8]) > 1000) { $passed_insz[1000]++; }
					else { $passed_insz[abs($read1[8])]++; }
				}
				else {
					$total_PE_failed_alignment++; $failed_mapq[$read1[4]]++; 
					if(abs($read1[8]) > 1000) { $failed_insz[1000]++; }
                                        else { $failed_insz[abs($read1[8])]++; }
				}

				if($read2[4] >= 29 && abs($read2[8]) >= 100 && abs($read2[8]) <= 1000) { 
					$total_PE_passed_alignment++; $passed_mapq[$read2[4]]++;
					if(abs($read2[8]) > 1000) { $passed_insz[1000]++; }
                                        else { $passed_insz[abs($read2[8])]++; }
					$chrom_pass_alignment_freq{$read2[2]}++;
				}
				else {
					$total_PE_passed_alignment++; $failed_mapq[$read2[4]]++;
					if(abs($read2[8]) > 1000) { $failed_insz[1000]++; }
                                        else { $failed_insz[abs($read1[8])]++; }
				}
			}
			else {
				$total_PE_crossmapped++;
				my $id = "$read1[2]:$read2[2]";
				if(!exists($crossmapped_hash{$id})) { $crossmapped_hash{$id} = 1; }
				else { $crossmapped_hash{$id}++; }
			}
		}
		else {
			if($read1[2] eq "*" && $read2[2] eq "*") { #unaligned both end
				$total_PE_unaligned++;
			}
			else {
				$total_SE_aligned++;
				if($read1[2] ne "*") {
					if($read1[4] >= 29) { $total_SE_passed_alignment++; $passed_mapq[$read1[4]]++; }
					else { $total_SE_failed_alignment++; $failed_mapq[$read1[4]]++; }
				}
				else { $total_SE_unaligned++; }

				if($read2[2] ne "*") {
					if($read2[4] >= 29) { $total_SE_passed_alignment++; $passed_mapq[$read2[4]]++; }
					else { $total_SE_failed_alignment++; $failed_mapq[$read2[4]]++; }
				}
				else { $total_SE_unaligned++; }
			}
		}

		return(0);
	}
	else {
		return(1);
	}
}
#######################################################################################################################################################################################################	
sub prepare_stats {
	my($sampleid,$output_path,$TOTAL_SUB_FASTQ,$type) = @_;
	
	my ($i,$r,@a,$input_file,$output_file);

	my @all_insertsize_count = (); my @pass_insertsize_count = (); my @fail_insertsize_count = ();
	for($i=0; $i<=1000; $i++) {	$all_insertsize_count[$i] = 0; $pass_insertsize_count[$i] = 0; $fail_insertsize_count[$i]=0; }
	for($i=1; $i<= $TOTAL_SUB_FASTQ; $i++) {
		$input_file = $output_path."/"."$sampleid.$i.$type.allinsertsize.distribution.txt";
		open(INFP,"<$input_file") || die "Could not open $input_file\n";
		while($r=<INFP>) {
			chomp($r);
			@a = split(/\t/,$r);
			$all_insertsize_count[$a[0]] += $a[1];
			$pass_insertsize_count[$a[0]] += $a[2];
			$fail_insertsize_count[$a[0]] += $a[3];
		}
		close(INFP);

		system("rm $input_file");
	}

	$output_file = $output_path."/"."$sampleid.$type.InsertSizeDistribution.txt";
	open(OUTFP,">$output_file") || die "Could not open $output_file\n";
	for($i=0; $i<=1000; $i++) { print OUTFP "$i\t$all_insertsize_count[$i]\t$pass_insertsize_count[$i]\t$fail_insertsize_count[$i]\n"; }
	close(OUTFP);

	my @mapqual_count = (); my @pass_mapqual_count = (); my @fail_mapqual_count = ();
	for($i=0; $i<=60; $i++) { $mapqual_count[$i] = 0; $pass_mapqual_count[$i] = 0; $fail_mapqual_count[$i] =0; }
	for($i=1; $i<= $TOTAL_SUB_FASTQ; $i++) {
		$input_file = $output_path."/"."$sampleid.$i.$type.allmappingquality.txt";
		open(INFP,"<$input_file") || die "Could not open $input_file\n";
		while($r=<INFP>) {
			chomp($r);
			@a = split(/\t/,$r);
			$mapqual_count[$a[0]] += $a[1];
			$pass_mapqual_count[$a[0]] += $a[2];
			$fail_mapqual_count[$a[0]] += $a[3];
		}
		close(INFP);

		system("rm $input_file");
	}

	$output_file = $output_path."/"."$sampleid.$type.MappingQualityDistribution.txt";
	open(OUTFP,">$output_file") || die "Could not open $output_file\n";
	for($i=0; $i<=60; $i++) { print OUTFP "$i\t$mapqual_count[$i]\t$pass_mapqual_count[$i]\t$fail_mapqual_count[$i]\n"; }
	close(OUTFP); 

	my @alignment_count = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	my ($total_ins,$total_ins_pass,$total_ins_fail,$total_mapq,$total_mapq_pass,$total_mapq_fail) = (0,0,0,0,0,0);
        for($i=1; $i<= $TOTAL_SUB_FASTQ; $i++) {
                $input_file = $output_path."/"."$sampleid.$i.$type.alignmentsummary.txt";
                open(INFP,"<$input_file") || die "Could not open $input_file\n";
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[0] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[1] += $a[1];
		$r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[2] += $a[1];
		$r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[3] += $a[1];
		$r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[4] += $a[1];
		$r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[5] += $a[1];
		$r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[6] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[7] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[8] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[9] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[10] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[11] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[12] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[13] += $a[1];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[14] += ($a[1]*$a[2]); $total_ins+= $a[2];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[15] += ($a[1]*$a[2]); $total_ins_pass += $a[2];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[16] += ($a[1]*$a[2]); $total_ins_fail += $a[2];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[17] += ($a[1]*$a[2]); $total_mapq+= $a[2];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[18] += ($a[1]*$a[2]); $total_mapq_pass += $a[2];
                $r=<INFP>; chomp($r); @a = split(/\t/,$r); $alignment_count[19] += ($a[1]*$a[2]); $total_mapq_fail += $a[2];
		close(INFP);

                system("rm $input_file");
	}

        $output_file = $output_path."/"."$sampleid.$type.OverallAlignmentSummary.txt";
        open(OUTFP,">$output_file") || die "Could not open $output_file\n";
	print OUTFP "TOTAL_READS\t$alignment_count[0]\n";
        if($alignment_count[0] > 0) { 
		print OUTFP "TOTAL_PE_READS\t$alignment_count[1]\t(",$alignment_count[1]/$alignment_count[0]*2*100,"\%)\n"; 
		print OUTFP "TOTAL_SE_READS\t$alignment_count[2]\t(",$alignment_count[2]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_ALIGNED\t$alignment_count[3]\t(",$alignment_count[3]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_SE_ALIGNMENT\t$alignment_count[4]\t(",$alignment_count[4]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_PE_ALIGNMENT\t$alignment_count[5]\t(",$alignment_count[5]/$alignment_count[0]*2*100,"\%)\n";
		print OUTFP "TOTAL_PASSED_ALIGNMENT\t$alignment_count[6]\t(",$alignment_count[6]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_FAILED_ALIGNMENT\t$alignment_count[7]\t(",$alignment_count[7]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_CROSSMAPPED\t$alignment_count[8]\t(",$alignment_count[8]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_UNALIGNED\t$alignment_count[9]\t(",$alignment_count[9]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_SE_UNALIGNED\t$alignment_count[10]\t(",$alignment_count[10]/$alignment_count[0]*100,"\%)\n";
		print OUTFP "TOTAL_PE_UNALIGNED\t$alignment_count[11]\t(",$alignment_count[11]/$alignment_count[0]*2*100,"\%)\n";
		print OUTFP "TOTAL_HARDCLIPPED_SET\t$alignment_count[12]\t(",$alignment_count[12]/$alignment_count[0]*100,"\%)\n";
	}
	else {
		print OUTFP "TOTAL_PE_READS\t0\t0\n";
                print OUTFP "TOTAL_SE_READS\t0\t0\n";
                print OUTFP "TOTAL_ALIGNED\t0\t0\n";
                print OUTFP "TOTAL_SE_ALIGNMENT\t0\t0\n";
                print OUTFP "TOTAL_PE_ALIGNMENT\t0\t0\n";
                print OUTFP "TOTAL_PASSED_ALIGNMENT\t0\t0\n";
                print OUTFP "TOTAL_FAILED_ALIGNMENT\t0\t0\n";
                print OUTFP "TOTAL_CROSSMAPPED\t0\t0\n";
                print OUTFP "TOTAL_UNALIGNED\t0\t0\n";
                print OUTFP "TOTAL_SE_UNALIGNED\t0\t0\n";
                print OUTFP "TOTAL_PE_UNALIGNED\t0\t0\n";
		print OUTFP "TOTAL_HARDCLIPPED_SET\t0\t0\n";
	}
		if($total_ins > 0) { print OUTFP "OVERALL_MEAN_INSERTSIZE\t",$alignment_count[14]/$total_ins,"\t$total_ins\n"; }
	else { print OUTFP "OVERALL_MEAN_INSERTSIZE\t0\t0\n"; }

	if($total_ins_pass > 0) { print OUTFP "PASSED_MEAN_INSERTSIZE\t",$alignment_count[15]/$total_ins_pass,"\t$total_ins_pass\n"; }
        else { print OUTFP "PASSED_MEAN_INSERTSIZE\t0\t0\n"; }

	if($total_ins_fail > 0) { print OUTFP "FAILED_MEAN_INSERTSIZE\t",$alignment_count[16]/$total_ins_fail,"\t$total_ins_fail\n"; }
        else { print OUTFP "FAILED_MEAN_INSERTSIZE\t0\t0\n"; }

	if($total_mapq > 0) { print OUTFP "OVERALL_MEAN_MAPQ\t",$alignment_count[17]/$total_mapq,"\t$total_mapq\n"; }
        else { print OUTFP "OVERALL_MEAN_INSERTSIZE\t0\t0\n"; }

        if($total_mapq_pass > 0) { print OUTFP "PASSED_MEAN_MAPQ\t",$alignment_count[18]/$total_mapq_pass,"\t$total_mapq_pass\n"; }
        else { print OUTFP "PASSED_MEAN_INSERTSIZE\t0\t0\n"; }

        if($total_mapq_fail > 0) { print OUTFP "FAILED_MEAN_MAPQ\t",$alignment_count[19]/$total_mapq_fail,"\t$total_mapq_fail\n"; }
        else { print OUTFP "FAILED_MEAN_INSERTSIZE\t0\t0\n"; }

	print OUTFP "ERROR_FLAG\t$alignment_count[13]\n";

	close(OUTFP);
}
################################################################################################################################################
