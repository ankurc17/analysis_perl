#!usr/perl -w

use strict;
use warnings;

$ENV{_JAVA_OPTIONS}="-Xmx8g";

print "Merging Sam Files\n";
merge_sam ($ARGV[0]);
print "Merging of Sam files finished\n";

####################################################################################################################

sub merge_sam {
	my ($config) = @_;
	my ($TOTALFILES, $picard, $sampleid, $working_dir);
	open (CFG, "$config") or die "Can't open $config\n";
	while (my $line = <CFG>) 
	{
		chomp ($line);
		my @records = split ("\t", $line);
		if ($records[0] eq "SAMPLEID") { $sampleid = $records[1]; }
		if ($records[0] eq "WORKING_DIR") { $working_dir = $records[1]; }
		if ($records[0] eq "PICARD_PATH") { $picard = $records[1]; }
	}

	open (INFP, "$working_dir/$sampleid.totalsubfastq.txt") or die "Can't open $working_dir/$sampleid.totalsubfastq.txt\n";
	while(my $r = <INFP>) { chomp($r); $TOTALFILES = $r; }
	close (INFP);
	
	print "File with sub-file count number opened successfully\n";

	print "Sample_Name:$sampleid\n";
	print "Working_dir:$working_dir\n";
	print "Picard:$picard\n";

	my $discordant_sam = "I=$working_dir/$sampleid.1.discordant.sam";
	my $splitters_sam = "I=$working_dir/$sampleid.1.splitters.sam";
	my $clean_sam = "I=$working_dir/$sampleid.1.clean.sam";
	my $unmapped_fq = "$working_dir/$sampleid.1.unmapped.fq";

	for(my $i=2; $i<=$TOTALFILES; $i++) 
	{ 
		$discordant_sam = "$discordant_sam I=$working_dir/$sampleid.$i.discordant.sam"; 
		$splitters_sam = "$splitters_sam I=$working_dir/$sampleid.$i.splitters.sam";
		$clean_sam = "$clean_sam I=$working_dir/$sampleid.$i.clean.sam";
		$unmapped_fq = "$unmapped_fq $working_dir/$sampleid.$i.unmapped.fq";
	}
	
	system("java -jar $picard/MergeSamFiles.jar $discordant_sam O=$working_dir/$sampleid.discordant.sam VALIDATION_STRINGENCY=SILENT MERGE_SEQUENCE_DICTIONARIES=true");
	system("java -jar $picard/MergeSamFiles.jar $splitters_sam O=$working_dir/$sampleid.splitters.sam VALIDATION_STRINGENCY=SILENT MERGE_SEQUENCE_DICTIONARIES=true");
	system("java -jar $picard/MergeSamFiles.jar $clean_sam O=$working_dir/$sampleid.clean.sam VALIDATION_STRINGENCY=SILENT MERGE_SEQUENCE_DICTIONARIES=true");
	system ("cat $unmapped_fq | gzip -c > $working_dir/$sampleid.unmapped.fq.gz");
	system ("rm $unmapped_fq");
}
