#!usr/perl -w
use strict;
use warnings;

print "Initiantng Samtools Process\n";

my ($sampleid,$working_dir,$thread,$samtools_path);
open (INFP, "$ARGV[0]");


while (my $l = <INFP>)
{
	chomp ($l);
	my @input = split("\t", $l);
	if ($input[0] eq "SAMPLEID") { $sampleid = $input[1]; }
	elsif ($input[0] eq "WORKING_DIR") { $working_dir = $input[1]; }
	elsif ($input[0] eq "TOTAL_THREADS") { $thread = $input[1]; }
	elsif ($input[0] eq "SAMTOOLS_PATH") { $samtools_path = $input[1]; }
}

my $samtools = "$samtools_path/samtools";
print "Converting Sam files to BAM files\n";

system ("$samtools view -bS -o $working_dir/$sampleid.ns.discordant.bam -@ $thread $working_dir/$sampleid.discordant.sam");
system ("$samtools view -bS -o $working_dir/$sampleid.ns.splitters.bam -@ $thread $working_dir/$sampleid.splitters.sam");
system ("$samtools view -bS -o $working_dir/$sampleid.ns.clean.bam -@ $thread $working_dir/$sampleid.clean.sam");

print "Sorting Bam Files\n";

system ("$samtools sort -@ $thread -m 8G -o $working_dir/$sampleid.discordant.sorted.bam -T $working_dir/samtools-tmp -O bam $working_dir/$sampleid.ns.discordant.bam");
system ("$samtools sort -@ $thread -m 8G -o $working_dir/$sampleid.splitters.sorted.bam -T $working_dir/samtools-tmp -O bam $working_dir/$sampleid.ns.splitters.bam");
system ("$samtools sort -@ $thread -m 8G -o $working_dir/$sampleid.clean.sorted.bam -T $working_dir/samtools-tmp -O bam $working_dir/$sampleid.ns.clean.bam");

print "Indexing sorted bam files\n";

system ("$samtools index $working_dir/$sampleid.discordant.sorted.bam");
system ("$samtools index $working_dir/$sampleid.splitters.sorted.bam");
system ("$samtools index $working_dir/$sampleid.clean.sorted.bam");

print "BAM files are sorted and indexed\n";

if (!-z "$working_dir/$sampleid.discordant.sorted.bam") { system ("rm $working_dir/$sampleid.ns.discordant.bam*"); } else { print "Discordant bam not sorted\n"; }
if (!-z "$working_dir/$sampleid.splitters.sorted.bam") { system ("rm $working_dir/$sampleid.ns.splitters.bam*"); } else { print "Splitters bam not sorted\n"; }
if (!-z "$working_dir/$sampleid.clean.sorted.bam") { system ("rm $working_dir/$sampleid.ns.clean.bam*"); } else { print "Clean bam not sorted\n"; }
