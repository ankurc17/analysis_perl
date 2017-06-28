#!usr/perl -w
use strict;
use warnings;

print "Initiantng Sambamba Process\n";

my ($sampleid,$working_dir,$thread,$sambamba);
open (INFP, "$ARGV[0]");


while (my $l = <INFP>)
{
	chomp ($l);
	my @input = split("\t", $l);
	if ($input[0] eq "SAMPLEID") { $sampleid = $input[1]; }
	elsif ($input[0] eq "WORKING_DIR") { $working_dir = $input[1]; }
	elsif ($input[0] eq "TOTAL_THREADS") { $thread = $input[1]; }
	elsif ($input[0] eq "SAMBAMBA") { $sambamba = $input[1]; }
}

print "Converting Sam files to BAM files\n";

system ("$sambamba view -f bam -S -o $working_dir/$sampleid.ns.discordant.bam -t $thread $working_dir/$sampleid.discordant.sam");
system ("$sambamba view -f bam -S -o $working_dir/$sampleid.ns.splitters.bam -t $thread $working_dir/$sampleid.splitters.sam");
system ("$sambamba view -f bam -S -o $working_dir/$sampleid.ns.clean.bam -t $thread $working_dir/$sampleid.clean.sam");

print "Sorting Bam Files\n";

system ("$sambamba sort -m 20G --tmpdir=$working_dir -t $thread -o $working_dir/$sampleid.discordant.sorted.bam $working_dir/$sampleid.ns.discordant.bam");
system ("$sambamba sort -m 20G --tmpdir=$working_dir -t $thread -o $working_dir/$sampleid.splitters.sorted.bam $working_dir/$sampleid.ns.splitters.bam");
system ("$sambamba sort -m 20G --tmpdir=$working_dir -t $thread -o $working_dir/$sampleid.clean.sorted.bam $working_dir/$sampleid.ns.clean.bam");

print "Indexing sorted bam files\n";

system ("$sambamba index -t $thread $working_dir/$sampleid.discordant.sorted.bam");
system ("$sambamba index -t $thread $working_dir/$sampleid.splitters.sorted.bam");
system ("$sambamba index -t $thread $working_dir/$sampleid.clean.sorted.bam");

print "BAM files are sorted and indexed\n";

if (!-z "$working_dir/$sampleid.discordant.sorted.bam") { system ("rm $working_dir/$sampleid.ns.discordant.bam*"); } else { print "Discordant bam not sorted\n"; }
if (!-z "$working_dir/$sampleid.splitters.sorted.bam") { system ("rm $working_dir/$sampleid.ns.splitters.bam*"); } else { print "Splitters bam not sorted\n"; }
if (!-z "$working_dir/$sampleid.clean.sorted.bam") { system ("rm $working_dir/$sampleid.ns.clean.bam*"); } else { print "Clean bam not sorted\n"; }
