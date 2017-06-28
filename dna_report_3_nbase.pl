#!/usr/bin/perl

use strict;
use warnings;

my ($Path,$Sample);

open (CFG, "$ARGV[0]") or die "Can't open config file\n";

while (my $l = <CFG>)
{
	chomp ($l);
	my @input = split ("\t", $l);
	if ($input[0] eq "WORKING_DIR") { $Path = $input[1]; }
	elsif ($input[0] eq "SAMPLEID") { $Sample = $input[1]; }
}

open (INFP, "$Path/$Sample.BaseComposition_R1.txt");
open (OUTFP, ">$Path/$Sample.TotalN_R1.txt");
while (my $line=<INFP>)
{
	chomp ($line);
	my @a = split(/\t/, $line);
	if ($a[0] eq "BasePosition") { next;}
	elsif (($a[1]==0) && ($a[2]==0) && ($a[3]==0) && ($a[4]==0)) { last; }
	else {
		my $N = 100 - ($a[1]+$a[2]+$a[3]+$a[4]);
		print OUTFP "$a[0]\t$N\n";
	}
}
close (INFP); close (OOUTFP);

open (INFP, "$Path/$Sample.BaseComposition_R2.txt");
open (OUTFP, ">$Path/$Sample.TotalN_R2.txt");
while (my $line=<INFP>)
{
	chomp ($line);
	my @a = split(/\t/, $line);
	if ($a[0] eq "BasePosition") { next;}
	elsif (($a[1]==0) && ($a[2]==0) && ($a[3]==0) && ($a[4]==0)) { last; }
	else {
		my $N = 100 - ($a[1]+$a[2]+$a[3]+$a[4]);
		print OUTFP "$a[0]\t$N\n";
	}
}
close (INFP); close (OOUTFP);

open (OUTFP, ">$Path/$Sample.TotalN.txt");
open (INFP, "$Path/$Sample.TotalN_R1.txt");
my %save4; my %save5; my %save6;

while(my $r = <INFP>) {
	chomp($r);
	my @a = split(/\t/,$r);
	next if ($r eq "BasePosition");
	$save4{$a[0]}="$a[1]";
	$save5{$a[0]}=$a[0];
}

open (INFP2,"$Path/$Sample.TotalN_R2.txt");

while (my $r = <INFP2>) {
	chomp ($r);
	my @a = split(/\t/,$r);
	next if ($r eq "BasePosition");
	if (defined ($save5{$a[0]}))
	{
		$save6{$a[0]}=$a[1];
	}
}
print OUTFP "BasePosition\tR1\tR2\n";
foreach my $key(sort {$a<=>$b} keys %save4)
{
	print OUTFP "$key\t$save4{$key}\t$save6{$key}\n";
}
close (INFP); close (INFP2); close (OUTFP);

open (INFP, "$Path/$Sample.AT.BaseComposition_R1.txt");
open (OUTFP, ">$Path/$Sample.AT.TotalN_R1.txt");
while (my $line=<INFP>)
{
        chomp ($line);
        my @a = split(/\t/, $line);
        if ($a[0] eq "BasePosition") { next;}
        elsif (($a[1]==0) && ($a[2]==0) && ($a[3]==0) && ($a[4]==0)) { last; }
        else {
                my $N = 100 - ($a[1]+$a[2]+$a[3]+$a[4]);
                print OUTFP "$a[0]\t$N\n";
        }
}
close (INFP); close (OUTFP);

open (INFP, "$Path/$Sample.AT.BaseComposition_R2.txt");
open (OUTFP, ">$Path/$Sample.AT.TotalN_R2.txt");
while (my $line=<INFP>)
{
        chomp ($line);
        my @a = split(/\t/, $line);
        if ($a[0] eq "BasePosition") { next;}
        elsif (($a[1]==0) && ($a[2]==0) && ($a[3]==0) && ($a[4]==0)) { last; }
        else {
                my $N = 100 - ($a[1]+$a[2]+$a[3]+$a[4]);
                print OUTFP "$a[0]\t$N\n";
        }
}
close (INFP); close (OOUTFP);

open (OUTFP, ">$Path/$Sample.AT.TotalN.txt");
open (INFP, "$Path/$Sample.AT.TotalN_R1.txt");
my %save; my %save2; my %save3;

while(my $r = <INFP>) {
        chomp($r);
        my @a = split(/\t/,$r);
        next if ($r eq "BasePosition");
        $save{$a[0]}="$a[1]";
        $save2{$a[0]}=$a[0];
}

open (INFP2,"$Path/$Sample.AT.TotalN_R2.txt");

while (my $r = <INFP2>) {
        chomp ($r);
        my @a = split(/\t/,$r);
        next if ($r eq "BasePosition");
        if (defined ($save2{$a[0]}))
        {
                $save3{$a[0]}=$a[1];
        }
}
print OUTFP "BasePosition\tR1\tR2\n";

foreach my $key(sort {$a<=>$b} keys %save)
{
        print OUTFP "$key\t$save{$key}\t$save3{$key}\n";
}
close (INFP); close (INFP2); close (OUTFP);
system ("rm $Path/$Sample.AT.TotalN_R2.txt $Path/$Sample.AT.TotalN_R1.txt $Path/$Sample.TotalN_R2.txt $Path/$Sample.TotalN_R1.txt");
