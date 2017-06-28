#!usr/perl -w
use strict;
use warnings;


$ENV{_JAVA_OPTIONS}="-Xmx8g";

print "Initiantng Variant generation process\n";

my ($sampleid,$working_dir,$thread,$gatk,$pedigree,$ref,$hapmap,$omni,$phase1snp,$dbsnp,$mills,$ploidy);

open (INFP, "$ARGV[0]");

while (my $l = <INFP>)
{
	chomp ($l);
	my @input = split("\t", $l);
	if ($input[0] eq "SAMPLEID") { $sampleid = $input[1]; }
	elsif ($input[0] eq "WORKING_DIR") { $working_dir = $input[1]; }
	elsif ($input[0] eq "TOTAL_THREADS") { $thread = $input[1]; }
	elsif ($input[0] eq "GATK") { $gatk = $input[1]; }
	elsif ($input[0] eq "PEDIGREE") { $pedigree = $input[1]; }
	elsif ($input[0] eq "REF") { $ref = $input[1]; }
	elsif ($input[0] eq "HAPMAP") { $hapmap = $input[1]; }
	elsif ($input[0] eq "OMNI") { $omni = $input[1]; }
	elsif ($input[0] eq "1000G") { $phase1snp = $input[1]; }
	elsif ($input[0] eq "DBSNP") { $dbsnp = $input[1]; }
	elsif ($input[0] eq "MILLSANDGOLDSTANDARD") { $mills = $input[1]; }
}

print "Mapping Pedigree Info\n";

open (PDG, "$pedigree") or die "Can't open pedigree file\n";

while (my $line = <PDG>)
{
	my @samples = split (" ", $line);
	if ($samples[1] eq $sampleid) { $ploidy = $samples[4]; }
}

#print "Generating gVCF\n";
#system ("java -jar -Xmx6g $gatk -T HaplotypeCaller -R $ref -I $working_dir/$sampleid.clean.sorted.bam -stand_emit_conf 10 -stand_call_conf 30 -o $working_dir/$sampleid.autosomes.g.vcf -nct $thread -ploidy $ploidy");

print "SNP VariantRecalibrating\n";

system ("java -jar -Xmx4g $gatk -T VariantRecalibrator -R $ref --input $working_dir/$sampleid.autosomes.g.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni -resource:1000G,known=false,training=true,truth=true,prior=10.0 $phase1snp -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP -recalFile $working_dir/$sampleid.recalibrate_SNP.recal -tranchesFile $working_dir/$sampleid.recalibrate_SNP.tranches -rscriptFile $working_dir/$sampleid.recalibrate_SNP_plots.R -nt 5");
system ("java -jar -Xmx6g $gatk -T ApplyRecalibration -R $ref --input $working_dir/$sampleid.autosomes.g.vcf -mode SNP --ts_filter_level 99.5 -recalFile $working_dir/$sampleid.recalibrate_SNP.recal -tranchesFile $working_dir/$sampleid.recalibrate_SNP.tranches -o $working_dir/$sampleid.recalibrated_SNPs.vcf -nt 5");

print "Indel VariantRecalibrating\n";

system ("java -jar -Xmx4g $gatk -T VariantRecalibrator -R $ref --input $working_dir/$sampleid.recalibrated_SNPs.vcf -resource:millsResource,known=false,training=true,truth=true,prior=12.0 $mills -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL -recalFile $working_dir/$sampleid.recalibrate_Indels.recal -tranchesFile $working_dir/$sampleid.recalibrate_Indels.tranches -rscriptFile $working_dir/$sampleid.recalibrate_SNP_plots.R -nt 5");
system ("java -jar -Xmx6g $gatk -T ApplyRecalibration -R $ref --input $working_dir/$sampleid.recalibrated_SNPs.vcf -mode INDEL -recalFile $working_dir/$sampleid.recalibrate_Indels.recal -tranchesFile $working_dir/$sampleid.recalibrate_Indels.tranches --ts_filter_level 99.0 -o $working_dir/$sampleid.recalibrated_snps_indels.vcf -nt 5");
