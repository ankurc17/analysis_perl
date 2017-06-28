######################################################################################################################
#		This program is to generate analysis summary file
#
######################################################################################################################
#use strict;
use warnings;
use Parallel::ForkManager;
use Statistics::R;
use File::stat;
use Time::localtime;
use FindBin qw($Bin);

my $CONFIG_FILE = $ARGV[0];

my ($lid, $sid, $genderdetails, $panel, $disease, $Pid, $Index, $ap) = ("","","","","","","",0,0,0,0,0,0,0,"","","");
my @info = ();
my $line = "";

open(INFP,"<$CONFIG_FILE") || die "Could not open config file\n";
print "\nInput parameter for the program read from config file....\n";
print "---------------------------------------------------------\n";
while($line = <INFP>) {
	chomp($line);
   	@info = split(/\t/,$line);
   	print "$line\n";
	if($info[0] eq "SAMPLEID")             { $sid = $info[1]; }
	elsif($info[0] eq "WORKING_DIR")	{$ap = $info[1]; }
}

my $fs = "$ap/$sid.FastqSummary.txt";
my $Trim = "$ap/$sid.AT.FastqSummary.txt";

my ($tpr, $tr, $td, $al, $abq, $tq30, $tq20, $tq10, $abqr1, $abqr2, $tdr1, $tdr2, $readtype, $arlr1, $arlr2,  $q30r1, $q20r1, $q10r1, $q30r2, $q20r2, $q10r2);

open (INFP, "<$fs") or die "Couldn't open FastQ summary\n";
while (my $in = <INFP>){
	chomp ($in);
	my @records = split(/\t/, $in);
	if ($records[0]=~/^TOTAL_PAIRED_READS/) {$tpr = $records[1];}
	if ($records[0]=~/^TOTAL_READS/) {$tr = $records[1];}
	if ($records[0]=~/^TOTAL_DATA$/) { if($records[1]=~/(\d+)\.+(\d+)\(+(\D+)\)+/) { $td = "$1.$2"; $td = (sprintf "%.4f", $td); $readtype = $3;} }
	if ($records[0]=~/^AVERAGE_READ_LENGTH/) {$al = $records[1];}
	if ($records[0]=~/^AVERAGE_BASE_QUALITY/) {$abq = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^TOTAL_DATA_MORE_THAN_Q30/) {$tq30 = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^TOTAL_DATA_MORE_THAN_Q20/) {$tq20 = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^TOTAL_DATA_MORE_THAN_Q10/) {$tq10 = sprintf "%.4f", $records[1];}
       	if ($records[0]=~/^READ1_AVG_BASE_QUALITY/) {$abqr1 = sprintf "%.4f", $records[1];}
       	if ($records[0]=~/^READ2_AVG_BASE_QUALITY/) {$abqr2 = sprintf "%.4f", $records[1];}
       	if ($records[0]=~/^READ1_TOTAL_DATA$/) {if($records[1]=~/(\d+)\.+(\d+)\(+(\D+)\)+/) { $tdr1= "$1.$2"; $tdr1 = (sprintf "%.4f", $tdr1); $readtype = $3;}}
	if ($records[0]=~/^READ2_TOTAL_DATA$/) {if($records[1]=~/(\d+)\.+(\d+)\(+(\D+)\)+/) { $tdr2= "$1.$2"; $tdr2 = (sprintf "%.4f", $tdr2); $readtype = $3;}}
        if ($records[0]=~/^READ1_AVG_READ_LENGTH/) {$arlr1 = $records[1];}
        if ($records[0]=~/^READ2_AVG_READ_LENGTH/) {$arlr2 = $records[1];}
        if ($records[0]=~/^READ1_TOTAL_DATA_MORE_THAN_Q30/) {@t1 = split(/\(/, $records[1]); $q30r1 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ1_TOTAL_DATA_MORE_THAN_Q20/) {@t1 = split(/\(/, $records[1]); $q20r1 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ1_TOTAL_DATA_MORE_THAN_Q10/) {@t1 = split(/\(/, $records[1]); $q10r1 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ2_TOTAL_DATA_MORE_THAN_Q30/) {@t1 = split(/\(/, $records[1]); $q30r2 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ2_TOTAL_DATA_MORE_THAN_Q20/) {@t1 = split(/\(/, $records[1]); $q20r2 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ2_TOTAL_DATA_MORE_THAN_Q10/) {@t1 = split(/\(/, $records[1]); $q10r2 = sprintf "%.4f", $t1[0];}
}
close (INFP);

my ($tapr,$ttr, $ttd, $tarl, $tabq, $ttq30, $ttq20, $ttq10, $tabqr1, $tabqr2, $ttdr1, $ttdr2, $adaptype, $tarlr1, $tarlr2, $tq30r1, $tq20r1, $tq10r1, $tq30r2, $tq20r2, $tq10r2);
open (INFP, "<$Trim") or die "Can't open $Trim\n";
while (my $in = <INFP>) {
	chomp ($in);
	my @records = split(/\t/, $in);
	if ($records[0]=~/^TOTAL_PAIRED_READS/) { $tapr  = $records[1];}
	if ($records[0]=~/^TOTAL_READS/) { $ttr = $records[1];}
	if ($records[0]=~/^TOTAL_DATA$/) { if($records[1]=~/(\d+)\.+(\d+)\(+(\D+)\)+/) { $ttd = "$1.$2"; $ttd = (sprintf "%.4f", $td); $adaptype = $3;} }
	if ($records[0]=~/^AVERAGE_READ_LENGTH/) {$tarl = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^AVERAGE_BASE_QUALITY/) {$tabq = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^TOTAL_DATA_MORE_THAN_Q30/) {$ttq30 = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^TOTAL_DATA_MORE_THAN_Q20/) {$ttq20 = sprintf "%.4f", $records[1];}
	if ($records[0]=~/^TOTAL_DATA_MORE_THAN_Q10/) {$ttq10 = sprintf "%.4f", $records[1];}
        if ($records[0]=~/^READ1_AVG_BASE_QUALITY/) {$tabqr1 = sprintf "%.4f", $records[1];}
        if ($records[0]=~/^READ2_AVG_BASE_QUALITY/) {$tabqr2 = sprintf "%.4f", $records[1];}
        if ($records[0]=~/^READ1_TOTAL_DATA$/) {if($records[1]=~/(\d+)\.+(\d+)\(+(\D+)\)+/) { $ttdr1= "$1.$2"; $tdr1 = (sprintf "%.4f", $tdr1); $adaptype = $3;}}
        if ($records[0]=~/^READ2_TOTAL_DATA$/) {if($records[1]=~/(\d+)\.+(\d+)\(+(\D+)\)+/) { $ttdr2= "$1.$2"; $tdr2 = (sprintf "%.4f", $tdr2); $adaptype = $3;}}
        if ($records[0]=~/^READ1_AVG_READ_LENGTH/) {$tarlr1 = $records[1];}
        if ($records[0]=~/^READ2_AVG_READ_LENGTH/) {$tarlr2 = $records[1];}
        if ($records[0]=~/^READ1_TOTAL_DATA_MORE_THAN_Q30/) {@t1 = split(/\(/, $records[1]); $tq30r1 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ1_TOTAL_DATA_MORE_THAN_Q20/) {@t1 = split(/\(/, $records[1]); $tq20r1 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ1_TOTAL_DATA_MORE_THAN_Q10/) {@t1 = split(/\(/, $records[1]); $tq10r1 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ2_TOTAL_DATA_MORE_THAN_Q30/) {@t1 = split(/\(/, $records[1]); $tq30r2 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ2_TOTAL_DATA_MORE_THAN_Q20/) {@t1 = split(/\(/, $records[1]); $tq20r2 = sprintf "%.4f", $t1[0];}
        if ($records[0]=~/^READ2_TOTAL_DATA_MORE_THAN_Q10/) {@t1 = split(/\(/, $records[1]); $tq10r2 = sprintf "%.4f", $t1[0];}
}
close (INFP);

my $clean_reads = `grep "^TOTAL_READS" $ap/$sid.clean.OverallAlignmentSummary.txt | cut -f2`;
my $discordant_reads = `grep "^TOTAL_READS" $ap/$sid.discordant.OverallAlignmentSummary.txt | cut -f2`;
my $splitters_read = `grep "^TOTAL_READS" $ap/$sid.splitters.OverallAlignmentSummary.txt | cut -f2`;

open (CFG, "$ARGV[0]");
my %hash = ();

while (my $l = <CFG>)
{
	chomp ($l);
	my @records = split ("\t", $l);
	$hash{$records[0]}=$records[2];
}

close CFG;

open (FH,">$ap/$sid.latext.tex");
print FH <<"END";
\\documentclass[a4paper,12pt]{article}
\\usepackage[a4paper,bindingoffset=0.2in,left=1in,right=1in,top=1in,bottom=1in,footskip=.15in]{geometry}
\\usepackage{titlesec}
\\titleformat*{\\section}{\\bfseries}
\\titleformat*{\\subsection}{\\bfseries}
\\usepackage{graphicx}
\\usepackage[utf8]{inputenc}
\\usepackage[table]{xcolor}
\\setlength{\\tabcolsep}{18pt}
\\renewcommand{\\arraystretch}{1.5}
\\usepackage{lipsum}
\\usepackage{float}
\\usepackage{marvosym}
\\usepackage{multirow}
\\usepackage{pbox}
\\usepackage{caption}
\\captionsetup{font=scriptsize,labelfont=scriptsize}
\\usepackage{fancyhdr}
\\pagestyle{fancy}
\\fancyhf{}
\\fancyhead[L]{\\includegraphics[height=.4in]{/home/ankur/logo.png}}
\\usepackage{hyperref}
\\hypersetup{
	colorlinks=true,
	linkcolor=blue,
	filecolor=magenta,
	urlcolor=blue,
}
END

print FH <<"END";
\\begin{document}

\\begin{titlepage}
\\begin{center}	
\\huge \\bfseries Bioinformatics Analysis Report \\bfseries \\\\[6 cm]

\\small \\bfseries Sample Name $sid \\\\[0.4 cm]
{\\small \\today} \\\\ [2 cm]
\\includegraphics[width=0.25\\textwidth]{/home/ankur/logo.png}~\\\\[2 cm]

\\end{center}
\\end{titlepage}

\\setcounter{tocdepth}{1}
\\setcounter{page}{1}
\\thispagestyle{empty}

\\listoftables
\\listoffigures
\\clearpage

\\begin{table}
\\centering
\\captionsetup{justification=centering,singlelinecheck=off}
\\caption {Version}
\\begin{tabular}{|c|c|}
\\hline

Tools & Version \\\\ \\hline
Trim Galore & $hash{'TRIM_GALORE'}\\\\ \\hline
BWA & $hash{'BWA'}\\\\ \\hline
Picard & $hash{'PICARD_PATH'}\\\\ \\hline
Samtools & $hash{'SAMTOOLS_PATH'} \\\\ \\hline
Samblaster & $hash{'SAMBLASTER'}\\\\ \\hline
GATK & $hash{'GATK'} \\\\ \\hline
SamBamba & $hash{'SAMBAMBA'} \\\\ \\hline
Human Reference & $hash{'REF'} \\\\ \\hline
HapMap & $hash{'HAPMAP'} \\\\ \\hline
OMNI & $hash{'OMNI'} \\\\ \\hline
1000 Genome & $hash{'1000G'} \\\\ \\hline
dbsNp & $hash{'DBSNP'} \\\\ \\hline
Mills and Gold Standard & $hash{'MILLSANDGOLDSTANDARD'} \\\\ \\hline
\\end{tabular}
\\end{table}
\\clearpage

\\begin{table}
\\centering
\\tiny
\\setlength\\tabcolsep{4pt}
\\begin{minipage}[b]{0.47\\linewidth}
\\centering
\\caption {FastQC Summary}
\\label{Table:1}
\\begin{tabular}{|p{4cm}|c|}
\\hline
\\multicolumn{2}{|>{\\columncolor[gray]{.7}}c|}{\\bfseries Overall QC Summary} \\\\
\\hline
Total Paired Reads             & $tpr\\\\ \\hline
Total Reads                    & $tr \\\\ \\hline
Total Data ($readtype)             & $td \\\\ \\hline
Average Read Length (bp)       & $al \\\\ \\hline
Average base quality (Phred)   & $abq \\\\ \\hline
Total Data \$>\$= Q30 (\\%)          & $tq30 \\\\ \\hline
Total Data \$>\$= Q20 (\\%)          & $tq20 \\\\ \\hline
Total Data \$>\$= Q10 (\\%)          & $tq10 \\\\ \\hline
\\end{tabular}
\\end{minipage}
\\hfill
\\begin{minipage}{0.47\\linewidth}
\\centering
\\tiny
\\caption {Read-Wise QC Summary}
\\label{Table:2}
\\begin{tabular}[t]{|p{4cm}|c|c|}
\\hline
\\multicolumn{3}{|>{\\columncolor[gray]{.7}}c|}{\\bfseries Read Summary} \\\\
\\hline
                               & Read1 & Read2\\\\ \\hline
Average base quality (Phred)   & $abqr1 & $abqr2 \\\\ \\hline 
Total Data ($readtype)         & $tdr1 & $tdr2 \\\\ \\hline 
Total Data \$>\$= Q30 (\\%)          & $q30r1& $q30r2 \\\\ \\hline
Total Data \$>\$= Q20 (\\%)          & $q20r1 & $q20r2 \\\\ \\hline
Total Data \$>\$= Q10 (\\%)          & $q10r1 & $q10r2 \\\\ \\hline
\\end{tabular}
\\end{minipage}
\\end{table}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/R1BaseComposition.png}
\\caption{Read 1 Base Composition}
\\label{Figure:minipage1}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/R2BaseComposition.png}
\\caption{Read 2 Base Composition}
\\label{Figure:minipage2}
\\end{minipage}
\\end{figure}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/R1BaseQualComposition.png}
\\caption{Read 1 Base Quality}
\\label{Figure:minipage3}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/R2BaseQualComposition.png}
\\caption{Read 2 Base Quality}
\\label{Figure:minipage4}
\\end{minipage}
\\end{figure}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/TotalN.png}
\\caption{N Distribution}
\\label{Figure:minipage5}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ReadGCContent.png}
\\caption{GC Content Distribution}
\\label{Figure:minipage6}
\\end{minipage}
\\end{figure}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ReadLength.png}
\\caption{Read Length Distribution}
\\label{Figure:minipage7}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ReadQualBoxplot.png}
\\caption{GC Content Distribution}
\\label{Figure:minipage8}
\\end{minipage}
\\end{figure}

\\clearpage
\\begin{table}
\\tiny
\\centering
\\setlength\\tabcolsep{4pt}
\\begin{minipage}[b]{0.48\\linewidth}
\\centering
\\caption {Trim Stats}
\\label{Table:3}
\\begin{tabular}{|p{4cm}|c|}
\\hline
\\multicolumn{2}{|>{\\columncolor[gray]{.7}}c|}{\\bfseries Adapter Trim Summary} \\\\
\\hline
Total Paired Reads             & $tapr\\\\ \\hline
Total Reads                    & $ttr \\\\ \\hline
Total Data ($adaptype)             & $ttd \\\\ \\hline
Average Read Length (bp)       & $tarl \\\\ \\hline
Average base quality (Phred)   & $tabq \\\\ \\hline
Total Data \$>\$= Q30 (\\%)          & $ttq30 \\\\ \\hline
Total Data \$>\$= Q20 (\\%)          & $ttq20 \\\\ \\hline
Total Data \$>\$= Q10 (\\%)          & $ttq10 \\\\ \\hline
\\end{tabular}
\\end{minipage}
\\hfill
\\begin{minipage}{0.47\\linewidth}
\\centering
\\tiny
\\caption {Read-Wise QC Summary}
\\label{Table:4}
\\begin{tabular}[t]{|p{4cm}|c|c|}
\\hline
\\multicolumn{3}{|>{\\columncolor[gray]{.7}}c|}{\\bfseries Apater Trim Read Summary} \\\\
\\hline
                               & Read1 & Read2\\\\ \\hline
Average base quality (Phred)   & $tabqr1 & $tabqr2 \\\\ \\hline
Total Data ($adaptype)         & $ttdr1 & $ttdr2 \\\\ \\hline
Total Data \$>\$= Q30 (\\%)          & $tq30r1& $tq30r2 \\\\ \\hline
Total Data \$>\$= Q20 (\\%)          & $tq20r1 & $tq20r2 \\\\ \\hline
Total Data \$>\$= Q10 (\\%)          & $tq10r1 & $tq10r2 \\\\ \\hline
\\end{tabular}
\\end{minipage}
\\end{table}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ATBaseCompR1.png}
\\caption{Adap Trim Read 1 Base Composition}
\\label{Figure:minipage9}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ATBaseCompR2.png}
\\caption{Adap Trim Read 2 Base Composition}
\\label{Figure:minipage10}
\\end{minipage}
\\end{figure}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ATBaseQualCompR1.png}
\\caption{Adap Trim Read 1 Base Quality}
\\label{Figure:minipage11}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ATBaseQualCompR2.png}
\\caption{Adap Trim Read 2 Base Quality}
\\label{Figure:minipage12}
\\end{minipage}
\\end{figure}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ATTotalN.png}
\\caption{Adap Trim N Distribution}
\\label{Figure:minipage13}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ATReadGCContent.png}
\\caption{Adap Trim GC Distribution}
\\label{Figure:minipage14}
\\end{minipage}
\\end{figure}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/ReadLength.png}
\\caption{Adap Trim Read Length Distribution}
\\label{Figure:minipage15}
\\end{minipage}
\\end{figure}

\\clearpage

\\begin{table}[ht!]
\\tiny
\\caption {Alignment Summary}
\\label{Table:5}
\\begin{tabular}{|p{4.5cm}|c|}
\\hline
\\multicolumn{2}{|>{\\columncolor[gray]{.7}}c|}{\\bfseries Alignment Summary} \\\\
\\hline
Read Type & Number of Reads \\\\ \\hline
Clean Reads & $clean_reads \\\\ \\hline
Discordant Reads & $discordant_reads \\\\ \\hline
Split Reads & $splitters_read \\\\ \\hline
\\end{tabular}
\\end{table}

\\begin{figure}
\\centering
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/Insert.png}
\\caption{Insert Size Distribution}
\\label{Figure:minipage16}
\\end{minipage}
\\quad
\\begin{minipage}[b]{0.45\\linewidth}
\\includegraphics[width=\\textwidth]{$ap/Mapping.png}
\\caption{Mapping Quality Distribution}
\\label{Figure:minipage17}
\\end{minipage}
\\end{figure}

END

print FH <<"END";
\\end{document}

END
	
close(FH);
print "Converting LaTex to PDF\n";
system ("pdflatex $ap/$sid.latext.tex");
system ("pdflatex $ap/$sid.latext.tex");
system ("pdflatex $ap/$sid.latext.tex");
system ("pdflatex $ap/$sid.latext.tex");
system ("mv $ap/$sid.latext.pdf $ap/$sid.AnalysisReport.pdf");
print "Conversion Sucessfull\n";
