#!usr/perl -w
use warnings;
use Statistics::R;

my $working_dir = $ARGV[1];
my $sampleid = $ARGV[0];

print "Generating Graphs for Raw Fastq Files\n";
make_read_quality_boxplot("$working_dir/$sampleid.Read_Qual_File.txt","$working_dir/ReadQualBoxplot.png");
make_base_composition_plot("$working_dir/$sampleid.BaseComposition_R1.txt","$working_dir/R1BaseComposition.png","Base composition distribution for R1");
make_base_composition_plot("$working_dir/$sampleid.BaseComposition_R2.txt","$working_dir/R2BaseComposition.png","Base composition distribution for R2");
make_base_quality_composition_plot("$working_dir/$sampleid.BaseQualComposition_R1.txt","$working_dir/R1BaseQualComposition.png","Base quality distribution for Read1");
make_base_quality_composition_plot("$working_dir/$sampleid.BaseQualComposition_R2.txt","$working_dir/R2BaseQualComposition.png","Base quality distribution for Read2");
make_gc_plot("$working_dir/$sampleid.ReadGCContent.txt","$working_dir/ReadGCContent.png");
make_length_plot("$working_dir/$sampleid.ReadLength.txt","$working_dir/ReadLength.png");
make_N_composition_plot ("$working_dir/$sampleid.TotalN.txt","$working_dir/TotalN.png","Read-Wise N Distribution");

print "Generating Graphs for Adap Trim Fastq Files\n";
make_base_composition_plot("$working_dir/$sampleid.AT.BaseComposition_R1.txt","$working_dir/ATBaseCompR1.png","Base composition distribution for R1 after AdapTrim");
make_base_composition_plot("$working_dir/$sampleid.AT.BaseComposition_R2.txt","$working_dir/ATBaseCompR2.png","Base composition distribution for R2 after AdapTrim");
make_base_quality_composition_plot("$working_dir/$sampleid.AT.BaseQualComposition_R1.txt","$working_dir/ATBaseQualCompR1.png","Base quality distribution for Read1 after AdapTrim");
make_base_quality_composition_plot("$working_dir/$sampleid.AT.BaseQualComposition_R2.txt","$working_dir/ATBaseQualCompR2.png","Base quality distribution for Read2 after AdapTrim");
make_length_plot("$working_dir/$sampleid.AT.ReadLength.txt","$working_dir/ATReadLength.png");
make_gc_plot("$working_dir/$sampleid.AT.ReadGCContent.txt","$working_dir/ATReadGCContent.png");
make_N_composition_plot ("$working_dir/$sampleid.AT.TotalN.txt","$working_dir/ATTotalN.png","Adap Trimmed Read-Wise N Distribution");

print "Generating Alignment Graphs\n";
make_insert_graph("$working_dir/$sampleid.insert.txt","$working_dir/Insert.png");
make_mapping_graph("$working_dir/$sampleid.mapping.txt","$working_dir/Mapping.png");
#############################################################################################################################################################################################
sub make_N_composition_plot {
	my($in_file,$out_file,$heading) = @_;
	my $R = Statistics::R->new();
	$R->startR;
	$R->set('Infilename',$in_file);
	$R->set('Outfilename',$out_file);
	$R->set('heading',$heading);

	$R->send(q`basecomposition_data <- read.table(file=Infilename, header=T, sep="\t")`);
	$R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
	$R->send(q`matplot(basecomposition_data$BasePosition,cbind(basecomposition_data$R1,basecomposition_data$R2),type="l",col=c("blue","red"),xlab="Position in read (bp)",ylab="Read percentage (%)",lty=c(1,1,1,1),lwd=c(2,2,2,2),main=heading)`);
	$R->send(q`legend("topright",c("Read 1","Read 2"),cex=1.0, bty="n", fill=c("blue","red"))`);
	$R->send(q`dev.off()`);
	$R->stopR() ;
}
#############################################################################################################################################################################################
sub make_mapping_graph {
	my($in_file,$out_file) = @_;
	my $R = Statistics::R->new();
	$R->startR; # R code starts here
	$R->set('Infilename',$in_file);
	$R->set('Outfilename1',$out_file);
	$R->send(q`data <- read.table(file=Infilename, header=TRUE, sep="\t")`);
	$R->send(q`png(filename=Outfilename1,height=4000,width=5000,res=600)`);
	$R->send(q`c1<-c(data$ReadPercentage_Clean,data$ReadPercentage_Splitters,data$ReadPercentage_Discordant)`);
	$R->send(q`color<-c("blue","red","green")`);
	$R->send(q`x<-barplot(c1, col=color, ylab="", xlab="", font.lab=2, cex.lab=0.75, cex.axis=0.75, cex.names=0.75, border="black", ylim=c(-2,100), main ="Mapping Quality Distribution")`);
	$R->send(q`mtext(side=2,"Reads (%)",line=3, cex=1, font=2)`);
	$R->send(q`mtext(side=1,"Quality Score",line=2, cex=1, font=2)`);
	$R->send(q`legend("topright",c("Clean","Splitters","Discordant"),cex=1.0, bty="n", fill=c("blue","red","green"))`);	
	$R->send(q`dev.off()`);
	$R->stopR() ;
}

#############################################################################################################################################################################################
sub make_insert_graph {
	my($in_file,$out_file) = @_;

	my $R = Statistics::R->new();
	$R->startR;
	$R->set('Infilename',$in_file);
	$R->set('Outfilename',$out_file);

	$R->send(q`insertsize_data <- read.table(file=Infilename, header=T, sep="\t")`);
	$R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
	$R->send(q`par(mai=c(1,1,1,1))`);
	$R->send(q`matplot(cbind(insertsize_data$ReadPercentage_Clean,insertsize_data$ReadPercentage_Splitters,insertsize_data$ReadPercentage_Discordant), lty=1:5, lwd=1, type="l", cex.axis=0.75, cex.lab=0.75, col=(c("blue","red","green")), xlab="", ylab="", main="Insert Size Distribution", font=2)`);
	$R->send(q`mtext(side=1,"Insert Size (bp)",line=2, cex=1, font=2)`);
	$R->send(q`mtext(side=2,"Total Reads (%)",line=3, cex=1, font=2)`);
	$R->send(q`legend("topright",c("Clean","Splitters","Discordant"),cex=1.0, bty="n", fill=c("blue","red","green"))`);
	$R->send(q`dev.off()`);
	$R->stopR() ; # R code ends here
}
#############################################################################################################################################################################################
sub make_gc_plot {
        my($in_file,$out_file) = @_;

        my $R = Statistics::R->new();
        $R->startR; # R code starts here

        $R->set('Infilename',$in_file);
        $R->set('Outfilename',$out_file);

        $R->send(q`gc_data <- read.table(file=Infilename, header=T, sep="\t")`);
        $R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
        $R->send(q`R1_data<-gc_data[1:10,2]`);
        $R->send(q`R2_data<-gc_data[1:10,3]`);
        $R->send(q`gc_data_matrix<-matrix(c(R1_data,R2_data),nrow=10,ncol=2)`);
        $R->send(q`barplot(t(gc_data_matrix),beside=TRUE,col=c("Blue","yellow"),names.arg=c("<10","<20","<30","<40","<50","<60","<70","<80","<90","<=100"), ylab="Read percentage (%)", xlab="Read GC (%)", main="Read GC Distribution",font.lab=2, cex.lab=1.0, cex.axis=1.0, cex.names=1.0, border="black", ylim=c(0,60))`);
        $R->send(q`legend("topright",c("Read1","Read2"),cex=1.0, bty="n", fill=c("blue","yellow"))`);
        $R->send(q`dev.off()`);

        $R->stopR() ; # R code ends here
}

##############################################################################################################################################################################################
sub make_read_quality_boxplot {
        my($in_file,$out_file) = @_;

        my $R = Statistics::R->new();
        $R->startR; # R code starts here

        $R->set('Infilename',$in_file);
        $R->set('Outfilename',$out_file);

        $R->send(q`readquality_data <- read.table(file=Infilename, header=T, sep="\t")`);
        $R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
        $R->send(q`boxplot(Percentage~Quality,data=readquality_data,notch=TRUE,col=(c("gold","red")),outline=FALSE,xlab="Quality Score (Phred)",ylab="Read percentage (%)",main="Read quality distribution",font.lab=4,cex.lab=1.0,cex.axis=0.6,cex.main=1.0,cex.sub=0.6)`);
        $R->send(q`dev.off()`);

        $R->stopR() ; # R code ends here
}

#############################################################################################################################################################################################
sub make_base_quality_composition_plot {
        my($in_file,$out_file,$heading) = @_;

        my $R = Statistics::R->new();
        $R->startR; # R code starts here

        $R->set('Infilename',$in_file);
        $R->set('Outfilename',$out_file);
        $R->set('heading',$heading);

        $R->send(q`basequality_data <- read.table(file=Infilename, header=T, sep="\t")`);
        $R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
        $R->send(q`v<-rbind(basequality_data$LessQ10,basequality_data$LessQ20,basequality_data$LessQ30,basequality_data$MoreEqQ30)`);
        $R->send(q`barplot(v,col=c("red","orange","brown","green"),xlab="Position in read (bp)",ylab="Read percentage (%)",names.arg=basequality_data$BasePosition,main=heading,ylim=c(0,110),axes="TRUE")`);
$R->send(q`legend("topleft",c("LessQ10","LessQ20","LessQ30","MoreEqQ30"),cex=1.0, bty="n", fill=c("red","orange","brown","green"),horiz=TRUE)`);
        $R->send(q`dev.off()`);

        $R->stopR() ; # R code ends here
}

#############################################################################################################################################################################################
sub make_base_composition_plot {
        my($in_file,$out_file,$heading) = @_;

        my $R = Statistics::R->new();
        $R->startR; # R code starts here

        $R->set('Infilename',$in_file);
        $R->set('Outfilename',$out_file);
        $R->set('heading',$heading);

        $R->send(q`basecomposition_data <- read.table(file=Infilename, header=T, sep="\t")`);
        $R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
        $R->send(q`matplot(basecomposition_data$BasePosition,cbind(basecomposition_data$PerA,basecomposition_data$PerC,basecomposition_data$PerG,basecomposition_data$PerT,basecomposition_data$PerN),type="l",col=c("green","blue","black","red","orange"),xlab="Position in read (bp)",ylab="Read percentage (%)",ylim=c(0,100),lty=c(1,1,1,1),lwd=c(2,2,2,2),main=heading)`);

        $R->send(q`legend("topright",c("A","C","G","T","N"),cex=1.0, bty="n", fill=c("green","blue","black","red","orange"))`);
        $R->send(q`dev.off()`);

        $R->stopR() ; # R code ends here
}

#############################################################################################################################################################################################
sub make_length_plot {
        my($in_file,$out_file) = @_;

        my $R = Statistics::R->new();
        $R->startR; # R code starts here

        $R->set('Infilename',$in_file);
        $R->set('Outfilename',$out_file);

        $R->send(q`length_data <- read.table(file=Infilename, header=T, sep="\t")`);
        $R->send(q`png(filename=Outfilename,height=3000,width=5000,res=600)`);
        $R->send(q`matplot(length_data$Readlength,cbind(length_data$Read1,length_data$Read2),type="l",col=c("blue","red"),xlab="Read length (bp)",ylab="Percentage of all reads (%)",ylim=c(0,100),lty=c(1,1),lwd=c(3,3),main="Read length distribution")`);

        $R->send(q`legend("topleft",c("Read1","Read2"),cex=1.0, bty="n", fill=c("blue","red"))`);
        $R->send(q`dev.off()`);

        $R->stopR() ; # R code ends here
}
