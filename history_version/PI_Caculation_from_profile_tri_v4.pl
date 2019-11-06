#!/usr/bin/perl
use warnings;
use strict;

#get option;
my $tableFlag = 0;
my $inFlag = 0;
my $outFlag = 0;
my $father_mutFlag = 0;
my $mother_mutFlag = 0;
my %opt;
my $father_mut;
my $mother_mut;

for (my $i=0; $i<@ARGV; $i++){
	if($ARGV[$i] eq '-table'){
		$tableFlag = 1;
		$opt{'table'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-in'){
		$inFlag = 1;
		$opt{'in'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-out'){
		$outFlag = 1;
		$opt{'out'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-father_mutation'){
		$father_mutFlag = 1;
		$opt{'father_mutation'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-mother_mutation'){
		$father_mutFlag = 1;
		$opt{'mother_mutation'} = $ARGV[++$i];
	}else{
		printCMD();
	}
}

if(@ARGV < 1){
	printCMD();
}

if(($tableFlag == 0) || ($inFlag == 0) || ($outFlag == 0)){
	printCMD();
}

if($father_mutFlag == 0){
	$father_mut = 0.002;
}else{
	$father_mut = $opt{'father_mutation'};
}
if($mother_mutFlag == 0){
	$mother_mut = 0.0005;
}else{
	$mother_mut = $opt{'mother_mutation'};
}

#open allelic frequency table and generate hash of hash;
my $af_table_name;
my %af_table;
my @panal;
open(AF, '< ', "$opt{'table'}") || die "can't open $ARGV[0] file:$!";
while (my $line = <AF>){
#	print "$line";
	chomp $line;
#	print "$line\n";
	if($line =~ /^#(\w+)/){
		$af_table_name = $1;
		next;
	}
	my @lines = split ("\t", $line);
#	print "@lines\n";
	my $locus = $lines[0];
#	print "$locus:\n";
	for my $i (1..$#lines){
		my ($bin, $frequency) = split("=", $lines[$i]);
		$af_table{$locus}{$bin}= "$frequency";
		my $output = "$bin".'='."$frequency";
#		print "$output\n";
	}	
}
close (AF);


##test:print each allele on each locus;
#print "\n\n\n";
#foreach my $locus (sort keys %af_table){
#	print "$locus: {\n";
#	for my $allele (sort keys %{ $af_table{$locus} }){
#		print "$allele=$af_table{$locus}{$allele}\n";
#	}
#	print "}\n";
#}

#get the minimal allele frequency of each locus, if allele of locus is not in the table, we use minimal allele frequency in that locus;
my $min;
foreach my $locus (sort keys %af_table){
	for my $allele (sort keys %{ $af_table{$locus} }){
		unless($min){
			$min = $af_table{$locus}{$allele};
			next;
		}
		if($af_table{$locus}{$allele} < $min){
			$min = $af_table{$locus}{$allele};
		}
	}
	$af_table{$locus}{min} = $min;
	$min = "";
}

##test:  weather the min allele frequency of each locus is in hash of hash %af_table;
#print "\n\n\n";
#foreach my $locus (sort keys %af_table){
#	print "$locus: {\n";
#	for my $allele (sort keys %{ $af_table{$locus} }){
#		print "$allele=$af_table{$locus}{$allele}\n";
#	}
#	print "}\n";
#}

#print allelic frequency table name and loci it contained;
@panal = keys %af_table;
print "the allelic frequency table is $af_table_name.\n";
print "the loci encluded in the allelic frequency table are: \n@panal\n\n";

#get locus, profiles from file, get PI of each locus and CPI;
my %pi_tri;
my %pi_mc;
my %pi_fc;
open(INPUT, '< ', "$opt{'in'}") || die "can't open $opt{'in'} file:$!";
open(OUTPUT, '> ', "$opt{'out'}") || die "can't open $opt{'out'} file:$!";
print OUTPUT "Locus\t","FatherGenotye\t","MotherGenotype\t","ChildGenotype\t","PI_MotherChild\t", "PI_FatherChild\t","PI\n";
while (my $input = <INPUT>){
	my $fa1; 
	my $fa2; 
	my $ma1;
	my $ma2;
	my $ca1; 
	my $ca2;
	my $tri_type;
	chomp $input;
	$input =~ s/[\r\n]//g;
	my @inputs = split "\t", $input;
	my $locus = $inputs[0];
	if($locus =~ /AMEL/){
		my @father_genotype = split "/", $inputs[3];
		$fa1 = $father_genotype[0];
		$fa2 = $father_genotype[1];
		my @mother_genotype = split "/", $inputs[1];
		$ma1 = $mother_genotype[0];
		$ma2 = $mother_genotype[1];
		my @child_genotype = split "/", $inputs[2];
		$ca1 = $child_genotype[0];
		$ca2 = $child_genotype[1];
	}else{
		my @father_genotype = sort {$a <=> $b} split "/", $inputs[3];
		$fa1 = $father_genotype[0];
		$fa2 = $father_genotype[1];
		my @mother_genotype = sort {$a <=> $b} split "/", $inputs[1];
		$ma1 = $mother_genotype[0];
		$ma2 = $mother_genotype[1];
		my @child_genotype = sort {$a <=> $b} split "/", $inputs[2];
		$ca1 = $child_genotype[0];
		$ca2 = $child_genotype[1];
	}
	if(exists $af_table{$locus}){
		my ($mc_pi, $fc_pi, $pi) = &get_pi_tri($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2);
		my $rounded_mc_pi = sprintf("%.4e", $mc_pi);
		$pi_mc{$locus} = $rounded_mc_pi;
		my $rounded_fc_pi = sprintf("%.4e", $fc_pi);
		$pi_fc{$locus} = $rounded_fc_pi;
		my $rounded_pi = sprintf("%.4e", $pi);
		$pi_tri{$locus} = $rounded_pi;
		print OUTPUT "$locus\t","$fa1/","$fa2\t", "$ma1/","$ma2\t","$ca1/","$ca2\t","$pi_mc{$locus}\t","$pi_fc{$locus}\t","$pi_tri{$locus}\n";
	}elsif($locus =~ /AMEL/){
		$pi_mc{$locus} = 1;
		$pi_fc{$locus} = 1;
 		$pi_tri{$locus} = 1;
		print OUTPUT "$locus\t","$fa1/","$fa2\t", "$ma1/","$ma2\t","$ca1/","$ca2\t","$pi_mc{$locus}\t","$pi_fc{$locus}\t","$pi_tri{$locus}\n";
	}else{
		print "the locus $locus is not in the $af_table_name allelic frequency table.\n";
		print OUTPUT "$locus\t","$fa1/","$fa2\t", "$ma1/","$ma2\t","$ca1/","$ca2\t","----\t", "----\t", "----\n";
	}
}

#caculate CPI;
my $cpi = 1;
for my $locus (keys %pi_tri){
	$cpi *= $pi_tri{$locus};
}
	
print OUTPUT 'CPI = ', "$cpi";
close(OUTPUT);
close(INPUT);

#define &get_pi subroutin
sub printCMD{
	print STDERR "\n\tUsage: perl PI_Caculation_from_profile.pl -table allele.frequency.table -input file -out put file [option]\n";
	print STDERR "\tThis program is used for caculation Paternal index in two couplet. For example:\n";
	print STDERR "\t\t1.without mutation\n";
	print STDERR "\t\tperl PI_Caculation_from_profile.pl -table /path/to/allele_frequency_table.txt -in case1.txt -out case1.pi.txt\n";
	print STDERR "\t\t2.with mutation\n";
	print STDERR "\t\tperl PI_Caculation_from_profile.pl -table /path/to/allele_frequency_table.txt -in case1.txt -out case1.pi.txt -mut 0.002\n";
	print STDERR "\toption:\n";
	print STDERR "\t\t-table Path to allele frequency table\n";
	print STDERR "\t\t-in input profile file (.txt)\n";
	print STDERR "\t\t-out output profile file (.pi.txt)\n";
	print STDERR "\t\t-mut mutation rate; in paternal mutation, we recommend 0.002; in maternal mutation, we recommend 0.0005\n";
	print STDERR "\n";
	exit;
}

sub get_pi_tri{
	#declare variable;
	my ($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2);
	my ($ms1_1, $ms1_2, $ms2_1, $ms2_2, $fs1_1, $fs1_2, $fs2_1, $fs2_2);
	my @ms;
	my @fs;
	my ($ms_min, $fs_min);
	my ($mProb1_1, $mProb1_2, $mProb2_1, $mProb2_2, $fProb1_1, $fProb1_2, $fProb2_1, $fProb2_2);
	my ($mc_x, $mc_y, $mc_pi, $fc_x, $fc_y, $fc_pi, $x, $y, $pi);
	#caculate mutation step of mother and father, and find the minimal setp of each
	($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2) = @_;
	$ms1_1 = abs($ma1 - $ca1);
	$ms1_2 = abs($ma1 - $ca2);
	$ms2_1 = abs($ma2 - $ca1);
	$ms2_2 = abs($ma2 - $ca2);
	@ms = sort {$a <=> $b} ($ms1_1, $ms1_2, $ms2_1,$ms2_2);
	$ms_min = $ms[0];
	

	$fs1_1 = abs($fa1 - $ca1);
	$fs1_2 = abs($fa1 - $ca2);
	$fs2_1 = abs($fa2 - $ca1);
	$fs2_2 = abs($fa2 - $ca2);	
	@fs = sort {$a <=> $b} ($fs1_1, $fs1_2, $fs2_1,$fs2_2);
	$fs_min = $fs[0];

	#caculate probability of each
	if($ms1_1 == 0){
		$mProb1_1 = 0.5;
	}else{
		if($ms1_1 == $ms_min){
			$mProb1_1 = 0.5*$mother_mut*0.5*0.1**($ms1_1 -1);
		}else{
			$mProb1_1 =0;
		}
	}
	if($ms1_2 == 0){
		$mProb1_2 = 0.5;
	}else{
		if($ms1_2 == $ms_min){
			$mProb1_2 = 0.5*$mother_mut*0.5*0.1**($ms1_2 -1);
		}else{
			$mProb1_2 =0;
		}
	}
	if($ms2_1 == 0){
		$mProb2_1 = 0.5;
	}else{
		if($ms2_1 == $ms_min){
			$mProb2_1 = 0.5*$mother_mut*0.5*0.1**($ms2_1 -1);
		}else{
			$mProb2_1 =0;
		}
	}
	if($ms2_2 == 0){
		$mProb2_2 = 0.5;
	}else{
		if($ms2_2 == $ms_min){
			$mProb2_2 = 0.5*$mother_mut*0.5*0.1**($ms2_2 -1);
		}else{
			$mProb2_2 =0;
		}
	}

	if($fs1_1 == 0){
		$fProb1_1 = 0.5;
	}else{
		if($fs1_1 == $fs_min){
			$fProb1_1 = 0.5*$father_mut*0.5*0.1**($fs1_1 -1);
		}else{
			$fProb1_1 =0;
		}
	}
	if($fs1_2 == 0){
		$fProb1_2 = 0.5;
	}else{
		if($fs1_2 == $fs_min){
			$fProb1_2 = 0.5*$father_mut*0.5*0.1**($fs1_2 -1);
		}else{
			$fProb1_2 =0;
		}
	}
	if($fs2_1 == 0){
		$fProb2_1 = 0.5;
	}else{
		if($fs2_1 == $fs_min){
			$fProb2_1 = 0.5*$father_mut*0.5*0.1**($fs2_1 -1);
		}else{
			$fProb2_1 =0;
		}
	}
	if($fs2_2 == 0){
		$fProb2_2 = 0.5;
	}else{
		if($fs2_2 == $fs_min){
			$fProb2_2 = 0.5*$father_mut*0.5*0.1**($fs2_2 -1);
		}else{
			$fProb2_2 =0;
		}
	}
	
	unless($af_table{$locus}{$ca1}){
		$af_table{$locus}{$ca1} = $af_table{$locus}{min};
	}
	unless($af_table{$locus}{$ca2}){
		$af_table{$locus}{$ca2} = $af_table{$locus}{min};
	}

	$mc_x = $mProb1_1*$af_table{$locus}{$ca2} + $mProb1_2*$af_table{$locus}{$ca1} + $mProb2_1*$af_table{$locus}{$ca2} + $mProb2_2*$af_table{$locus}{$ca1};
	$mc_y = 2*$af_table{$locus}{$ca1}*$af_table{$locus}{$ca2};
	$mc_pi = $mc_x/$mc_y;

	$fc_x = $fProb1_1*$af_table{$locus}{$ca2} + $fProb1_2*$af_table{$locus}{$ca1} + $fProb2_1*$af_table{$locus}{$ca2} + $fProb2_2*$af_table{$locus}{$ca1};
	$fc_y = 2*$af_table{$locus}{$ca1}*$af_table{$locus}{$ca2};
	$fc_pi = $fc_x/$fc_y;
	
	$x = $mProb1_1 * ($fProb1_2 + $fProb2_2) + $mProb1_2*($fProb1_1 + $fProb2_1) + $mProb2_1*($fProb1_2 + $fProb2_2) + $mProb2_2*($fProb1_1 + $fProb2_1);
	$y =  $mProb1_1*$af_table{$locus}{$ca2} + $mProb1_2*$af_table{$locus}{$ca1} +  $mProb2_1*$af_table{$locus}{$ca2} + $mProb2_2*$af_table{$locus}{$ca1};
	$pi = $x/$y;
	return ($mc_pi, $fc_pi, $pi);
} 	

