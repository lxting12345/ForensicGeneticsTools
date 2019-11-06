#!/usr/bin/perl
use warnings;
use strict;

#get option;
my $table_flag = 0;
my $in_flag = 0;
my $out_flag = 0;
my $father_mutation_flag = 0;
my $mother_mutation_flag = 0;
my $precision_flag = 0;
my $float_flag = 0;
my $scientific_flag = 0;
my %opt;
my $father_mut;
my $mother_mut;
my $mut;
my $precision;
my $format;

for (my $i=0; $i<@ARGV; $i++){
	if($ARGV[$i] eq '-table'){
		$table_flag = 1;
		$opt{'table'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-in'){
		$in_flag = 1;
		$opt{'in'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-out'){
		$out_flag = 1;
		$opt{'out'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-father_mutation'){
		$father_mutation_flag = 1;
		$opt{'father_mutation'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-mother_mutation'){
		$mother_mutation_flag = 1;
		$opt{'mother_mutation'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-precision'){
		$precision_flag = 1;
		$opt{'precision'} = $ARGV[++$i];
	}elsif($ARGV[$i] eq '-float'){
		$float_flag = 1;
	}elsif($ARGV[$i] eq '-scientifc'){
		$scientific_flag = 1;
	}else{
		printCMD();
	}
}

if(@ARGV < 1){
	printCMD();
}

if(($table_flag == 0) || ($in_flag == 0) || ($out_flag == 0)){
	printCMD();
}

if($father_mutation_flag == 0){
	$father_mut = 0.002;		#default is 0.002 in trio
}else{
	$father_mut = $opt{'father_mutation'};
}
if($mother_mutation_flag == 0){
	$mother_mut = 0.0005;		#default is 0.0005 in trio
}else{
	$mother_mut = $opt{'mother_mutation'};
}

if($precision_flag == 0){
	$precision = 4;			#default is 4 digital
}else{
	$precision = $opt{'precision'};
}

$format = '%' . '.' . "$precision";

if( ($float_flag == 0) && ($scientific_flag == 0) ){
	$format .= 'e';
}elsif( ($float_flag ==1) && ($scientific_flag == 0) ){
	$format .= 'f';
}elsif( ($float_flag == 0) && ($scientific_flag == 1) ){
	$format .= 'e';
}else{
	printCMD();
}

#open allelic frequency table and generate hash of hash;
my $af_table_name;
my %af_table;
open(AF, '< ', "$opt{'table'}") || die "can't open $ARGV[0] file:$!";
while (my $line = <AF>){
	chomp $line;
	if($line =~ /^#(\w+)/){
		$af_table_name = $1;
		next;
	}
	my @lines = split ("\t", $line);
	my $locus = $lines[0];
	$locus = uc $locus;
	for my $i (1..$#lines){
		my ($bin, $frequency) = split("=", $lines[$i]);
		$af_table{$locus}{$bin}= "$frequency";
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

#output allele frequency table to see wether min is in the table;
print "\n\n\n";
foreach my $locus (sort keys %af_table){
	print "$locus: \n";
	for my $allele (sort keys %{ $af_table{$locus} }){
		print "$allele=$af_table{$locus}{$allele}\n";
	}
	print "\n";
}

#print allelic frequency table name and loci it contained;
my @panal = keys %af_table;
print "the allelic frequency table is $af_table_name.\n";
print "the loci encluded in the allelic frequency table are: \n@panal\n\n";

#get locus, profiles from file, get PI of each locus;
my %pi_mc;
my %pi_fc;
my %pi_tri;
open(INPUT, '< ', "$opt{'in'}") || die "can't open $opt{'in'} file:$!";
open(OUTPUT, '> ', "$opt{'out'}") || die "can't open $opt{'out'} file:$!";
print OUTPUT "Locus\t","MatherGenotye\t","ChildGenotype\t","FatherGenotype\t","PI_MotherChild\t", "PI_FatherChild\t","PI_tri\n";
while (my $input = <INPUT>){
	my ($locus, $ma1, $ma2, $ca1, $ca2, $fa1, $fa2);
	chomp $input;
	$input =~ s/[\r\n]//g;
	my @inputs = split "\t", $input;
	$locus = $inputs[0];
	$locus = uc $locus;
	if($locus =~ /AMEL/){
		if(defined $inputs[1]){
			my @mother_genotype = split "/", $inputs[1];
			$ma1 = $mother_genotype[0];
			$ma2 = $mother_genotype[1];
		}else{
			$ma1 = "";
			$ma2 = "";
		}

		my @child_genotype = split "/", $inputs[2];
		$ca1 = $child_genotype[0];
		$ca2 = $child_genotype[1];

		if(defined $inputs[3]){
			my @father_genotype = split "/", $inputs[3]; 
			$fa1 = $father_genotype[0];
			$fa2 = $father_genotype[1];
		}else{
			$fa1 = "";
			$fa2 = "";
		}

	}else{
		if(defined $inputs[1]){
			my @mother_genotype = sort {$a <=> $b} split "/", $inputs[1];
			$ma1 = $mother_genotype[0];
			$ma2 = $mother_genotype[1];
		}else{
			$ma1 = "";
			$ma2 = "";
		}

		my @child_genotype = sort {$a <=> $b} split "/", $inputs[2];
		$ca1 = $child_genotype[0];
		$ca2 = $child_genotype[1];

		if(defined $inputs[3]){
			my @father_genotype = sort {$a <=> $b} split "/", $inputs[3];
			$fa1 = $father_genotype[0];
			$fa2 = $father_genotype[1];
		}else{
			$fa1 = "";
			$fa2 = "";
		}	
	}

	if(exists $af_table{$locus}){
		unless($ca1 && $ca2){
			print "ERROR: please input genotype in the file!\n";
		}
		unless(($ma1 && $ma2) && ($fa1 && $fa2)){
			print "ERROR: please input genotype in the file!\n";
		}
		if($ma1 && $ma2){
			if($mother_mutation_flag == 0){
				$mut = 0.002;		#default is 0.002 in duo
			}else{
				$mut = $opt{'mother_mutation'};
			}
			
			my $mc_pi = &get_pi_duo($locus, $ma1, $ma2, $ca1, $ca2);
			my $rounded_mc_pi = sprintf("$format", $mc_pi);
			$pi_mc{$locus} = $rounded_mc_pi;
		}
		if($fa1 && $fa2){
			if($father_mutation_flag == 0){
				$mut = 0.002;		#default is 0.002 in duo
			}else{
				$mut = $opt{'father_mutation'};
			}
			my $fc_pi = &get_pi_duo($locus, $fa1, $fa2, $ca1, $ca2);
			my $rounded_fc_pi = sprintf("$format", $fc_pi);
			$pi_fc{$locus} = $rounded_fc_pi;
		}
		if(  ($ma1 && $ma2) && ($fa1 && $fa2)  ){
			my $pi = &get_pi_tro($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2);
			my $rounded_pi = sprintf("$format", $pi);
			$pi_tri{$locus} = $rounded_pi;
		}


		if($inputs[1] && (!(defined $inputs[3]))){
			print OUTPUT "$locus\t","$ma1/","$ma2\t","$ca1/","$ca2\t", "--/","--\t", "$pi_mc{$locus}\t","----\t","----\n";
		}elsif(($inputs[1] eq "") && $inputs[3]){
			print OUTPUT "$locus\t","--/","--\t","$ca1/","$ca2\t", "$fa1/","$fa2\t", "----\t","$pi_fc{$locus}\t","----\n";
		}elsif($inputs[1] && $inputs[3]){
			print OUTPUT "$locus\t","$ma1/","$ma2\t","$ca1/","$ca2\t", "$fa1/","$fa2\t", "$pi_mc{$locus}\t","$pi_fc{$locus}\t","$pi_tri{$locus}\n";
		}else{
			print "ERROR:please input genotype!\n";
		}
	}elsif($locus =~ /AMEL/){
		$pi_mc{$locus} = 1;
		$pi_fc{$locus} = 1;
 		$pi_tri{$locus} = 1;
		if( defined $inputs[1] && ( !(defined $inputs[3] ) )){
			print OUTPUT "$locus\t","$ma1/","$ma2\t","$ca1/","$ca2\t", "--/","--\t", "$pi_mc{$locus}\t","----\t","----\n";
		}elsif(($inputs[1] eq "") && $inputs[3]){
			print OUTPUT "$locus\t","--/","--\t","$ca1/","$ca2\t", "$fa1/","$fa2\t", "----\t","$pi_fc{$locus}\t","----\n";
		}elsif($inputs[1] && $inputs[3]){
			print OUTPUT "$locus\t","$ma1/","$ma2\t","$ca1/","$ca2\t", "$fa1/","$fa2\t", "$pi_mc{$locus}\t","$pi_fc{$locus}\t","$pi_tri{$locus}\n";
		}else{
			print "ERROR:please input genotype!\n";
		}
	}else{
		print "the locus $locus is not in the $af_table_name allelic frequency table.\n";
		print OUTPUT "$locus\t","$fa1/","$fa2\t", "$ma1/","$ma2\t","$ca1/","$ca2\t","----\t", "----\t", "----\n";
	}
}



#caculate CPI and RCP;
my $cpi_tri = 1;
my $cpi_mc = 1;
my $cpi_fc = 1;
my ($RCP_mc, $RCP_fc, $RCP_tri);
for my $locus (keys %pi_mc){
	$cpi_mc *= $pi_mc{$locus};
}
$cpi_mc = sprintf("$format", $cpi_mc);
$RCP_mc = $cpi_mc / ($cpi_mc + 1);
$RCP_mc = sprintf("$format", $RCP_mc);

for my $locus (keys %pi_fc){
	$cpi_fc *= $pi_fc{$locus};
}
$cpi_fc = sprintf("$format", $cpi_fc);
$RCP_fc = $cpi_fc / ($cpi_fc + 1);
$RCP_fc = sprintf("$format", $RCP_fc);

for my $locus (keys %pi_tri){
	$cpi_tri *= $pi_tri{$locus};
}
$cpi_tri = sprintf("$format", $cpi_tri);
$RCP_tri = $cpi_tri / ($cpi_tri + 1);
$RCP_tri = sprintf("$format", $RCP_tri);

print OUTPUT 'CPI_MC = ', "$cpi_mc\t", 'RCP_MC = ', "$RCP_mc\n";
print OUTPUT 'CPI_FC = ', "$cpi_fc\t", 'RCP_FC = ', "$RCP_fc\n";
print OUTPUT 'CPI_TRI = ', "$cpi_tri\t", 'RCP_TRI = ', "$RCP_tri\n";
close(OUTPUT);
close(INPUT);


#define &get_pi subroutin
sub printCMD{
	print STDERR "\n\tUsage: perl PI_Caculation_from_profile.pl -table allele.frequency.table -input file -out put file [option]\n";
	print STDERR "\tThis program is used for caculation Paternal index in two couplet or triplet. For example:\n";
	print STDERR "\t\t1.without mutation\n";
	print STDERR "\t\tperl PI_Caculation_from_profile.pl -table /path/to/allele_frequency_table.txt -in case1.txt -out case1.pi.txt\n";
	print STDERR "\t\t2.with mutation and setting mutation value\n";
	print STDERR "\t\tperl PI_Caculation_from_profile.pl -table /path/to/allele_frequency_table.txt -in case1.txt -out case1.pi.txt -father_mutation 0.003\n";
	print STDERR "\toption:\n";
	print STDERR "\t\t-table: Path to allele frequency table\n";
	print STDERR "\t\t-in: input profile file (.txt)\n";
	print STDERR "\t\t-out: output profile file (.pi.txt)\n";
	print STDERR "\t\t-father_mutation[float]:  in paternal mutation, we recommend 0.002(default)\n";
	print STDERR "\t\t-mother_mutation[float]:  in maternal mutation, we recommend 0.0005(default)\n";
	print STDERR "\t\t-precision[int]: specify a precision\n";
	print STDERR "\t\t-float: specify the display format in decimal float, if with mutation, the value is small, we recommend use -scientific \n";
	print STDERR "\t\t-scientific: specify the display format in scientific notation\n";
	print STDERR "\n";
	exit;
}

sub get_pi_duo{
		#declare variable;
	my ($locus, $ma1, $ma2, $ca1, $ca2);
	my ($ms1_1, $ms1_2, $ms2_1, $ms2_2,);
	my $match_flag;
	my @ms;
	my $ms_min;
	my ($mProb1_1, $mProb1_2, $mProb2_1, $mProb2_2);
	my ($mc_x, $mc_y, $mc_pi);
	#judge wether match
	($locus, $ma1, $ma2, $ca1, $ca2) = @_;
	if(($ca1 == $ma1)||($ca1 == $ma2)){
		$match_flag = 1;
	}elsif(($ca2 == $ma1) || ($ca2 == $ma2)){
		$match_flag = 1;
	}else{
		$match_flag = 0;
	}
	#caculate mutation steps, and find the smallest value but not 0
	$ms1_1 = abs($ma1 - $ca1);
	$ms1_2 = abs($ma1 - $ca2);
	$ms2_1 = abs($ma2 - $ca1);
	$ms2_2 = abs($ma2 - $ca2);
	@ms = sort {$a <=> $b} ($ms1_1, $ms1_2, $ms2_1,$ms2_2);
	for my $i (@ms){
		if($i == 0){
			next;
		}else{
			$ms_min = $i;
			last;
		}
	}	
	#caculate probability of each
	if($ms1_1 == 0){
		$mProb1_1 = 0.5;
	}else{
		if($match_flag){
			$mProb1_1 =0;
		}else{
			if($ms1_1 == $ms_min){
				$mProb1_1 = 0.5*$mut*0.5*0.1**($ms1_1 -1);
			}else{
				$mProb1_1 = 0;
			}
		}
	}

	if($ms1_2 == 0){
		$mProb1_2 = 0.5;
	}else{
		if($match_flag){
			$mProb1_2 =0;
		}else{
			if($ms1_2 == $ms_min){
				$mProb1_2 = 0.5*$mut*0.5*0.1**($ms1_2 -1);
			}else{
				$mProb1_2 = 0;
			}
		}
	}


	if($ms2_1 == 0){
		$mProb2_1 = 0.5;
	}else{
		if($match_flag){
			$mProb2_1 =0;
		}else{
			if($ms2_1 == $ms_min){
				$mProb2_1 = 0.5*$mut*0.5*0.1**($ms2_1 -1);
			}else{
				$mProb2_1 = 0;
			}
		}
	}

	if($ms2_2 == 0){
		$mProb2_2 = 0.5;
	}else{
		if($match_flag){
			$mProb2_2 =0;
		}else{
			if($ms2_2 == $ms_min){
				$mProb2_2 = 0.5*$mut*0.5*0.1**($ms2_2 -1);
			}else{
				$mProb2_2 = 0;
			}
		}
	}

	unless($af_table{$locus}{$ca1}){
		$af_table{$locus}{$ca1} = $af_table{$locus}{min};
	}
	unless($af_table{$locus}{$ca2}){
		$af_table{$locus}{$ca2} = $af_table{$locus}{min};
	}

	if($ma1 && $ma2){
		$mc_x = $mProb1_1*$af_table{$locus}{$ca2} + $mProb1_2*$af_table{$locus}{$ca1} + $mProb2_1*$af_table{$locus}{$ca2} + $mProb2_2*$af_table{$locus}{$ca1};
		$mc_y = 2*$af_table{$locus}{$ca1}*$af_table{$locus}{$ca2};
		$mc_pi = $mc_x/$mc_y;
	}else{
		$mc_pi = "";
	}

	return ($mc_pi);
} 
	

sub get_pi_tro{
	#declare variable;
	my ($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2);
	my ($ms1_1, $ms1_2, $ms2_1, $ms2_2, $fs1_1, $fs1_2, $fs2_1, $fs2_2);
	my ($mother_match_flag, $father_match_flag, $match_flag);
	my @ms;
	my @fs;
	my ($ms_min, $fs_min, $mfs_min);
	my ($mProb1_1, $mProb1_2, $mProb2_1, $mProb2_2, $fProb1_1, $fProb1_2, $fProb2_1, $fProb2_2);
	my ($mc_x, $mc_y, $mc_pi, $fc_x, $fc_y, $fc_pi, $x, $y, $pi);
	#judge wether match
	($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2) = @_;
	if(($ca1 == $ma1)||($ca1 == $ma2)){
		$mother_match_flag = 1;
	}elsif(($ca2 == $ma1) || ($ca2 == $ma2)){
		$mother_match_flag = 1;
	}else{
		$mother_match_flag = 0;
	}
	if(($ca1 == $fa1)||($ca1 == $fa2)){
		$father_match_flag = 1;
	}elsif(($ca2 == $fa1) || ($ca2 == $fa2)){
		$father_match_flag = 1;
	}else{
		$father_match_flag = 0;
	}
	if( (($ca1 == $ma1)||($ca1 == $ma2))  &&  (($ca2 == $fa1) || ($ca2 == $fa2)) ){
		$match_flag = 1;
	}elsif( (($ca2 == $ma1) || ($ca2 == $ma2))  &&  (($ca1 == $fa1) || ($ca1 == $fa2))  ){
		$match_flag = 1;
	}else{
		$match_flag = 0;
	}
	#caculate mutation steps, and find the smallest value but not 0
	$ms1_1 = abs($ma1 - $ca1);
	$ms1_2 = abs($ma1 - $ca2);
	$ms2_1 = abs($ma2 - $ca1);
	$ms2_2 = abs($ma2 - $ca2);
	@ms = sort {$a <=> $b} ($ms1_1, $ms1_2, $ms2_1,$ms2_2);
	for my $i (@ms){
		if($i == 0){
			next;
		}else{
			$ms_min = $i;
			last;
		}
	}	
		
	
	$fs1_1 = abs($fa1 - $ca1);
	$fs1_2 = abs($fa1 - $ca2);
	$fs2_1 = abs($fa2 - $ca1);
	$fs2_2 = abs($fa2 - $ca2);
	@fs = sort {$a <=> $b} ($fs1_1, $fs1_2, $fs2_1,$fs2_2);
	for my $i (@fs){
		if($i == 0){
			next;
		}else{
			$fs_min = $i;
			last;
		}
	}
	my $s_min;
	if(!($ms_min)){
		$s_min = $fs_min;
	}elsif(!($fs_min)){
		$s_min = $ms_min;
	}elsif($ms_min < $fs_min){
		$s_min = $ms_min;
	}else{
		$s_min = $fs_min;
	}
		

	#caculate probability of each
	if($match_flag){
		if($ms1_1 == 0){
			$mProb1_1 = 0.5;
		}else{
			$mProb1_1 =0;
		}
	}else{
		if($ms1_1 == 0){ 
			$mProb1_1 = 0.5; 
		}elsif($ms1_1 == $s_min){
			$mProb1_1 = 0.5*$mother_mut*0.5*0.1**($ms1_1 -1);
		}else{
			$mProb1_1 =0;
		}
	}
	
	if($match_flag){
		if($ms1_2 == 0){
			$mProb1_2 = 0.5;
		}else{
			$mProb1_2 =0;
		}
	}else{
		if($ms1_2 == 0){ 
			$mProb1_2 = 0.5; 
		}elsif($ms1_2 == $s_min){
			$mProb1_2 = 0.5*$mother_mut*0.5*0.1**($ms1_2 -1);
		}else{
			$mProb1_2 =0;
		}
	}

	if($match_flag){
		if($ms2_1 == 0){
			$mProb2_1 = 0.5;
		}else{
			$mProb2_1 =0;
		}
	}else{
		if($ms2_1 == 0){ 
			$mProb2_1 = 0.5; 
		}elsif($ms2_1 == $s_min){
			$mProb2_1 = 0.5*$mother_mut*0.5*0.1**($ms2_1 -1);
		}else{
			$mProb2_1 =0;
		}
	}
	
	if($match_flag){
		if($ms2_2 == 0){
			$mProb2_2 = 0.5;
		}else{
			$mProb2_2 =0;
		}
	}else{
		if($ms2_2 == 0){ 
			$mProb2_2 = 0.5; 
		}elsif($ms2_2 == $s_min){
			$mProb2_2 = 0.5*$mother_mut*0.5*0.1**($ms2_2 -1);
		}else{
			$mProb2_2 =0;
		}
	}

	if($match_flag){
		if($fs1_1 == 0){
			$fProb1_1 = 0.5;
		}else{
			$fProb1_1 = 0;
		}
	}else{
		if($fs1_1 == 0){ 
			$fProb1_1 = 0.5;
		}elsif($fs1_1 == $s_min){
			$fProb1_1 = 0.5*$father_mut*0.5*0.1**($fs1_1 -1);
		}else{
			$fProb1_1 = 0;
		}
	}

	if($match_flag){
		if($fs1_2 == 0){
			$fProb1_2 = 0.5;
		}else{
			$fProb1_2 = 0;
		}
	}else{
		if($fs1_2 == 0){ 
			$fProb1_2 = 0.5;
		}elsif($fs1_2 == $s_min){
			$fProb1_2 = 0.5*$father_mut*0.5*0.1**($fs1_2 -1);
		}else{
			$fProb1_2 = 0;
		}
	}
	
	if($match_flag){
		if($fs2_1 == 0){
			$fProb2_1 = 0.5;
		}else{
			$fProb2_1 = 0;
		}
	}else{
		if($fs2_1 == 0){ 
			$fProb2_1 = 0.5;
		}elsif($fs2_1 == $s_min){
			$fProb2_1 = 0.5*$father_mut*0.5*0.1**($fs2_1 -1);
		}else{
			$fProb2_1 = 0;
		}
	}

	if($match_flag){
		if($fs2_2 == 0){
			$fProb2_2 = 0.5;
		}else{
			$fProb2_2 = 0;
		}
	}else{
		if($fs2_2 == 0){ 
			$fProb2_2 = 0.5;
		}elsif($fs2_2 == $s_min){
			$fProb2_2 = 0.5*$father_mut*0.5*0.1**($fs2_2 -1);
		}else{
			$fProb2_2 = 0;
		}
	}

	unless($af_table{$locus}{$ca1}){
		$af_table{$locus}{$ca1} = $af_table{$locus}{min};
	}
	unless($af_table{$locus}{$ca2}){
		$af_table{$locus}{$ca2} = $af_table{$locus}{min};
	}
	
	$x = $mProb1_1 * ($fProb1_2 + $fProb2_2) + $mProb1_2*($fProb1_1 + $fProb2_1) + $mProb2_1*($fProb1_2 + $fProb2_2) + $mProb2_2*($fProb1_1 + $fProb2_1);
	$y = 0;
	if($mother_match_flag){
		if($ms1_1 == 0){
			$y =  $mProb1_1*$af_table{$locus}{$ca2} ;
		}
		if($ms1_2 == 0){
			$y = $y + $mProb1_2*$af_table{$locus}{$ca1};
		}
		if($ms2_1 == 0){
			$y = $y + $mProb2_1*$af_table{$locus}{$ca2} ;
		}
		if($ms2_2 == 0){
			$y = $y + $mProb2_2*$af_table{$locus}{$ca1};
		}
	}elsif( (!($mother_match_flag)) && ($father_match_flag) ){
		if($fs1_1 == 0){
			$y = $fProb1_1*$af_table{$locus}{$ca2};
		}
		if($fs1_2 == 0){
			$y = $y + $fProb1_2*$af_table{$locus}{$ca1};
		}
		if($fs2_1 == 0){
			$y = $y + $fProb2_1*$af_table{$locus}{$ca2};
		}
		if($fs2_2 == 0){
			$y = $y + $fProb2_2*$af_table{$locus}{$ca1};
		}
	}else{
		if($ms1_1 == $s_min){
			$y =  $mProb1_1*$af_table{$locus}{$ca2} ;
		}
		if($ms1_2 == $s_min){
			$y = $y + $mProb1_2*$af_table{$locus}{$ca1};
		}
		if($ms2_1 == $s_min){
			$y = $y + $mProb2_1*$af_table{$locus}{$ca2} ;
		}
		if($ms2_2 == $s_min){
			$y = $y + $mProb2_2*$af_table{$locus}{$ca1};
		}
		if($fs1_1 == $s_min){
			$y = $fProb1_1*$af_table{$locus}{$ca2};
		}
		if($fs1_2 == $s_min){
			$y = $y + $fProb1_2*$af_table{$locus}{$ca1};
		}
		if($fs2_1 == $s_min){
			$y = $y + $fProb2_1*$af_table{$locus}{$ca2};
		}
		if($fs2_2 == $s_min){
			$y = $y + $fProb2_2*$af_table{$locus}{$ca1};
		}
		
	}
	$pi = $x/$y;
	return ($pi);
} 

