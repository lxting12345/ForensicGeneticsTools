#!/usr/bin/perl
use warnings;
use strict;

#get option;
my $table_flag = 0;
my $in_flag = 0;
my $out_flag = 0;
my $father_mutation_flag = 0;
my $mother_mutation_flag = 0;
my %opt;
my $father_mut;
my $mother_mut;

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
	$father_mut = 0.002;		#default is 0.002
}else{
	$father_mut = $opt{'father_mutation'};
}
if($mother_mutation_flag == 0){
	$mother_mut = 0.0005;		#default is 0.0005
}else{
	$mother_mut = $opt{'mother_mutation'};
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
	my ($ma1, $ma2, $ca1, $ca2, $fa1, $fa2);
	chomp $input;
	$input =~ s/[\r\n]//g;
	my @inputs = split "\t", $input;
	my $locus = $inputs[0];
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
		
		my ($mc_pi, $fc_pi, $pi) = &get_pi_tri($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2);
		if($mc_pi ne ""){
			my $rounded_mc_pi = sprintf("%.4e", $mc_pi);
			$pi_mc{$locus} = $rounded_mc_pi;	
		}
		if($fc_pi ne ""){
			my $rounded_fc_pi = sprintf("%.4e", $fc_pi);
			$pi_fc{$locus} = $rounded_fc_pi;
		}
		if(defined $pi){
			my $rounded_pi = sprintf("%.4e", $pi);
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
$cpi_mc = sprintf("%.4e", $cpi_mc);
$RCP_mc = $cpi_mc / ($cpi_mc +1) ;
for my $locus (keys %pi_fc){
	$cpi_fc *= $pi_fc{$locus};
}
$cpi_fc = sprintf("%.4e", $cpi_fc);
$RCP_fc = $cpi_fc / ($cpi_fc + 1);
for my $locus (keys %pi_tri){
	$cpi_tri *= $pi_tri{$locus};
}
$cpi_tri = sprintf("%.4e", $cpi_tri);
$RCP_tri = $cpi_tri / ($cpi_tri + 1);
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
	print STDERR "\t\t2.with mutation\n";
	print STDERR "\t\tperl PI_Caculation_from_profile.pl -table /path/to/allele_frequency_table.txt -in case1.txt -out case1.pi.txt -father_mutation 0.002\n";
	print STDERR "\toption:\n";
	print STDERR "\t\t-table: Path to allele frequency table\n";
	print STDERR "\t\t-in: input profile file (.txt)\n";
	print STDERR "\t\t-out: output profile file (.pi.txt)\n";
	print STDERR "\t\t-father_mutation/-mother_mutation:  in paternal mutation, we recommend 0.002; in maternal mutation, we recommend 0.0005\n";
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
	if($ma1 && $ma2){
		$ms1_1 = abs($ma1 - $ca1);
		$ms1_2 = abs($ma1 - $ca2);
		$ms2_1 = abs($ma2 - $ca1);
		$ms2_2 = abs($ma2 - $ca2);
		@ms = sort {$a <=> $b} ($ms1_1, $ms1_2, $ms2_1,$ms2_2);
		$ms_min = $ms[0];
	}	
	if($fa1 && $fa2){
		$fs1_1 = abs($fa1 - $ca1);
		$fs1_2 = abs($fa1 - $ca2);
		$fs2_1 = abs($fa2 - $ca1);
		$fs2_2 = abs($fa2 - $ca2);	
		@fs = sort {$a <=> $b} ($fs1_1, $fs1_2, $fs2_1,$fs2_2);
		$fs_min = $fs[0];
	}
	#caculate probability of each
	if(defined $ms1_1){
		if($ms1_1 == 0){
			$mProb1_1 = 0.5;
		}else{
			if($ms1_1 == $ms_min){
				$mProb1_1 = 0.5*$mother_mut*0.5*0.1**($ms1_1 -1);
			}else{
				$mProb1_1 =0;
			}
		}
	}
	if(defined $ms1_2){
		if($ms1_2 == 0){
			$mProb1_2 = 0.5;
		}else{
			if($ms1_2 == $ms_min){
				$mProb1_2 = 0.5*$mother_mut*0.5*0.1**($ms1_2 -1);
			}else{
				$mProb1_2 =0;
			}
		}
	}
	if(defined $ms2_1){
		if($ms2_1 == 0){
			$mProb2_1 = 0.5;
		}else{
			if($ms2_1 == $ms_min){
				$mProb2_1 = 0.5*$mother_mut*0.5*0.1**($ms2_1 -1);
			}else{
				$mProb2_1 =0;
			}
		}
	}
	if(defined $ms2_2){
		if($ms2_2 == 0){
			$mProb2_2 = 0.5;
		}else{
			if($ms2_2 == $ms_min){
				$mProb2_2 = 0.5*$mother_mut*0.5*0.1**($ms2_2 -1);
			}else{
				$mProb2_2 =0;
			}
		}
	}

	if(defined $fs1_1){
		if($fs1_1 == 0){
			$fProb1_1 = 0.5;
		}else{
			if($fs1_1 == $fs_min){
				$fProb1_1 = 0.5*$father_mut*0.5*0.1**($fs1_1 -1);
			}else{
				$fProb1_1 =0;
			}
		}
	}

	if(defined $fs1_2){
		if($fs1_2 == 0){
			$fProb1_2 = 0.5;
		}else{
			if($fs1_2 == $fs_min){
				$fProb1_2 = 0.5*$father_mut*0.5*0.1**($fs1_2 -1);
			}else{
				$fProb1_2 =0;
			}
		}
	}

	if(defined $fs2_1){
		if($fs2_1 == 0){
			$fProb2_1 = 0.5;
		}else{
			if($fs2_1 == $fs_min){
				$fProb2_1 = 0.5*$father_mut*0.5*0.1**($fs2_1 -1);
			}else{
				$fProb2_1 =0;
			}
		}
	}

	if(defined $fs2_2){
		if($fs2_2 == 0){
			$fProb2_2 = 0.5;
		}else{
			if($fs2_2 == $fs_min){
				$fProb2_2 = 0.5*$father_mut*0.5*0.1**($fs2_2 -1);
			}else{
				$fProb2_2 =0;
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

	if($fa1 && $fa2){
		$fc_x = $fProb1_1*$af_table{$locus}{$ca2} + $fProb1_2*$af_table{$locus}{$ca1} + $fProb2_1*$af_table{$locus}{$ca2} + $fProb2_2*$af_table{$locus}{$ca1};
		$fc_y = 2*$af_table{$locus}{$ca1}*$af_table{$locus}{$ca2};
		$fc_pi = $fc_x/$fc_y;
	}else{
		$fc_pi = "";
	}

	if(($ma1 && $ma2) && ($fa1 && $fa2)){
		$x = $mProb1_1 * ($fProb1_2 + $fProb2_2) + $mProb1_2*($fProb1_1 + $fProb2_1) + $mProb2_1*($fProb1_2 + $fProb2_2) + $mProb2_2*($fProb1_1 + $fProb2_1);
		$y =  $mProb1_1*$af_table{$locus}{$ca2} + $mProb1_2*$af_table{$locus}{$ca1} +  $mProb2_1*$af_table{$locus}{$ca2} + $mProb2_2*$af_table{$locus}{$ca1};
		$pi = $x/$y;
	}
	return ($mc_pi, $fc_pi, $pi);
} 

