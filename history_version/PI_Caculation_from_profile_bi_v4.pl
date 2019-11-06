#!/usr/bin/perl
use warnings;
use strict;

#get option;
my $tableFlag = 0;
my $inFlag = 0;
my $outFlag = 0;
my $mutFlag = 0;
my %opt;

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
	}elsif($ARGV[$i] eq '-mut'){
		$mutFlag = 1;
		$opt{'mut'} = $ARGV[++$i];
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

#get locus, C_profile and F_profile from profile file, get PI of each locus and CPI;
my %pi_locus;
open(INPUT, '< ', "$opt{'in'}") || die "can't open $opt{'in'} file:$!";
open(OUTPUT, '> ', "$opt{'out'}") || die "can't open $opt{'out'} file:$!";
while (my $input = <INPUT>){
	my @Child_genotype; 
	my $ca1; 
	my $ca2; 
	my @Father_genotype; 
	my $fa1; 
	my $fa2; 
	chomp $input;
	my @inputs = split "\t", $input;
	my $locus = $inputs[0];
	unless($locus =~ /AMEL/){
		@Child_genotype = sort {$a <=> $b} split "/", $inputs[1];
		$ca1 = $Child_genotype[0];
		$ca2 = $Child_genotype[1];
		@Father_genotype = sort {$a <=> $b} split "/", $inputs[2];
		$fa1 = $Father_genotype[0];
		$fa2 = $Father_genotype[1];
	}
	if(exists $af_table{$locus}){
		if(($fa1 == $ca1) || ($fa1 == $ca2) || ($fa2 == $ca1) || ($fa2 == $ca2)){
			my ($fomula, $pi) = &get_pi($inputs[0], $inputs[1], $inputs[2]);
			my $rounded_pi = sprintf("%.4f", $pi);
			$pi_locus{$locus} = $rounded_pi;
			print OUTPUT "$locus\t", "$fomula\t", "$rounded_pi\n";
		}else{
			my ($fomula, $pi) = &get_mutation_pi($inputs[0], $inputs[1], $inputs[2]);
			$pi_locus{$locus} = $pi;
			print OUTPUT "$locus\t", "$fomula\t", "$pi\n";
		}
	}elsif($locus =~ /AMEL/){
 		print OUTPUT "$locus\t", "1\t", "1\n";
	}else{
		print "the locus $locus is not in the $af_table_name allelic frequency table.\n";
		print OUTPUT "$locus\t", "----\t", "----\n";
	}
}

#caculate CPI;
my $cpi = 1;
for my $locus (keys %pi_locus){
	$cpi *= $pi_locus{$locus};
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

sub get_pi {

	#declare variable;
	my $pi;
	my $pog1;
	my $pog2;
	my $fomula;

	#get profile of mother, child and father;
	my ($locus, $C_profile, $F_profile) = @_;
	my @Child_genotype = sort {$a <=> $b} split "/", $C_profile;
	my $ca1 = $Child_genotype[0];
	my $ca2 = $Child_genotype[1];
	my @Father_genotype = sort {$a <=> $b} split "/", $F_profile;
	my $fa1 = $Father_genotype[0];
	my $fa2 = $Father_genotype[1];

	#according profile, caculate Paternal Index; 
	if($ca1 == $ca2){
		$pog1 = $ca1;
		if($fa1 == $fa2){
			if($fa1 == $ca1){
				if($af_table{$locus}{$pog1}){
					$pi = 1/$af_table{$locus}->{$pog1};
					$fomula = '1/'."$af_table{$locus}->{$pog1}";
				}else{
					$pi = 1/$af_table{$locus}->{min};
					$fomula = '1/'."$af_table{$locus}->{min}";
				}
				return($fomula, $pi);
			}else{print "on the locus $locus, the father can't provide POG\n";}
		}
		if($fa1 != $fa2){
			if(($fa1 == $pog1) || ($fa2 == $pog1)){
				if($af_table{$locus}{$pog1}){
					$pi = 1/(2*$af_table{$locus}->{$pog1});
					$fomula = '1/(2*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = 1/(2*$af_table{$locus}->{min});
					$fomula = '1/(2*'."$af_table{$locus}->{min}".')';
				}
				return($fomula, $pi);
			}else{print "on the locus $locus, the father can't provide POG\n";}
		}
	}
	if($ca1 != $ca2){
		if($fa1 == $fa2){
			if($fa1 == $ca1){
				$pog1 = $ca1;
				if($af_table{$locus}{$pog1}){
					$pi = 1/(2*$af_table{$locus}->{$pog1});
					$fomula = '1/(2*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = 1/(2*$af_table{$locus}->{min});
					$fomula = '1/(2*'."$af_table{$locus}->{min}".')';
				}
				return($fomula, $pi);
			}elsif($fa1 == $ca2){
				$pog1 = $ca2;
				if($af_table{$locus}{$pog1}){
					$pi = 1/(2*$af_table{$locus}->{$pog1});
					$fomula = '1/(2*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = 1/(2*$af_table{$locus}->{min});
					$fomula = '1/(2*'."$af_table{$locus}->{min}".')';
				}
				return($fomula, $pi);
			}else{print "on the locus $locus, the father can't provide POG\n";}
		}
		if($fa1 != $fa2){
			if(($fa1 == $ca1) && ($fa2 == $ca2)){
				$pog1 = $ca1;
				$pog2 = $ca2;
				if(($af_table{$locus}->{$pog1}) && ($af_table{$locus}->{$pog2})){
					$pi = ($af_table{$locus}->{$pog1} + $af_table{$locus}->{$pog2})/(4*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
					$fomula = '('."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{$pog2}".')/(4*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
				}elsif($af_table{$locus}->{$pog1}){
					$pi = ($af_table{$locus}->{$pog1} + $af_table{$locus}->{min})/(4*$af_table{$locus}->{$pog1}*$af_table{$locus}->{min});
					$fomula = '('."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{min}".')/(4*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{min}".')';
				}elsif($af_table{$locus}->{$pog2}){
					$pi = ($af_table{$locus}->{min} + $af_table{$locus}->{$pog2})/(4*$af_table{$locus}->{min}*$af_table{$locus}->{$pog2});
					$fomula = '('."$af_table{$locus}->{min}".'+'."$af_table{$locus}->{$pog2}".')/(4*'."$af_table{$locus}->{min}".'*'."$af_table{$locus}->{$pog2}".')';
				}else{
					$pi = ($af_table{$locus}->{min} + $af_table{$locus}->{min})/(4*$af_table{$locus}->{min}*$af_table{$locus}->{min});
					$fomula = '('."$af_table{$locus}->{min}".'+'."$af_table{$locus}->{min}".')/(4*'."$af_table{$locus}->{min}".'*'."$af_table{$locus}->{min}".')';
				}
				return($fomula, $pi);
			}elsif(($fa1 == $ca1) || ($fa2 == $ca1)){
				$pog1 = $ca1;
				if($af_table{$locus}->{$pog1}){
					$pi = 1/(4*$af_table{$locus}->{$pog1});
					$fomula = '1/(4*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = 1/(4*$af_table{$locus}->{min});
					$fomula = '1/(4*'."$af_table{$locus}->{min}".')';
				}
				return($fomula, $pi);
			}elsif(($fa1 == $ca2) || ($fa2 == $ca2)){
				$pog1 = $ca2;
				if($af_table{$locus}->{$pog1}){
					$pi = 1/(4*$af_table{$locus}->{$pog1});
					$fomula = '1/(4*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = 1/(4*$af_table{$locus}->{min});
					$fomula = '1/(4*'."$af_table{$locus}->{min}".')';
				}
				return($fomula, $pi);
			}else{print "on the locus $locus, the father can't provide POG\n";}
		}
	}
}

#define &get_mutation_pi subroutin
sub get_mutation_pi{
	#declare variable;
	my $pi;
	my $pog1;
	my $pog2;
	my $fomula;
	my $u;
	#get profile of mother, child and father;
	my ($locus, $C_profile, $F_profile) = @_;
	my @Child_genotype = sort {$a <=> $b} split "/", $C_profile;
	my $ca1 = $Child_genotype[0];
	my $ca2 = $Child_genotype[1];
	my @Father_genotype = sort {$a <=> $b} split "/", $F_profile;
	my $fa1 = $Father_genotype[0];
	my $fa2 = $Father_genotype[1];

	#according profile, caculate Paternal Index;
	if(exists $opt{'mut'}){
		$u = $opt{'mut'};
	}else{
		print "Error: this case exists mutation, please with optin -mut\n";
		printCMD();
		exit;
	}
	my $s; 
	if($ca1 == $ca2){
		$pog1 = $ca1;
		if($fa1 == $fa2){
			$s = abs($fa1-$pog1);
			if($af_table{$locus}->{$pog1}){
				$pi = ($u*0.1**($s-1))/(2*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'*0.1^('."$s".'-1)]/(2*'."$af_table{$locus}->{$pog1}".')';
			}else{
				$pi = ($u*0.1**($s-1))/(2*$af_table{$locus}->{min});
				$fomula = '['."$u".'*0.1^('."$s".'-1)]/(2*'."$af_table{$locus}->{min}".')';
			}
		}
		if($fa1 != $fa2){
			my $s1 = abs($fa1 - $pog1);
			my $s2 = abs($fa2 - $pog1);
			if($s1 < $s2){
				$s = $s1;
				if($af_table{$locus}->{$pog1}){
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{$pog1});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{min});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{min}".')';
				}
			}elsif($s1 > $s2){
				$s = $s2;
				if($af_table{$locus}->{$pog1}){
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{$pog1});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{min});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{min}".')';
				}
			}else{
				$s = $s1;
				if($af_table{$locus}->{$pog1}){
					$pi = ($u*0.1**($s-1))/(2*$af_table{$locus}->{$pog1});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(2*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = ($u*0.1**($s-1))/(2*$af_table{$locus}->{min});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(2*'."$af_table{$locus}->{min}".')';
				}
			}
		}
	}
	if($ca1 != $ca2){
		if($fa1 == $fa2){
			$pog1 = $ca1;
			$pog2 = $ca2;
			my $s1 = abs($fa1 - $pog1);
			my $s2 = abs($fa1 - $pog2);
			if($s1 < $s2){
				$s = $s1;
				if($af_table{$locus}->{$pog1}){
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{$pog1});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{$pog1}".')';
				}else{
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{min});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{min}".')';
				}
			}elsif($s1 > $s2){
				$s = $s2;
				if($af_table{$locus}->{$pog1}){
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{$pog2});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{$pog2}".')';
				}else{
					$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{min});
					$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{min}".')';
				}
			}else{
				$s = $s1;
				$pi = ($u*0.1**($s-1)*($af_table{$locus}->{$pog1}+$af_table{$locus}->{$pog2}))/(4*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
				$fomula = '['."$u".'*0.1^('."$s".'-1)*('."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{$pog2}".')]/4*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
			}
		}
		if($fa1 != $fa2){
			my $s1 = abs($fa1 - $ca1);
			my $s2 = abs($fa1 - $ca2);
			my $s3 = abs($fa2 - $ca1);
			my $s4 = abs($fa2 - $ca2);
			my @steps = ($s1,$s2,$s3,$s4);
			@steps = sort {$a <=> $b} @steps;
			my $min = $steps[0];
			$s = $min;
			if(($s1 == $min)&&($s2 != $min)&&($s3 != $min)&&($s4 != $min)){	#s1(7/11; 6/9) and (7/11; 6/13)(7/11; 8/9)(7/11; 8/13)
				$pog1 = $ca1;
				unless($af_table{$locus}->{$pog1}){
					$pog1 = 'min';
				}
				$pi = ($u*0.1**($s-1))/(8*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'0.1^('."$s".'-1)]/(8*'."$af_table{$locus}->{$pog1}".')';
			}elsif(($s1 != $min)&&($s2 == $min)&&($s3 != $min)&&($s4 != $min)){	#s2(7/11; 10/13) and (7/11; 12/13)
				$pog1 = $ca2;
				unless($af_table{$locus}->{$pog1}){
					$pog1 = 'min';
				}
				$pi = ($u*0.1**($s-1))/(8*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'0.1^('."$s".'-1)]/(8*'."$af_table{$locus}->{$pog1}".')';
			}elsif(($s1 != $min)&&($s2 != $min)&&($s3 == $min)&&($s4 != $min)){	#s3(7/11,5/6) and (7,11; 5,8)
				$pog1 = $ca1;
				unless($af_table{$locus}->{$pog1}){
					$pog1 = 'min';
				}
				$pi = ($u*0.1**($s-1))/(8*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'0.1^('."$s".'-1)]/(8*'."$af_table{$locus}->{$pog1}".')';
			}elsif(($s1 != $min)&&($s2 != $min)&&($s3 != $min)&&($s4 == $min)){	#s4(7/11; 5/10) and (7/11; 5/12)(7/11; 9/10)(7/11; 9/12)
				$pog1 = $ca2;
				unless($af_table{$locus}->{$pog1}){
					$pog1 = 'min';
				}
				$pi = ($u*0.1**($s-1))/(8*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'0.1^('."$s".'-1)]/(8*'."$af_table{$locus}->{$pog1}".')';
			}elsif(($s1 == $min)&&($s2 == $min)&&($s3 != $min)&&($s4 != $min)){	#s1s2(7/9; 8/11)
				$pog1 = $ca1;
				$pog2 = $ca2;
				if((! $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
					$pog2 = 'min';
				}elsif((! $af_table{$locus}->{$pog1}) && ( $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
				}elsif(( $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog2 = 'min';
				}
				$pi = ($u*0.1**($s-1)*($af_table{$locus}->{$pog1} + $af_table{$locus}->{$pog2}))/(8*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
				$fomula = '['."$u".'*0.1^('."$s".'-1)*('."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{$pog2}".')]/(8*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
			}elsif(($s1 == $min)&&($s2 != $min)&&($s3 == $min)&&($s4 != $min)){	#s1/s3(7,11; 6,8)
				$pog1 = $ca1;
				unless($af_table{$locus}->{$pog1}){
					$pog1 = 'min';
				}
				$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{$pog1}".')';
			}elsif(($s1 == $min)&&($s2 != $min)&&($s3 != $min)&&($s4 == $min)){	#s1s4(7/11; 6/10) and (7/11; 6/12)(7/11; 8/10)(7/11; 8/12)
				$pog1 = $ca1;
				$pog2 = $ca2;
				if((! $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
					$pog2 = 'min';
				}elsif((! $af_table{$locus}->{$pog1}) && ( $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
				}elsif(( $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog2 = 'min';
				}
				$pi = ($u*0.1**($s-1)*($af_table{$locus}->{$pog1} + $af_table{$locus}->{$pog2}))/(8*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
				$fomula = '['."$u".'*0.1^('."$s".'-1)*('."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{$pog2}".')]/(8*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
			}elsif(($s1 != $min)&&($s2 == $min)&&($s3 != $min)&&($s4 == $min)){	#s2s4(7/11; 10/12)
				$pog1 = $ca2;
				unless($af_table{$locus}->{$pog1}){
					$pog1 = 'min';
				}
				$pi = ($u*0.1**($s-1))/(4*$af_table{$locus}->{$pog1});
				$fomula = '['."$u".'*0.1^('."$s".'-1)]/(4*'."$af_table{$locus}->{$pog1}".')';
			}elsif(($s1 != $min)&&($s2 != $min)&&($s3 == $min)&&($s4 == $min)){	#s3s4(7/9; 5/8)
				$pog1 = $ca1;
				$pog2 = $ca2;
				if((! $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
					$pog2 = 'min';
				}elsif((! $af_table{$locus}->{$pog1}) && ( $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
				}elsif(( $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog2 = 'min';
				}
				$pi = ($u*0.1**($s-1)*($af_table{$locus}->{$pog1} + $af_table{$locus}->{$pog2}))/(8*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
				$fomula = '['."$u".'*0.1^('."$s".'-1)*('."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{$pog2}".')]/(8*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
			}elsif(($s1 == $min)&&($s2 == $min)&&($s3 != $min)&&($s4 == $min)){	#s1s2s4(7/9; 8/10)
				$pog1 = $ca1;
				$pog2 = $ca2;
				if((! $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
					$pog2 = 'min';
				}elsif((! $af_table{$locus}->{$pog1}) && ( $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
				}elsif(( $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog2 = 'min';
				}
				$pi = ($u*0.1**($s-1)*(2*$af_table{$locus}->{$pog1} + $af_table{$locus}->{$pog2}))/(8*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
				$fomula = '['."$u".'*0.1^('."$s".'-1)*(2*'."$af_table{$locus}->{$pog1}".'+'."$af_table{$locus}->{$pog2}".')]/(8*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
			}elsif(($s1 == $min)&&($s2 != $min)&&($s3 == $min)&&($s4 == $min)){	#s1s3s4(7/9; 6/8)
				$pog1 = $ca1;
				$pog2 = $ca2;
				if((! $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
					$pog2 = 'min';
				}elsif((! $af_table{$locus}->{$pog1}) && ( $af_table{$locus}->{$pog2})){
					$pog1 = 'min';
				}elsif(( $af_table{$locus}->{$pog1}) && (! $af_table{$locus}->{$pog2})){
					$pog2 = 'min';
				}
				$pi = ($u*0.1**($s-1)*($af_table{$locus}->{$pog1} + 2*$af_table{$locus}->{$pog2}))/(8*$af_table{$locus}->{$pog1}*$af_table{$locus}->{$pog2});
				$fomula = '['."$u".'*0.1^('."$s".'-1)*('."$af_table{$locus}->{$pog1}".'+2*'."$af_table{$locus}->{$pog2}".')]/(8*'."$af_table{$locus}->{$pog1}".'*'."$af_table{$locus}->{$pog2}".')';
			
			}
		}
	}
	return($fomula, $pi);
}

