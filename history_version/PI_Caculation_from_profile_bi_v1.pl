#!/usr/bin/perl
use warnings;
use strict;

#open allelic frequency table and generate hash of hash;
my $af_table_name;
my %af_table;
my @panal;
open(AF, '< ', "$ARGV[0]") || die "can't open $ARGV[0] file:$!";
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


#test:print each allele on each locus;
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

#test:  weather the min allele frequency of each locus is in hash of hash %af_table;
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

#get locus, child_profile and father_profile from profile file, get PI of each locus and CPI;
my %pi_locus;
open(INPUT, '< ', "$ARGV[1]") || die "can't open $ARGV[1] file:$!";
open(OUTPUT, '> ', "$ARGV[2]") || die "can't open $ARGV[2] file:$!";
while (my $input = <INPUT>){
	chomp $input;
	my @inputs = split "\t", $input;
	my $locus = $inputs[0];
	if(exists $af_table{$locus}){
		my ($fomula, $pi) = &get_pi($inputs[0], $inputs[1], $inputs[2]);
		$pi_locus{$locus} = $pi;
		print OUTPUT "$locus\t", "$fomula\t", "$pi\n";
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

#define subroutin
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


