#!/usr/bin/perl
use warnings;
use strict;

#usage: perl PI_Caculation.pl ./AF/CNAS_AF_39.txt

my $af_table_name;
my %af_table;
my @panal;
open(AF, '< ', "$ARGV[0]") || die "can't open $ARGV[0] file:$!";
while (my $line = <AF>){
	chomp $line;
	if($line =~ /^#(\w+)/){
		$af_table_name = uc $1;
		next;
	}
	my @lines = split ("\t", $line);
	my $locus = $lines[0];
	for my $i (1..$#lines){
		my ($bin, $frenquency) = split("=", $lines[$i]);
		$af_table{$locus}{$bin}= "$frenquency";
	}
}
@panal = keys %af_table;

print "the allelic frequent table is $af_table_name.\n";
print "the loci encluded in the allelic table are: \n@panal\n\n";
print $af_table{"D19S433"}->{"16.2"}, "\n";
close (AF);

my $input;
my @inputs;
my $pi_value;

while ($input = <STDIN>){
	chomp $input;
	@inputs = split " ", $input;
	my $locus = $inputs[0];
	unless($af_table{$locus}){
		print "the locus is not in the $af_table_name allelic frequency table.\n";
		exit;
	}else{
		&get_pi($inputs[0], $inputs[1], $inputs[2], $inputs[3]);
	}
}




sub get_pi {

#declare variable;
	my $pi;
	my $pog1;
	my $pog2;

#get profile of mother, child and father;
	my ($locus, $M_profile, $C_profile, $F_profile) = @_;
	my @Mother_genotype = sort {$a <=> $b} split '/', $M_profile;
	my $ma1 = $Mother_genotype[0];
	my $ma2 = $Mother_genotype[1];
	my @Child_genotype = sort {$a <=> $b} split "/", $C_profile;
	my $ca1 = $Child_genotype[0];
	my $ca2 = $Child_genotype[1];
	my @Father_genotype = sort {$a <=> $b} split "/", $F_profile;
	my $fa1 = $Father_genotype[0];
	my $fa2 = $Father_genotype[1];
	print "$ma1,$ma2\t", "$ca1,$ca2\t", "$fa1,$fa2\t", "\n";
#according profile, caculate PI;         
	if($ma1 == $ma2){								#1 mother is homozygote(pp,?,?);      
		if(($ca1 == $ma1) || ($ca2 == $ma1)){					#weather child contain one allele of mother;
			if($ca1 == $ca2){						# 1.1 child is homozygote(pp,pp,?);
				$pog1 = $ca1;						#determin pog
				if(($fa1 == $pog1) || ($fa2 == $pog1)){			#weather father can provide POG;
					if($fa1 == $fa2){				#1.1.1 father is homzygote(pp,pp,pp);
						$pi = 1/$af_table{$locus}->{$pog1};
						print 'PI_',"$locus",' = 1/',"$af_table{$locus}->{$pog1}",' = ',"$pi", "\n";
					}
					if($fa1 != $fa2){				#1.1.2 father is heterzygote(pp,pp,pq);
						$pi = 1/(2*$af_table{$locus}->{$pog1});
						print 'PI_',"$locus",' = 1/(2*',"$af_table{$locus}->{$pog1}",') = ',"$pi", "\n";
					}
				}else{
					print "on the locus $locus, the father can't provide POG.\n";
				}
			}
			if($ca1 != $ca2){						#1.2 child is heterozygote(pp,pq,?);
				if($ca1 == $ma1){					#determin pog
					$pog1 = $ca2;
				}
				if($ca2 == $ma1){
					$pog1 = $ca1;
				}
				if(($fa1 == $pog1) || ($fa2 == $pog1)){			#weather father can provide POG;
					if($fa1 == $fa2){				#1.2.1 father is homozygote and can provide pog (pp pq qq);
						$pi = 1/$af_table{$locus}-> {$pog1};
						print 'PI_',"$locus",' = 1/',"$af_table{$locus}->{$pog1}",' = ',"$pi", "\n";
					}
					if($fa1 != $fa2){				#1.2.2 father is heterozygote(pp pq qr, pp pq pq);
						$pi = 1/(2*$af_table{$locus}-> {$pog1});
						print 'PI_',"$locus",' = 1/(2*',"$af_table{$locus}->{$pog1}",') = ',"$pi", "\n";
					}
				}else{
					print "on the locus $locus, the father can't provide POG.\n";
				}
			}
		}else{
			print "on the locus $locus, the mother can't provide MOG.\n";
		}
	}

	if($ma1 != $ma2){									#2.mother is heterozygote(pq,?,?);
		if((($ca1 == $ma1) || ($ca2 == $ma1)) || (($ca1 == $ma2) || ($ca2 == $ma2))){	#weather child contain one allele of mother;
			if($ca1 == $ca2){							#2.1child is homozygote(pq,pp,?; pq,qq,?);
				$pog1 = $ca1;							#determin POG;
				if(($fa1 == $pog1) || ($fa2 == $pog1)){				#wether the father can provide pog;
					if($fa1 == $fa2){					#2.1.1 father is homozygote(pq,pp,pp; pq,qq,qq);
						$pi = 1/$af_table{$locus}-> {$pog1};
						print 'PI_',"$locus",' = 1/',"$af_table{$locus}->{$pog1}",' = ',"$pi", "\n";
					}
					if($fa1 != $fa2){					#2.1.2 father is heterozygote(pq,pp,pr/pq; pq,qq,pq/qr);
						$pi = 1/(2*$af_table{$locus}-> {$pog1});
						print 'PI_',"$locus",' = 1/(2*',"$af_table{$locus}->{$pog1}",') = ',"$pi", "\n";
					}
				}else{
					print "on the locus $locus, the father can't provide POG.\n";
				}
			}
			if($ca1 != $ca2){							#child is heterozygote(pq,?,?);
				if(($ca1 == $ma1) && ($ca2 == $ma2)){				#2.2 child is heterozygote and same to the mather(pq,pq,?)	
					my $pog1 = $ca1;					#determin pog
					my $pog2 = $ca2;
					if((($fa1 == $pog1) || ($fa2 == $pog1)) || (($fa1 == $pog2) || ($fa2 == $pog2))){	#weather father can provide pog;
						if($fa1 == $fa2){				#2.2.1 father is homozygote(pq, pq, pp; pq,pq,qq)
							if(($fa1 == $pog1) || ($fa1 == $pog2)){
								$pi = 1/($af_table{$locus}-> {$pog1} + $af_table{$locus}-> {$pog2});
								print 'PI_',"$locus",' = 1/(',"$af_table{$locus}->{$pog1}",'+', "$af_table{$locus}->{$pog2}",') = ',"$pi", "\n";
							}
						}
						if($fa1 != $fa2){				#father is heterozygote(pq,pq,?);
							if(($fa1 == $pog1) && ($fa2 == $pog2)){	#2.2.2 father is same to child(pq,pq,pq)
								$pi = 1/($af_table{$locus}-> {$pog1} + $af_table{$locus}-> {$pog2});
								print 'PI_',"$locus",' = 1/(',"$af_table{$locus}->{$pog1}",'+', "$af_table{$locus}->{$pog2}",') = ',"$pi", "\n";
							}
							if((($fa1 == $pog1) && ($fa2 != $pog2)) || ($fa2 == $pog1) || ($fa1 == $pog2) || (($fa2 == $pog2) && ($fa1 != $pog1))){
												#2.2.3 fathe is hterozygote and not same to child(pq,pq,pr;pq,pq,op;  pq,pq,qr;pq,pq,oq);	
								$pi = 1/(2*($af_table{$locus}-> {$pog1} + $af_table{$locus}-> {$pog2}));
								print 'PI_',"$locus",' = 1/[2*(',"$af_table{$locus}->{$pog1}",'+', "$af_table{$locus}->{$pog2}",')] = ',"$pi", "\n";
							}
						}
					}else{
						print "on the locus $locus, the fater can't provide POG.\n";
					}
				}
				if((($ca1 == $ma1) && ($ca2 != $ma2)) || ($ca1 == $ma2)){	#2.3 child is heterzygote but different to mather(pq,pr,?; pq,qr,?)
					$pog1 = $ca2;
					if(($fa1 == $pog1) || ($fa2 == $pog1)){
						if($fa1 == $fa2){				#2.3.1 father is homo(pq,pr,rr)
							$pi = 1/$af_table{$locus}-> {$pog1};
							print 'PI_',"$locus",' = 1/',"$af_table{$locus}->{$pog1}",' = ',"$pi", "\n";
						}
						if($fa1 != $fa2){				#2.3.2 father is heter(pq,pr,rs)
							$pi = 1/(2*$af_table{$locus}-> {$pog1});
							print 'PI_',"$locus",' = 1/(2*',"$af_table{$locus}->{$pog1}", ') = ', "$pi", "\n";
						}
					}else{
						print "on the locus $locus, the father can't provide POG.\n";
					}
				}
				if(($ca2 == $ma1) || (($ca2 == $ma2) && ($ca1 != $ma1))){	#2.3 child is hereo but diff to mather(pq,rp,?)
					$pog1 = $ca1;
					if(($fa1 == $pog1) || ($fa2 == $pog1)){
						if($fa1 == $fa2){				#2.3.1 father is homo(pq,pr,rr)
							$pi = 1/$af_table{$locus}-> {$pog1};
							print 'PI_',"$locus",' = 1/',"$af_table{$locus}->{$pog1}",' = ',"$pi", "\n";
						}
						if($fa1 != $fa2){				#2.3.2 father is heter(pq,pr,rs)
							$pi = 1/(2*$af_table{$locus}-> {$pog1});
							print 'PI_',"$locus",' = 1/(2*',"$af_table{$locus}->{$pog1}",') = ',"$pi", "\n";
						}
					}else{
						print "on the locus $locus, the father can't provide POG.\n";
					}
				}
			}
				
		}else{
			print "on the locus $locus, the mother can't provide MOG.\n";
		}
	}
}	


