#!/usr/bin/perl -w
 
use Tk;
use Tk::Table;
use strict;

#   Type names		Extension(s)	Mac File Type(s)
#---------------------------------------------------------
my @types =(
	["Text files",           [qw/.txt .csv/]],
	["Text files",           '',             'TEXT'],
	["Perl Scripts",         '.pl',		'TEXT'],
	["C Source Files",	['.c', '.h']],
	["All Source Files",     [qw/.tcl .c .h/]],
	["Image Files",		'.gif'],
	["Image Files",		['.jpeg', '.jpg']],
	["Image Files",   	'',		[qw/GIFF JPEG/]],
	["All files",		'*']
);

my $mw = MainWindow->new;
$mw->title("PI_Caculator");

#deal with frame1
my $af_file;
my $Fram1 = $mw->Frame()->pack(-side => 'top');
my $Lab1Fram1 = $Fram1->Label(-text => "1.please upload Allelic Frequency table:")->pack(-side => "left");
my $Ent1Fram1 = $Fram1->Entry()->pack(-side => 'left', -expand => 'yes', -fill => 'x');
my $But1Frame1 = $Fram1->Button(-text => "browse...",
                                -command => sub{GetAllelicFrequencyFile($Ent1Fram1)})->pack(-side => 'left');

sub GetAllelicFrequencyFile {
	my $ent = shift;	
	$af_file = $mw->getOpenFile(-filetypes => \@types); 
	if (defined $af_file and $af_file ne '') {
		$ent->delete(0, 'end');
		$ent->insert(0, $af_file);
		$ent->xview('end');
	}
	print $af_file,"\n";
}

#deal with frame2
my $geno_file;
my $Fram2 = $mw->Frame()->pack(-side => 'top');
my $lab1Fram2 = $Fram2->Label(-text => "2.please upload Genotype table:")->pack(-side => "left");
my $ent1Fram2 = $Fram2->Entry()->pack(-side => 'left', -expand => 'yes', -fill => 'x');
my $but1Frame2 = $Fram2->Button(-text => "browse...",
                                -command => sub{GetGenotypeFile($ent1Fram2)})->pack(-side => 'left');

sub GetGenotypeFile {
	my $ent = shift;	
	$geno_file = $mw->getOpenFile(-filetypes => \@types); 
	if (defined $geno_file and $geno_file ne '') {
		$ent->delete(0, 'end');
		$ent->insert(0, $geno_file);
		$ent->xview('end');
	}
	print "$geno_file\n";
}

#deal with frame3
my $father_mut;
my $mother_mut;
my $Fram3 = $mw->Frame()->pack(-side => 'top');
my $lab1Fram3 = $Fram3->Label(-text => "3.please set following parameter if exists mutation: father_mutation:")->pack(-side => 'left');
my $ent1Fram3 = $Fram3->Entry()->pack(-side => 'left', -expand => 'yes', -fill => 'x');
my $lab2Fram3 = $Fram3->Label(-text => "mother_mutation:")->pack(-side => 'left');
my $ent2Fram3 = $Fram3->Entry()->pack(-side => 'left', -expand => 'yes', -fill => 'x');

#test wether get value
#my $but1Frame3 = $Fram3->Button(-text => "test",
#                                -command => sub{$father_mut = $ent1Fram3->get();
#						$mother_mut = $ent2Fram3->get();
#						print "$father_mut $mother_mut\n";})->pack(-side => 'left');

#deal with frame4
my $format;
my $display_format = "e";
my $Fram4 = $mw->Frame()->pack(-side => 'top');
my $lab1Fram4 = $Fram4->Label(-text => "4.display format:")->pack(-side => 'left');
my $radio1Fram4 = $Fram4->Radiobutton(-text => "float", -value => "f", -variable => \$display_format)->pack(-side => 'left');
my $radio2Fram4 = $Fram4->Radiobutton(-text => "scientific  ", -value => "e", -variable => \$display_format)->pack(-side => 'left');

my $digit_format;
my $ent1Fram4 = $Fram4->Entry()->pack(-side => 'left', -expand => 'yes', -fill => 'x');
my $lab2Fram4 = $Fram4->Label(-text => "digit")->pack(-side => 'left');
#test
#my $but1Fram4 = $Fram4->Button(-text => "test",
#                                -command => sub{$digit_format = $ent1Fram4->get();
#						print "$display_format\t$digit_format\n";})->pack(-side => 'left');


#deal with frame5
my $af_table_name;
my @panal;
my %af_table;
my @kit;
my %result_table;
my %pi_mc;
my %pi_fc;
my %pi_tri;

my $Fram5 = $mw->Frame()->pack(-side => 'top');
my $but1Fram5 = $Fram5->Button(-text => "Caculate",
                            -command => \&Caculator)->pack(-side => 'left');
my $but4Fram5 = $Fram5->Button(-text => "exit",
                            -command => sub {exit})->pack(-side => 'left');

sub Caculator{
	&OpenAFtable; 
	&checkMut;
	&checkFormat;
	&openGenofile;
#	&getCPI;
	
	my $subwin1 = $mw->Toplevel;
	$subwin1->title('result');
	my $Frame1Subwin1 = $subwin1->Frame()->pack();
	my $row_number = $#kit+1;
	my $table = $Frame1Subwin1->Table(-columns => 7,
					-rows => "$row_number",  #need change row number
					-fixedrows => 1,
					-scrollbars => 'oe',
					-relief => 'raised');
	my @col_name = ('Locus', 'MotherGenotype', 'ChildGenotype', 'FatherGenotype', 'PI_MC', 'PI_FC', 'PI_Tri');
	foreach my $col (1 .. 7){
		my $tmp_label = $table->Label(-text => "$col_name[$col-1]", -relief =>'raised');
		$table->put(0, $col, $tmp_label);
	}
	foreach my $row (1 .. $row_number){
		foreach my $col (1 .. 7){
			
			my $tmp_label = $table->Label(-text => "$result_table{$kit[$row-1]}->[$col-1]",
                                  -padx => 2,
                                  -anchor => 'w',
                                  -background => 'white',
                                  -relief => "groove");
			$table->put($row, $col, $tmp_label);
		}
	}
	$table->pack();
	my $Frame2Subwin1 = $subwin1->Frame()->pack();
	my $but1 = $Frame2Subwin1->Button(-text => "save", -command => sub{&saveResult})->pack(-side => "left");
	my $but2 = $Frame2Subwin1->Button(-text => "close", -command => [$subwin1 => 'destroy'])->pack(-side => "left");
	undef @kit;
}


sub saveResult{
	my $save_FH;	
	my $save = $mw->getSaveFile(-filetypes => \@types,
				-initialfile => 'Untitled',
				-defaultextension => '.txt');
	open($save_FH, '>', $save) || die "can't open $save file:$!";
	for my $locus (@kit){
		for my $i (1..7){
			if($i < 7){
				print $save_FH "$result_table{$locus}->[$i-1]\t";
			}else{
				print $save_FH "$result_table{$locus}->[$i-1]\n";
			}
		}
	}
	close($save_FH);
	undef @kit;
}


sub openGenofile{
	my $genoFH;
	open($genoFH, '<', $geno_file) || die "can't open $geno_file file: $!\n";
	while (my $input = <$genoFH>){
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
		my $locus = uc($inputs[0]);
		push(@kit,$locus);
		$result_table{$locus} = \@inputs;
		my @mother_genotype = split "/", $inputs[1];
		$ma1 = $mother_genotype[0];
		$ma2 = $mother_genotype[1];
		my @child_genotype = split "/", $inputs[2];
		$ca1 = $child_genotype[0];
		$ca2 = $child_genotype[1];
		my @father_genotype = split "/", $inputs[3];
		$fa1 = $father_genotype[0];
		$fa2 = $father_genotype[1];
		if(exists $af_table{$locus}){
			my ($mc_pi, $fc_pi, $pi) = &get_pi_tri($locus, $fa1, $fa2, $ma1, $ma2, $ca1, $ca2);
			if(defined $mc_pi){
				my $rounded_mc_pi = sprintf("$format", $mc_pi);
				$pi_mc{$locus} = $rounded_mc_pi;
				push( @{$result_table{$locus}}, $pi_mc{$locus})	;
			}
			if(defined $fc_pi){
				my $rounded_fc_pi = sprintf("$format", $fc_pi);
				$pi_fc{$locus} = $rounded_fc_pi;
				push( @{$result_table{$locus}}, $pi_fc{$locus})	;
			}
			if(defined $pi){
				my $rounded_pi = sprintf("$format", $pi);
				$pi_tri{$locus} = $rounded_pi;
				push( @{$result_table{$locus}}, $pi_tri{$locus});
			}	
		}elsif($locus =~ /AMEL/){
			$pi_mc{$locus} = 1;
			push( @{$result_table{$locus}}, $pi_mc{$locus})	;
			$pi_fc{$locus} = 1;
			push( @{$result_table{$locus}}, $pi_fc{$locus})	;
 			$pi_tri{$locus} = 1;
			push( @{$result_table{$locus}}, $pi_tri{$locus});
		}
	}

	#test wether get correct value
	for my $locus (@kit){
		print "$locus\t","$pi_mc{$locus}\t","$pi_fc{$locus}\t$pi_tri{$locus}\n";
	}
	
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
	}

	if($fa1 && $fa2){
		$fc_x = $fProb1_1*$af_table{$locus}{$ca2} + $fProb1_2*$af_table{$locus}{$ca1} + $fProb2_1*$af_table{$locus}{$ca2} + $fProb2_2*$af_table{$locus}{$ca1};
		$fc_y = 2*$af_table{$locus}{$ca1}*$af_table{$locus}{$ca2};
		$fc_pi = $fc_x/$fc_y;
	}

	if(($ma1 && $ma2) && ($fa1 && $fa2)){
		$x = $mProb1_1 * ($fProb1_2 + $fProb2_2) + $mProb1_2*($fProb1_1 + $fProb2_1) + $mProb2_1*($fProb1_2 + $fProb2_2) + $mProb2_2*($fProb1_1 + $fProb2_1);
		$y =  $mProb1_1*$af_table{$locus}{$ca2} + $mProb1_2*$af_table{$locus}{$ca1} +  $mProb2_1*$af_table{$locus}{$ca2} + $mProb2_2*$af_table{$locus}{$ca1};
		$pi = $x/$y;
	}
	return ($mc_pi, $fc_pi, $pi);
} 

sub checkFormat{
	my $digit = $ent1Fram4->get();
	if(defined $digit and $digit ne ''){
		$digit_format = '.'."$digit";
	}else{
		$digit_format = '';
	}
	$format = '%'."$digit_format"."$display_format";
	#test
	#print "$display_format\t$digit_format\t$format\n";
}

sub checkMut{
	$father_mut = $ent1Fram3->get();
	unless(defined $father_mut and $father_mut ne ''){
		$father_mut = 0.002;
	}
	$mother_mut = $ent2Fram3->get();
	unless(defined $mother_mut and $mother_mut ne ''){
		$mother_mut = 0.0005;
	}
	#test wether get the correct value
	#print "$father_mut\t$mother_mut\n";
}

sub OpenAFtable{
	#open allelic frequency table and generate hash of hash
	#and get minimal alleleic frequency
	my $AFfileHandle;	
	open($AFfileHandle, '< ', "$af_file") || die "can't open $af_file file:$!";
	while (my $line = <$AFfileHandle>){
		chomp $line;
		if($line =~ /^#(\w+)/){
			$af_table_name = $1;
			next;
		}
		my @lines = split ("\t", $line);
		my $locus = uc($lines[0]);
		for my $i (1..$#lines){
			my ($bin, $frequency) = split("=", $lines[$i]);
			$af_table{$locus}{$bin}= "$frequency";
		}
	}	
	close ($AFfileHandle);
	@panal = keys %af_table;
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
}


=pod
#sub getCPI{
#	my ($cpi_mc,$cpi_fc,$cpi_tri) = (1, 1, 1);
#	for my $locus (@panal){
#		if(defined $pi_mc{$locus}){
#			$cpi_mc *= $pi_mc{$locus};
#		}
#		if(defined $pi_fc{$locus}){
#			$cpi_fc *= $pi_fc{$locus};
#		}
#		if(defined $pi_tri{$locus}){
#			$cpi_tri *= $pi_tri{$locus};
#		}
#	}
#	return ($cpi_mc, $cpi_fc, $cpi_tri);
#}

=cut
MainLoop;
