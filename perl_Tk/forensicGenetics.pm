



sub checkAF{
	my $subwin1 = $mw->Toplevel;
	$subwin1->title("check allelic frequency table");
	my $win1Fram1 = $subwin1->Fram()->pack(-side => 'top');
	my $table1Win1Fram1 = $win1Fram1->Table(-columns => "$colums",
                                                -rows => "$numberOfLocus",
                                                -fixedrows => 1,
                                                -scrollbars => 'oe',
                                                -relief => 'raised');

	
	
}
