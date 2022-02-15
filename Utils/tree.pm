package YIW::tree;

#use strict;
#use warnings;

BEGIN {
	require Exporter;
# set the version for version checking
	our $VERSION = 1.00;
# Inherit from Exporter to export functions and variables
	our @ISA = qw(Exporter);
# Functions and variables which are exported by default
	our @EXPORT = qw(read_newick process_newick parse_newick write_newick write_tree_rec clean_tree label_tree copy_tree list_tree_nodes process_treenames tree_leaf_sets tree_collect_tbl tree_ultra tree_rescale tree_height tree_weight tree_hortonorder);
# Functions and variables which can be optionally exported
	our @EXPORT_OK = qw();
}

my $ndelim = "[/|\]";

return 1;

############################################################
#	read_newick($ftre,$EPSILON)
#		assumes single Newick tree in a file; allocates space; returns pointer to the root; makes ultra-short branches random ~$EPSILON
#	process_newick($rnod,$line,$EPSILON)
#		assumes single Newick line; allocates space; builds from $rnod; makes ultra-short branches random ~$EPSILON
#	parse_newick($rnod,$line,$EPSILON)
#		assumes clean Newick line; allocates space; builds from $rnod; makes ultra-short branches random ~$EPSILON
#	write_newick($root,$name,$flag,$pref)
#	write_tree_rec($root,$flag);
#	clean_tree($rnod,$tag1,tag2,...)
#	label_tree($rnod,$modl)
#		to "lab"
#		0 - PHYLIP, 1 - PAML, 2 - COUNT, 3 - HIGGS
#	copy_tree($rfro)
#	list_tree_nodes($rnod,$rlin,$rlle)
#	process_treenames($rnod,$word)
#	tree_leaf_sets($rnod)
#		to "snod"
#	tree_collect_tbl($rnod,$wadd)
#		to "tbl","add"
#	tree_ultra($rnod)
#		to "l", "tbl" and "deep"
#	tree_rescale($rnod,$quot)
#		to "l", "tbl" and "deep"
#	tree_height($rnod)
#		to "high"
#	tree_weight($rnod,$wsum,$wadd)
#		to "high", "tbl", "add", "wt"
#	tree_hortonorder($rnod)
#		to "hord"
############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Basic Fields
# n	node name (must be defined at leaves) [string]
# d	list of descendants [reference to array of references]
# p	parent [reference]
# l	in-branch length [real number]
# b	in-branch bootstrap [real number]
# i	in-branch info [string] (as in :[info]length)
# head	tree header ("[xxx]")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################################
#	read_newick($name,$EPSILON)
#	August 12 2015
############################################################
sub read_newick
{
 my $name = shift;
 my $EPSILON = shift;
 my $line = "";
 
 open HAND,"<$name" or die "Can't read \"$name\"";
 while(<HAND>){
  chomp;
  $line .= $_;
 }
 close HAND;
 
 my %tmp = ();
 process_newick(\%tmp,$line,$EPSILON);
 return \%tmp;
}

############################################################
#	process_newick($rnod,$line,$EPSILON)
#	February 19 2020
############################################################
sub process_newick
{
 my $rnod = shift;
 my $line = shift;
 my $EPSILON = shift;

 my ($head) = ($line=~m/^(\[.+\])/);	# get header
 $$rnod{"head"} = $head if($head ne "");
 
 $line =~ s/^\[.+\]\s*//;		# clean header
 $line =~ s/^[^(]*// if($line=~m/\(/);	# clean leader
 $line =~ s/;.*$//;			# clean trailer
 $line =~ tr/ //d;			# clean spaces
 #$line =~ s/\):[0-9.eE+-]+$/\)/;	# clean terminal branch length
 $line =~ s/\)[^)]+$/\)/;		# clean terminal labels and branch length #@@# NB: temporary solution!

 parse_newick($rnod,$line,$EPSILON);
 $$rnod{"p"} = "";
}

############################################################
#	parse_newick($rnod,$line,$EPSILON)
#	August 12 2015
############################################################
sub parse_newick
{
 my $rnod = shift;
 my $ilin = shift;
 my $EPSILON = shift;
 
 my $line = unwrap_line($ilin);		# strip (...)
 my @desc = split_line($line);		# parse descendants
 	#print "{$ilin}\n";
 	#print "[$line]\n";
 	#for(my $i=0;$i<@desc;$i++){
 	# printf "\t[ %s ]\n",$desc[$i];
 	#}
 if(@desc==1){				# none - terminal node
  $$rnod{"n"} = $desc[0];
  return;
 }
 
 my @tmp = ();				# internal node - set up
 $$rnod{"d"} = \@tmp;			# an array of descendants
 my $rdes = $$rnod{"d"};
 
 for(my $i=0;$i<@desc;$i++){		# descendants - one by one...
  my %tmpn = ();
  push @$rdes,\%tmpn;				# record in parent
  #pure style
  #my ($trex,$lenn) = ($desc[$i] =~ m/^(.*):(\d+\.*\d*e*[-+]*\d*)$/);	# extract branch length
  #BEAST style
  	#printf "***%s***\n",$desc[$i];
  my ($trex,$lenx) = ($desc[$i] =~ m/^(.*):([^:]+)$/);			# before last ":", after
  $trex = $desc[$i] unless($trex ne "");				# no branch length
  	#printf "[%s|%s]\n",$trex,$lenx;
  my ($tren,$boot) = ($trex =~ m/^(\(.+\))([^)]+)$/);			# extract bootstrap or label
  $tren = $trex unless($tren ne "");					# no bootstrap or label
  my ($iadd,$lenn) = ($lenx =~ m/\[([^\]]+)\](\d+\.*\d*e*[-+]*\d*)/);	# extract info and branch length
  $lenn = $lenx unless($lenn ne "");					# no info
  $tmpn{"l"} = $lenn;				# record branch length
  $tmpn{"l"} = $EPSILON*(0.5+rand()/2) if($lenn<$EPSILON/2);
  $tmpn{"b"} = $boot;				# record bootstrap or label
  $tmpn{"i"} = $iadd;				# record additional info
  $tmpn{"p"} = $rnod;				# record parent node - just in case
  parse_newick(\%tmpn,$tren,$EPSILON);		# go further down the tree
 }
}

############################################################
#	unwrap_line($line);
############################################################
sub unwrap_line
{
 my $line = shift;
 my $pbal = find_balance($line);
 if($pbal<0){ die "Unbalanced expression [$line]";}
 
 return $line if($pbal+1<length $line);		# not wrapped
 $line =~ s/^\s*\(//;
 $line =~ s/\)\s*$//;
 return $line;
}

############################################################
#	split_line($line);
############################################################
sub split_line
{
 my $line = shift;
 my @lout = ();
 while(length $line){
  my $pbal = find_balance($line);
  if($pbal<0){ die "Unbalanced expression [$line]";}
  push @lout,(substr $line,0,$pbal);	# left part
  $line = substr $line,$pbal+1;		# right part
 }
 return @lout;
}

############################################################
#	find_balance($line);
############################################################
sub find_balance
{
 my $line = shift;
 my $bal = 0;my $pos = 0;
 while($line=~m/(.)/g){
  if($1 eq ","){
   last if($bal==0);
  }elsif($1 eq "("){
   $bal++;
  }elsif($1 eq ")"){
   $bal--;
  }
  return -1 if($bal<0);		# -1 - unbalanced
  $pos++;
 }
 return -1 if($bal>0);		# -1 - unbalanced
 return $pos;
}

############################################################
#	write_newick($root,$name,$flag,$pref)
#	0x02	- include lengths
#	0x04	- include bootstraps
#	0x08	- include labels
#	returns prepared line (without ";")
############################################################
sub write_newick
{
 my $rroo = shift;
 my $name = shift;
 my $flag = shift;
 my $pref = shift;
 
 my $line = write_tree_rec($rroo,$flag);
 $line = "[".$pref."] ".$line if($pref ne "");
 open XNAH,">$name" or die "Can't write \"$name\"";
 print XNAH "$line;\n";
 close XNAH;
 return $line;
}

############################################################
#	write_tree_rec($root,$flag)
############################################################
sub write_tree_rec
{
 my $rnod = shift;
 my $flag = shift;
 my $rdes = $$rnod{"d"};
 my $line = "";
 
 unless(@$rdes){			# terminal node
  return $$rnod{"n"};
 }

 for(my $i=0;$i<@$rdes;$i++){
  next if($$rdes[$i] eq "");
  $line .= "," if($line ne "");
  $line .= write_tree_rec($$rdes[$i],$flag);
  $line .= $$rdes[$i]{"b"} if($$rdes[$i]{"b"} ne "" and $flag&0x04);
  $line .= "[".$$rdes[$i]{"lab"}."]" if($$rdes[$i]{"lab"} ne "" and $flag&0x08 and $$rdes[$i]{"n"} eq "");
  $line .= ":".$$rdes[$i]{"l"} if($$rdes[$i]{"l"} ne "" and $flag&0x02);
 }
 $line = "(".$line.")";
 return $line;
}

############################################################
#	print_tree($rnod,$level)
############################################################
# not used - debugging procedure
sub print_tree
{
 my $rnod = shift;
 my $level = shift;
 
 my $blen = $$rnod{"l"};
 my $name = $$rnod{"n"};
 my $rdes = $$rnod{"d"};
 if($name ne ""){
  for(my $i=0;$i<$level;$i++){ print "\t";} print "$name\n";
 }else{
  for(my $i=0;$i<$level;$i++){ print "\t";} print "INTERNAL\n";
 }
 for(my $i=0;$i<$level;$i++){ print "\t";} printf "L= %.4f\n",$blen;
 my $rlst = $$rnod{"dis"};
 foreach my $name (sort keys %$rlst){
  my $deep = $$rlst{$name};
  for(my $i=0;$i<$level;$i++){ print "\t";} printf "d= %.4f\t%s\n",$deep,$name;
 }
 
 	#for(my $i=0;$i<$level;$i++){ print "\t";} printf "SL= %s\n",$rrsl;
 	#for(my $i=0;$i<$level;$i++){ print "\t";} printf "SR= %s\n",$rrsr;
 
 for(my $i=0;$i<@$rdes;$i++){ print_tree($$rdes[$i],$level+1);}
}

############################################################
#	clean_tree($rnod,$tag1,tag2,...)
############################################################
sub clean_tree
{
 my $rnod = shift;
 my @ltag = @_;

 my $rdes = $$rnod{"d"};			# descendants

 foreach my $xtag (@ltag){ delete $$rnod{$xtag};}

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  clean_tree($$rdes[$i],@ltag);				# continue
 }
}

############################################################
#	label_tree($rnod,$modl)
############################################################
sub label_tree
{
 my $rnod = shift;
 my $modl = shift;

 my $ncnt = 0;
 my $lcnt = 0;
 
 if($modl==1){			# PAML mode
  label_leaves_rec($rnod,\$lcnt,1);
  $ncnt = $lcnt;
  label_nodes_rec($rnod,\$ncnt);
 }elsif($modl==2){		# COUNT mode
  label_nodre_rec($rnod,\$ncnt);
  label_leaves_rec($rnod,\$lcnt,1);
 }elsif($modl==3){		# HIGGS mode
  label_leaves_rec($rnod,\$lcnt,1);
  $ncnt = $lcnt;
  label_nodre_rec($rnod,\$ncnt);
 }else{				# PHYLIP mode
  label_nodes_rec($rnod,\$ncnt);
  label_leaves_rec($rnod,\$lcnt,1);
 }

 return $lcnt;
}

############################################################
#	label_leaves_rec($rnod,$rcnt,$flab)
############################################################
sub label_leaves_rec
{
 my $rnod = shift;
 my $rcnt = shift;
 my $flab = shift;

 my $rdes = $$rnod{"d"};			# descendants

 unless(@$rdes){				# terminal node - label and go
  $$rcnt++;
  if($flab==1){ $$rnod{"lab"} = $$rcnt;}
  else{ $$rnod{"lab"} = $$rnod{"n"};}
  return;
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  label_leaves_rec($$rdes[$i],$rcnt,$flab);			# continue
 }
}

############################################################
#	label_nodes_rec($rnod,$rcnt)
############################################################
sub label_nodes_rec
{
 my $rnod = shift;
 my $rcnt = shift;

 my $rdes = $$rnod{"d"};			# descendants

 unless(@$rdes){				# terminal node - do nothing
  return;
 }

 $$rcnt++;
 $$rnod{"lab"} = $$rcnt;			# label by number

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  label_nodes_rec($$rdes[$i],$rcnt);			# continue
 }
}

############################################################
#	label_nodre_rec($rnod,$rcnt)
############################################################
sub label_nodre_rec
{
 my $rnod = shift;
 my $rcnt = shift;

 my $rdes = $$rnod{"d"};			# descendants

 unless(@$rdes){				# terminal node - do nothing
  return;
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  label_nodre_rec($$rdes[$i],$rcnt);			# continue
 }

 $$rcnt++;
 $$rnod{"lab"} = $$rcnt;			# label by number
}

############################################################
#	copy_tree($rfro)
############################################################
sub copy_tree
{
 my $rfro = shift;
 my %tmpn = %$rfro;
 copy_tree_rec($rfro,\%tmpn);
 return \%tmpn;
}

############################################################
#	copy_tree_rec($rfro,$rtoo)
############################################################
sub copy_tree_rec
{
 my $rfro = shift;
 my $rtoo = shift;

 my $rdfr = $$rfro{"d"};			# descendants in the original
 return unless(@$rdfr);				# nothing to do
 
 my @tmpd = ();
 $$rtoo{"d"} = \@tmpd;				# descendants in the copy

 for(my $i=0;$i<@$rdfr;$i++){			# scan descendants
  my $rdes = $$rdfr[$i];				# i-th descendant
  my %tmpn = %$rdes;					# copy bulk
  $tmpn{"p"} = $rtoo;					# descendant->current
  $tmpd[$i] = \%tmpn;					# current->descendant
  copy_tree_rec($rdes,\%tmpn);				# proceed
 }
}

############################################################
#	list_tree_nodes($rnod,$rlin,$rlle)
############################################################
sub list_tree_nodes
{
 my $rnod = shift;
 my $rlin = shift;
 my $rlle = shift;

 my $rdes = $$rnod{"d"};			# descendants

 if(@$rdes){ push @$rlin,$rnod;}		# internal node
 else{ push @$rlle,$rnod;}			# terminal node

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  list_tree_nodes($$rdes[$i],$rlin,$rlle);		# continue
 }
}

############################################################
#	process_treenames($rnod,$word)
############################################################
sub process_treenames
{
 my $rnod = shift;
 my $word = shift;

 return if($word<=0);

 my $rdes = $$rnod{"d"};			# descendants

 unless(@$rdes){				# internal node
  my $name = (split/$ndelim/,$$rnod{"n"})[$word-1];
  $$rnod{"n"} = $name if($name ne "");
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  process_treenames($$rdes[$i],$word);			# continue
 }
}

############################################################
#	tree_leaf_sets($rnod)
############################################################
sub tree_leaf_sets
{
 my $rnod = shift;

 my %lset = ();

 my $rdes = $$rnod{"d"};
 
 unless(@$rdes){			# terminal node
  $lset{$$rnod{"n"}} = 1;
 }

 for(my $i=0;$i<@$rdes;$i++){
  my $rset = tree_leaf_sets($$rdes[$i]);
  foreach my $nn (keys %$rset){ $lset{$nn} = 1;}  
 }
 $$rnod{"snod"} = \%lset;
 return \%lset;
}

############################################################
#	tree_collect_tbl($rnod,$wadd)
############################################################
sub tree_collect_tbl
{
 my $rnod = shift;
 my $wadd = shift;
 
 $$rnod{"tbl"} = 0;				# init
 my $rdes = $$rnod{"d"};			# descendants
 if(@$rdes==0){					# bonus to leaves
  $$rnod{"tbl"} += $wadd*$$rnod{"high"};
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  $$rnod{"tbl"} += tree_collect_tbl($$rdes[$i],$wadd);	# continue
 }
 
 return $$rnod{"tbl"} + $$rnod{"l"};		# add branch
}

############################################################
#	tree_ultra($rnod)
############################################################
sub tree_ultra
{
 my $rnod = shift;

 $$rnod{"tbl"} = 0;				# init
 $$rnod{"deep"} = 0;				# init

 my $rdes = $$rnod{"d"};			# descendants

 if(@$rdes==0){					# leaf
  return;
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  tree_ultra($$rdes[$i]);				# continue forward
 }

 my $sumtl = 0;
 my $sumdr = 0;
 for(my $i=0;$i<@$rdes;$i++){			# scan descendants 2
  my $rden = $$rdes[$i];
  $sumtl += $$rden{"l"} + $$rden{"tbl"};
  $sumdr += ($$rden{"l"} + $$rden{"tbl"})/($$rden{"l"} + $$rden{"deep"});
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants 3
  my $rden = $$rdes[$i];
  my $quot = $sumtl/$sumdr/($$rden{"l"} + $$rden{"deep"});
  tree_rescale($rden,$quot);
  $$rnod{"tbl"} += $$rden{"l"} + $$rden{"tbl"};
  $$rnod{"deep"} = $$rden{"l"} + $$rden{"deep"} if($i==0);
 }
 return;
}

############################################################
#	tree_rescale($rnod,$quot)
############################################################
sub tree_rescale
{
 my $rnod = shift;
 my $quot = shift;
 
 $$rnod{"tbl"} = 0;				# init
 $$rnod{"deep"} = 0;				# init
 $$rnod{"l"} *= $quot;				# rescaled

 my $rdes = $$rnod{"d"};			# descendants

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  my $rden = $$rdes[$i];
  tree_rescale($rden,$quot);
  $$rnod{"tbl"} += $$rden{"l"} + $$rden{"tbl"};
  $$rnod{"deep"} = $$rden{"l"} + $$rden{"deep"} if($i==0);
 }
}

############################################################
#	tree_height($rnod)
############################################################
sub tree_height
{
 my $rnod = shift;

 
 my $high = $$rnod{"l"};			# init
 my $rpar =  $$rnod{"p"};			# parent
 $high += $$rpar{"high"} if($rpar ne "");		# + parent
 $$rnod{"high"} = $high;			# this node

 my $rdes = $$rnod{"d"};			# descendants

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  tree_height($$rdes[$i]);				# continue
 }
}

############################################################
#	tree_weight($rnod,$wsum,$wadd)
############################################################
sub tree_weight
{
 my $rnod = shift;
 my $wsum = shift;
 my $wadd = shift;
 
 
 tree_height($rnod);
 
 my $tble = tree_collect_tbl($rnod,$wadd);

 tree_collect_wt($rnod,$wsum);
}

############################################################
#	tree_collect_wt($rnod,$nodw)
############################################################
sub tree_collect_wt
{
 my $rnod = shift;
 my $nodw = shift;
 
 $$rnod{"wt"} = $nodw;				# init

 my $rdes = $$rnod{"d"};			# descendants

 my $sumx = 0;
 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  my $rden = $$rdes[$i];				# i-th descendant
  $sumx += $$rnod{"l"}+$$rden{"l"}+$$rden{"tbl"};
 }

 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  my $rden = $$rdes[$i];				# i-th descendant
  my $desw = $nodw*($$rnod{"l"}+$$rden{"l"}+$$rden{"tbl"})/$sumx;
  tree_collect_wt($$rdes[$i],$desw);
 }
}

############################################################
#	tree_hortonorder($rnod)
############################################################
sub tree_hortonorder
{
 my $rnod = shift;

 
 my $rdes = $$rnod{"d"};			# descendants
 
 if(@$rdes==0){					# leaf
  $$rnod{"hord"} = 1;
  return;
 }

 my @hode = ();
 
 for(my $i=0;$i<@$rdes;$i++){			# scan descendants
  my $rden = $$rdes[$i];
  tree_hortonorder($rden);					# continue
  $hode[$$rden{"hord"}]++;
 }

 my $horc = @hode;
 $horc-- if($hode[$horc-1]<2);
 $$rnod{"hord"} = $horc;
}

