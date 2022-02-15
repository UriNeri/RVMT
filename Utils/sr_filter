#!/usr/bin/perl5.16

############################################################
#	Template:	0.4, February 02, 2009
#	Parent Path:	/home/wolf/bin on frosty
############################################################

############################################################
#	System etc
############################################################
$Pver = "0.8";
$Pdat = "May 10 2020";

($Pnam) = ($0 =~ m/([^\/]+)$/);
($Path) = ($0 =~ m/^(.+)\/[^\/]+$/);$Path = "." unless($Path);

$Ppid = $$;	# Process ID

############################################################
#	Definitions
############################################################
$gncut = 0x7fffffff;
$gnbar = 0;
$grcut = 1.0;
$grbar = 0.0;
$iccut = -100.0;
$icbar = 100.0;
$hocut = -100.0;
$hobar = 100.0;
$pseudo = 0;
$pseudq = 0.0;
$igno = 0;
$prof = 0;
$prxx = 0;
$cons = 0;
$hpro = 0;
$grcon = 0.499;
$hocon = 0.333;
$grswe = 0.499;
$nacon = "CONSENSUS";
$EPSILON = 1e-6;

############################################################
#	Data
############################################################
%aafreq = (
"A" => 0.07422,
"C" => 0.02469,
"D" => 0.05363,
"E" => 0.05431,
"F" => 0.04742,
"G" => 0.07415,
"H" => 0.02621,
"I" => 0.06792,
"K" => 0.05816,
"L" => 0.09891,
"M" => 0.02499,
"N" => 0.04465,
"P" => 0.03854,
"Q" => 0.03426,
"R" => 0.05161,
"S" => 0.05723,
"T" => 0.05089,
"V" => 0.07292,
"W" => 0.01303,
"Y" => 0.03228,
);

%simmat = (
 "AA" =>   4 , "AR" =>  -1 , "AN" =>  -1 , "AD" =>  -2 , "AC" =>   0 , "AQ" =>  -1 , "AE" =>  -1 , "AG" =>   0 , "AH" =>  -2 , "AI" =>  -1 , "AL" =>  -1 , "AK" =>  -1 , "AM" =>  -1 , "AF" =>  -2 , "AP" =>  -1 , "AS" =>   1 , "AT" =>   0 , "AW" =>  -3 , "AY" =>  -2 , "AV" =>   0 ,
 "RA" =>  -1 , "RR" =>   5 , "RN" =>   0 , "RD" =>  -1 , "RC" =>  -3 , "RQ" =>   1 , "RE" =>   0 , "RG" =>  -2 , "RH" =>   0 , "RI" =>  -3 , "RL" =>  -2 , "RK" =>   2 , "RM" =>  -1 , "RF" =>  -3 , "RP" =>  -2 , "RS" =>  -1 , "RT" =>  -1 , "RW" =>  -3 , "RY" =>  -2 , "RV" =>  -2 ,
 "NA" =>  -1 , "NR" =>   0 , "NN" =>   6 , "ND" =>   1 , "NC" =>  -2 , "NQ" =>   0 , "NE" =>   0 , "NG" =>   0 , "NH" =>   1 , "NI" =>  -3 , "NL" =>  -3 , "NK" =>   0 , "NM" =>  -2 , "NF" =>  -3 , "NP" =>  -2 , "NS" =>   1 , "NT" =>   0 , "NW" =>  -4 , "NY" =>  -2 , "NV" =>  -3 ,
 "DA" =>  -2 , "DR" =>  -1 , "DN" =>   1 , "DD" =>   6 , "DC" =>  -3 , "DQ" =>   0 , "DE" =>   2 , "DG" =>  -1 , "DH" =>  -1 , "DI" =>  -3 , "DL" =>  -3 , "DK" =>  -1 , "DM" =>  -3 , "DF" =>  -3 , "DP" =>  -1 , "DS" =>   0 , "DT" =>  -1 , "DW" =>  -4 , "DY" =>  -3 , "DV" =>  -3 ,
 "CA" =>   0 , "CR" =>  -3 , "CN" =>  -2 , "CD" =>  -3 , "CC" =>   9 , "CQ" =>  -3 , "CE" =>  -3 , "CG" =>  -2 , "CH" =>  -3 , "CI" =>  -1 , "CL" =>  -1 , "CK" =>  -3 , "CM" =>  -1 , "CF" =>  -2 , "CP" =>  -3 , "CS" =>  -1 , "CT" =>  -1 , "CW" =>  -2 , "CY" =>  -2 , "CV" =>  -1 ,
 "QA" =>  -1 , "QR" =>   1 , "QN" =>   0 , "QD" =>   0 , "QC" =>  -3 , "QQ" =>   5 , "QE" =>   2 , "QG" =>  -2 , "QH" =>   1 , "QI" =>  -3 , "QL" =>  -2 , "QK" =>   1 , "QM" =>   0 , "QF" =>  -3 , "QP" =>  -1 , "QS" =>   0 , "QT" =>  -1 , "QW" =>  -2 , "QY" =>  -1 , "QV" =>  -2 ,
 "EA" =>  -1 , "ER" =>   0 , "EN" =>   0 , "ED" =>   2 , "EC" =>  -3 , "EQ" =>   2 , "EE" =>   5 , "EG" =>  -2 , "EH" =>   0 , "EI" =>  -3 , "EL" =>  -3 , "EK" =>   1 , "EM" =>  -2 , "EF" =>  -3 , "EP" =>  -1 , "ES" =>   0 , "ET" =>  -1 , "EW" =>  -3 , "EY" =>  -2 , "EV" =>  -2 ,
 "GA" =>   0 , "GR" =>  -2 , "GN" =>   0 , "GD" =>  -1 , "GC" =>  -2 , "GQ" =>  -2 , "GE" =>  -2 , "GG" =>   6 , "GH" =>  -2 , "GI" =>  -3 , "GL" =>  -4 , "GK" =>  -1 , "GM" =>  -2 , "GF" =>  -3 , "GP" =>  -2 , "GS" =>   0 , "GT" =>  -2 , "GW" =>  -2 , "GY" =>  -3 , "GV" =>  -3 ,
 "HA" =>  -2 , "HR" =>   0 , "HN" =>   1 , "HD" =>  -1 , "HC" =>  -3 , "HQ" =>   1 , "HE" =>   0 , "HG" =>  -2 , "HH" =>   7 , "HI" =>  -3 , "HL" =>  -3 , "HK" =>  -1 , "HM" =>  -1 , "HF" =>  -1 , "HP" =>  -2 , "HS" =>  -1 , "HT" =>  -2 , "HW" =>  -2 , "HY" =>   2 , "HV" =>  -3 ,
 "IA" =>  -1 , "IR" =>  -3 , "IN" =>  -3 , "ID" =>  -3 , "IC" =>  -1 , "IQ" =>  -3 , "IE" =>  -3 , "IG" =>  -3 , "IH" =>  -3 , "II" =>   4 , "IL" =>   2 , "IK" =>  -3 , "IM" =>   1 , "IF" =>   0 , "IP" =>  -3 , "IS" =>  -2 , "IT" =>  -1 , "IW" =>  -2 , "IY" =>  -1 , "IV" =>   3 ,
 "LA" =>  -1 , "LR" =>  -2 , "LN" =>  -3 , "LD" =>  -3 , "LC" =>  -1 , "LQ" =>  -2 , "LE" =>  -3 , "LG" =>  -4 , "LH" =>  -3 , "LI" =>   2 , "LL" =>   4 , "LK" =>  -2 , "LM" =>   2 , "LF" =>   0 , "LP" =>  -3 , "LS" =>  -2 , "LT" =>  -1 , "LW" =>  -2 , "LY" =>  -1 , "LV" =>   1 ,
 "KA" =>  -1 , "KR" =>   2 , "KN" =>   0 , "KD" =>  -1 , "KC" =>  -3 , "KQ" =>   1 , "KE" =>   1 , "KG" =>  -1 , "KH" =>  -1 , "KI" =>  -3 , "KL" =>  -2 , "KK" =>   4 , "KM" =>  -1 , "KF" =>  -3 , "KP" =>  -1 , "KS" =>   0 , "KT" =>  -1 , "KW" =>  -3 , "KY" =>  -2 , "KV" =>  -2 ,
 "MA" =>  -1 , "MR" =>  -1 , "MN" =>  -2 , "MD" =>  -3 , "MC" =>  -1 , "MQ" =>   0 , "ME" =>  -2 , "MG" =>  -2 , "MH" =>  -1 , "MI" =>   1 , "ML" =>   2 , "MK" =>  -1 , "MM" =>   5 , "MF" =>   0 , "MP" =>  -2 , "MS" =>  -1 , "MT" =>  -1 , "MW" =>  -1 , "MY" =>  -1 , "MV" =>   1 ,
 "FA" =>  -2 , "FR" =>  -3 , "FN" =>  -3 , "FD" =>  -3 , "FC" =>  -2 , "FQ" =>  -3 , "FE" =>  -3 , "FG" =>  -3 , "FH" =>  -1 , "FI" =>   0 , "FL" =>   0 , "FK" =>  -3 , "FM" =>   0 , "FF" =>   6 , "FP" =>  -4 , "FS" =>  -2 , "FT" =>  -2 , "FW" =>   1 , "FY" =>   3 , "FV" =>  -1 ,
 "PA" =>  -1 , "PR" =>  -2 , "PN" =>  -2 , "PD" =>  -1 , "PC" =>  -3 , "PQ" =>  -1 , "PE" =>  -1 , "PG" =>  -2 , "PH" =>  -2 , "PI" =>  -3 , "PL" =>  -3 , "PK" =>  -1 , "PM" =>  -2 , "PF" =>  -4 , "PP" =>   7 , "PS" =>  -1 , "PT" =>  -1 , "PW" =>  -4 , "PY" =>  -3 , "PV" =>  -2 ,
 "SA" =>   1 , "SR" =>  -1 , "SN" =>   1 , "SD" =>   0 , "SC" =>  -1 , "SQ" =>   0 , "SE" =>   0 , "SG" =>   0 , "SH" =>  -1 , "SI" =>  -2 , "SL" =>  -2 , "SK" =>   0 , "SM" =>  -1 , "SF" =>  -2 , "SP" =>  -1 , "SS" =>   4 , "ST" =>   1 , "SW" =>  -3 , "SY" =>  -2 , "SV" =>  -2 ,
 "TA" =>   0 , "TR" =>  -1 , "TN" =>   0 , "TD" =>  -1 , "TC" =>  -1 , "TQ" =>  -1 , "TE" =>  -1 , "TG" =>  -2 , "TH" =>  -2 , "TI" =>  -1 , "TL" =>  -1 , "TK" =>  -1 , "TM" =>  -1 , "TF" =>  -2 , "TP" =>  -1 , "TS" =>   1 , "TT" =>   4 , "TW" =>  -2 , "TY" =>  -2 , "TV" =>   0 ,
 "WA" =>  -3 , "WR" =>  -3 , "WN" =>  -4 , "WD" =>  -4 , "WC" =>  -2 , "WQ" =>  -2 , "WE" =>  -3 , "WG" =>  -2 , "WH" =>  -2 , "WI" =>  -2 , "WL" =>  -2 , "WK" =>  -3 , "WM" =>  -1 , "WF" =>   1 , "WP" =>  -4 , "WS" =>  -3 , "WT" =>  -2 , "WW" =>  10 , "WY" =>   2 , "WV" =>  -3 ,
 "YA" =>  -2 , "YR" =>  -2 , "YN" =>  -2 , "YD" =>  -3 , "YC" =>  -2 , "YQ" =>  -1 , "YE" =>  -2 , "YG" =>  -3 , "YH" =>   2 , "YI" =>  -1 , "YL" =>  -1 , "YK" =>  -2 , "YM" =>  -1 , "YF" =>   3 , "YP" =>  -3 , "YS" =>  -2 , "YT" =>  -2 , "YW" =>   2 , "YY" =>   6 , "YV" =>  -1 ,
 "VA" =>   0 , "VR" =>  -2 , "VN" =>  -3 , "VD" =>  -3 , "VC" =>  -1 , "VQ" =>  -2 , "VE" =>  -2 , "VG" =>  -3 , "VH" =>  -3 , "VI" =>   3 , "VL" =>   1 , "VK" =>  -2 , "VM" =>   1 , "VF" =>  -1 , "VP" =>  -2 , "VS" =>  -2 , "VT" =>   0 , "VW" =>  -3 , "VY" =>  -1 , "VV" =>   4 ,
);

$xchar = "x";

$log2 = log(2.0);

%ntfreq = (
"T" => 0.28,
"C" => 0.22,
"A" => 0.28,
"G" => 0.22,
);

%sntmat = (
 "TT" =>   4 , "TC" =>  -1 , "TA" =>  -4 , "TG" =>  -4 ,
 "CT" =>  -1 , "CC" =>   4 , "CA" =>  -4 , "CG" =>  -4 ,
 "AT" =>  -4 , "AC" =>  -4 , "AA" =>   4 , "AG" =>  -1 ,
 "GT" =>  -4 , "GC" =>  -4 , "GA" =>  -1 , "GG" =>   4 ,
);

############################################################
#	Global variables
############################################################
%aafrel; %expect;

############################################################
#	Instructions etc
############################################################
$Instructions = <<EOINPUT;

$Path/$Pnam $Pver, $Pdat

Use: $Pnam srfile [options]

Options:

	-nt		nucleotides (default aa)

	-gncut=n	no more than n gaps (default $gncut, <0 count from all)

	-gnbar=n	no less than n gaps (default $gnbar, <0 count from all)

	-grcut=r	no more than r fraction of gaps (default $grcut)

	-grbar=r	no less than r fraction of gaps (default $grbar)

	-iccut=i	no less than i bits of info in a column (default $iccut)

	-icbar=i	no more than i bits of info in a column (default $icbar)

	-hocut=h	homogeneity no less than h (default $hocut)

	-hobar=h	homogeneity no more than h (default $hobar)

	-pse=count	set pseudocount to count (default $pseudq * numseq)

	-psq=quot	use quot * numseq for pseudocount (default $pseudq)

	-aaonly		ignore gap/ambiguity for homogeneity (default use)

	-grswe=r	maximum fraction of gaps for weights (default $grswe, <0 if no weighting)

	-profile	print profile in seqcolumns

	-profake	print profile in seqcolumns using 1st sequence as consensus

	-homogen	print homogeneity values (0-10)

	-consens	print consensus only (implies -hocut)

	-conplus	print sequences plus consensus (implies -hocut)

	-gcon=r		maximum fraction of gaps for consensus (default $grcon; implies -con*)

	-hcon=h		minimum homogeneity for consensus (default $hocon; implies -con*)

	-ncon=x		consensus sequence name (default $nacon; implies -con*)
EOINPUT

############################################################
#	code start
############################################################

#---	get and process arguments --------------------------
my_args(@ARGV);
-t STDIN and !@myGlobList and print $Instructions and exit 0;
push @myGlobList,"-" unless(-t STDIN);

$gncut = $myOptList{"gncut"}+0 if($myOptList{"gncut"} ne "");
$gnbar = $myOptList{"gnbar"}+0 if($myOptList{"gnbar"} ne "");
$grcut = $myOptList{"grcut"}+0 if($myOptList{"grcut"} ne "");
$grbar = $myOptList{"grbar"}+0 if($myOptList{"grbar"} ne "");
$iccut = $myOptList{"iccut"}+0 if($myOptList{"iccut"} ne "");
$icbar = $myOptList{"icbar"}+0 if($myOptList{"icbar"} ne "");
$hocut = $myOptList{"hocut"}+0 if($myOptList{"hocut"} ne "");
$hobar = $myOptList{"hobar"}+0 if($myOptList{"hobar"} ne "");
$pseudo = $myOptList{"pse"}+0 if($myOptList{"pse"} ne "");
$pseudq = $myOptList{"psq"}+0 if($myOptList{"psq"} ne "");
$grswe = $myOptList{"grswe"}+0 if($myOptList{"grswe"} ne "");
$prof = 1 if(exists $myOptList{"profile"});
$prxx = 1 if(exists $myOptList{"profake"});
$hpro = 1 if(exists $myOptList{"homogen"});
$cons = 2 if(exists $myOptList{"consens"});
$cons = 1 if(exists $myOptList{"conplus"});
$igno = 1 if(exists $myOptList{"aaonly"});
$grcon = $myOptList{"gcon"} if(exists $myOptList{"gcon"});
$hocon = $myOptList{"hcon"} if(exists $myOptList{"hcon"});
$nacon = $myOptList{"ncon"} if($myOptList{"ncon"} ne "");

if($prxx>0){			# fake profile, no filtering, no consensus
 $prof = 1;
 $cons = 0;
 $gncut = 0x7fffffff;
 $gnbar = 0;
 $grcut = 1.0;
 $grbar = 0.0;
 $iccut = -100.0;
 $icbar = 100.0;
 $hocut = -100.0;
 $hobar = 100.0;
}

if(exists $myOptList{"nt"}){
 %aafreq = %ntfreq;
 %simmat = %sntmat;
 $xchar = "n";
}

$hocut = 0 if(($cons>0 or $hpro>0) and $hocut<0);

$grswe = -1 if($myOptList{"gncut"} ne "" or $myOptList{"gnbar"} ne "");

if($iccut>=0.0){
 foreach my $aa (keys %aafreq){ $aafrel{$aa} = log($aafreq{$aa});}
}
if($hocut>=0.0){
 foreach my $aa (keys %aafreq){
  $expect{$aa} = 0;
  foreach my $bb (keys %aafreq){ $expect{$aa} += $simmat{$aa.$bb}*$aafreq{$bb};}
 }
}

#---	rest of the code -----------------------------------
foreach my $name (@myGlobList){
#---	read alignment -------------------------------------
 my @seq = (); my @sid = ();
 my $len = -1;
 	#printf STDERR "Reading [%s]\n",$name;
 open HAND,"<$name" or die "Can't read \"$name\"";
 while(<HAND>){
  chomp;
  my ($id,$sq) = split/\s+/;
  my $ll = length $sq;
  next unless($ll);
  if($len<0){
   $len = $ll;
  }else{
   die "Unequal sequence lengths at \"$id\": $ll vs $len" if($ll!=$len);
  }
  push @sid,$id;
  my @ss = split//,(uc $sq);
  for(my $i=0;$i<$len;$i++){ $seq[$i] .= $ss[$i];}
 }
 close HAND;
 
 unless(@seq>0 and @sid>0){ print STDERR "No sequences in \"$name\"\n"; next;}

 $gncut = @sid + $gncut if($gncut<0);
 $gnbar = @sid + $gnbar if($gnbar<0);
 
#---	compute alignment statistics -----------------------
 my @pro = (); my @swe = ();
 ali_profile(\@seq,\@pro,\@swe);
 
#---	print and filter -----------------------------------
 ali_print(\@seq,\@sid,\@pro,\@swe);
} 

############################################################
#	ali_profile($rseq,$rpro,$rswe)
############################################################
#	$rseq	- seqcolumns!
sub ali_profile
{
 my $rseq = shift;
 my $rpro = shift;
 my $rswe = shift;

 my $alen = @$rseq;
 my $nseq = length $$rseq[0];
 
 weight_sequences($rseq,$rswe,$alen,$nseq);
 
 	#for(my $i=0;$i<$nseq;$i++){
 	# printf "%d\t%.4f\n",$i+1,$seqw[$i];
 	#}
 my $pse = 0;						# pseudocounts
 if($pseudo>0){ $pse = $pseudo;}else{ $pse = $pseudq*$nseq;}
 
 my $lns = log($nseq+$pse);				# log effective number of sequences

 for(my $i=0;$i<$alen;$i++){				# by position
  my @seq = split //,$$rseq[$i];
  my %pro = ();
  my %fre = ();
  my $inf = 0;
  my $hom = 0;
  my $nga = 0;							# compute gaps
  for(my $j=0;$j<$nseq;$j++){ $nga += $$rswe[$j] if($seq[$j] eq "-");}
  my $caa = "-";
  if($iccut>=0 or $hocut>=0){					# compute frequencies
   foreach my $aa (keys %aafreq){ $fre{$aa} = $pse*$aafreq{$aa};}
   for(my $j=0;$j<$nseq;$j++){
    $ch = $seq[$j];
    if($aafreq{$ch}>0){ $fre{$ch} += $$rswe[$j];}
    elsif($igno==0){ foreach my $aa (keys %aafreq){ $fre{$aa} += $aafreq{$aa}*$$rswe[$j];}}
   }
  }
  	#printf "Pos %d\tG= %.3f\n",$i+1,$nga;
  	#for(my $j=0;$j<$nseq;$j++){ printf "%s\n",$seq[$j];}
  	#foreach my $aa (sort keys %aafreq){
  	# printf "\t%s\t%.4f\n",$aa,$fre{$aa};
  	#}
  	#print "\n";
  if($iccut>=0){						# compute information
   foreach my $aa (keys %fre){ $inf += $fre{$aa}*(log($fre{$aa}) - $lns - $aafrel{$aa}) if($fre{$aa}>0);}
   $inf /= ($nseq+$pse)*$log2;
  }
  if($hocut>=0){						# compute homogeneity & consensus
   my $bestaa = "";my $bestwt = -100.0*($nseq+$pse);
   my $pcnt = $nseq;
   if($igno>0){
    $pcnt = 0;
    foreach my $aa (keys %aafreq){ $pcnt += $fre{$aa};}
    $pcnt = $nseq if($pcnt<=0);
    	#printf "%d\t%d\t%.3f\n",$i,$nseq,$pcnt;
    	#exit;
   }
   foreach my $aa (keys %aafreq){
    my $wt = 0;
    foreach my $bb (keys %fre){ $wt += $simmat{$aa.$bb}*$fre{$bb};}
    if($wt>$bestwt){ $bestwt = $wt;$bestaa = $aa;}
   }
   $hom = ($bestwt/($pcnt+$pse) - $expect{$bestaa})/($simmat{$bestaa.$bestaa} - $expect{$bestaa}) if($bestaa ne "" and $bestwt/($nseq+$pse)>$expect{$bestaa});
   if($nga<=$nseq*$grcon){
    if($hom>=$hocon){ $caa = $bestaa;}else{ $caa = $xchar;}
   }
  }
  $pro{"g"} = $nga;
  $pro{"i"} = $inf;
  $pro{"h"} = $hom;
  $pro{"c"} = $caa;
  $$rpro[$i] = \%pro;
 }
}

############################################################
#	weight_sequences($rseq,$rswe,$alen,$nseq);
############################################################
#	$rseq	- seqcolumns!
#
# 	Implements:
# Henikoff S, Henikoff JG.
# Position-based sequence weights.
# J Mol Biol. 243, 574-578, 1994.
# PMID: 7966282
# 	Uses:
# $grswe - maximum fraction of gaps to include a column into calculation (<0 if no weighting)
# %aafreq - valid vs. invalid amino acids
sub weight_sequences
{
 my $rseq = shift;
 my $rswe = shift;
 my $alen = shift;
 my $nseq = shift;

 if($grswe<0){
  for(my $i=0;$i<$nseq;$i++){ $$rswe[$i] = 1;}
  return;
 }else{
  for(my $i=0;$i<$nseq;$i++){ $$rswe[$i] = 1/$nseq;}
 }
 
 foreach my $sq (@$rseq){
  my %aaw = (); my $naa = 0; my $npo = 0;
  my @aco = split//,$sq;
  for(my $i=0;$i<@aco;$i++){
   my $aa = $aco[$i];
   next unless($aafreq{$aa});
   $npo++; $aaw{$aa}++;
  }
  $naa = scalar keys %aaw;
  next unless($naa>0);
  next unless($npo>=$nseq*$grswe);
  foreach my $aa (keys %aaw){
   $aaw{$aa} = 1/$naa/$aaw{$aa};
  }
  for(my $i=0;$i<@aco;$i++){
   my $aa = $aco[$i];
   $$rswe[$i] += $aaw{$aa};
  }
 }
 
 my $wwt = 0;
 for(my $i=0;$i<$nseq;$i++){ $wwt += $$rswe[$i];}
 $wwt /= $nseq; $wwt = 1 if($wwt<=0);
 
 for(my $i=0;$i<$nseq;$i++){ $$rswe[$i] /= $wwt;}
}

############################################################
#	ali_print($rseq,$rsid,$rpro,$rswe)
############################################################
#	$rseq	- seqcolumns!
sub ali_print
{
 my $rseq = shift;
 my $rsid = shift;
 my $rpro = shift;
 my $rswe = shift;

 my $alen = @$rseq;
 my $nseq = length $$rseq[0];
 my $ncon = ($cons>0)?1:0;

 if($prof){					# print seqcolumns with data
  printf ">$nacon\t# 1.0\n" if($cons>0);
  for(my $j=0;$j<$nseq;$j++){ printf ">%s\t# %.4f\n",$$rsid[$j],$$rswe[$j];}
  for(my $i=0;$i<$alen;$i++){
   my $rppp = $$rpro[$i];
   next if($$rppp{"g"}>$gncut+$EPSILON or $$rppp{"g"}<$gnbar-$EPSILON);
   next if($$rppp{"g"}>$grcut*$nseq+$EPSILON or $$rppp{"g"}<$grbar*$nseq-$EPSILON);
   next if($$rppp{"i"}<$iccut-$EPSILON or $$rppp{"i"}>$icbar+$EPSILON);
   next if($$rppp{"h"}<$hocut-$EPSILON or $$rppp{"h"}>$hobar+$EPSILON);
   if($prxx>0){
    $$rppp{"c"} = substr($$rseq[$i],0,1);
    if($$rppp{"c"} eq "-"){
     $$rppp{"g"} = $nseq; $$rppp{"h"} = 0;
    }else{
     $$rppp{"g"} = 0; $$rppp{"h"} = 1;
    }
   }
   printf "%s",$$rppp{"c"} if($cons>0);
   printf "%s\t%d\t%.3f\t%.3f\t%.3f\t%s\n",$$rseq[$i],$nseq+$ncon,$$rppp{"g"},$$rppp{"i"},$$rppp{"h"},$$rppp{"c"};
  }
  return;
 }
						# print seqrows without data
 if($cons>0){						# print consensus
  my $sout = "";
  for(my $i=0;$i<$alen;$i++){
   my $rppp = $$rpro[$i];
   next if($$rppp{"g"}>$gncut+$EPSILON or $$rppp{"g"}<$gnbar-$EPSILON);
   next if($$rppp{"g"}>$grcut*$nseq+$EPSILON or $$rppp{"g"}<$grbar*$nseq-$EPSILON);
   next if($$rppp{"i"}<$iccut-$EPSILON or $$rppp{"i"}>$icbar+$EPSILON);
   next if($$rppp{"h"}<$hocut-$EPSILON or $$rppp{"h"}>$hobar+$EPSILON);
   $sout .= $$rppp{"c"};
  }
  print "$nacon\t$sout\n";
 }
 if($cons<2){
  for(my $j=0;$j<$nseq;$j++){				# print the rest
   my $sout = "";
   for(my $i=0;$i<$alen;$i++){
    my $rppp = $$rpro[$i];
    next if($$rppp{"g"}>$gncut+$EPSILON or $$rppp{"g"}<$gnbar-$EPSILON);
    next if($$rppp{"g"}>$grcut*$nseq+$EPSILON or $$rppp{"g"}<$grbar*$nseq-$EPSILON);
    next if($$rppp{"i"}<$iccut-$EPSILON or $$rppp{"i"}>$icbar+$EPSILON);
    next if($$rppp{"h"}<$hocut-$EPSILON or $$rppp{"h"}>$hobar+$EPSILON);
    $sout .= substr $$rseq[$i],$j,1;
   }
   printf "%s\t%s\n",$$rsid[$j],$sout;
  }
 }
 if($hpro>0){						# print homogeneity
  my $sout = "";
  for(my $i=0;$i<$alen;$i++){
   my $rppp = $$rpro[$i];
   next if($$rppp{"g"}>$gncut+$EPSILON or $$rppp{"g"}<$gnbar-$EPSILON);
   next if($$rppp{"g"}>$grcut*$nseq+$EPSILON or $$rppp{"g"}<$grbar*$nseq-$EPSILON);
   next if($$rppp{"i"}<$iccut-$EPSILON or $$rppp{"i"}>$icbar+$EPSILON);
   next if($$rppp{"h"}<$hocut-$EPSILON or $$rppp{"h"}>$hobar+$EPSILON);
   my $hval = int($$rppp{"h"}*10+0.5);
   my $hcha = chr(ord(0)+$hval);
   $sout .= ($$rppp{"c"} eq "-")?"-":$hcha;
  }
  print "$nacon|H\t$sout\n";
 }
}

############################################################
#	max min
############################################################
sub max
{
 my $max = shift;
 while(@_){
  my $x = shift;
  $max = $x if($x>$max);
 }
 return $max;
}

sub min
{
 my $min = shift;
 while(@_){
  my $x = shift;
  $min = $x if($x<$min);
 }
 return $min;
}


############################################################
#	my_args
############################################################
sub my_args
{
 my $nop = scalar @_;
 for($i=0;$i<$nop;$i++){
  if($_[$i]=~m/^-[^=]+=$/){
   my $new = $_[$i].$_[$i+1];
   splice @_,$i,2,$new;
  }
 }
 foreach my $arg (@_){
  my ($opt) = ($arg =~ m/^-(\w+)/);
  my ($val) = ($arg =~ m/^-\w+=(.*)/);
  if($opt){
   $myOptList{$opt} = $val;
  }else{
   push @myGlobList,$arg;
  }
 }
}
