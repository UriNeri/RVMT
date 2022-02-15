package YIW::basic;

use strict;
#use warnings;

BEGIN {
	require Exporter;
# set the version for version checking
	our $VERSION = 1.00;
# Inherit from Exporter to export functions and variables
	our @ISA = qw(Exporter);
# Functions and variables which are exported by default
	our @EXPORT = qw(%myOptList @myGlobList $myCmdLine max min maxr minr int_commify log_this_point);
# Functions and variables which can be optionally exported
	our @EXPORT_OK = qw();
}

our %myOptList;
our @myGlobList;
our $myCmdLine;

return 1;

############################################################
#	int_commify($inum)
#	max(@arra) min(@arra)
#	maxr($rarr) minr($rarr)
#	log_this_point($pnam,$ppid,$mess)
#	YIW::basic::my_args($flag,@ARGV)
#	YIW::basic::tempdir_make_safe($dtag)
############################################################

############################################################
#	int_commify($inum)
############################################################
# @$rpro is a sorted numerical array; last element >=max({$xxxx})
# returns index, such that $$rpro[$i]<=$xxxx
sub int_commify
{
 my $inum = shift;
 1 while($inum =~s/^(-?\d+)(\d{3})/$1,$2/);
 return $inum;
}

############################################################
#	max() min()
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
#	maxr($rarr) minr($rarr)
############################################################
sub maxr
{
 my $rarr = shift;
 my $max = $$rarr[0];
 for(my $i=1;$i<@$rarr;$i++){
  $max = $$rarr[$i] if($$rarr[$i]>$max);
 }
 return $max;
}

sub minr
{
 my $rarr = shift;
 my $min = $$rarr[0];
 for(my $i=1;$i<@$rarr;$i++){
  $min = $$rarr[$i] if($$rarr[$i]<$min);
 }
 return $min;
}

############################################################
#	log_this_point($flog,$prog,$mess)
############################################################
sub log_this_point
{
 my $flog = shift;
 my $prog = shift;
 my $mess = shift;

 open HANZ,">>log.$flog.txt" or print STDERR "$prog:\tCan't append \"log.$flog.txt\"\n" and return;
 printf HANZ "%s\t%s\t%s\n",(scalar localtime()),$prog,$mess;
 close HANZ;
}

############################################################
#	my_args($rarg,$flag)
############################################################
sub my_args
{
 my $rarg = shift;
 my $flag = shift;

 my $cmdl = join "\t",@$rarg;
 
 $myCmdLine = $cmdl; $myCmdLine =~ tr/\t/ /;
 
 $cmdl =~ s/=\t+/=/g;
 
 foreach my $arg (split/\t/,$cmdl){
  my ($opt) = ($arg =~ m/^-(\w+)/);
  my ($val) = ($arg =~ m/^-\S+=(.*)/);
  if($opt ne ""){
   $myOptList{$opt} = $val;
  }else{
   push @myGlobList,$arg;
  }
 }
 push @myGlobList,"-" if($flag and not -t STDIN);
}

############################################################
#	tempdir_make_safe($dtag)
############################################################
# makes directory d.dtag.nnnnn/
sub tempdir_make_safe
{
 my $dtag = shift;
 if($dtag eq ""){
  ($dtag) = ($0 =~ m/([^\/]+)$/);
 }
 
 my $dd = "";
 
 while(1){
  my $r1 = 1 + int(rand(0x7ffffffd));
  $dd = join ".",("d",$dtag,$r1);
  unless(-e $dd){
   mkdir $dd or die "Can't make \"$dd\"";
   last;
  }
 }
 return $dd;
}
