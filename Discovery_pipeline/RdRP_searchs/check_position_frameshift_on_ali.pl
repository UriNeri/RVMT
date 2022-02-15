#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $cmd='';
my $out='';
GetOptions ('help' => \$h, 'h' => \$h);
if ($h==1 || $ARGV[0] eq ""){ # If asked for help or did not set up any argument
	print "# Script to get the position of the frameshift on an alignment
# Arguments :
# toto
";
	die "\n";
}

my $in_faa_ali="Clean_frameshifts_seq_mafft_auto.faa";
my $c_c="";
my %store;
open my $faa,"<",$in_faa_ali;
while(<$faa>){
      chomp($_);
      if ($_=~/^>(\S+)/){$c_c=$1;}
      else{$store{$c_c}.=$_;}
}
close $faa;

my $out_file="Position_frameshift_on_ali.tsv";
open my $s1,">",$out_file;
print $s1 "Sequence\tPosition_on_sequence\tPosition_in_ali\n";
foreach my $seq (sort keys %store){
      my @t=split(/\|/,$seq);
      print "looking at $t[0] = $t[1]\n";
      my $seq_c=$store{$seq};
      my @tab_residues=split("",$seq_c);
      my $n_total=0;
      my $n_nongap=0;
      my $tag=0;
      foreach my $pos (@tab_residues){
            $n_total++;
            if ($pos eq "-"){}
            else{$n_nongap++;}
            if ($n_nongap == $t[1]){
                  print "found $n_nongap -> $n_total\n";
                  print $s1 $seq."\t".$n_nongap."\t".$n_total."\n";
                  $tag=1;
                  last;
            }
      }
      if ($tag==0){
            print "### PBLM WITH THIS SEQUENCE\n";
            <STDIN>;
      }
}
close $s1;
