#!/usr/bin/perl

# John M. Gaspar
# Oct. 2015

# After annotation by VEP using --refseq, and reordering by
#   VEPtoGemini.pl, add gene "SYMBOL" and HGVS fields to
#   primary annotation.

use strict;
use warnings;

sub usage {
  print q(Usage: perl addHGVS.pl  <infile1>  <infile2>  <outfile>
  Required:
    <infile1>    Input VCF file containing VEP annotations
                   (including SYMBOL, Codons, CDS_position,
                     Amino_acids, and Protein_position)
    <infile2>    BED file listing locations of primers
    <outfile>    Output VCF file
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

# open files
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(BED, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]");

# load gene symbols by position
my %reg;
my %pos;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      print "Warning! Skipping amplicon $spl[3] -- ",
      "located at chromosomes $spl[0] and $div[0]!?\n";
    } else {
      my $st; my $end;
      if ($spl[1] < $div[1]) {
        $st = $spl[2] + 1;
        $end = $div[1] + 1;
      } else {
        $st = $div[2] + 1;
        $end = $spl[1] + 1;
      }
      my @brk = split('_', $spl[3]);
      for (my $y = $st-1; $y < $end; $y++) {
        $reg{"$spl[0]\t$y"} = $brk[0];
      }
    }
    delete $pos{$spl[3]};
  } else {
    $pos{$spl[3]} = "$spl[0]\t$spl[1]\t$spl[2]";
  }
}
close BED;

# indexes of CSQ annotation fields
my $sym = -1; my $cod; my $cpos; my $ami; my $ppos;

# parse variants in VCF
while (my $line = <IN>) {

  if (substr($line, 0, 1) eq '#') {
    # load positions of fields for CSQ annotation
    if ($line =~ m/CSQ.*Format: /) {
      my @spl = split(/\|/, substr($', 0, -3));
      $sym = getIdx("SYMBOL", @spl);
      $cod = getIdx("Codons", @spl);
      $cpos = getIdx("CDS_position", @spl);
      $ami = getIdx("Amino_acids", @spl);
      $ppos = getIdx("Protein_position", @spl);
      $line = substr($line, 0, -3) . "|HGVSc|HGVSp\">\n";
    }
    print OUT $line;
    next;
  }

  die "Error! Cannot find CSQ format header in $ARGV[0]\n"
    if ($sym == -1);

  # check variants
  chomp $line;
  my @spl = split("\t", $line);
  my @div = split(';', $spl[7]);
  my $ins = 0; my $del = 0;
  for (my $x = 0; $x < scalar @div; $x++) {
    my @cut = split('=', $div[$x]);

    # load variant type
    if ($cut[0] eq "TYPE") {
      if ($cut[1] eq "ins") {
        $ins = 1;
      } elsif ($cut[1] eq "del") {
        $del = 1;
      }

    } elsif ($cut[0] eq "CSQ") {

      # get gene name
      my $gene = "";
      if (! exists $reg{"$spl[0]\t$spl[1]"}) {
        print "Warning! Unknown gene for variant at $spl[0], $spl[1]\n";
      }
      $gene = $reg{"$spl[0]\t$spl[1]"};

      # adjust annotations
      my @zap = split(',', $cut[1]);
      for (my $z = 0; $z < scalar @zap; $z++) {
        my @arr = split(/\|/, $zap[$z]);

        # fill in gene name
        if (scalar @arr == $sym || $arr[$sym] eq "") {
          $arr[$sym] = $gene;
        }

        # construct HGVSc
        my $hgc = "";
        if ($arr[$cod]) {
          $hgc = "c.";
          my @br1 = split(/\//, $arr[$cpos]);
          my $cp = $arr[$cod];
          $cp =~ tr/a-z//d;
          my @br2 = split(/\//, $cp);
          die "Error! Cannot parse nt change in variant:\n$line\n"
            if (scalar @br1 < 2 || (scalar @br2 < 2 && (! $del)));

          # insertion
          if ($ins) {
            $br1[0] =~ tr/-/_/;
            $hgc .= "$br1[0]ins$br2[1]";

          # deletion
          } elsif ($del) {
            $br1[0] =~ tr/-/_/;
            $hgc .= "$br1[0]del$br2[0]";

          # substitution
          } else {
            $hgc .= "$br1[0]$br2[0]>$br2[1]";
          }
        }

        # construct HGVSp
        my $hgp = "";
        my @ab1 = split(/\//, $arr[$ami]);
        if (scalar @ab1 > 1) {
          $hgp = "p.";
          my @ab2 = split(/\//, $arr[$ppos]);

          # insertion
          if ($ins) {
            my $ch = "ins";
            my @zz = split('-', $ab2[0]);
            if ($ab1[1] =~ m/X/) {
              # frameshift
              my $base = substr($ab1[0], 0, 1);
              $base = '?' if ($base eq '-');
              $hgp .= $base . "$zz[0]fs";
            } else {
              if ($ab1[0] ne '-') {
                # extra residues in annotation
                if ($ab1[0] ne substr($ab1[1], 0, length $ab1[0])) {
                  if ($ab1[0] ne substr($ab1[1], -length $ab1[0])) {
                    $ch = "delins";  # mixed deletion-insertion
                  } else {
                    $ab1[1] = substr($ab1[1], 0, -length $ab1[0]);
                  }
                } else {
                  $ab1[1] = substr($ab1[1], length $ab1[0]);
                }
              } else {
                $ab1[0] = '?';
              }
              if (scalar @zz > 1) {
                $hgp .= "$ab1[0]$zz[0]_?$zz[1]$ch$ab1[1]";
              } else {
                $hgp .= "$ab1[0]$zz[0]$ch$ab1[1]";
              }
            }

          # deletion
          } elsif ($del) {
            my $ch = "del";
            my @zz = split('-', $ab2[0]);
            if ($ab1[1] =~ m/X/) {
              # frameshift
              my $base = substr($ab1[0], 0, 1);
              $base = '?' if ($base eq '-');
              $hgp .= $base . "$zz[0]fs";
            } else {
              if ($ab1[1] ne '-') {
                # extra residues in annotation
                if ($ab1[1] ne substr($ab1[0], 0, length $ab1[1])) {
                  if ($ab1[1] ne substr($ab1[0], -length $ab1[1])) {
                    $ch = "delins$ab1[1]";  # mixed deletion-insertion
                  } else {
                    $ab1[0] = substr($ab1[0], 0, -length $ab1[1]);
                    $zz[1] -= length $ab1[1];
                  }
                } else {
                  $ab1[0] = substr($ab1[0], length $ab1[1]);
                  $zz[0] += length $ab1[1];
                }
              }
              if (length $ab1[0] > 1) {
                # multiple-aa deletion
                $hgp .= substr($ab1[0], 0, 1) . "$zz[0]_" .
                  substr($ab1[0], -1) . "$zz[1]$ch";
              } else {
                $hgp .= "$ab1[0]$zz[0]$ch";
              }
            }

          # substitution
          } else {
            $hgp .= "$ab1[0]$ab2[0]$ab1[1]";
            $hgp .= "ext*?" if ($ab1[0] eq '*');  # stop loss
          }
        }

        # adjust annotation
        $zap[$z] = join('|', @arr);
        $zap[$z] .= "|$hgc|$hgp";
      }

      # print output
      $cut[1] = join(',', @zap);
      $div[$x] = join('=', @cut);
      $spl[7] = join(';', @div);
      print OUT join("\t", @spl), "\n";

      last;
    }
  }

}
close IN;
close OUT;

# return array index for item
sub getIdx {
  my $tag = shift @_;
  my @arr = @_;
  my $x;
  for ($x = 0; $x < scalar @arr; $x++) {
    last if ($tag eq $arr[$x]);
  }
  die "Error! Cannot find $tag in CSQ field\n"
    if ($x == scalar @arr);
  return $x;
}
