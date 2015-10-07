#!/usr/bin/perl

# John M. Gaspar
# Aug. 2015

# Filter VEP annotations in a VCF: choose the most severe
#   based on IMPACT, ranked by PolyPhen and SIFT scores.
# Also, change 'AF' label to 'GQ' so that gemini will
#   load them (under 'gt_quals').

use strict;
use warnings;

sub usage {
  print q(Usage: perl VEPtoGemini.pl  <infile>  <outfile>  [option]
  Required:
    <infile>     Input VCF file containing VEP annotations
                   (including 'IMPACT', 'PolyPhen', and 'SIFT')
    <outfile>    Output VCF file
  Optional:
      -gq        Option to convert 'AF' genotype labels to 'GQ'
                   (loaded by gemini as 'gt_quals')
);
  exit;
}

usage() if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

# open files
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(OUT, ">$ARGV[1]");

# option to convert 'AF' to 'GQ'
my $gq = 0;
$gq = 1 if (scalar @ARGV > 2 && $ARGV[2] eq "-gq");

# rank of IMPACT levels
my %rank = ("MODIFIER" => 0, "LOW" => 1,
  "MODERATE" => 2, "HIGH" => 3);

# indexes of CSQ annotation fields
my $imp = -1; my $pph; my $sft;

# parse VCF
while (my $line = <IN>) {

  if (substr($line, 0, 1) eq '#') {

    # change 'AF' to 'GQ'
    if ($gq && $line =~ m/ID=AF/) {
      $line =~ s/AF/GQ/;
    }

    # load positions of fields for CSQ annotation
    if ($line =~ m/CSQ.*Format: /) {
      my @spl = split(/\|/, substr($', 0, -3));
      $imp = getIdx("IMPACT", @spl);
      $pph = getIdx("PolyPhen", @spl);
      $sft = getIdx("SIFT", @spl);
    }

    print OUT $line;
    next;
  }

  die "Error! Cannot find CSQ format header in $ARGV[0]\n"
    if ($imp == -1);

  chomp $line;
  my @spl = split("\t", $line);
  my @div = split(';', $spl[7]);
  for (my $x = 0; $x < scalar @div; $x++) {
    my @cut = split('=', $div[$x]);
    if ($cut[0] eq "CSQ") {

      my @res = ();    # best match(es)
      my @sum = ();    # best match(es): sums of PolyPhen and SIFT
      my $Mcq = -1;    # best match(es): rank of IMPACT

      # compare each annotation
      my @com = split(',', $cut[1]);
      for (my $y = 0; $y < scalar @com; $y++) {
        my @brk = split(/\|/, $com[$y]);
        if (! exists $rank{$brk[$imp]}) {
          die "Error! Unknown variant impact: $brk[$imp]\n";
        }
        my $cq = $rank{$brk[$imp]};
        my $pp = ($brk[$pph] =~ m/\(([0-9.]+)\)/ ? $1 : -2);
        my $sf = ($brk[$sft] =~ m/\(([0-9.]+)\)/ ? 1 - $1 : -2);

        # replace if IMPACT is greater
        if ($cq > $Mcq) {
          @res = ();
          $res[0] = $com[$y];
          $sum[0] = $pp + $sf;
          $Mcq = $cq;
        } elsif ($cq == $Mcq) {
          # same IMPACT: add to @res array
          my $z;
          for ($z = 0; $z < scalar @res; $z++) {
            last if ($pp + $sf > $sum[$z]);
          }
          splice(@res, $z, 0, $com[$y]);
          splice(@sum, $z, 0, $pp + $sf);
        }

      }
      die "Error! Cannot load annotation in VCF record:\n$line\n"
        if (! @res);
      $div[$x] = "CSQ=" . join(',', @res);
    }
  }

  # print output
  $spl[7] = join(';', @div);
  $spl[8] =~ s/AF/GQ/ if ($gq);
  print OUT join("\t", @spl), "\n";
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
