#!/usr/bin/env perl
use strict;
use warnings;
my $p=1;
my $h=0;
while(<>){
  if(/```\{r inv-\d+,.*/){
    $p=0;next
  }
  if(/```\{r head-\d+,*/){
    $h=1;next
  }
  if($h==1){
    my $line = $_; 
    $line =~ /\s*h(\d)\(\"(.*)\"\)/;
    my $n = $1; my $t = $2;
    print "# $t\n" if ($n =~ /1/);
    print "## $t\n" if ($n =~ /2/);
    print "### $t\n" if ($n =~ /3/);
    print "#### $t\n" if ($n =~ /4/);
    print "##### $t\n" if ($n =~ /5/);
    print "###### $t\n" if ($n =~ /6/);
    $h=0;$p=0;next
  }
  if(/```/ && $p == 0){
      $p=1;next
  }
  print if($p)
}
