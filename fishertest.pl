#!/usr/bin/perl -w

use strict;
#use Text::NSP::Measures::2D::Fisher::left;
use Text::NSP::Measures::2D::Fisher::right;

my $group1 = $ARGV[0]; #significant position
my $group2 = $ARGV[1]; #no significant position

my $command1 = 'awk '."'".'{for (i=5;i<=8;i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}'."'".' '.$group1;
my $command2 = 'awk '."'".'{for (i=5;i<=8;i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}'."'".' '.$group2;

my $tmp1 = qx($command1);
my $tmp2 = qx($command2);
my $errorCode;

my @group1 = split /\n/, $tmp1;
my @group2 = split /\n/, $tmp2;

foreach (@group1)
{
	print "$_\n";
}
foreach (@group2)
{
	print "$_\n";
}

#for A fisher's exact test
my $n11 = $group1[0];# count mismatch with A group1
my $np1 = $n11+$group2[0];  # mismatch with A group1 + mismatch with A group2
my $n1p = $group1[0]+$group1[1]+$group1[2]+$group1[3];  #mismatch all group1
my $npp = $n1p+$group2[0]+$group2[1]+$group2[2]+$group2[3];  #mismatch all group1 + mismatch all group2
 
my $left_value = calculateStatistic( n11=>$n11,
                                    n1p=>$n1p,
                                    np1=>$np1,
                                    npp=>$npp);


 
if( ($errorCode = getErrorCode()))
{
  print STDERR $errorCode." - ".getErrorMessage();
}
else
{
  print "$n11\t$np1\t$n1p\t$npp\tvalue for A is $left_value\n";
}

#for T fisher's exact test
$n11 = $group1[1];# count mismatch with T group1
$np1 = $n11+$group2[1];  # mismatch with T group1 + mismatch with T group2
$n1p = $group1[0]+$group1[1]+$group1[2]+$group1[3];  #mismatch all group1
$npp = $n1p+$group2[0]+$group2[1]+$group2[2]+$group2[3];  #mismatch all group1 + mismatch all group2
 
$left_value = calculateStatistic( n11=>$n11,
                                    n1p=>$n1p,
                                    np1=>$np1,
                                    npp=>$npp);


 
if( ($errorCode = getErrorCode()))
{
  print STDERR $errorCode." - ".getErrorMessage();
}
else
{
  print "$n11\t$np1\t$n1p\t$npp\tvalue for T is $left_value\n";
}

# #for C fisher's exact test
$n11 = $group1[2];# count mismatch with C group1
$np1 = $n11+$group2[2];  # mismatch with C group1 + mismatch with C group2
$n1p = $group1[0]+$group1[1]+$group1[2]+$group1[3];  #mismatch all group1
$npp = $n1p+$group2[0]+$group2[1]+$group2[2]+$group2[3];  #mismatch all group1 + mismatch all group2
 
$left_value = calculateStatistic( n11=>$n11,
                                    n1p=>$n1p,
                                    np1=>$np1,
                                    npp=>$npp);


 
if( ($errorCode = getErrorCode()))
{
  print STDERR $errorCode." - ".getErrorMessage();
}
else
{
  print "$n11\t$np1\t$n1p\t$npp\tvalue for C is $left_value\n";
}

# #for G fisher's exact test
$n11 = $group1[3];# count mismatch with G group1
$np1 = $n11+$group2[3];  # mismatch with G group1 + mismatch with G group2
$n1p = $group1[0]+$group1[1]+$group1[2]+$group1[3];  #mismatch all group1
$npp = $n1p+$group2[0]+$group2[1]+$group2[2]+$group2[3];  #mismatch all group1 + mismatch all group2
 
$left_value = calculateStatistic( n11=>$n11,
                                    n1p=>$n1p,
                                    np1=>$np1,
                                    npp=>$npp);


 
if( ($errorCode = getErrorCode()))
{
  print STDERR $errorCode." - ".getErrorMessage();
}
else
{
  print "$n11\t$np1\t$n1p\t$npp\tvalue for G is $left_value\n";
}
