#!/usr/bin/perl -w


##################################################
# input hash of real and hash of all
# get position real count mismatch and count all read + count position>1
# create 2 file - real and random
#
##################################################

use strict;
use Storable;



    my $file = $ARGV[0]; #firstbase
    my $name = $ARGV[1];

my %sigposition = %{retrieve('ref_hash')}; #from retrive.pl $ref{$chr}{$pos}{$strand} = $base;
#my %reference = %{retrieve('ref_hash')};
my %query1 = %{retrieve("$file")}; #from sam.pl $sam{$chro}{$position}{$strand}{$firstbase} = 1;
my %refgenome = %{retrieve('file_ref')}; #my $base = $reference{$chr}{$pos};


open(my $fh, '>', "$name".'_firsebase_mm_sig.txt');
 open(my $fh5, '>', "$name".'_firsebase_mm_nosig.txt');



foreach my $chromosome (sort keys %sigposition) 
{
    foreach my $position (sort keys %{ $sigposition{$chromosome} }) 
    {
        foreach my $strand (sort keys %{ $sigposition{$chromosome}{$position} }) 
        {
            my $base = $sigposition{$chromosome}{$position}{$strand};
            my $positionbeforesig = 0;
    
            
            my $numA = 0;
            my $numT = 0;
             my $numC = 0;
             my $numG = 0;
             my $numN = 0;


            

            if ($strand eq "+") 
            {
                $positionbeforesig = $position-1;

            }
            elsif($strand eq "-")
            {
                $positionbeforesig = $position+1;

            }

                $numA = $query1{$chromosome}{$positionbeforesig}{$strand}{"A"}; #count beforesig that first
                $numT = $query1{$chromosome}{$positionbeforesig}{$strand}{"T"};
                $numC = $query1{$chromosome}{$positionbeforesig}{$strand}{"C"};
                $numG = $query1{$chromosome}{$positionbeforesig}{$strand}{"G"};
                $numN = $query1{$chromosome}{$positionbeforesig}{$strand}{"N"};

             
                delete $query1{$chromosome}{$positionbeforesig}{$strand}{"A"};
                delete $query1{$chromosome}{$positionbeforesig}{$strand}{"T"};
                delete $query1{$chromosome}{$positionbeforesig}{$strand}{"C"};
                delete $query1{$chromosome}{$positionbeforesig}{$strand}{"G"};
                delete $query1{$chromosome}{$positionbeforesig}{$strand}{"N"};

                if ((!$numA)) { $numA = 0;}
                if ((!$numT)) { $numT = 0;}
                if ((!$numC)) { $numC = 0;}
                if ((!$numG)) { $numG = 0;}
                if ((!$numN)) { $numN = 0;}
            
     
                print $fh "$chromosome\t$position\t$strand\t$base\t$numA\t$numT\t$numC\t$numG\t$numN\n";

        }
        
    }
}
close $fh;

my $baseref;
 
 foreach my $chr (sort keys %query1) 
{
    foreach my $pos (sort keys %{ $query1{$chr} }) 
    {
        foreach my $str(sort keys %{ $query1{$chr}{$pos} }) 
        {
            if ($str eq "+") 
            {
                
                $baseref = $refgenome{$chr}{$pos};

            }
            elsif($str eq "-")
            {
                $baseref = $refgenome{$chr}{$pos};
                $baseref =~ tr/ATGCatgc/TACGtacg/;
            }

            my $matchA = 0;
            my $matchT = 0;
             my $matchC = 0;
             my $matchG = 0;
             my $matchN = 0;

                $matchA = $query1{$chr}{$pos}{$str}{"A"}; #first base but no sig
                $matchT = $query1{$chr}{$pos}{$str}{"T"};
                $matchC = $query1{$chr}{$pos}{$str}{"C"};
                $matchG = $query1{$chr}{$pos}{$str}{"G"};
                $matchN = $query1{$chr}{$pos}{$str}{"N"};

                if ((!$matchA)) { $matchA = 0;}
                if ((!$matchT)) { $matchT = 0;}
                if ((!$matchC)) { $matchC = 0;}
                if ((!$matchG)) { $matchG = 0;}
                if ((!$matchN)) { $matchN = 0;}

                print $fh5 "$chr\t$pos\t$str\t$baseref\t$matchA\t$matchT\t$matchC\t$matchG\t$matchN\n";
        }
    }
}

close $fh5;

    
