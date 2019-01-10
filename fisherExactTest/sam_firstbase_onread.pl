#!/usr/bin/perl -w

use strict;
use Storable;

    my $file = $ARGV[0];  #name of sample ./uploads/Test_27/27-1-PR-allLab.fasta
    open(INFO, $file);		# Open the file
    my @seq = <INFO>;		# Read it into an array
    close(INFO);                # Close the file
    
    my %sam =();
    my @chars = ();
    my $firstbase;
    my $length;
    my $strand;
    my %refgenome = %{retrieve('file_ref')}; #my $base = $reference{$chr}{$pos};

    foreach my $s (@seq)
    {
        chomp($s);
       	my $letter = substr($s, 0, 1);
       	if ($letter ne "@") 
       	{
       		my @tmp = split /\t/, $s;
       		my $name = $tmp[0];
       		my $flag = $tmp[1];
          my $sequence = $tmp[9];
          my $chro = $tmp[2];
          my $position = $tmp[3];
          my $base;

   		    if ($flag & 16) 
   		    { 
				    $strand = "-";
            $length = length($sequence);
            $position = $position+$length-1;
    
            my $revcomp = reverse $sequence;
            $revcomp =~ tr/ATGCatgc/TACGtacg/;
           # $sequence = $revcomp;
            $firstbase = substr($revcomp,0,1);
            $base = $refgenome{$chro}{$position};
            $base =~ tr/ATGCatgc/TACGtacg/;
            if ($base ne $firstbase) #count mismatch
            {
              if (exists $sam{$chro}{$position}{$strand}{$firstbase}) 
                {
                  $sam{$chro}{$position}{$strand}{$firstbase} = $sam{$chro}{$position}{$strand}{$firstbase}+1;
                }
                else
                {
                  $sam{$chro}{$position}{$strand}{$firstbase} = 1;
                }
            }
				  } 
          else
          {
            $strand = "+";            
            $firstbase = substr($sequence, 0, 1);
            $base = $refgenome{$chro}{$position};
            if ($base ne $firstbase) #count mismatch
            {
                if (exists $sam{$chro}{$position}{$strand}{$firstbase}) 
                {
                  $sam{$chro}{$position}{$strand}{$firstbase} = $sam{$chro}{$position}{$strand}{$firstbase}+1;
                }
                else
                {
                  $sam{$chro}{$position}{$strand}{$firstbase} = 1;
                }
            }

          }
     
    
       		


       	}
    }
     
 	store \%sam, "$file"."_hash_rev_first";
# $hashref = retrieve('file');


# foreach my $chromosome (sort keys %sam) 
# {
#     foreach my $position (sort keys %{ $sam{$chromosome} }) 
#     {
#         foreach my $strand (sort keys %{ $sam{$chromosome}{$position} }) 
#         {
#             foreach my $base (sort keys %{ $sam{$chromosome}{$position}{$strand} }) 
#             {
#                 print "$chromosome\t$position\t$strand\t$base\t$sam{$chromosome}{$position}{$strand}{$base}\n";

#             }
#         }
#     }
# }
    
    

    
    
    
    
    
    
    
