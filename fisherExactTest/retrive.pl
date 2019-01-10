#!/usr/bin/perl -w

use strict;
use Storable;
    


    my $file = $ARGV[0];  #sigposition
    open(INFO, $file);      # Open the file
    my @seq = <INFO>;       # Read it into an array
    close(INFO);                # Close the file
    

my %reference = %{retrieve('file_ref')};
my %ref = ();


foreach my $p (@seq)
    {
        chomp ($p);
        my @tmp = split /\t/, $p;
        my $chr = $tmp[0];
        my $pos = $tmp[1];
        my $strand = $tmp[2];
        my $base = $reference{$chr}{$pos};
        #print "$chr\t$pos\t$strand\t$base\n";
        $ref{$chr}{$pos}{$strand} = $base;

    }

    store \%ref, 'ref_hash';

# foreach my $chromosome (sort keys %reference) 
# {
#     foreach my $position (sort keys %{ $reference{$chromosome} }) 
#     {
#         print "$name, $subject: $grades{$name}{$subject}\n";
#     }
# }

 
    
    
    
    
    
    
    
    
    
    
