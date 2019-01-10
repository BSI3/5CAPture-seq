#!/usr/bin/perl -w

use strict;
use Storable;

    my $file = $ARGV[0];  #name of sample ./uploads/Test_27/27-1-PR-allLab.fasta
    open(INFO, $file);		# Open the file
    my @seq = <INFO>;		# Read it into an array
    close(INFO);                # Close the file
    

    my $allSeq;
    my $length_seq;
    my $name;
    my $s;
    my %read =();
    my @chars = ();
    my $position;

    foreach $s (@seq)
    {
        chomp($s);
        if($s =~ />/)
        {
             if ($allSeq)
             {
                
                @chars = split //, $allSeq;
                $position = 0;
                foreach (@chars)
                {
                    $position++;
                    $read{$name}{$position} = $_;
                }
            
                $allSeq = "";   
            }
       
            my $tmpname =  substr $s, 1;
            my @words = split / /, $tmpname;
            $name = $words[0];
        }
        else
        {
            $allSeq = $allSeq.$s;
        }
    }
     
    if ($allSeq)
    {
        
        @chars = split //, $allSeq;
        $position = 0;
        foreach (@chars)
        {
            chomp($_);
            $position++;
            $read{$name}{$position} = $_;
        }
        $allSeq = "";
    }
    
  store \%read, 'file_ref';
# $hashref = retrieve('file');

    # while ((my $key, my $value) = each(%read))
    # {

    #     print "$key\t";
    #     my $length_seq = length ($value);
    #     print "$length_seq\n";

    # }
    
    
    
    
    
    
    
    
    
    
    
