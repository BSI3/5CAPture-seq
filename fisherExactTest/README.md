# 5CAPture-seq
calculate Fisher exact test
1. fasta.pl - convert reference genome to file_ref{chromosome}{position} = base
2. retrive.pl - convert 20,502 significant positions to ref_hash{chromosome}{position}{strand} = base
3. sam_firstbase_onread.pl - extract first base of all read from sam file to file.sam_hash_rev 
4. retrive_mismatchbeforesig.pl - get base of mismatch significant position
5. fishertest.pl - calculate fisher exact test of 2 groups
