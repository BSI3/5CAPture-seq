## Determining the Start Site of 5'Capped Transcripts
Maintainer: chadapohn.chaosrikul@gmail.com

### Analysis of Dinucleotide Frequencies
CAP protein contains specific dinucleotide frequencies and it locates at the start site of transcripts. We'd like to test whether our lab protocal able to capture the CAP protein. So we performed Fisher Exact Test to compare 5'capped and random dinucleotide positions.

The main files are:
* `dint_freq/`
  - `random_nt.ipynb`
  - `fisher_exact_test.ipynb`

### Exploratory of Data Analysis - Transcriptomic Factors
We'd like to find out whether the occupancies of transcriptomic factors including Nucleosome, Elongating RNA Polymerase II, Ribosome can be used to determine the transcriptions start site

The main file is:
* `occupancy_plots.ipynb`


## Fisher Exact Test
Maintainer: pavita.kae@biotec.or.th
The main files are:
* `fisherExactTest/`
  - `fasta.pl`
  - `retrive.pl`
  - `sam_firstbase_onread.pl`
  - `retrive_mismatchbeforesig.pl`
  - `fishertest.pl`

