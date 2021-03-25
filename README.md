# DnoisE
## _Distance denoise by Entropy_
### _An open source parallelizable alternative to Unoise_
by Adrià Antich (CEAB-CSIC, Center of Advanced Studies of Blanes)

Here we present a new program to denoise sequence data sets from Illumina using D (distance) corrected by Entropy.
DnoisE is a denoising software that use Unoise equation algorithm to detect incorrect sequences from PCR and sequencing errors. For coding barcodes where Entropy of nucleotides is depending on the position in the codon, correction is needed to avoid removal of correct sequences that differ from spected correct ones by a change in position 3 of the codon (highly variable in nature).
DnoisE has been tested with COI barcode region in Antich et al. (2021).

Pros versus Unoise:

1 - DnoisE can correct the heigh of differences depending on the position of nucleotides in the codon based on Entropy values of each position.

2 - DnoisE algorithm is parallelizable leading computational speed depending on computational hardware and a very good option if multicore computer is available. 

3 - DnoisE is written in Python3 and Open access which makes it user changeable.

4 - It also accepts both .csv and .fasta inputs and can return both too. 

5 - DnoisE allows to user to choose the joining method. As derived from Edgar's equation (beta(d)=.5^(alpha\*d+1)), Unoise algorithm joins incorrect sequences to the most abundant mother sequence that satisfies the equation (ratio criteria). From our point of view this can lead to an overjoining of sequences to the most abundant one. We have developed an algorithm that returns two extra types of joining criteria outputs. The first one (d criteria) joins sequences to the most similar and most abundant sequence. The later, joins sequences to that for which the ratio abundance divided per beta(d) is the lowest (ratio_d criteria).


*INPUT FILES*

Input files can be as both .csv format and .fasta

If input file is a fasta, sequence must be in a single row and both id and size must have ";" closing.


      >Seq_000000012;size=433081;
      TTTGAGTTCAATACAAAGTCATTCAGGAGCTGCTATTGACTTAGCTATCTTCAGTTTACATCTTTCAGGAGCTTCTTCGATTCTAGGAGCAATTAATTTTATTTCTACCATTATAAATATGCGAAATCCTGGACAAACATTTTATCGCATTCCTTTATTTGTTTGATCGATTTTCGTAACTGCTTTACTACTATTATTAGCAGTACCAGTTTTAGCAGGAGCTATTACCATGTTACTAACTGATCGTAATTTTAATACAGCCTTTTTTGACCCTTCTGGAGGTGGTGATCCTGTACTTTATCAACATTTATTT

If input file is a .csv, separation between columns can be specified usind -p parameter (see help)


