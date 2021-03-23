# DnoisE
## _Distance denoise by Entropy_
### _An open source parallelizable alternative to Unoise_
by Adrià Antich (CEAB-CSIC, Center of Advanced Studies of Blanes)

Here we present a new program to denoise sequence data sets from Illumina using D (distance) corrected by Entropy.
DnoisE is a denoising software that use Unoise equation algorithm to detect incorrect sequences from PCR and sequencing errors. For coding barcodes where Entropy of nucleotides is depending on the position in the codon, correction is needed to avoid removal of correct sequences that differ from spected correct ones by a change in position 3 of the codon (highly variable in nature).
DnoisE has been tested with COI barcode region in Antich et al. (2021).

Pros versus Unoise:
1 - DnoisE algorithm is parallelizable leading computational speed depending on computational hardware and a very good option if multicore computer is available. 

2- DnoisE is written all in Python3 which makes it user changeable.

3- It also accepts both .csv and .fasta inputs and can return both too. 

4- 
