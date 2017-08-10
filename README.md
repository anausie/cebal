#### preProcessed454.pl

Requirements:
-sffinfo
-sfffile 
-seqclean

This script combines the usage of these tools in order to preprocess the 454 sequencing data which is provided in SFF format. The sffinfo and sfffile tools used the features represented by upper and lowercase of the sequences of nucleotides, associated with the 454 consideration of which is of good or bad quality. This script, trims the sequence adaptor, converts all the nucleotides to uppercase in order to set all the nucleotides as a good quality, to further execute seqclean. The final output is further preprocessed by the Mothur tool in order to perform a controlled quality trimming.
