# tANI_Matrix
This perl script allows one to reconstruct the results of https://doi.org/10.1101/2020.01.15.908137.

You will need to have the latest version of blast+ installed.* 

To run, one must include the script in the same working directory as the genomes one wishes to make a phylogeny of. Genome files should be in fasta format.

Furthermore, one may use genomes with any degree of genome completion; however, low quality assemblies with contigs smaller than 1020 nt will result in large loss of meaningfull phylogenetic information. Therefore one should always strive to use genome assemblies with the lowest number of these small contigs as possible. In general, the higher the level of completion, the better.

The script includes a number of features...



*However if this is not an option, there are commented lines in the FormatDatabase and BlastGenome subroutines which can be uncommented into use if needed.
