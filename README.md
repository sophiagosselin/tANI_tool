# tANI_Matrix
This perl script allows one to reconstruct the results of https://doi.org/10.1101/2020.01.15.908137.

You will need to have the latest version of blast+ installed.* Additionally you will need a version of Perl that supports Perl threads. 

To run, include the script in the same working directory as the genomes you wish to make a phylogeny of. Genome files should be in fasta format.

Furthermore, one may use genomes with any degree of genome completion; however, low quality assemblies with contigs smaller than 1020 nt will result in large loss of meaningfull phylogenetic information. Therefore one should always strive to use genome assemblies with the lowest number of these small contigs as possible. In general, the higher the level of completion, the better, and always be critical of phylogenies built from lower quality assemblies.


The script includes a number of features...
NOTE: default ID, CV, and EVAL values are set to match those of Gosselin et al. 2020.


-H -> provides a list of the avilable input options with a small amount of information.

-ID -> Indentity cutoff value. This can be on a 0-1 or 0-100 scale. This will control the lowest possible percent identity value that a hit within the blast output may have without being filtered out. The default value is 0.7.

-CV -> The coverage cutoff value. This can be on a 0-1 or 0-100 scale. This controls the lowest possible coverage value that a hit within the blast output may have without being filtered out. The default value is 0.7.

-BT -> This controls the number of bootstraps to be made. The default value is 0.

-EVAL -> This allows you to adjust the e-value cutoff for blastn. The default value is 1E-4

-TASK -> Allows one to change blastn's task if you so desire. Refer to the blast manual for specific options. Input should be entered as it would be for blast if the -task flag was invoked (i.e. blastn, megablast, etc).


*However if this is not an option, there are commented lines in the FormatDatabase and BlastGenome subroutines which can be uncommented into use if needed.
