# tANI_tool

## As of 03/30/2023 the new version of the code should be functional; however the backup functionality is still broken. If you need to restart a run, do so in a clearn directory. See tANI_tool.pl Legacy versions are located in a seperate directory. 

These perl scripts allow one to reconstruct the results of https://doi.org/10.1101/2020.01.15.908137.

I would suggest that you use low_mem which has a number of imporvements. Prinicpally, it decreases the amount of memory needed to run and provides more user friendly directory structures. This will increase the storage space needed while running an analysis, but not by much. Additionally the low-mem version provides additional checkpointing processes to help reduce loss of time due to outside factors.


You will need to have the latest version of blast+ installed.* Additionally you will need a version of Perl that supports Perl threads. 
*However if this is not an option, there are commented lines in the FormatDatabase and BlastGenome subroutines which can be uncommented into use if needed.

To run, include the script in the same working directory as the genomes you wish to make a phylogeny of. Genome files should be in fasta format.

Furthermore, one may use genomes with any degree of genome completion; however, low quality assemblies with contigs smaller than 1020 nt will result in large loss of meaningfull phylogenetic information. Therefore one should always strive to use genome assemblies with the lowest number of these small contigs as possible. In general, the higher the level of completion, the better, and always be critical of phylogenies built from lower quality assemblies.

NOTE: default ID, CV, and EVAL values are set to match those of Gosselin et al. 2020.

Usage: perl tANI.pl -id percent ID cutoff -cv coverage cutoff -boot bootstrap #

Required Inputs:
	[id]: Percent identity cutoff for inclusion of BLAST hit in tANI calculation. Default: .7
	[cv]: Percent coverage cutoff for inclusion of BLAST hit in tANI calculation. Default: .7

	Optional Inputs:
	[e]: Evalue cutoff for inclusion. Default: 1e-4
	[task]: Setting BLAST uses for its search criteria (see -task in BLAST).
	[boot OR bt]: Number of non-parametric tANI bootstraps. Default: 0
	[v]: Verbosity level. 1 for key checkpoints only. 2 for all messages. Default: 0
	[t]: Thread count. Default will use half of available cores.\n";
