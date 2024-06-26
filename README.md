# tANI_tool version 1.4.0

#### As of 04/03/2023 the new version of the code (tANI_tool.pl) is functional. I would highly reccomend to use it instead of the original legacy versions. If you do need or desire one of the legacy packages, they are still available (see tANI_original for the version used in the original manuscript, or tANI_low_mem, which was released alongside it). 

## General Information

Pairwise whole genome comparison via total average nucleotid identity (tANI). Non-parametric bootstrap capabilities included.
Please cite
	
	Improving Phylogenies Based on Average Nucleotide Identity, Incorporating Saturation Correction and Nonparametric Bootstrap Support
	
	Sophia Gosselin, Matthew S Fullmer, Yutian Feng, Johann Peter Gogarten
	
	DOI: https://doi.org/10.1093/sysbio/syab060


tANI_tool.pl takes a set of whole genomes (or incomplete assemblies) and allows a user to compute genome-genome distance values per: https://doi.org/10.1093/sysbio/syab060. The script will output matrices for the tANI (total average nucleotide identity) metric as per the paper, as well as three other whole genome metrics (AF (alignment fraction), jANI (jspecies-ANI), and a modified gANI (whole genome-ANI)). It will do so for the original genomes, and if requested, use a non parametric bootstrap approach to create bootstrapped matrices for these metrics.


Note that the gANI and AF metrics do not apply to only ORF's as per the original implementation (see https://doi.org/10.1093/nar/gkv657), but instead are applied to the whole genome broken into 1020nt fragments. 


## Important Considerations


One may use genomes with any degree of genome completion; however, low quality assemblies with contigs smaller than 1020 nt will result in large losses of meaningfull phylogenetic information. Therefore one should always strive to use genome assemblies with the lowest number of these small contigs as possible. In general, the higher the level of completion, the better, and always be critical of phylogenies built from lower quality assemblies. For the purposes of the original manuscript, anything below 80% completion was discarded, but better results are obtained when using a more stringent cutoff (90% and up).


Be warry of comparisons between genomes with wildly different sizes (e.g. a 2MB genome vs a 6MB genome). Such differences will lead to potentially inflated AF results, and hence an inflated tANI distance between these taxa.


If the output matrix reports a value of exactly 13 between two taxa, then there were no BLAST hits between the two genomes (for that specific query-subject pair) that passed one of the following cutoffs: coverage, percent identity, or e-value. If this is true for a single query across all subjects then that genome may simply be too divergent from the rest of your samples to get an accurate estimation of distance. 

You may also see values of 0 between sufficiently similar taxa within the bootstrapped outputs. These two genomes are not necessarily identical. A value of 0 could indicate that the randomly selected sample of fragments from the query genome were sufficiently similar to the subject genome AND that the randomly selected fragements collectively comprised a length longer than the original query genome. A tANI calculation of this type of sample would result in a negative tANI score; hence we instead output a 0 value. Check the original non-bootstrapped matrix for an accurate assesement of the tANI between those two genomes.  


## Usage and Help Text.


To run, include the script in the same working directory as the genomes you wish to compute pair-wise comparisons for. Your genome files should be in fasta format, and have one of the following extentions: .fna, .fasta, .contig, .contigs. Be aware that input files within the home directory will be edited to remove special characters; however, the original unedited inputs can be found in "intermediates/unchanged_inputs" after initial setup.

#### IMPORTANT: tANI tool has a checkpointing system. If your run is interupted simply rerun your original command in the starting directory, and the code will backup from the logs file therin.

Note that the checkpointing system can cause issues if you attempt to rerun the script from a directory where tANI ran to completion (even if there were errors along the way). So if you do encounter a bug please try to rerun the program after removing the files created by tANI_tool.pl.

### Dependencies:

	perl v5.36.0 and up
	
	BLAST v2.11.0 and up

### Usage:

	perl tANI_tool.pl -id your_%_here -cv your_%_here -boot bootstrap_#_here


Optional Inputs:

	[id]: Percent identity cutoff for inclusion of BLAST hit in tANI calculation. Default: .7

	[cv]: Percent coverage cutoff for inclusion of BLAST hit in tANI calculation. Default: .7

	[e]: Evalue cutoff for inclusion. Default: 1e-4

	[task]: Setting BLAST uses for its search criteria (see -task in BLAST).

	[boot OR bt]: Number of non-parametric tANI bootstraps. Default: 0

	[v]: Verbosity level. 1 for key checkpoints only. 2 for all messages. Default: 0

	[log OR l]: Name of file to print logs to. If none is provided program prints messages to screen only. Default: None

	[t]: Thread count. Default will use half of available cores.
	


#### NOTE: Default ID, CV, and E values are set to match those of Gosselin et al. 2020.
