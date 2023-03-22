#!/usr/bin/perl -w
use warnings;
use strict;
use Thread::Semaphore;
use threads;
use threads::shared;
use Getopt::Long;
use Scalar::Util;
use File::Copy;

#globals
my @input_genomes = glob "*.fna *.fasta *.contig *.contigs";
my $evalue = "1E-4";
my $task = "";
my $bootnum = my $verbosity = 0;
my ($identity_cutoff,$coverage_cutoff,$help,$thread_limit);
#is shared for multithreading
my (%matrix,%all_lengths,%bootmatrix,$semaphore):shared;

#get_inputs
GetOptions('t=i' => \$thread_limit, 'v=i' =>\$verbosity, 'task=s' =>\$task, 'id=s' => \$identity_cutoff, 'ev=s' => \$evalue, 'cv=s' => \$coverage_cutoff, 'boot=s' => \$bootnum, 'help+' => \$help, 'h+' => \$help);

#check for help call
if($help==1){
	die
	"\ntANI tool v1.2.5 Updated from tANI_low_mem.pl\n
	Pairwise whole genome comparison via total average nucleotid identity (tANI). Non-parametric bootstrap capabilities included.\n
	Please cite\: \"Improving Phylogenies Based on Average Nucleotide Identity, Incorporating Saturation Correction and Nonparametric Bootstrap Support\"\n
	Sophia Gosselin, Matthew S Fullmer, Yutian Feng, Johann Peter Gogarten\n
	DOI: https\:\/\/doi\.org\/10\.1093\/sysbio\/syab060\n

	Input genomes should be placed in the directory you are executing this command from.
	Only files with the extensions: fasta, fna, contig, or contigs, will be recognized.
	Furthermore ensure these sequence files are in FASTA format.
	The original input files can be found in intermediates/unchanged_inputs after initial setup is complete.

	Please report any problems to the GitHub: https\:\/\/github.com\/sophiagosselin\/tANI_Matrix
	OR to sophia.gosselin\@uconn.edu


	Usage: perl tANI.pl -id percent ID cutoff -cv coverage cutoff -boot bootstrap #

	IMPORTANT: tANI tool has a checkpointing system.
	 If your run is interupted simply rerun your original command in the starting directory.

	Required Inputs:
	[id]: Percent identity cutoff for inclusion of BLAST hit in tANI calculation. Default: .7
	[cv]: Percent coverage cutoff for inclusion of BLAST hit in tANI calculation. Default: .7

	Optional Inputs:
	[e]: Evalue cutoff for inclusion. Default: 1e-4
	[task]: Setting BLAST uses for its search criteria (see -task in BLAST).
	[boot]: Number of non-parametric tANI bootstraps. Default: 0
	[v]: Verbosity level. 1 for key checkpoints only. 2 for all messages. Default: 0
	[t]: Thread count. Default will use half of available cores.\n";
}

VERBOSEPRINT(1, "Initializing\n");
my($split_genome_ref) = SETUP(@input_genomes);
MAIN($split_genome_ref,@input_genomes);


sub SETUP{
	my @input_files = @_;

	#check ID and CV inputs, reformat if needed.
	$identity_cutoff = SETUP_CV_ID($identity_cutoff,"No identity threshold specified. Defaulting to 0.7\n");
	$coverage_cutoff = SETUP_CV_ID($coverage_cutoff,"No coverage threshold specified. Defaulting to 0.7\n")

	#check for directory presence, creates if not already
	DIRECTORY_CHECK("intermediates","intermediates/splits","intermediates/calc_backup","intermediates/blast_output","intermediates/unchanged_inputs","intermediates/blastdb","Outputs","Outputs/AF","Outputs/Distance","Outputs/gANI""Outputs/jANI",)

	#finds core count if no thread limit was specified
	if($thread_limit = ""){
		$thread_limit = ((my $cores = CORE_COUNT())/2);
		VERBOSEPRINT(1,"No thread limit specified. Using half ($thread_limit) of number of cores ($cores) as thread limit.\n");
	}
	#sets up thread limit
	$semaphore = Thread::Semaphore->new($thread_limit);

	#check if backup logs exist. Load into memory if present.
	my($split_genome_backup_ref) = RECOVER_SETUP_INFORMATION("setup.log");
	my(%split_genomes_backup) = %{$split_genome_backup_ref};

	#removes any already prepped genomes from the SETUP GENOME subroutine
	foreach my $backup (keys %split_genomes_backup){
		my ($previously_ran_input) = ($backup=~/intermediates\/splits\/(.*?)\.split/);
		my $index = 0;
		$index ++ until $input_files[$index] eq $previously_ran_input;
		splice(@input_files, $index, 1);
	}

	#sends input files to be prepared with the setup genomes subroutine. Reutns array of fragmented genomes
	my(@array_of_hashes) = @{(my($split_genome_ref)=THREAD_MANAGER("SETUP_GENOME", \@input_files))};
	my(%split_genomes_and_lengths)=MERGE_HASHES(@array_of_hashes);

	#merge backup with new calcs if nessecary.
	%split_genomes_and_lengths = MERGE_HASHES(\%split_genomes_backup,\%split_genomes_and_lengths);

	return(\%split_genomes_and_lengths);
}

sub SETUP_GENOME{
	#takes 1 genome as input. Returns split genome location and whole genome length as hash
	my $input_file = shift;
	my %return_info;
	#checks if file has already been processed via backup subroutine
	next if($split_genomes_backup{$input_file});
	#prepare genome files for downstream use and check R/W privelages
	FILE_I_O_CHECK($input_file);
	STANDARDIZE_FASTA($input_file);

	#split genomes into fragments. Return length and fragment count
	my($input_split,$genome_length)=SPLIT_FASTA($input_file);
	$return_info{$input_split}=$genome_length;

	#make BLAST database from input file
	MAKE_BLAST_DATABASE($input_file,"intermediates/blastdb/$input_file");

	#in case of crash, saves related genome information
	PRINT_TO_FILE("$input_file\t$genome_length\n","setup.log");

	#for return
	return(\%return_info);
}

sub SETUP_CV_ID{
	my $variable = shift;
	my $error_message = shift;
	if($variable = ""){
		VERBOSEPRINT(1,$error_message);
		$variable = 0.7;
	}
	else{
		if($variable > 1){
			$variable=$variable/100;
		}
	}
	return($variable);
}

sub MAIN{
	my %split_genomes_and_lengths = %{shift}:shared;
	my @input_files = @_;
	my (@split_genomes,@names_for_output);

	#creates input array from hash
	foreach my $input (keys %split_genomes_and_lengths){
		push(@split_genomes,$input);
	}

	#comprehensive all vs. all BLAST searches and database creation
	my(@blast_outputs) = @{(my($blast_outputs_reference)=THREAD_MANAGER("BLAST_GENOMES",\@split_genomes,\@input_files))};

	#clear memory
	$blast_outputs_reference = "";

	#load fragment lengths for all genomes into memory
	my(@array_of_hashes) = @{(my($fragment_l_ref)=THREAD_MANAGER("GET_ACESSION_LENGTHS",\@split_genomes))};
	my(%fragment_lengths)=MERGE_HASHES(@array_of_hashes):shared;

	#load contig lengths for all genomes into memory
	my(@array_of_hashes2) = @{(my($conting_l_ref)=THREAD_MANAGER("GET_ACESSION_LENGTHS",\@input_files))};
	my(%contig_lengths)=MERGE_HASHES(@array_of_hashes2):shared;

	#clear memory
	(@array_of_hashes = @array_of_hashes2 = $fragment_l_ref = $conting_l_ref = ());

	#Calculating distance, ANI, AF, and gANI
	my(@calculations)=@{(my($calc_reference)=THREAD_MANAGER("CALCULATE_METRICS",\@blast_outputs))};

	#clear memory
	($calc_reference)=();

	#creates array of genome names without file extensions
	#and pushes result for self v self comp to @calculations for output
	foreach my $entry (@input_files){
		my($genome_name)=($entry=~/(.*?)\..*/);
		push(@names_for_output,$genome_name);
		push(@calculations,"$genome_name\&$genome_name\t0\t0\t0\t13");
	}

	#send calculations to output files
	OUTPUT("original",\@calculations,\@names_for_output);

	#if bootstrapping is needed, begin here
	if($bootnum eq 0){
		VERBOSEPRINT(0,"No bootstraped matrices requested. Terminating run.\n");
		return();
	}
	else{
		my @best_hit_logs = glob "intermediates/calc_backup/*.log";
		for (my $boot_counter=0; $boot_counter<$bootnum; $boot_counter++){
			my(@bootstrap_calculations)=@{(my $calc_reference)=THREAD_MANAGER("BOOTSTRAP",\@best_hit_logs))};
			#need to add a method for pushing self-self to array
			OUTPUT($boot_counter,\@bootstrap_calculations,\@names_for_output);
		}
	}
}

sub THREAD_MANAGER{
	#takes a subroutine, the array of values to thread over and any values specific to subroutine as inputs
	#creates new threads up to the thread limit
	#returns an array of outputs
	#IMPORTANT!!!!
	#every sub called by this function must end with $semaphore->up;
	my $subroutine = shift;
	my @thread_input = @{shift};
	my @sub_specific_arguments = @_;
	my @threads;
	foreach my $threadable_input (@thread_input){
		# request a thread slot, waiting if none are available:
		$semaphore->down();
		my ($thread) = threads->create({'context' => 'list'}, \&$subroutine,$threadable_input,@sub_specific_arguments);
		push(@threads,$thread);
	}
	my @return_values = $_->join() for @threads;
	return(\@return_values);
}

sub BLAST_GENOMES {
	#takes one genome, and an array of all genomes as input
	#BLASTS the first genome against all other genomes. Returns list of BLAST output files
	#NOTES: Needs a mechanism to backup and check if a search has already been completed.
	my $query_genome = shift;
	my @database_genomes = @{shift};
	my ($query_file_handle) = ($query_genome=~/.*?\/(.*?)\..*/);
	my @return_files;
	foreach my $database (@database_genomes){
		next if($database eq $query_genome);
		my ($database_file_handle) = ($database=~/(.*?)\..*/);
		system("blastn -db intermediates/blastdb/$database -query $query_genome -evalue $evalue -outfmt 6 -out 'intermediates/blast_output/$query_file_handle\&$database_file_handle\.blast' $task");
		push(@return_files,"intermediates/blast_output/$query_file_handle\&$database_file_handle\.blast");
	}
	$semaphore->up;
	return(@return_files);
}

sub GET_ACESSION_LENGTHS{
	#takes a genome file as input
	#returns list of acessions and their associated nucleotide length
	my $input_file = shift;
	my $nucl_length = 0;
	my $acession ="";
	my %nucl_lengths
	open(IN, "< $input_file");
	while(<IN>){
		chomp;
		if($_=~/\>/){
			if(!$acession = ""){
				$nucl_lengths{$acession}=$nucl_length;
			}
			$acession = $_;
			$nucl_length = 0;
		}
		else{
			$nucl_length += length($_);
		}
	}
	$nucl_lengths{$acession}=$nucl_length;
	close IN;
	return(\%nucl_lengths);
}

sub CALCULATE_METRICS{
	#takes a blast output file as input
	#calculates jANI, gANI, AF, and tANI for the genome - genome comparison
	#saves information to backup files, and returns calculated values as a tab seperated string.
	my $blast_input = shift;
	my (%best_hits,$shorter_gene);
	my $gANI_numerator = my $gANI_denominator = my $jANI_numerator = 0;

	#get query and database genome names
	my ($query_handle,$database_handle) = ($blast_input =~/(.*?)\&(.*?)\..*/);

	open(IN, "< $blast_input");
	while(<IN>){
		chomp;
		my @blast_results = split(/\t/,@_);
		#skips match if a better hit was already found for this a given query
		next if($best_hits{$blast_results[0]});
		#skips match if it does not meet identity cutoff
		next if($blast_results[2] < $identity_cutoff);
		#determines the "shorter gene" (the fragment, or the contig it is matching to). Most likely to be fragment.
		if ($contig_lengths{$blast_results[1]} <= $fragment_lengths{$blast_results[0]}){
			$shorter_gene = $contig_lengths{$blast_results[1]};
		}
		else{
			$shorter_gene = $fragment_lengths{$blast_results[0]};
		}
		#skips match if it does not meet the coverage cutoff.
		next if (($blast_results[3]/$shorter_gene) < $coverage_cutoff);
		#stores match to memory to prevent matching the same query fragment again, and for backup downstream
		$best_hits{$blast_results[0]} = "$blast_results[2]\t$blast_results[3]\n";
		#stores values for whole genome calculations
		$gANI_numerator += ($blast_results[3]*($blast_results[2]/100));
		$gANI_denominator += $shorter_gene;
		$jANI_numerator += $Blast_Lines[2];
		#backs up best hits for bootstrapping
		PRINT_TO_FILE("$blast_results[0]\t$blast_results[2]\t$blast_results[3]\n","intermediates/calc_backup/$blast_input.log");
	}
	close IN;

	#calculation time
	my $total_matches = keys %best_hits;
	my $jANI = ($jANI_numerator/$total_matches);
	my $gANI = ($gANI_numerator/$gANI_denominator);
	my $AF = ($gANI_denominator/$split_genomes_and_lengths{"intermediates/splits/$query_handle.split"});
	my $tANI = -log($gANINumerator/$split_genomes_and_lengths{"intermediates/splits/$query_handle.split"});

	#return calculations with the following format:
	#Query_Genome&Database_Genome	jANI	gANI	AF	tANI
	my $return_info = "$query_handle\&$database_handle\t$jANI\t$gANI\t$AF\t$tANI";

	#for backup
	PRINT_TO_FILE("$return_info\n","original_calculations.log");

	return($return_info);
}

sub BOOTSTRAP{
	#takes a best hit log file as input
	#choses at random with replacement a set of best hits (it does a non-parametric bootstrap)
	#returns the bootstrapped metrics
	my $best_hits_file = shift;
	my (%best_hits,%random_sample);

	#reads in log, gets number of best hits, and data
	open(IN, "< $best_hits_file");
	while(<IN>){
		my @best_hit_data = split(/\t/,@_);
		$best_hits{$best_hits_data[0]}="$best_hits_data[1]\t$best_bits_data[2]";
	}

	#takes a random sample with replacement of best hits
	my @best_hits_query_names = keys %best_hits;
	my $number_of_hits = keys %best_hits;
	for(my $rng_counter=0; $rng_counter<$number_of_hits; $rng_counter++){
		$random_sample{$rng_counter} = $best_hits_query_names[rand @best_hits_query_names];
	}

	#calculates metrics from the random sample
	

}

sub OUTPUT{
	#takes an inidcator string to append to the end of file names before the extension
	#and an array of the four calculations to print out (jANI,gANI,AF,tANI).
	#each entry of the array must be formated as: Query_Genome&Database_Genome	jANI	gANI	AF	tANI
	#and an array of headers
	#prints resulting matrices to files.
	my $string_to_append = shift;
	my @unhashed_information = @{shift};
	my @names = @{shift};
	my (%jANI,%gANI,%AF,%tANI,@names);

	#prep data for printing
	foreach my $entry (@unhashed_information){
		my @split = split(/\t/,@_);
		$jANI{$split[0]}=$split[1];
		$gANI{$split[0]}=$split[2];
		$AF{$split[0]}=$split[3];
		$tANI{$split[0]}=$split[4];
	}
	#clear memory
	@unhashed_information=();

	#create matrices for printing
	my($jANI_matrix)=MATRIX_FROM_HASH(\%jANI,@names);
	my($gANI_matrix)=MATRIX_FROM_HASH(\%gANI,@names);
	my($AF_matrix)=MATRIX_FROM_HASH(\%AF,@names);
	my($tANI_matrix)=MATRIX_FROM_HASH(\%tANI,@names);

	#print matrices to file
	PRINT_TO_FILE($jANI_matrix,"jANI_$string_to_append.matrix");
	PRINT_TO_FILE($gANI_matrix,"gANI_$string_to_append.matrix");
	PRINT_TO_FILE($AF_matrix,"AF_$string_to_append.matrix");
	PRINT_TO_FILE($tANI_matrix,"tANI_$string_to_append.matrix");
}

sub MATRIX_FROM_HASH{
	#takes hash of data to convert into a 2D matrix, as well as an array of headers
	#returns the 2D matrix as a string
	my $hashref = shift;
	my %matrix_data = %{$hashref};
	my @matrix_header = @_;
	my $matrix_string ="";

	#create matrix header
	foreach my $header (sort @matrix_header){
		$matrix_string.="$header\t";
	}

	#convert hash to matrix
	my $current_query_name = "";
	foreach my $entry (sort keys %matrix_data){
		my($query_name)=($entry=~/(.*?)\&.*/);
		if(!$query_name eq $current_query_name){
			$matrix_string.="\n$query_name\t";
			$current_query_name = $query_name;
		}
		$matrix_string.="$matrix_data{$entry}\t";
	}
	return($matrix_string);
}


sub DIRECTORY_CHECK{
	#checks if directories exists. If not, the sub creates it.
	foreach my $directory (@_){
		unless(-d $directory){
			mkdir($directory);
		}
	}
}

sub FILE_I_O_CHECK{
	#checks a given file path for existance, and R/W privelages
	my ($path) = shift;
	if(!-e $path){
		die VERBOSEPRINT(0,"PATH: $path does not exist. Check file name for errors. If no name is present, make sure sequences are in run directory.\nAlso make sure to run in a directory that has not previously completed a search.\n");
	}
	my($read_privelage,$write_privelage) = (-r $path, -w _);
	if(defined($read_privelage & $write_privelage)){
		VERBOSEPRINT(2,"All privelages present for $path.\n");
	}
	else{
		die VERBOSEPRINT(0,"Missing privelages for $path. Check read and write privelages. Read $read_privelage, Write $write_privelage.\n");
	}
}

sub VERBOSEPRINT{
	(my $verblevel, my $message) = @_;
	if($verblevel <= $verbosity){
		print "$message\n";
	}
}

sub STANDARDIZE_FASTA {
	#removes most unique characters from annotation lines and simplifies the line
	#makes later searches and moving of files much easier.
	#additionally moves original file into a directory while leaving changed file behind
	my $fastafile = shift;
	open(IN, "< $fastafile");
	my ($annotation) = ($fastafile=~/(.*?)\..*/);
	$annotation=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?]/\_/g;
	open(OUT, "+> temp.fasta");
	while(<IN>){
		if($_=~/\>/){
			print OUT "\>$annotation\n";
		}
		else{
			print OUT $_;
		}
	}
	close IN;
	close OUT;
	move($fastafile,"intermediates/unchanged_inputs/");
	rename "temp.fasta", $fastafile;
}

sub SPLIT_FASTA{
	#splits input genome into  seperately annotated 1020nt long fragments
	#prints these fragments to a new file and returns that file name
	my $fastafile = shift;
	my $fragment_count = 0;
	my $whole_genome_length = 0;
	my($fragment_asc) = ($fastafile=~/(.*?)\..*/);
	open(IN, "< $fastafile");
	open(OUT, "+> intermediates/splits/$fastafile.split")
	while(<IN>){
		chomp;
		if($_=~/\>/){
			next;
		}
		else{
			my @sequence_fragments = ($_ =~ /(.{1,1020})/g);
			$whole_genome_length += lenght($_);

			#remove any fragment under 100nt in size (field standard)
			if($sequence_fragments[$sequence_fragments] < 100){
				#check this code. Not sure what I was thinking, and whether this actually works as intended or not!
				pop(@sequence_fragments);
			}

			foreach my $fragment (@sequence_fragments){
				print OUT ">$fragment_asc\n$fragment\n";
				$fragment_count++;
			}

		}
	}
	close IN;
	close OUT;
	return("intermediates/splits/$fastafile.split",$whole_genome_length);
}

sub MAKE_BLAST_DATABASE{
	#takes a fasta file, and an output location as input
	#creates a BLAST searchable database from input
	my $fasta_file = shift;
	my $output_location = shift;
	system("makeblastdb -dbtype nucl -in $fastafile -out $output_location");
}

sub PRINT_TO_FILE{
	#takes text and file names as input
	#prints text to file
	my $text_to_print = shift;
	my $file_handle = shift;
	open(BACKUP, ">> $file_handle");
	print BACKUP "$text_to_print";
	close BACKUP;
}

sub RECOVER_SETUP_INFORMATION{
	#recovers split genome locations and whole genome lengths.
	VERBOSEPRINT(1, "Recovering setup data from previous run.\n");
	my $infile = shift;
	my %genome_setup_backup;
	open(BACKUP, "< $infile");
	while(<BACKUP>){
		chomp;
		my @recovered = split(/\t/,@_);
		$genome_setup_backup{$recovered[0]}=$recovered[1];
	}
	close BACKUP;
	return(\%genome_setup_backup);
}

sub CORE_COUNT{
	use Sys::Info;
	my $count = Sys::CPU::cpu_count;
	return($count);
}

sub MERGE_HASHES{
	my @hashes_to_merge = @_;
	my %new_hash;
	foreach my $hashref (@hashes_to_merge){
		$new_hash = {%$new_hash,%$hashref};
	}
	return(%new_hash);
}
