#!/usr/bin/perl -w
use warnings;
use strict;
use threads;
use threads::shared;
use Thread::Semaphore;
use Getopt::Long;
use Scalar::Util;
use File::Copy;
use Cwd;

#globals
my @input_genomes = glob "*.fna *.fasta *.contig *.contigs";
my $evalue = "1E-4";
my $task = "";
my $bootnum = my $verbosity = my $help = my $printlog = 0;
my ($identity_cutoff,$coverage_cutoff,$thread_limit);
#is shared for multithreading
my ($semaphore,%fragment_lengths,%contig_lengths,%split_genomes_and_lengths):shared;

#get_inputs
GetOptions('t=i' => \$thread_limit, 'v=i' =>\$verbosity, 'task=s' =>\$task, 'id=s' => \$identity_cutoff, 'ev=s' => \$evalue, 'cv=s' => \$coverage_cutoff, 'boot=s' => \$bootnum, 'bt=s' => \$bootnum, 'help+' => \$help, 'h+' => \$help, 'log=s' => \$printlog, 'l=s' => \$printlog);

#check for help call
if($help==1){
	die
	"\ntANI tool v1.3.0 Updated from tANI_low_mem.pl\n
	Pairwise whole genome comparison via total average nucleotid identity (tANI). Non-parametric bootstrap capabilities included.\n
	Please cite\: \"Improving Phylogenies Based on Average Nucleotide Identity, Incorporating Saturation Correction and Nonparametric Bootstrap Support\"\n
	Sophia Gosselin, Matthew S Fullmer, Yutian Feng, Johann Peter Gogarten\n
	DOI: https\:\/\/doi\.org\/10\.1093\/sysbio\/syab060\n

	Input genomes should be placed in the directory you are executing this command from.
	Only files with the extensions: fasta, fna, contig, or contigs, will be recognized.
	Furthermore ensure these sequence files are in FASTA format.
	Be aware that input files within the home directory will be edited to remove special characters;
	However, the original unedited inputs can be found in intermediates/unchanged_inputs after initial setup.

	Please report any problems to the GitHub: https\:\/\/github.com\/sophiagosselin\/tANI_Matrix
	OR to sophia.gosselin\@uconn.edu


	Usage: perl tANI_tool.pl -id percent ID cutoff -cv coverage cutoff -boot bootstrap #

	IMPORTANT: tANI tool has a checkpointing system.
	 If your run is interupted simply rerun your original command in the starting directory.

	Required Inputs:
	[id]: Percent identity cutoff for inclusion of BLAST hit in tANI calculation. Default: .7
	[cv]: Percent coverage cutoff for inclusion of BLAST hit in tANI calculation. Default: .7

	Optional Inputs:
	[e]: Evalue cutoff for inclusion. Default: 1e-4
	[task]: Setting BLAST uses for its search criteria (see -task in BLAST).
	[boot OR bt]: Number of non-parametric tANI bootstraps. Default: 0
	[v]: Verbosity level. 1 for key checkpoints only. 2 for all messages. Default: 0
	[log OR l]: Name of file to print logs to. If none is provided program prints messages to screen only. Default: None
	[t]: Thread count. Default will use half of available cores.\n";
}

#NOTES FOR FUTURE DEV.
#Replace your current hash system with the more complicated hash sub key system. It should simplify things

VERBOSEPRINT(1, "Initializing\n");

my @renamed_input_genomes = SETUP(@input_genomes);
MAIN(@renamed_input_genomes);
VERBOSEPRINT(0,"All processes complete. See output directory for results.\n");

sub SETUP{
	my @input_files = @_;
	my %processed_genomes_information_from_backup;

	#check ID and CV inputs, reformat if needed.
	$identity_cutoff = SETUP_CV_ID($identity_cutoff,"No identity threshold specified. Defaulting to 0.7\n");
	$coverage_cutoff = SETUP_CV_ID($coverage_cutoff,"No coverage threshold specified. Defaulting to 0.7\n");

	#check for directory presence, creates if not already
	DIRECTORY_CHECK("intermediates","intermediates/splits","intermediates/calc_backup","intermediates/blast_output","intermediates/unchanged_inputs","intermediates/blastdb","outputs","outputs/AF","outputs/tANI","outputs/gANI","outputs/jANI");

	#renames inputs for ease of use. Saves old versions elsewhere
	my(@renamed_input_files) = RENAME_INPUT_FILES(@input_files);
	my(@files_for_setup) = @renamed_input_files;

	#finds core count if no thread limit was specified
	if(!defined $thread_limit){
		$thread_limit = 4;
		VERBOSEPRINT(1,"No thread limit specified. Using 4 threads as limit. Please speficy if you have more cores available\n");
		#$thread_limit = ((my $cores = CORE_COUNT())/2);
	}
	#sets up thread limit
	$semaphore = Thread::Semaphore->new($thread_limit);

	#check if backup logs exist. Load into memory if present.
	if(-e "setup.log"){
		(%processed_genomes_information_from_backup)=RECOVER_INFORMATION_HASH("setup.log");

		#removes any already prepped genomes from the SETUP GENOME subroutine
		foreach my $backup (keys %processed_genomes_information_from_backup){
			my($lookup)=($backup=~/.*\/(.*?)\.split/);
			my $index = 0;
			$index ++ until $files_for_setup[$index] eq $lookup;
			splice(@files_for_setup, $index, 1);
		}
	}

	#sends input files to be prepared with the setup genomes subroutine. Returns array of fragmented genomes then dereferences it
	my($array_of_hashes_ref)=THREAD_MANAGER("SETUP_GENOME", \@files_for_setup);
	my @array_of_hashes = @{$array_of_hashes_ref};
	my(%processed_genomes_information_new) = MERGE_HASHES(@array_of_hashes);

	#merge backup with new calcs if nessecary
	if(%processed_genomes_information_from_backup){
		(%processed_genomes_information_new) = MERGE_HASHES(\%processed_genomes_information_from_backup,\%processed_genomes_information_new);
	}
	#set info to shared hash
	%split_genomes_and_lengths = %processed_genomes_information_new;
	return(@renamed_input_files);
}

sub SETUP_GENOME{
	#takes 1 genome as input. Returns split genome location and whole genome length as hash
	my $input_file = shift;
	my %return_info;
	#prepare genome files for downstream use and check R/W privelages
	FILE_I_O_CHECK($input_file);
	STANDARDIZE_FASTA($input_file);
	VERBOSEPRINT(2,"$input_file has been standardized. Moving to splitting.\n");

	#split genomes into fragments. Return length and fragment count
	my($input_split,$genome_length)=SPLIT_FASTA($input_file);
	$return_info{$input_split}=$genome_length;
	VERBOSEPRINT(2,"$input_file has been fragmented. 1020nt fragment file located at $input_split\n");

	#make BLAST database from file
	VERBOSEPRINT(2,"Creating blastn database for $input_file\n");
	MAKE_BLAST_DATABASE($input_file,"intermediates/blastdb/$input_file");

	#in case of crash, saves related genome information
	PRINT_TO_FILE("intermediates/splits/$input_file.split\t$genome_length\n","setup.log");

	#free up thread
	$semaphore->up;

	#for return
	return(\%return_info);
}

sub SETUP_CV_ID{
	my $variable = shift;
	my $error_message = shift;
	if(!defined $variable){
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
	my (@input_files) = @_;
	my (@split_genomes,@names_for_output,@blast_input_backup);

	#creates input array from hash
	foreach my $input (keys %split_genomes_and_lengths){
		push(@split_genomes,$input);
	}

	#comprehensive all vs. all BLAST searches and database creation
	#make BLAST database from input file
	if(-e "blast_database.log"){
		@blast_input_backup = RECOVER_INFORMATION_ARRAY("blast_database.log");
	}

	VERBOSEPRINT(1,"Running all vs all BLAST.\n");
	my($blast_outputs_reference)=THREAD_MANAGER("BLAST_GENOMES",\@split_genomes,\@input_files,\@blast_input_backup);
	my(@blast_outputs_refs) = @{$blast_outputs_reference};
	#dereference array
	my @blast_outputs;
	foreach my $ref (@blast_outputs_refs){
		push(@blast_outputs,@{$ref});
	}

	VERBOSEPRINT(1,"Loading fragment and genome lengths into memory.\n");
	#load fragment lengths for all genomes into memory
	my($fragment_l_ref)=THREAD_MANAGER("GET_ACESSION_LENGTHS",\@split_genomes);
	my(@array_of_hashes) = @{$fragment_l_ref};
	%fragment_lengths = MERGE_HASHES(@array_of_hashes);

	#load contig lengths for all genomes into memory
	my($conting_l_ref)=THREAD_MANAGER("GET_ACESSION_LENGTHS",\@input_files);
	my(@array_of_hashes2) = @{$conting_l_ref};
	%contig_lengths = MERGE_HASHES(@array_of_hashes2);

	VERBOSEPRINT(1,"Calculating metrics for all comparisons.\n");
	#Calculating distance, ANI, AF, and gANI
	my($calc_reference)=THREAD_MANAGER("CALCULATE_METRICS",\@blast_outputs);
	my(@calculations)=@{$calc_reference};

	#clear memory
	(@array_of_hashes = @blast_outputs_refs = @array_of_hashes2 = $fragment_l_ref = $blast_outputs_reference = $conting_l_ref = $calc_reference = ());

	#creates array of genome names without file extensions
	#and pushes result for self v self comp to @calculations for output
	VERBOSEPRINT(1,"Printing original all vs all matrix to output files.\n");
	foreach my $entry (@input_files){
		my($genome_name)=($entry=~/(.*?)\..*/);
		push(@names_for_output,$genome_name);
		push(@calculations,"$genome_name\&$genome_name\t100\t100\t1\t0");
	}

	#send calculations to output files
	OUTPUT("original",\@calculations,\@names_for_output);

	#if bootstrapping is needed, begin here
	if($bootnum eq 0){
		VERBOSEPRINT(1,"No bootstraped matrices requested. Terminating run.\n");
		return();
	}
	else{
		VERBOSEPRINT(1,"Beginning non-parametric bootstrap process.\n");
		my @best_hit_logs = glob "intermediates/calc_backup/*.log";
		for (my $boot_counter=0; $boot_counter<$bootnum; $boot_counter++){
			my $for_verbose = ($boot_counter+1);
			VERBOSEPRINT(1,"Starting bootstap $for_verbose.\n");
			my $calc_reference=THREAD_MANAGER("BOOTSTRAP",\@best_hit_logs);
			my(@bootstrap_calculations)=@{$calc_reference};
			foreach my $entry (@input_files){
				my($genome_name)=($entry=~/(.*?)\..*/);
				push(@bootstrap_calculations,"$genome_name\&$genome_name\t100\t100\t1\t0");
			}
			#need to add a method for pushing self-self to array
			VERBOSEPRINT(1,"Printing bootstap $for_verbose to file.\n");
			OUTPUT($boot_counter,\@bootstrap_calculations,\@names_for_output);
		}
		return();
	}
}

sub THREAD_MANAGER{
	#takes a subroutine, the array of values to thread over and any values specific to subroutine as inputs
	#creates new threads up to the thread limit
	#returns an array of outputs
	#IMPORTANT!!!!
	#every sub called by this function must end with $semaphore->up;
	my $subroutine = shift;
	my $array_ref = shift;
	my @thread_input = @{$array_ref};
	my @sub_specific_arguments = @_;
	my (@threads,@return_values);
	foreach my $threadable_input (@thread_input){
		# request a thread slot, waiting if none are available:
		$semaphore->down();
		my($thread)=threads->create(\&$subroutine,$threadable_input,@sub_specific_arguments);
		push(@threads,$thread);
	}
	foreach (@threads){
		my $return_value = $_->join();
		push(@return_values,$return_value);
	}
	return(\@return_values);
}

sub BLAST_GENOMES {
	#takes one genome, and an array of all genomes as input
	#BLASTS the first genome against all other genomes. Returns list of BLAST output files
	#NOTES: Needs a mechanism to backup and check if a search has already been completed.
	my $query_genome = shift;
	my $array_ref = shift;
	my $array_ref2 = shift;
	my @database_genomes = @{$array_ref};
	my @completed_searches = @{$array_ref2};
	my ($query_file_handle) = ($query_genome=~/.*\/(.*?)\..*/);
	#print "$query_file_handle\n" or die VERBOSEPRINT(0, "Could not properly retrieve file handle for $query_genome.\n");
	my @return_files;
	foreach my $database (@database_genomes){
		my ($database_file_handle) = ($database=~/(.*?)\..*/);
		next if($database_file_handle eq $query_file_handle);
		if(!grep(/^$query_file_handle\&$database_file_handle$/, @completed_searches)){
			system("blastn -db intermediates/blastdb/$database -query $query_genome -evalue $evalue -outfmt 6 -out 'intermediates/blast_output/$query_file_handle\&$database_file_handle\.blast' $task");
			PRINT_TO_FILE("$query_file_handle&$database_file_handle\n","blast_database.log");
		}
		push(@return_files,"intermediates/blast_output/$query_file_handle\&$database_file_handle\.blast");
	}
	$semaphore->up;
	return(\@return_files);
}

sub GET_ACESSION_LENGTHS{
	#takes a genome file as input
	#returns list of acessions and their associated nucleotide length
	my $input_file = shift;
	my $dna_length = 0;
	my $acession ="";
	my %nucl_lengths;
	open(my $fh, "< $input_file");
	while(<$fh>){
		chomp;
		if($_=~/\>/){
			if(!$acession eq ""){
				$nucl_lengths{$acession}=$dna_length;
			}
			($acession) = ($_=~/\>(.*)/);
			$dna_length = 0;
		}
		else{
			$dna_length += length($_);
		}
	}
	$nucl_lengths{$acession}=$dna_length;
	close $fh;
	$semaphore->up;
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
	my ($query_handle,$database_handle) = ($blast_input =~/.*\/(.*?)\&(.*?)\.blast/);

	open(my $fh, "< $blast_input");
	while(<$fh>){
		chomp;
		my @blast_results = split(/\t/,$_);
		#skips match if a better hit was already found for this a given query
		next if($best_hits{$blast_results[0]});
		#skips match if it does not meet identity cutoff
		next if(($blast_results[2]/100) < $identity_cutoff);
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
		$jANI_numerator += $blast_results[2];
		#backs up best hits for bootstrapping
		#PRINT_TO_FILE("$blast_results[0]\t$blast_results[2]\t$blast_results[3]\t$shorter_gene\n","intermediates/calc_backup/$query_handle\&$database_handle.log");
	}
	close $fh;

	#calculation time
	my $total_matches = keys %best_hits;
	my $jANI = ($jANI_numerator/$total_matches);
	my $gANI = ($gANI_numerator/$gANI_denominator);
	my $AF = ($gANI_denominator/$split_genomes_and_lengths{"intermediates/splits/$query_handle.fasta.split"});
	my $tANI = -log($gANI_numerator/$split_genomes_and_lengths{"intermediates/splits/$query_handle.fasta.split"});

	#return calculations with the following format:
	#Query_Genome&Database_Genome	jANI	gANI	AF	tANI
	my $return_info = "$query_handle\&$database_handle\t$jANI\t$gANI\t$AF\t$tANI";

	#for backup
	PRINT_TO_FILE("$return_info\n","original_calculations.log");

	$semaphore->up;
	return($return_info);
}

sub BOOTSTRAP{
	#takes a best hit log file as input
	#choses at random with replacement a set of best hits (it does a non-parametric bootstrap)
	#returns the bootstrapped metrics
	my $best_hits_file = shift;
	my (%best_hits,%random_sample);
	my $gANI_numerator = my $gANI_denominator = my $jANI_numerator = 0;

	#reads in log, gets number of best hits, and data
	open(my $fh, "< $best_hits_file");
	while(<$fh>){
		chomp;
		my @best_hits_data = split(/\t/,$_);
		$best_hits{$best_hits_data[0]}="$best_hits_data[1]\t$best_hits_data[2]\t$best_hits_data[3]";
	}
	close $fh;

	#gets handles for output
	my($query_handle,$database_handle)=($best_hits_file=~/.*\/(.*?)\&(.*?)\.log/);

	#takes a random sample with replacement of best hits
	my @best_hits_query_names = keys %best_hits;
	my $number_of_hits = keys %best_hits;
	for(my $rng_counter=0; $rng_counter<$number_of_hits; $rng_counter++){
		$random_sample{$rng_counter} = $best_hits_query_names[rand @best_hits_query_names];
	}

	#calculates metrics from the random sample
	foreach my $r_sample (values %random_sample){
		my @random_sample_data = split(/\t/,$best_hits{$r_sample});
		$gANI_numerator += ($random_sample_data[1]*($random_sample_data[0]/100));
		$gANI_denominator += $random_sample_data[2];
		$jANI_numerator += $random_sample_data[0];
	}

	my $jANI = ($jANI_numerator/$number_of_hits);
	my $gANI = ($gANI_numerator/$gANI_denominator);
	my $AF = ($gANI_denominator/$split_genomes_and_lengths{"intermediates/splits/$query_handle.fasta.split"});
	my $tANI = -log($gANI_numerator/$split_genomes_and_lengths{"intermediates/splits/$query_handle.fasta.split"});

	my $return_info = "$query_handle\&$database_handle\t$jANI\t$gANI\t$AF\t$tANI";

	$semaphore->up;
	return($return_info);
}

sub OUTPUT{
	#takes an inidcator string to append to the end of file names before the extension
	#and an array of the four calculations to print out (jANI,gANI,AF,tANI).
	#each entry of the array must be formated as: Query_Genome&Database_Genome	jANI	gANI	AF	tANI
	#and an array of headers
	#prints resulting matrices to files.
	my $string_to_append = shift;
	my $array_ref = shift;
	my $array_ref2 = shift;
	my @unhashed_information = @{$array_ref};
	my @names = @{$array_ref2};
	my (%jANI,%gANI,%AF,%tANI);

	#prep data for printing
	foreach my $entry (@unhashed_information){
		my @split = split(/\t/,$entry);
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
	PRINT_TO_FILE($jANI_matrix,"outputs/jANI/jANI_$string_to_append.matrix");
	PRINT_TO_FILE($gANI_matrix,"outputs/gANI/gANI_$string_to_append.matrix");
	PRINT_TO_FILE($AF_matrix,"outputs/AF/AF_$string_to_append.matrix");
	PRINT_TO_FILE($tANI_matrix,"outputs/tANI/tANI_$string_to_append.matrix");
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
		if($query_name ne $current_query_name){
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
	#prints a message to screen if it's verbosity level is below the users desired level.
	#additionally if a log file name has been provided, prints message to file
	(my $verblevel, my $message) = @_;
	if($verblevel <= $verbosity){
		print "$message\n";
	}
	if(!$printlog == 0){
		PRINT_TO_FILE($message,$printlog);
	}
}

sub STANDARDIZE_FASTA {
	#removes most unique characters from annotation lines and simplifies the line
	#makes later searches and moving of files much easier.
	#removes all line breaks from sequences, and gives every contig a unique ID
	#additionally moves original file into a directory while leaving changed file behind
	my $fastafile = shift;
	my $counter = 0;
	my ($annotation) = ($fastafile=~/(.*?)\..*/);
	$annotation=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?]/\_/g;
	open(my $fh, "< $fastafile");
	open(my $OUT, "+> temp.fasta.$fastafile");
	while(<$fh>){
		if($_=~/\>/){
			if($counter == 0){
			}
			else{
				print $OUT "\n";
			}
			print $OUT "\>$annotation\_contig\_$counter\n";
			$counter++;
		}
		else{
			chomp;
			print $OUT $_;
		}
	}
	close $fh;
	close $OUT;
	unlink $fastafile;
	rename("temp.fasta.$fastafile",$fastafile);
}

sub SPLIT_FASTA{
	#splits input genome into  seperately annotated 1020nt long fragments
	#prints these fragments to a new file and returns that file name
	my $fastafile = shift;
	my $fragment_count = 0;
	my $whole_genome_length = 0;
	my($fragment_asc) = ($fastafile=~/(.*?)\..*/);
	open(my $fh, "< $fastafile") or die VERBOSEPRINT(0,"Could not open $fastafile. Try providing the software with more RAM.\n If this does not help, please contact the dev.\n");
	open(my $OUT, "+> intermediates/splits/$fastafile.split");
	while(<$fh>){
		chomp;
		if($_=~/\>/){
			next;
		}
		else{
			my @sequence_fragments = ($_ =~ /(.{1,1020})/g);
			$whole_genome_length += length($_);

			#remove any fragment under 100nt in size (field standard)
			for (my $index=0; $index<$#sequence_fragments; $index++){
				if(length($sequence_fragments[$index]) < 100){
					splice(@sequence_fragments, $index, 1);
				}
			}

			foreach my $fragment (@sequence_fragments){
				print $OUT ">$fragment_asc\_fragment\_$fragment_count\n$fragment\n";
				$fragment_count++;
			}

		}
	}
	close $fh;
	close $OUT;

	return("intermediates/splits/$fastafile.split",$whole_genome_length);
}

sub MAKE_BLAST_DATABASE{
	#takes a fasta file, and an output location as input
	#creates a BLAST searchable database from input
	my $fasta_file = shift;
	my $output_location = shift;
	system("makeblastdb -dbtype nucl -in $fasta_file -out $output_location");
}

sub PRINT_TO_FILE{
	#takes text and file names as input
	#prints text to file
	my $text_to_print = shift;
	my $file_handle = shift;
	open(my $fh, ">> $file_handle");
	print $fh "$text_to_print";
	close $fh;
}

sub RECOVER_INFORMATION_HASH{
	#recovers key and value data from log file
	VERBOSEPRINT(1, "Recovering setup data from previous run.\n");
	my $infile = shift;
	my %backup;
	open(my $BACKUP, "< $infile");
	while(<$BACKUP>){
		chomp;
		my @recovered = split(/\t/,$_);
		$backup{$recovered[0]}=$recovered[1];
	}
	close $BACKUP;
	return(%backup);
}

sub RECOVER_INFORMATION_ARRAY{
	#recovers array from log file.
	VERBOSEPRINT(1, "Recovering BLAST search data from previous run.\n");
	my $infile = shift;
	my @backup;
	open(my $BACKUP, "< $infile");
	while(<$BACKUP>){
		chomp;
		my @recovered = split(/\t/,$_);
		push(@backup,$recovered[0]);
	}
	close $BACKUP;
	return(@backup);
}

#May remake in future. For now, scrapping feature
#sub CORE_COUNT{
#	use Sys::Info;
#	my $count = Sys::CPU::cpu_count;
#	return($count);
#}

sub MERGE_HASHES{
	my @hashes_to_merge = @_;
	my %new_hash;
	foreach my $hashref (@hashes_to_merge){
		my %temp_hash = %{$hashref};
		%new_hash = (%new_hash,%temp_hash);
	}
	my $test = (\%new_hash);
	return(%new_hash);
}

sub RENAME_INPUT_FILES{
	my @files_to_rename = @_;
	my (@renamed_files,$rename);
	foreach my $file (@files_to_rename){
		my $prehandle = $file;
		$prehandle=~s/[\ \[\]\(\)\:\;\/\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?]/\_/g;
		my($handle)=($prehandle=~/(.*)\..*/);
		$handle=~s/\./\_/g;
		my $new_file_name = ("$handle".".fasta");
		if($new_file_name eq $file){
			copy("$file","intermediates/unchanged_inputs/$file") or die VERBOSEPRINT(0,"Copying of files failed. Check your permissions.\n");
		}
		else{
			copy("$file","$new_file_name") or die VERBOSEPRINT(0,"Copying of files failed. Check your permissions.\n");
			copy("$file","intermediates/unchanged_inputs/$file") or die VERBOSEPRINT(0,"Copying of files failed. Check your permissions.\n");
			unlink($file);
		}
		push(@renamed_files,$new_file_name);
	}
	return @renamed_files;
}

sub TEST_HASH{
	#holdover testing tool incase needed for debugging.
	#prints contents of hash to screen as well as name of the hash provided
	my $hash_name = shift;
	my $hash_ref = shift;
	my %hash = %{$hash_ref};
	foreach my $key (keys %hash){
		print "$hash_name\t$key\t$hash{$key}\n";
	}
}
