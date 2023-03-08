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
my ($identity_cutoff,$coverage_cutoff,$help,$thread_limit,$semaphore);
#is shared for multithreading
my (%matrix,%all_lengths,%bootmatrix):shared;

#get_inputs
GetOptions('t=i' => \$thread_limit, 'v=i' =>\$verbosity, 'task=s' =>\$task, 'id=s' => \$identity_cutoff, 'ev=s' => \$evalue, 'cv=s' => \$coverage_cutoff, 'boot=s' => \$bootnum, 'help+' => \$help, 'h+' => \$help);

#check for help call
if($help==1){
	die
	"\ntANI tool v1.2.1 Updated from tANI_low_mem.pl\n
	Pairwise whole genome comparison via total average nucleotid identity (tANI). Non-parametric bootstrap capabilities included.\n
	Please cite\: \"Improving Phylogenies Based on Average Nucleotide Identity, Incorporating Saturation Correction and Nonparametric Bootstrap Support\"\n
	Sophia Gosselin, Matthew S Fullmer, Yutian Feng, Johann Peter Gogarten\n
	DOI: https\:\/\/doi\.org\/10\.1093\/sysbio\/syab060\n

	Input genomes should be placed in the directory you are executing this command from.
	Only files with the extensions: fasta, fna, contig, or contigs, will be recognized.
	Furthermore ensure these sequence files are in FASTA format.

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

&VERBOSEPRINT(1, "Initializing\n");

my($fragment_counts_ref,$genome_length_ref,$split_genome_ref) = SETUP(@input_genomes);
MAIN($fragment_counts_ref,$genome_length_ref,$split_genome_ref,@input_genomes);


sub SETUP{
	my @input_files = @_;
	my @split_genomes;

	#check ID inputs
	if($identity_cutoff = ""){
		VERBOSEPRINT(1,"No identity threshold specified. Defaulting to 0.7\n")
		$identity_cutoff = 0.7;
	}
	else{
		if($identity_cutoff > 1){
			$identity_cutoff=$identity_cutoff/100;
		}
	}

	#check CV input
	if($coverage_cutoff = ""){
		VERBOSEPRINT(1,"No coverage threshold specified. Defaulting to 0.7\n")
		$coverage_cutoff = 0.7;
	}
	else{
		if($coverage_cutoff > 1){
			$coverage_cutoff=$coverage_cutoff/100;
		}
	}

	#check for directory presence, creates if not already
	DIRECTORY_CHECK("intermediates","intermediates/splits","intermediates/blast_output","intermediates/blastdb","Outputs","Outputs/AF","Outputs/Distance","Outputs/gANI""Outputs/jANI",)

	#finds core count if no thread limit was specified
	if($thread_limit = ""){
		$thread_limit = ((my $cores = CORE_COUNT())/2);
		VERBOSEPRINT(1,"No thread limit specified. Using half ($thread_limit) of number of cores ($cores) as thread limit.\n");
	}
	#sets up thread limit
	$semaphore = Thread::Semaphore->new($thread_limit);

	#check if backup logs exist. Load into memory if present.
	my ($fragment_counts_hashref,$genome_lengths_hashref) = RECOVER_GENOME_INFORMATION("genome_information_backup.log");
	my %fragment_counts = %{$fragment_counts_hashref};
	my %genome_lengths = %{$genome_lengths_hashref};

	foreach my $input_file (@input_files){
		#checks if file has already been processed via backup subroutine
		next if($fragment_counts{$input_file} && $genome_lengths{$input_file});
		#prepare genome files for downstream use and check R/W privelages
		FILE_I_O_CHECK($input_file);
		STANDARDIZE_FASTA($input_file);

		#split genomes into fragments. Return length and fragment count
		my($input_genome_length,$input_fragment_count,$input_split)=SPLIT_FASTA($input_file);
		push(@split_genomes,$input_split);

		#make BLAST database from fragment file
		MAKE_BLAST_DATABASE($input_split,"intermediates/blastdb/$input_file");

		#in case of crash, saves related genome information
		BACKUP_TO_FILE("$input_file\t$input_genome_length\t$input_fragment_count\n","genome_information_backup.log");

		#for return
		$fragment_counts{$input_file}=$input_fragment_count;
		$genome_lengths{$input_file}=$input_genome_length;
	}
	return(\%fragment_counts,\%genome_lengths,\@split_genomes);
}

sub MAIN{
	my %fragments_per_genome = %{shift};
	my %lengths_of_genomes = %{shift};
	my @split_genomes = @{shift};
	my @input_files = @_;
	#comprehensive all vs. all BLAST searches and database creation
	#if this syntax works I will be gobsmacked
	my(@blast_outputs) = @{(my($blast_outputs_reference)=THREAD_MANAGER("BlastGenome",\@split_genomes,\@input_files))};

	#ANYTHING BELOW THIS POINT STILL DOES NOT WORK.
	#Calculating distance, ANI, AF, and gANI
	THREAD_MANAGER("Calculate_ANI",\@input_genomes);
	#outs to file. Sends file type
	THREAD_MANAGER("orig","orig",%matrix);
	#bootstrapping
	for (my $bc=0; $bc<$bootnum; $bc++){
		threaded("Bootstrap_ANI",\@input_genomes);
		Outfile($bc,"matrix",%bootmatrix);
		%bootmatrix=();
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
		my ($thread) = threads->create(\&$subroutine,$threadable_input,@sub_specific_arguments);
		push(@threads,$thread);
	}
	my @return_values = $_->join() for @threads;
	return(\@return_values);
}

sub BlastGenome {
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
		system("blastn -db intermediates/blastdb/$database -query $query_genome -evalue $evalue -outfmt 6 -out 'intermediates/blast_output/$query_file_handle\_$database_file_handle\.blast' $task");
		push(@return_files,"intermediates/blast_output/$query_file_handle\_$database_file_handle\.blast");
	}
	return(@return_files);
}

#All calculations are done here; name misleading
sub Calculate_ANI{
	#readin query information
	open(SPLIT, "< intermediates/splits/$_[0].split");
	my(%qsplitsize);
	while(<SPLIT>){
		if($_=~/\>/){
			my($annotation,$splitsize)=($_=~/\>(.*\_.*?\_.*?\_(.*))/);
			$qsplitsize{$annotation}=$splitsize;
		}
		else{}
	}
	close SPLIT;
	foreach my $database (@input_genomes){
		if($database eq $_[0]){
			$matrix{"$_[0]\+$_[0]"} = "1.0\t1.0\t1.0\t0.0";
			next;
		}
		else{}
		#readin database information
		open(DB, "< intermediates/splits/$database.split");
		my(%dbsplitsize,%dbcontigsize);
		while(<DB>){
			if($_=~/\>/){
				my($annotation,$consize)=($_=~/\>(.*)\_.*?\_(.*?)\_.*/);
				next if exists $dbcontigsize{$annotation};
				$dbcontigsize{$annotation}=$consize;
			}
			else{}
		}
		close DB;
		#start filtering BLAST output
		my(%BestHits,$gANINumerator,$TotalShort,$TotalID,$shortergene);
		my $HitCounter = 0;
		open(BLAST_FILE, "< intermediates/blast_output/$_[0]+$database.blast");
		while(<BLAST_FILE>){
			my @Blast_Lines = split;
			next if exists $BestHits{$Blast_Lines[0]};
			next if $Blast_Lines[2] < $identity_cutoff;
			if ($dbcontigsize{$Blast_Lines[1]} <= $qsplitsize{$Blast_Lines[0]}){
				$shortergene = $dbcontigsize{$Blast_Lines[1]};
			}
			else{
				$shortergene = $qsplitsize{$Blast_Lines[0]};
			}
			next if ($shortergene == 0);
			next if (($Blast_Lines[3]/$shortergene) < $coverage_cutoff);
			my $gANIindv = ($Blast_Lines[3]*($Blast_Lines[2]/100));
			$BestHits{$Blast_Lines[0]} = 1;
			$gANINumerator += $gANIindv;
			$TotalShort += $shortergene;
			$HitCounter ++;
			$TotalID += $Blast_Lines[2];
		}
		close BLAST_FILE;
		if($HitCounter == 0){
			$matrix{"$_[0]\+$database"} = "0.00\t0.00\t0.00\t13.00";
		}
		else{
			my $jANI = ($TotalID/$HitCounter);
			my $gANI = ($gANINumerator/$TotalShort);
			my $AF = ($TotalShort/$all_lengths{$_[0]});
			my $gDistance = -log($gANINumerator/$all_lengths{$_[0]});
			$matrix{"$_[0]\+$database"} = "$jANI\t$gANI\t$AF\t$gDistance";
		}
	}
}

sub Bootstrap_ANI{
	#readin query information
	open(SPLIT, "< intermediates/splits/$_[0].split");
	my(%qsplitsize,%randomsample,$bootsize);
	my $counter = 0;
	while(<SPLIT>){
		if($_=~/\>/){
			my($annotation,$splitsize)=($_=~/\>(.*\_.*?\_.*?\_(.*))/);
			$qsplitsize{$annotation}=$splitsize;
			$counter++;
		}
		else{}
	}
	close SPLIT;
	my @q_keys = keys %qsplitsize;
	for(my $rngcount=0; $rngcount<$counter; $rngcount++){
		$randomsample{$rngcount} = $q_keys[rand @q_keys];
	}
	foreach my $database (@input_genomes){
		if($database eq $_[0]){
			$bootmatrix{"$_[0]\+$_[0]"} = "1.0\t1.0\t1.0\t0.0";
			next;
		}
		else{}
		#readin database information
		open(DB, "< intermediates/splits/$database.split");
		my(%dbsplitsize,%dbcontigsize);
		while(<DB>){
			if($_=~/\>/){
				my($annotation,$consize)=($_=~/\>(.*)\_.*?\_(.*?)\_.*/);
				next if exists $dbcontigsize{$annotation};
				$dbcontigsize{$annotation}=$consize;
			}
			else{}
		}
		close DB;
		#start filtering BLAST output
		my(%BestHits,$gANINumerator,$TotalShort,$TotalID,$shortergene);
		my $HitCounter = 0;
		open(BLAST_FILE, "< intermediates/blast_output/$_[0]+$database.blast");
		while(<BLAST_FILE>){
			my @Blast_Lines = split;
			next if exists $BestHits{$Blast_Lines[0]};
			next if $Blast_Lines[2] < $identity_cutoff;
			if ($dbcontigsize{$Blast_Lines[1]} <= $qsplitsize{$Blast_Lines[0]}){
				$shortergene = $dbcontigsize{$Blast_Lines[1]};
			}
			else{
				$shortergene = $qsplitsize{$Blast_Lines[0]};
			}
			next if ($shortergene == 0);
			next if (($Blast_Lines[3]/$shortergene) < $coverage_cutoff);
			my $gANIindv = ($Blast_Lines[3]*($Blast_Lines[2]/100));
			$BestHits{$Blast_Lines[0]} = "$gANIindv\t$Blast_Lines[2]\t$shortergene";
		}
		close BLAST_FILE;
		foreach my $sample (values %randomsample){
			if(exists $BestHits{$sample}){
				my @split = split(/	/,$BestHits{$sample});
				$gANINumerator += $split[0];
				$TotalShort += $split[2];
				$HitCounter++;
				$TotalID += $split[1];
			}
			else{
			}
		}
		if($HitCounter == 0){
			$bootmatrix{"$_[0]\+$database"} = "0.00\t0.00\t0.00\t13.00";
		}
		else{
			my $jANI = ($TotalID/$HitCounter);
			my $gANI = ($gANINumerator/$TotalShort);
			my $AF = ($TotalShort/$all_lengths{$_[0]});
			my $gDistance = -log($gANINumerator/$all_lengths{$_[0]});
			$bootmatrix{"$_[0]\+$database"} = "$jANI\t$gANI\t$AF\t$gDistance";
		}
	}
}

sub Outfile{
	#number or indicator to go after the _ but before file extension.
	my $u1 = shift;
	#file extension
	my $u2 = shift;
	#content of the matrices
	my %toprint = @_;
	#checker tool
	my $currentline = "hi";
	my %handles;
	my @outfiles = ("Outputs/jANI/jANI","Outputs/gANI/gANI","Outputs/AF/AF","Outputs/Distance/Distance");
	foreach my $out1 (@outfiles){
		open (my $fh, "+>", "$out1\_$u1\.$u2");
		print {$fh} join("\t",sort @input_genomes);
		#try a print join approach when you have time
		$handles{$out1}=$fh;
	}
	foreach my $value (sort keys %toprint){
		my ($check) = ($value =~ /(.*?)\+.*/);
		if ($check eq $currentline){
			my $num = 0;
			my @outputs = split(/	/,$toprint{$value});
			foreach my $out (@outfiles){
				print {$handles{$out}} "\t$outputs[$num]";
				$num++;
			}
		}
		else{
			$currentline = $check;
			my $num = 0;
			my @outputs = split(/	/,$toprint{$value});
			foreach my $out (@outfiles){
				print {$handles{$out}} "\n$check\t$outputs[$num]";
				$num++;
			}
		}
	}
	foreach my $out2 (@outfiles){
		close $handles{$out2};
	}
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
	#removes most unique characters from annotation lines
	#makes later searches and moving of files much easier.
	my $fastafile = shift;
	open(IN, "< $fastafile");
	open(OUT, "+> temp.fasta");
	while(<IN>){
		if($_=~/\>/){
			$_=~s/[\ \[\]\(\)\:\;\/\.\-\~\`\!\@\#\$\%\^\&\*\=\+\{\}\?]/\_/g;
			print OUT $_;
		}
		else{
			print OUT $_;
		}
	}
	close IN;
	close OUT;
	unlink $fastafile;
	rename "temp.fasta", $fastafile;
}

sub SPLIT_FASTA{
	#splits input genome into  seperately annotated 1020nt long fragments
	#prints these fragments to a new file
	#returns total length of genome (sans any tiny fragments that were filtered outs) and fragment count
	my $fastafile = shift;
	my $genome_length = my $fragment_count = 0;
	open(IN, "< $fastafile");
	open(OUT, "+> intermediates/splits/$fastafile.split")
	while(<IN>){
		chomp;
		if($_=~/\>/){
			next;
		}
		else{
			my $sublength = length($_);
			my @sequence_fragments = ($_ =~ /(.{1,1020})/g);

			#remove any fragment under 100nt in size (field standard)
			if($sequence_fragments[$#sequence_fragments] < 100){
				my $lost_length = pop(@sequence_fragments);
				$sublength-= $lost_length;
			}

			$genome_length+= $sublength;
			foreach my $fragment (@sequence_fragments){
				print OUT ">$fragment_asc\n$fragment\n";
				$fragment_count++;
			}

		}
	}
	close IN;
	close OUT;
	return($genome_length,$fragment_count);
}

sub MAKE_BLAST_DATABASE{
	#takes a fasta file, and an output location as input
	#creates a BLAST searchable database from input
	my $fasta_file = shift;
	my $output_location = shift;
	system("makeblastdb -dbtype nucl -in $fastafile -out $output_location");
}

sub BACKUP_TO_FILE{
	#takes text and file names as input
	#prints text to file
	my $text_to_backup = shift;
	my $file_handle = shift;
	open(BACKUP, ">> $file_handle");
	print BACKUP "$text_to_backup";
	close BACKUP;
}

sub RECOVER_GENOME_INFORMATION{
	#recovers array values for searched and unsearched queries
	VERBOSEPRINT(1, "Recovering from previous run.\n");
	my $infile = shift;
	my %length_backup;
	my %fragment_backup;
	open(BACKUP, "< $infile");
	while(<BACKUP>){
		chomp;
		my @split_recover = split(/\t/,$_);
		$length_backup{$split[0]}=$split[1];
		$fragment_backup{$split[0]}=$split[2];
	}
	close BACKUP;
	return(\%fragment_backup,\%length_backup);
}

sub CORE_COUNT{
	use Sys::Info;
	my $count = Sys::CPU::cpu_count;
	return($count);
}
