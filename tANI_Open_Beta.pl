#!/usr/bin/perl -w
use warnings;
use strict;
use threads;
use threads::shared;

mkdir "Outputs";
mkdir "Outputs/AF";
mkdir "Outputs/Distance";
mkdir "Outputs/gANI";
mkdir "Outputs/jANI";
#this program uses FASTA formated genomes/contigs. This should catch most of the extensions for those and related files.
my @genomefiles = glob "*.fna *.fasta *.faa *.contig *.fa *.contigs";
my ($coverage,$identity,$bootnum);
my $evalue = "1E-4";
#variable for the blast search such that the user can switch between different default configurations
my $task = "";
my (%Bootout,%matrix,%AllFileContent,%GenomeLength,%ContigLength,%HitsforBoot,%Results,%BootResults):shared;


#Parse Input
if (!exists $ARGV[0]){
	die "\nNo inputs detected, please use Default, or type -H for assistance.\n\n";
}
elsif ($ARGV[0] eq "-H"){
	die "\nThe following options are available:\n\n\t-ID:# - The identity cutoff value (0-100).\n\t-CV:# - The coverage cutoff value (0-100).\n\t-BT:# - The number of nonparametric bootstraps for tree building (0-n).\n\t-EVAL:# - Change E-Value cutoff of BLAST search.\n\t-TASK:'task name' If you wish to change which default settings BLAST uses for its search criteria.\n\tDefault - Uses default search criteria as described in Gosselin et al. 2019.\n\nIf you encounter negative values in your matrix, or critical errors send an email to seangosselinofficial\@gmail.com. I will try to answer in a reasonable fassion.\n\n\n";
}
elsif ($ARGV[0] eq "Default"){
	$identity = .7;
	$coverage = .7;
	$bootnum = 100;
}
else{
	foreach my $Inputs (@ARGV){
		if ($Inputs =~ /\-ID\:(.*)/){
			if ($1 < 1){
				$identity = ($1*100);
			}
			else{
				$identity = $1;
			}
		}
		elsif  ($Inputs =~ /\-CV\:(.*)/){
			if ($1 > 1){
				$coverage = ($1/100);
			}
			else{
				$coverage = $1;
			}
		}
		elsif  ($Inputs =~ /\-BT\:(.*)/){
			$bootnum = $1;
		}
		elsif  ($Inputs =~ /\-EVAL\:(.*)/){
			$evalue = $1;
		}
		elsif  ($Inputs =~ /\-TASK\:(.*)/){
			$task = "-task $1";
		}
	}
}
#check for genomes, if none, see else
if (@genomefiles){
	my (@thrs,@thrs2,@thrs3);
	#splits genomes into fragmented files, and creates blast databases
	foreach my $filename (@genomefiles){
		my $thr = threads ->create(\&FirstThreadedTask,$filename);
		push @thrs, $thr;
	}
	$_->join() for @thrs;
	#comprehensive all vs. all BLAST searches and database creation
	foreach my $filename (@genomefiles){
		my $thr = threads ->create(\&BlastGenome,$filename);
		push @thrs2, $thr;
	}
	$_->join() for @thrs2;
	#Calculating distance, ANI, AF, and gANI
	foreach my $filename (@genomefiles){
		my $thr = threads ->create(\&Calculate_ANI,$filename);
		push @thrs3, $thr;
	}
	$_->join() for @thrs3;
	my $orig = "orig";
	#outs to file. Sends file type
	Outfile($orig,$orig,%matrix);
	#bootstrapping
	if ($bootnum == 0){
	}
	else{
		for (my $bc=0; $bc<$bootnum; $bc++){
			BootThread($bc);
		}
	}
}


else{
	print "\nNo genomes present in directory that are recognized.\n\nAcceptable formats - .contig .fna .faa .fasta\n\n";
}


sub FirstThreadedTask{
	#splits genome contigs into 1020 fragments
	my %FileContent = SplitGenomes($_[0]);
	open (SPLIT, "+> $_[0].split");
	foreach my $contigname (sort keys %FileContent){
		print SPLIT ">$contigname\n$FileContent{$contigname}\n";
	}
	close SPLIT;
	#blasts each genome against one another
	FormatDatabase($_[0]);
}

sub BootThread{
	my @thrs4;
	%Bootout = ();
	foreach my $filename (@genomefiles){
		my $thr = threads ->create(\&Sample_and_Calc,$filename);
		push @thrs4, $thr;
	}
	$_->join() for @thrs4;
	#outs to file
	Outfile($_[0],"matrix",%Bootout);
}

#splits genomes into $chop_length BP long chunks and returns these chuncks to %AllFileContent
sub SplitGenomes {
	open (FILE_TO_SPLIT, "< $_[0]");
	my ($sequenceinput,$wholegenome,%TempHash,$Origcontigname);
	my $counter = 0;
	while(<FILE_TO_SPLIT>){
		if ($_=~/\>/){
			if ($counter == 0){
				#stores contig name in convient format
				($Origcontigname) = ($_=~/\>(.*?)[\n\s]/);
				chomp $Origcontigname;
				$Origcontigname =~ s/\s//g;
				$counter++;
				next;
			}
			else{
				#splits genome contig into 1020 long chuncks and places them in Temp hash for output
				$ContigLength{$Origcontigname}=length($sequenceinput);
				my @cutcontig =($sequenceinput =~ /(.{1,1020})/g);
				$sequenceinput = undef;
				foreach my $cuts (@cutcontig){
					next if length($cuts) < 100;
					$wholegenome += (length($cuts));
					my $contigname = "$Origcontigname"."\_$counter";
					$TempHash{$contigname}=$cuts;
					$AllFileContent{$contigname}=(length($cuts));
					$counter++;
				}
				($Origcontigname) = ($_=~/\>(.*?)[\n\s]/);
				chomp $Origcontigname;
			}
		}
		else{
			chomp $_;
			$sequenceinput .= $_;
		}
	}
	#stores length of each contig for later calculations
	$ContigLength{$Origcontigname}=length($sequenceinput);
	my @cutcontig =($sequenceinput =~ /(.{1,1020})/g);
	#stores the final contig
	foreach my $cuts (@cutcontig){
		next if length($cuts) < 100;
		$wholegenome += (length($cuts));
		my $contigname = "$Origcontigname"."\_$counter";
		$TempHash{$contigname}=$cuts;
		$AllFileContent{$contigname}=(length($cuts));
		$counter++;
	}
	close FILE_TO_SPLIT;
	$GenomeLength{$_[0]}=$wholegenome;
	return(%TempHash);
}


sub FormatDatabase{
	if (-e "$_[0].nhr"){
	}
	else{
		#command for legacy blast, nessecary to replicate original results EXCATLY
		#system("formatdb -p F -i $_[0]");
		#commanf for nonlegacy
		#system("makeblastdb -dbtype nucl -in $_[0] -input_type fasta -max_file_sz 2GB");
		#linux command
		system("makeblastdb -dbtype nucl -in $_[0] -input_type fasta -max_file_sz 2GB");
	}
}

sub BlastGenome {
	#self explanitory
	foreach my $database (@genomefiles){
		if (-e "$_[0]+$database.blast"){
			my $size = -s "$_[0]+$database.blast";
			if($size==0){
			}
			else{
				next;
			}
		}
		elsif ($database eq $_[0]){
			next;
		}
		else{
		}
		#command for legacy blast, nessecary to replicate original results EXCATLY
		#system("blastall -p blastn -d $database -i $_[0].split -e 1E-4 -m 8 -K 5 -v 5 -o $_[0]+$database.blast");
		#command for nonlegacy
		#system("blastn -db $database -query $_[0].split -max_target_seqs 9999999 -evalue 1E-4 -outfmt 6 -out $_[0]+$database.blast");
		#Command for more distant sequences WIP/Testing
		system("blastn -db $database -query $_[0].split -evalue $evalue -outfmt 6 -out '$_[0]+$database.blast' $task");
	}
}

#All calculations are done here; name misleading
sub Calculate_ANI{
	my @Blastfiles = glob "$_[0]+*.blast";
	foreach my $blst (@Blastfiles){
		my (%BestHits,$shortergene,$gANINumerator,$TotalID,$TotalShort);
		open (BLAST_FILE, "< $blst");
		my ($name) = ($blst=~/.*?\+(.*?)\..*/);
		my ($HitCounter,$LinkerCode) = 0;
		while (<BLAST_FILE>){
			my @Blast_Lines = split;
			#skip non-best hits
			next if exists $BestHits{$Blast_Lines[0]};
			#skip below ID threshold
			next if $Blast_Lines[2] < $identity;
			#finds shorter gene
			my $gene1 = $ContigLength{$Blast_Lines[1]};
			my $gene2 = $AllFileContent{$Blast_Lines[0]};
			if ($gene1 <= $gene2){
				$shortergene = $gene1;
			}
			else{
				$shortergene = $gene2;
			}
			next if ($shortergene == 0);
			#skips below coverage threshold
			next if (($Blast_Lines[3]/$shortergene) < $coverage);
			#prevents hits against the same "gene" contributing to final calc
			my $booty = ($Blast_Lines[3]*($Blast_Lines[2]/100));
			$BestHits{$Blast_Lines[0]} = 1;
			$gANINumerator += $booty;
			$TotalShort += $shortergene;
			$HitCounter ++;
			$TotalID += $Blast_Lines[2];
			#stores hit for later use in bootstrapping
			$HitsforBoot{"$Blast_Lines[0]\t$name"}="$Blast_Lines[2]\t$shortergene\t$booty";
			#print "DATABASE $Blast_Lines[0]\t$name\n";
		}
		close BLAST_FILE;
		if($HitCounter==0){
			$matrix{"$_[0]\+$name"} = "0.00\t0.00\t0.00\t13.00";
		}
		else{
			#print "Orig Hitc$HitCounter\tNum:$gANINumerator\n";
			my $jANI = ($TotalID/$HitCounter);
			my $gANI = ($gANINumerator/$TotalShort);
			my $AF = ($TotalShort/$GenomeLength{$_[0]});
			my $gDistance = -log($gANINumerator/$GenomeLength{$_[0]});
			$matrix{"$_[0]\+$name"} = "$jANI\t$gANI\t$AF\t$gDistance";
			#print "Query: $_[0], Database: $name\njANInum: $TotalID\ngANInum: $gANINumerator\nGene_Count: $HitCounter\ngDistance: $gDistance\nQuery Length: $GenomeLength{$_[0]}\nTotal Short: $TotalShort\nJANI:\t$jANI\nGANI:\t$gANI\nAF:\t$AF\nDIST:\t$gDistance\n\n";
		}
	}
	$matrix{"$_[0]\+$_[0]"} = "1.0\t1.0\t1.0\t0.0";
}

sub Sample_and_Calc{
	#Creates bootstrapped results for a single query seq
	my $splitfile = "$_[0].split";
	my @Blastfiles = glob "$_[0]+*.blast";
	my %BootSample = Boot_Sample($splitfile);
	#gets blast data for each sample against all genomes
	foreach my $blastfile (sort @Blastfiles){
		my ($gANINumerator,$TotalID,$TotalShort);
		my $HitCounter = 0;
		my ($append) = ($blastfile =~ /.*?\+(.*?)\..*/);
		foreach my $sample (values %BootSample){
			#print "APPEND = $append\:YOUFUCK\n";
			#print "PRESAM = $sample\n";
			my $newsam = $sample;
			$newsam.="\t$append";
			#print "SAMPLE = $newsam\n";
			if(exists $HitsforBoot{$newsam}){
				my $temp =$HitsforBoot{$newsam};
				my @single = split(/	/,$temp);
				$gANINumerator += $single[2];
				$TotalShort += $single[1];
				$HitCounter ++;
				$TotalID += $single[0];
			}
			else{
				next;
			}
		}
		if($HitCounter==0){
			$Bootout{"$_[0]\+$append"} = "0.00\t0.00\t0.00\t13.00";
		}
		else{
			#print "Num:$gANINumerator \tL:$TotalShort \tHits:$HitCounter\n";
			my $jANI = ($TotalID/$HitCounter);
			my $gANI = ($gANINumerator/$TotalShort);
			my $AF = ($TotalShort/$GenomeLength{$_[0]});
			my $gDistance = -log($gANINumerator/$GenomeLength{$_[0]});
			if($gDistance<0){
				$gDistance = 0;
			}
			else{}
			$Bootout{"$_[0]\+$append"} = "$jANI\t$gANI\t$AF\t$gDistance";
		}
	}
	$Bootout{"$_[0]\+$_[0]"} = "1.0\t1.0\t1.0\t0.0";
}


sub Boot_Sample{
	#creates a random sample of the query (with replacement)
	my (%Distribution,%RandomSample);
	my ($counter) = 0;
	open (SPLIT, "< $_[0]");
	while (<SPLIT>){
		if ($_ =~ />.*/){
			my ($matchingname) = ($_=~/\>(.*?)[\n\s]/);
			$Distribution{$counter}=$matchingname;
			$counter++;
		}
		else{
			next;
		}
	}
	close SPLIT;
	for(my $rngcount=0; $rngcount<$counter; $rngcount++){
	$RandomSample{$rngcount} = $Distribution{(keys %Distribution)[rand keys %Distribution ]};
	#print "RNSAMP = $RandomSample{$rngcount}\n";
	}
	return (%RandomSample);
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
		print {$fh} join("\t",sort @genomefiles);
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
