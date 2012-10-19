#!/usr/bin/perl

# convert ucsc refseq table file to gff3


use strict;
use Getopt::Long;
use Pod::Usage;
use Net::FTP;
use Bio::SeqFeature::Lite;
use FindBin qw($Bin);
use tim_data_helper qw(
	format_with_commas
);
use tim_file_helper qw(
	open_to_read_fh
	open_to_write_fh
);
my $VERSION = '1.8.5';

print "\n A script to convert UCSC tables to GFF3 files\n\n";


### Quick help
unless (@ARGV) { 
	# when no command line options are present
	# print SYNOPSIS
	pod2usage( {
		'-verbose' => 0, 
		'-exitval' => 1,
	} );
}



### Command line options
my (
	$ftp_file,
	$database,
	$host,
	$do_chromo,
	$refseqstatusf,
	$refseqsumf,
	$ensemblnamef,
	$ensemblsourcef,
	$kgxreff,
	$chromof,
	$user_source,
	$do_gene,
	$do_cds,
	$do_utr,
	$do_codon,
	$gz,
	$help,
	$print_version,
);
my @genetables;
GetOptions( 
	'ftp=s'      => \$ftp_file, # which database table to retrieve
	'db=s'       => \$database, # which ucsc genome to use
	'host=s'     => \$host, # the ftp server to connect to
	'chr!'       => \$do_chromo, # include the chromosome file from ftp
	'table=s'    => \@genetables, # the input gene table files
	'status=s'   => \$refseqstatusf, # the refseqstatus file
	'sum=s'      => \$refseqsumf, # the refseqsummary file
	'kgxref=s'   => \$kgxreff, # the kgXref info file
	'ensname=s'  => \$ensemblnamef, # the ensemblToGeneName file
	'enssrc=s'   => \$ensemblsourcef, # the ensemblSource file
	'chromo=s'   => \$chromof, # a chromosome file
	'source=s'   => \$user_source, # user provided source
	'gene!'      => \$do_gene, # include genes in output
	'cds!'       => \$do_cds, # include CDS in output
	'utr!'       => \$do_utr, # include UTRs in output
	'codon!'     => \$do_codon, # include start & stop codons in output
	'gz!'        => \$gz, # compress file
	'help'       => \$help, # request help
	'version'    => \$print_version, # print the version
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

# Print help
if ($help) {
	# print entire POD
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

# Print version
if ($print_version) {
	print " Biotoolbox script ucsc_table2gff3.pl, version $VERSION\n\n";
	exit;
}





### Check requirements and defaults
unless (@genetables or $ftp_file) {
	die " Specify either an input table file or a FTP table!\n";
}
if ($ftp_file) {
	unless ($ftp_file =~ m/^refgene|ensgene|xenorefgene|known|all$/i) {
		die " requested table '$ftp_file' by FTP not supported! see help\n";
	}
	unless (defined $database) {
		die " a UCSC genome database must be provided! see help\n";
	}
	unless (defined $do_chromo) {
		$do_chromo = 1;
	}
	unless (defined $host) {
		$host = 'hgdownload.cse.ucsc.edu';
	}
}
unless (defined $do_gene) {
	$do_gene = 1;
}
unless (defined $do_utr) {
	$do_utr = 1;
}
unless (defined $do_cds) {
	$do_cds = 1;
}
my $start_time = time;




### Fetch files if requested
if ($ftp_file) {
	
	# collect the requested files by ftp
	my @files = fetch_files_by_ftp();
	
	# push file names into appropriate variables
	foreach my $file (@files) {
		if ($file =~ /refgene|ensgene|knowngene/i) {
			push @genetables, $file;
		}
		elsif ($file =~ /summary/i) {
			$refseqsumf = $file;
		}
		elsif ($file =~ /status/i) {
			$refseqstatusf = $file;
		}
		elsif ($file =~ /ensembltogene/i) {
			$ensemblnamef = $file;
		}
		elsif ($file =~ /ensemblsource/i) {
			$ensemblsourcef = $file;
		}
		elsif ($file =~ /kgxref/i) {
			$kgxreff = $file;
		}
		elsif ($file =~ /chrom/i) {
			$chromof = $file;
		}
	}
}




### Process the gene tables

# input accessory files
	# the load_extra_data() will read the appropriate file, if available,
	# and return the hash of the data
	# it is generic for handling multiple data types
	# pass the type of table we're working with
my $refseqsum   = load_extra_data('summary');

my $refseqstat  = load_extra_data('status');

my $kgxref      = load_extra_data('kgxref');

my $ensembldata = load_extra_ensembl_data();


# initialize globals
my $chromosome_done = 0; # boolean indicating chromosomes are written
my $source; # a re-usable global, may change depending on input table

# walk through the input tables
foreach my $file (@genetables) {
	
	# open output file
	my ($outfile, $gff_fh) = open_output_gff($file);
	
	# process chromosome
	if ($chromof and !$chromosome_done) {
		# if there is only one genetable, we will prepend the chromosomes 
		# to that output file, otherwise we'll make a separate gff file
		# I'm making this assumption because the chromosomes only need to be 
		# defined once when loading Bio::DB::SeqFeature::Store database
		# If user is collecting multiple gene tables, then separate files 
		# are ok, probably preferable, than a gigantic one
		print " Writing chromosome features....\n";
		
		if (scalar @genetables > 1) {
			# let's write a separate chromosome gff file
			
			# open new filehandle
			my ($chromo_outfile, $chromo_gff_fh) = open_output_gff($chromof);
			
			# convert the chromosomes
			print_chromosomes($chromo_gff_fh);
			
			# done
			$chromo_gff_fh->close;
			print " Wrote chromosome GFF file '$chromo_outfile'\n"; 
			$chromosome_done = 1;
		}
		else {
			# let's write to one gff file
			print_chromosomes($gff_fh);
			$chromosome_done = 1;
		}
	}	
	
	# set the source 
	if (defined $user_source) {
		$source = $user_source;
	}
	else {
		# determine from the input filename
		if ($file =~ /xenorefgene/i) {
			$source = 'xenoRefGene';
		}
		elsif ($file =~ /refgene/i) {
			$source = 'refGene';
		}
		elsif ($file =~ /ensgene/i) {
			$source = 'ensGene';
		}
		elsif ($file =~ /knowngene/i) {
			$source = 'knownGene';
		}
		else {
			$source = 'UCSC';
		}
	}
	
	# open the input gene table
	my $table_fh = open_to_read_fh($file) or
		die " unable to open gene table file '$file'!\n";
	
	# convert the table depending on what it is
	print " Converting gene table '$file' features....\n";
	my $count = process_gene_table($table_fh, $gff_fh);
	
	# report outcomes
	print "  converted ", format_with_commas($count->{gene}), 
		" gene features\n" if $count->{gene} > 0;
	print "  converted ", format_with_commas($count->{mrna}), 
		" mRNA transcripts\n" if $count->{mrna} > 0;
	print "  converted ", format_with_commas($count->{pseudogene}), 
		" pseudogene transcripts\n" if $count->{pseudogene} > 0;
	print "  converted ", format_with_commas($count->{ncrna}), 
		" ncRNA transcripts\n" if $count->{ncrna} > 0;
	print "  converted ", format_with_commas($count->{mirna}), 
		" miRNA transcripts\n" if $count->{mirna} > 0;
	print "  converted ", format_with_commas($count->{snrna}), 
		" snRNA transcripts\n" if $count->{snrna} > 0;
	print "  converted ", format_with_commas($count->{snorna}), 
		" snoRNA transcripts\n" if $count->{snorna} > 0;
	print "  converted ", format_with_commas($count->{trna}), 
		" tRNA transcripts\n" if $count->{trna} > 0;
	print "  converted ", format_with_commas($count->{rrna}), 
		" rRNA transcripts\n" if $count->{rrna} > 0;
	print "  converted ", format_with_commas($count->{other}), 
		" other transcripts\n" if $count->{other} > 0;
	
	# Finished
	printf "  wrote file '$outfile' in %.1f minutes\n", 
		(time - $start_time)/60;
	
}



### Finish
exit;





#########################  Subroutines  #######################################

sub fetch_files_by_ftp {
	
	
	# generate ftp request list
	my @ftp_files;
	if ($ftp_file eq 'all') {
		@ftp_files = qw(
			refgene
			ensgene
			xenorefgene
			known
		);
	}
	elsif ($ftp_file =~ /,/) {
		@ftp_files = split /,/, $ftp_file;
	}
	else {
		push @ftp_files, $ftp_file;
	}
	
	# generate list of files
	my @files;
	foreach my $item (@ftp_files) {
		if ($item =~ m/^xeno/i) {
			push @files, qw(
				xenoRefGene.txt.gz 
				refSeqStatus.txt.gz 
				refSeqSummary.txt.gz
			);
		}
		elsif ($item =~ m/refgene/i) {
			push @files, qw(
				refGene.txt.gz 
				refSeqStatus.txt.gz 
				refSeqSummary.txt.gz
			);
		}
		elsif ($item =~ m/ensgene/i) {
			push @files, qw(
				ensGene.txt.gz 
				ensemblToGeneName.txt.gz
				ensemblSource.txt.gz
			);
		}
		elsif ($item =~ m/known/i) {
			push @files, qw(
				knownGene.txt.gz 
				kgXref.txt.gz 
			);
		}
	}
	# this might seem convulated....
	# but we're putting all the file names in a single array
	# instead of specific global variables
	# to make retrieving through FTP a little easier
	# plus, not all files may be available for each species, e.g. knownGene
	# we also rename the files after downloading them
	
	# we will sort out the list of downloaded files later and assign them 
	# to specific global filename variables
	
	# add chromosome file if requested
	if ($do_chromo) {
		push @files, 'chromInfo.txt.gz';
	}
	
	# set the path based on user provided database
	my $path = 'goldenPath/' . $database . '/database/';
	
	# initiate connection
	print " Connecting to $host....\n";
	my $ftp = Net::FTP->new($host) or die "Cannot connect! $@";
	$ftp->login or die "Cannot login! " . $ftp->message;
	
	# prepare for download
	$ftp->cwd($path) or 
		die "Cannot change working directory to '$path'! " . $ftp->message;
	$ftp->binary;
	
	# download requested files
	my @fetched_files;
	foreach my $file (@files) {
		print "  fetching $file....\n";
		# prepend the local file name with the database
		my $new_file = $database . '_' . $file;
		
		# fetch
		if ($ftp->get($file, $new_file) ) { 
			push @fetched_files, $new_file;
		}
		else {	
			warn "Cannot get file $file: " . $ftp->message;
		}
	}
	$ftp->quit;
	
	print " Finished\n";
	return @fetched_files;
}




sub load_extra_data {
	# this will load extra tables of information into hash tables
	# this includes tables of descriptive information, such as 
	# status, summaries, common gene names, etc.
	
	# this sub is designed to be generic to be reusable for multiple 
	# data files
	
	# identify the appropriate file to use based on the type of 
	# table information being loaded
	# the file name should've been provided by command line or ftp
	my $type = shift;
	my %data;
	my $file;
	if ($type eq 'summary') {
		$file = $refseqsumf;
	}
	elsif ($type eq 'status') {
		$file = $refseqstatusf;
	}
	elsif ($type eq 'kgxref') {
		$file = $kgxreff;
	}
	
	# the appropriate data table file wasn't provided for the requested 
	# data table, return an empty hash
	return \%data unless defined $file;
	
	# load file
	my $fh = open_to_read_fh($file) or 
		die " unable to open $type file '$file'!\n";
	
	# load into hash
	while (my $line = $fh->getline) {
		
		# process line
		chomp $line;
		next if ($line =~ /^#/);
		my @line_data = split /\t/, $line;
		
		# the unique id should be the first element in the array
		# take it off the array, since it doesn't need to be stored there too
		my $id = shift @line_data;
		
		# check for duplicate lines
		if (exists $data{$id} ) {
			warn "  $type line for identifier $id exists twice!\n";
			next;
		}
		
		# store data into hash
		$data{$id} = [@line_data];
	}
	
	# finish
	print " Loaded ", format_with_commas( scalar(keys %data) ), 
		" transcripts from $type file '$file'\n";
	$fh->close;
	return \%data;

	### refSeqStatus table
	# 0	mrnaAcc	RefSeq gene accession name
	# 1	status	Status ('Unknown', 'Reviewed', 'Validated', 'Provisional', 'Predicted', 'Inferred')
	# 2	molecule type ('DNA', 'RNA', 'ds-RNA', 'ds-mRNA', 'ds-rRNA', 'mRNA', 'ms-DNA', 'ms-RNA', 'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'ss-DNA', 'ss-RNA', 'ss-snoRNA', 'tRNA', 'cRNA', 'ss-cRNA', 'ds-cRNA', 'ms-rRNA')	values	molecule type
	
	### refSeqSummary table
	# 0	RefSeq mRNA accession
	# 1	completeness	FullLength ('Unknown', 'Complete5End', 'Complete3End', 'FullLength', 'IncompleteBothEnds', 'Incomplete5End', 'Incomplete3End', 'Partial')	
	# 1	summary	 	text	values	Summary comments
	
	### kgXref table
	# 0	kgID	Known Gene ID
	# 1	mRNA	mRNA ID
	# 2	spID	SWISS-PROT protein Accession number
	# 3	spDisplayID	 SWISS-PROT display ID
	# 4	geneSymbol	Gene Symbol
	# 5	refseq	 RefSeq ID
	# 6	protAcc	 NCBI protein Accession number
	# 7	description	Description
	
	### ensemblToGeneName table
	# 0 Ensembl transcript ID
	# 1 gene name
}


sub load_extra_ensembl_data {
	# we will combine the ensemblToGeneName and ensemblSource data tables
	# into a single hash keyed by the ensGene transcript ID
	# Both tables are very simple two columns, so just trying to conserve
	# memory by combining them
	
	# initialize
	my %data;
		# key will be the ensembl transcript id
		# value will be anonymous array of [name, source]
	
	# load ensemblToGeneName first
	if ($ensemblnamef) {
		
		# open file
		my $fh = open_to_read_fh($ensemblnamef) or 
			die " unable to open file '$ensemblsourcef'!\n";
		
		# load into hash
		my $count = 0;
		while (my $line = $fh->getline) {
			
			# process line
			chomp $line;
			next if ($line =~ /^#/);
			my @line_data = split /\t/, $line;
			if (scalar @line_data != 2) {
				die " file $ensemblnamef doesn't seem right!? Line has " .
					scalar @line_data . " elements!\n";
			}
			
			# store data into hash
			$data{ $line_data[0] }->[0] = $line_data[1];
			$count++;
		}
		
		# finish
		print " Loaded ", format_with_commas($count), 
			" names from file '$ensemblnamef'\n";
		$fh->close;
	}
	
	# load ensemblSource second
	if ($ensemblsourcef) {
	
		# open file
		my $fh = open_to_read_fh($ensemblsourcef) or 
			die " unable to open file '$ensemblsourcef'!\n";
		
		# load into hash
		my $count = 0;
		while (my $line = $fh->getline) {
			
			# process line
			chomp $line;
			next if ($line =~ /^#/);
			my @line_data = split /\t/, $line;
			if (scalar @line_data != 2) {
				die " file $ensemblsourcef doesn't seem right!? Line has " .
					scalar @line_data . " elements!\n";
			}
			
			# store data into hash
			$data{ $line_data[0] }->[1] = $line_data[1];
			$count++;
		}
		
		# finish
		print " Loaded ", format_with_commas($count), 
			" transcript types from file '$ensemblsourcef'\n";
		$fh->close;
	}
	
	# done
	return \%data;
}



sub open_output_gff {
	
	# prepare output file name
	my $file = shift;
	my $outfile = $file;
	$outfile =~ s/\.txt(?:\.gz)?$//i; # remove the extension
	$outfile .= '.gff3';
	if ($gz) {
		$outfile .= '.gz';
	}
	
	# open file handle
	my $fh = open_to_write_fh($outfile, $gz) or
		die " unable to open file '$outfile' for writing!\n";
	
	# print comments
	$fh->print( "##gff-version 3\n");
	$fh->print( "# Generated on " . localtime(time) . "\n");
	$fh->print( "# UCSC table file $file\n");
	
	# finish
	return ($outfile, $fh);
}


sub process_line_data {
	
	my $line = shift;
	my %data;
	
	# load the relevant data from the table line into the hash
	# using the identified column indices
	chomp $line;
	my @linedata = split /\t/, $line;
	
	# we're identifying the type of table based on the number of columns
	# maybe not the best or accurate, but it works for now
	
	# don't forget to convert start from 0 to 1-based coordinates
	
	if (scalar @linedata == 16) {
		# a gene prediction table, e.g. refGene, ensGene, xenoRefGene
		
		# 0  bin
		# 1  name
		# 2  chrom
		# 3  strand
		# 4  txStart
		# 5  txEnd
		# 6  cdsStart
		# 7  cdsEnd
		# 8  exonCount
		# 9  exonStarts
		# 10 exonEnds
		# 11 score
		# 12 name2
		# 13 cdsStartStat
		# 14 cdsEndStat
		# 15 exonFrames
		
		$data{name}        = $linedata[1];
		$data{chrom}       = $linedata[2];
		$data{strand}      = $linedata[3];
		$data{txStart}     = $linedata[4] + 1;
		$data{txEnd}       = $linedata[5];
		$data{cdsStart}    = $linedata[6] + 1;
		$data{cdsEnd}      = $linedata[7];
		$data{exonCount}   = $linedata[8];
		$data{exonStarts}  = [ map {$_ += 1} ( split ",", $linedata[9] ) ];
		$data{exonEnds}    = [ ( split ",", $linedata[10] ) ];
		$data{name2}       = $linedata[12];
#		$data{exonFrames}  = [ ( split ",", $linedata[15] ) ];
		$data{note}        = $refseqsum->{ $linedata[1] }->[1] || undef;
		$data{status}      = $refseqstat->{ $linedata[1] }->[0] || undef;
		$data{completeness} = $refseqsum->{ $linedata[1] }->[0] || undef;
	}
	elsif (scalar @linedata == 12) {
		# a knownGene table
		
		# 0 name	known gene identifier
		# 1 chrom	Reference sequence chromosome or scaffold
		# 2 strand	+ or - for strand
		# 3 txStart	Transcription start position
		# 4 txEnd	Transcription end position
		# 5 cdsStart	Coding region start
		# 6 cdsEnd	Coding region end
		# 7 exonCount	Number of exons
		# 8 exonStarts	Exon start positions
		# 9 exonEnds	Exon end positions
		# 10 proteinID	UniProt display ID for Known Genes, UniProt accession or RefSeq protein ID for UCSC Genes
		# 11 alignID	Unique identifier for each (known gene, alignment position) pair
		
		$data{name}       = $kgxref->{ $linedata[0] }->[0] ||
							$linedata[0];
		$data{chrom}      = $linedata[1];
		$data{strand}     = $linedata[2];
		$data{txStart}    = $linedata[3] + 1;
		$data{txEnd}      = $linedata[4];
		$data{cdsStart}   = $linedata[5] + 1;
		$data{cdsEnd}     = $linedata[6];
		$data{exonCount}  = $linedata[7];
		$data{exonStarts} = [ map {$_ += 1} ( split ",", $linedata[8] ) ];
		$data{exonEnds}    = [ ( split ",", $linedata[9] ) ];
		$data{name2}       = $kgxref->{ $linedata[0] }->[3] || # geneSymbol
							 $kgxref->{ $linedata[0] }->[0] || # mRNA id
							 $kgxref->{ $linedata[0] }->[4] || # refSeq id
							 $linedata[0]; # ugly default
#		$data{exonFrames}  = undef; # not present in this table
		$data{note}        = $kgxref->{ $linedata[0] }->[6] || undef;
		$data{refseq}      = $kgxref->{ $linedata[0] }->[4] || undef;
		$data{status}      = $refseqstat->{ $data{refseq} }->[0] || undef;
		$data{completeness} = $refseqsum->{ $data{refseq} }->[0] || undef;
		$data{spid}        = $kgxref->{ $linedata[0] }->[1] || undef; # SwissProt ID
		$data{spdid}       = $kgxref->{ $linedata[0] }->[2] || undef; # SwissProt display ID
		$data{protacc}     = $kgxref->{ $linedata[0] }->[5] || undef; # NCBI protein accession
	}
	
	else {
		# unrecognized
		warn " Unrecognized line! There are " .  scalar(@linedata) . 
			" columns for this line! skipping\n";
	}

	return \%data;
}




sub process_gene_table {
	
	my ($table_fh, $gff_fh) = @_;
	
	# initialize 
#	my $current_chrom;
	my %gene2seqf; # hash to assemble genes and/or transcripts for this chromosome
	my %id2count; # hash to aid in generating unique primary IDs
	my %counts = (
		'gene'       => 0,
		'mrna'       => 0,
		'pseudogene' => 0,
		'ncrna'      => 0,
		'mirna'      => 0,
		'snrna'      => 0,
		'snorna'     => 0,
		'trna'       => 0,
		'rrna'       => 0,
		'other'      => 0,
	);
	
	
	#### Main Loop
	while (my $line = $table_fh->getline) {
		
		## process the row from the gene table
		next if $line =~ /^#/;
		my $linedata = process_line_data($line);
		
		## generate the transcript
		my $transcript = generate_new_transcript($linedata, \%id2count);
				
		## count the transcript type
		my $type = $transcript->primary_tag;
		if ($type eq 'mRNA') {
			$counts{mrna}++;
		}
		elsif ($type eq 'pseudogene') {
			$counts{pseudogene}++;
		}
		elsif ($type eq 'ncRNA') {
			$counts{ncrna}++;
		}
		elsif ($type eq 'miRNA') {
			$counts{mirna}++;
		}
		elsif ($type eq 'snRNA') {
			$counts{snrna}++;
		}
		elsif ($type eq 'snoRNA') {
			$counts{snorna}++;
		}
		elsif ($type eq 'tRNA') {
			$counts{trna}++;
		}
		elsif ($type eq 'rRNA') {
			$counts{rrna}++;
		}
		else {
			$counts{other}++;
		}
		
		
		## add transcript to gene if requested
		if ($do_gene) {
			# assemble transcripts into genes
			# multiple transcripts may be associated with a single gene
			# genes are store in the gene2seqf hash
			# there may be more than one gene with each gene name, each with 
			# a non-overlapping transcript (how complicated!)
			
			my $gene;
			if (exists $gene2seqf{ lc $linedata->{name2} }) {
				# we already have a gene for this transcript
				
				# pull out the gene seqfeature(s)
				my $genes = $gene2seqf{ lc $linedata->{name2} };
				
				# check that the current transcript intersects with the gene
				# sometimes we can have two separate transcripts with the 
				# same gene name, but located on opposite ends of the chromosome
				# part of a gene family, but unlikely the same gene 200 Mb in 
				# length
				foreach my $g (@$genes) {
					if ( $g->overlaps($transcript, 'strong') ) {
						# gene and transcript overlap on the same strand
						# we found the intersecting gene
						$gene = $g;
						last;
					}
				}
				
				# we have a gene for our transcript
				if ($gene) {
					# update the gene coordinates if necessary
					if ( ($linedata->{txStart}) < $gene->start) {
						# update the transcription start position
						$gene->start( $linedata->{txStart} );
					}
					if ($linedata->{txEnd} > $gene->end) {
						# update the transcription stop position
						$gene->end( $linedata->{txEnd} );
					}
				}
				
				# NONE of the genes and our transcript overlap
				else {
					# must make a new gene
					$gene = generate_new_gene($linedata, \%id2count);
					$counts{gene}++;
					
					# store the new gene oject into the gene hash
					push @{ $genes }, $gene;
				}
			}
			
			else {
				# generate new gene SeqFeature object
				$gene = generate_new_gene($linedata, \%id2count);
				$counts{gene}++;
				
				# store the gene oject into the gene hash
				$gene2seqf{ lc $linedata->{name2} } = [ $gene ];
			} 
			
			
			# associate our transcript with the gene
			$gene->add_SeqFeature($transcript);
			
			# Update gene note if necessary
			unless ($gene->has_tag('Note')) {
				# unless it somehow wasn't available for a previous transcript, 
				# but is now, we'll add it now
				# we won't check if transcripts for the same gene have the 
				# same note or not, why wouldn't they????
				if (defined $linedata->{note}) {
					# make note if it exists
					$gene->add_tag_value('Note', $linedata->{note});
				}
			} 
			
		}
		else {
			# do not assemble transcripts into genes
			# we will still use the gene2seqf hash, just organized by 
			# transcript name
			# there may be more than one transcript with the same name
			
			if (exists $gene2seqf{ lc $linedata->{name} }) {
				push @{ $gene2seqf{ lc $linedata->{name} } }, $transcript;
			}
			else {
				$gene2seqf{ lc $linedata->{name} } = [ $transcript ];
			}
		}
		
		
	} # Finished working through the table
	
	
	
	#### Finished
	
	# print remaining current genes and transcripts
	print_current_gene_list(\%gene2seqf, $gff_fh);
	
	# return the counts
	return \%counts;
}



sub generate_new_gene {
	my ($linedata, $id2counts) = @_;
	
	# Set the gene name
	# in most cases this will be the name2 item from the gene table
	# except for some ncRNA and ensGene transcripts
	my $name;
	if ($linedata->{name} =~ /^ENS/i) {
		# an ensGene transcript, look up the common name if possible
		if (defined $ensembldata->{ $linedata->{name} }->[0] ) {
			
			# use the common name as both the gene name
			$name  = $ensembldata->{ $linedata->{name} }->[0];
		}
		else {
			# use the name2 value
			$name  = $linedata->{name2};
		}
	}
	elsif (!defined $linedata->{name2}) {
		# some genes, notably some ncRNA genes, have no gene or name2 entry
		# we'll fake it and assign the transcript name
		# change it in linedata hash to propagate it in downstream code
		$linedata->{name2} = $linedata->{name};
		$name = $linedata->{name};
	}
	else {
		# default for everything else
		$name = $linedata->{name2};
	}
	
	# Uniqueify the gene ID and name
	my $alias;
	if (exists $id2counts->{ lc $name }) {
		# we've encountered this transcript ID before
		
		# the original name should become the alias
		$alias = $name;
		
		# then make name unique by appending a number
		$name = $name . '.' . $id2counts->{ lc $name };
		
		# remember this one
		# using alias value because that was the original name
		$id2counts->{ lc $alias } += 1;
	}
	else {
		# this is the first transcript with this id
		# set the id counter
		$id2counts->{lc $name} = 1;
	}
	
	
	# generate the gene SeqFeature object
	my $gene = Bio::SeqFeature::Lite->new(
		-seq_id        => $linedata->{chrom},
		-source        => $source,
		-primary_tag   => 'gene',
		-start         => $linedata->{txStart},
		-end           => $linedata->{txEnd},
		-strand        => $linedata->{strand} eq '+' ? 1 : -1,
		-phase         => '.',
		-display_name  => $name,
		-primary_id    => $name,
	);
	
	
	# add original gene name as an alias
	if (defined $alias) {
		$gene->add_tag_value('Alias', $alias);
	}
	
	# add the original ENSDARG identifier as an Alias in addition to ID
	# for ensGene transcripts
	if ($linedata->{name} =~ /^ENS/i) {
		$gene->add_tag_value('Alias', $linedata->{name2});
	}
	
	# add status if possible
	if (defined $linedata->{status} ) {
		$gene->add_tag_value( 'status', $linedata->{status} );
	}
	
	# add Note if possible
	if (defined $linedata->{note} ) {
		$gene->add_tag_value( 'Note', $linedata->{note} );
	}
	
	# add refSeq identifier if possible
	if (defined $linedata->{refseq}) {
		$gene->add_tag_value('refSeq', $linedata->{refseq});
	}
	
	# add SwissProt identifier if possible
	if (defined $linedata->{spid}) {
		$gene->add_tag_value('swiss_prot', $linedata->{spid});
	}
	
	# add SwissProt display identifier if possible
	if (defined $linedata->{spdid}) {
		$gene->add_tag_value('swiss_prot_display_id', $linedata->{spdid});
	}
	
	# add NCBI protein access identifier if possible
	if (defined $linedata->{protacc}) {
		$gene->add_tag_value('ncbi_protein_access', $linedata->{protacc});
	}
	
	# finished
	return $gene;
}



sub generate_new_transcript {
	my ($linedata, $id2counts) = @_;
	
	# Uniqueify the transcript ID and name
	my ($id, $name, $alias);
	if (exists $id2counts->{ lc $linedata->{name} } ) {
		# we've encountered this transcript ID before
		
		# now need to make ID unique by appending a number
		$id    = $linedata->{name} . '.' . $id2counts->{ lc $linedata->{name} };
		$name  = $id;
		$alias = $linedata->{name};
		
		# remember this one
		$id2counts->{ lc $linedata->{name} } += 1;
	}
	else {
		# this is the first transcript with this id
		$id   = $linedata->{name};
		$name = $linedata->{name};
		$id2counts->{lc $id} = 1;
	}
	
	# Generate the transcript SeqFeature object
	my $transcript = Bio::SeqFeature::Lite->new(
		-seq_id        => $linedata->{chrom},
		-source        => $source,
		-start         => $linedata->{txStart},
		-end           => $linedata->{txEnd},
		-strand        => $linedata->{strand} eq '+' ? 1 : -1,
		-phase         => '.',
		-display_name  => $name,
		-primary_id    => $id,
	);
	
	# Attempt to identify the transcript type
	if (
		$linedata->{cdsStart} - 1 == $linedata->{txEnd} and 
		$linedata->{cdsEnd} == $linedata->{txEnd}
		# we need to subtract 1 to the cdsStart to compensate for converting 
		# to 1-based coordinates
	) {
		# there appears to be no coding potential when 
		# txEnd = cdsStart = cdsEnd
		# if you'll look, all of the exon phases should also be -1
		
		# check if we have a ensGene transcript, we may have the type
		if (
			$linedata->{name} =~ /^ENS/i and 
			defined $ensembldata->{ $linedata->{name} }->[1]
		) {
			# this is a ensGene transcript
			
			# these should be fairly typical standards
			# snRNA, rRNA, pseudogene, etc
			$transcript->primary_tag( 
				$ensembldata->{ $linedata->{name} }->[1] );
		}
		
		# otherwise, we may be able to infer some certain 
		# types from the gene name
		
		elsif ($linedata->{name2} =~ /^LOC\d+/) {
			# empirical testing seems to suggest that all the noncoding 
			# genes with a name like LOC123456 are pseudogenes
			# well, at least with hg18, it may not be true for others
			$transcript->primary_tag('pseudogene');
		}
		elsif ($linedata->{name2} =~ /^mir/i) {
			# a noncoding gene whose name begins with mir is likely a 
			# a micro RNA
			$transcript->primary_tag('miRNA');
		}
		elsif ($linedata->{name2} =~ /^snr/i) {
			# a noncoding gene whose name begins with snr is likely a 
			# a snRNA
			$transcript->primary_tag('snRNA');
		}
		elsif ($linedata->{name2} =~ /^sno/i) {
			# a noncoding gene whose name begins with sno is likely a 
			# a snoRNA
			$transcript->primary_tag('snoRNA');
		}
		else {
			# a generic ncRNA
			$transcript->primary_tag('ncRNA');
		}
	}
	else {
		# the transcript has an identifiable CDS
		$transcript->primary_tag('mRNA');
	}
	
	
	# add the Ensembl Gene name if it is an ensGene transcript
	if ($linedata->{name} =~ /^ENS/i) {
		# if we have loaded the EnsemblGeneName data hash
		# we should be able to find the real gene name
		if (defined $ensembldata->{ $linedata->{name} }->[0] ) {
			# we will put the common gene name as an alias
			$transcript->add_tag_value('Alias', 
				$ensembldata->{ $linedata->{name} }->[0] );
		}
	}
	
	# add original transcript name as an alias
	if (defined $alias) {
		$transcript->add_tag_value('Alias', $alias);
	}
	
	# add gene name as an alias
	if (defined $linedata->{name2}) {
		$transcript->add_tag_value('Alias', $linedata->{name2});
	}
	
	# add a status for the transcript
	if (defined $linedata->{status} ) {
		$transcript->add_tag_value( 'status', $linedata->{status} );
	}
	
	# add the completeness value for the tag
	if (defined $linedata->{completeness} ) {
		$transcript->add_tag_value( 'completeness', $linedata->{completeness} );
	}
	
	# add Note if possible
	if (defined $linedata->{note} ) {
		$transcript->add_tag_value( 'Note', $linedata->{note} );
	}
	
	# add refSeq identifier if possible
	if (defined $linedata->{refseq}) {
		$transcript->add_tag_value('refSeq', $linedata->{refseq});
	}
	
	# add SwissProt identifier if possible
	if (defined $linedata->{spid}) {
		$transcript->add_tag_value('swiss_prot', $linedata->{spid});
	}
	
	# add SwissProt display identifier if possible
	if (defined $linedata->{spdid}) {
		$transcript->add_tag_value('swiss_prot_display_id', $linedata->{spdid});
	}
	
	# add NCBI protein access identifier if possible
	if (defined $linedata->{protacc}) {
		$transcript->add_tag_value('ncbi_protein_access', $linedata->{protacc});
	}
	
	# add the exons
	add_exons($transcript, $linedata);
	
	# add CDS, UTRs, and codons if necessary
	if ($transcript->primary_tag eq 'mRNA') {
		
		if ($do_codon) {
			add_codons($transcript, $linedata);
		}
		
		if ($do_utr) {
			add_utrs($transcript, $linedata);
		}
		
		if ($do_cds) {
			add_cds($transcript, $linedata);
		}
	}
	
	# transcript is complete
	return $transcript;
}



sub add_exons {
	my ($transcript, $linedata) = @_;
	
	
	# Add the exons
	for (my $i = 0; $i < $linedata->{exonCount}; $i++) {
			
		# transform index for reverse strands
		# this will allow numbering from 5'->3'
		my $number; 
		if ($transcript->strand == 1) {
			# forward strand
			$number = $i;
		}
		else {
			# reverse strand
			$number = abs( $i - $linedata->{exonCount} + 1);
		}
		
		# build the exon seqfeature
		my $exon = Bio::SeqFeature::Lite->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source,
			-primary_tag   => 'exon',
			-start         => $linedata->{exonStarts}->[$i],
			-end           => $linedata->{exonEnds}->[$i],
			-strand        => $transcript->strand,
			-primary_id    => $transcript->primary_id . ".exon$number",
			-display_name  => $transcript->display_name . ".exon$number",
		);
		
		# associate with transcript
		$transcript->add_SeqFeature($exon);
	}
}



sub add_utrs {
	my ($transcript, $linedata) = @_;
	
	# we will scan each exon and look for a potential utr and build it
	my @utrs;
	for (my $i = 0; $i < $linedata->{exonCount}; $i++) {
		
		# transform index for reverse strands
		# this will allow numbering from 5'->3'
		my $number; 
		if ($transcript->strand == 1) {
			# forward strand
			$number = $i;
		}
		else {
			# reverse strand
			$number = abs( $i - $linedata->{exonCount} + 1);
		}
		
		# identify UTRs
		# we will identify by comparing the cdsStart and cdsStop relative
		# to the exon coordinates
		# the primary tag is determined by the exon strand orientation
		my ($start, $stop, $tag);
		
		# 5'UTR forward, 3'UTR reverse
		if (
			$linedata->{exonStarts}->[$i] < $linedata->{cdsStart}
			and
			$linedata->{exonEnds}->[$i] < $linedata->{cdsStart}
		) {
			# the exon start/end is entirely before the cdsStart
			$start = $linedata->{exonStarts}->[$i];
			$stop  = $linedata->{exonEnds}->[$i];
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';
		}
		
		# Split 5'UTR forward, 3'UTR reverse
		elsif (
			$linedata->{exonStarts}->[$i] < $linedata->{cdsStart}
			and
			$linedata->{exonEnds}->[$i] >= $linedata->{cdsStart}
		) {
			# the start/stop codon is in this exon
			# we need to make the UTR out of a portion of this exon 
			$start = $linedata->{exonStarts}->[$i];
			$stop  = $linedata->{cdsStart} - 1;
			$tag   = $transcript->strand == 1 ? 'five_prime_UTR' : 'three_prime_UTR';
		}
		
		# CDS only
		elsif (
			$linedata->{exonStarts}->[$i] >= $linedata->{cdsStart}
			and
			$linedata->{exonEnds}->[$i] < $linedata->{cdsEnd}
		) {
			# CDS only exon
			next;
		}
		
		# Split 3'UTR forward, 5'UTR reverse
		elsif (
			$linedata->{exonStarts}->[$i] <= $linedata->{cdsEnd}
			and
			$linedata->{exonEnds}->[$i] > $linedata->{cdsEnd}
		) {
			# the stop/start codon is in this exon
			# we need to make the UTR out of a portion of this exon 
			$start = $linedata->{cdsEnd} + 1;
			$stop  = $linedata->{exonEnds}->[$i];
			$tag   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}
	
		# 3'UTR forward, 5'UTR reverse
		elsif (
			$linedata->{exonStarts}->[$i] > $linedata->{cdsEnd}
			and
			$linedata->{exonEnds}->[$i] > $linedata->{cdsEnd}
		) {
			# the exon start/end is entirely after the cdsStop
			# we have a 3'UTR
			$start = $linedata->{exonStarts}->[$i];
			$stop  = $linedata->{exonEnds}->[$i];
			$tag   = $transcript->strand == 1 ? 'three_prime_UTR' : 'five_prime_UTR';
		}
		
		else {
			# something else?
			next;
		}
			
		# build the utr object
		my $utr = Bio::SeqFeature::Lite->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source,
			-start         => $start,
			-end           => $stop,
			-strand        => $transcript->strand,
			-phase         => '.',
			-primary_tag   => $tag,
			-primary_id    => $transcript->primary_id . ".utr$number",
			-display_name  => $transcript->display_name . ".utr$number",
		);
		
		# store this utr seqfeature in a temporary array
		push @utrs, $utr;
	}
	
	# associate found UTRs with the transcript
	foreach my $utr (@utrs) {
		$transcript->add_SeqFeature($utr);
	}
}



sub add_cds {
	my ($transcript, $linedata) = @_;
	
	# we will scan each exon and look for a potential CDS and build it
	my @cdss;
	my $phase = 0; # initialize CDS phase and keep track as we process CDSs 
	for (my $i = 0; $i < $linedata->{exonCount}; $i++) {
		
		# transform index for reverse strands
		my $j;
		if ($transcript->strand == 1) {
			# forward strand
			$j = $i;
		}
		else {
			# reverse strand
			# flip the index for exon starts and stops so that we 
			# always progress 5' -> 3' 
			# this ensures the phase is accurate from the start codon
			$j = abs( $i - $linedata->{exonCount} + 1);
		}
		
		# identify CDSs
		# we will identify by comparing the cdsStart and cdsStop relative
		# to the exon coordinates
		my ($start, $stop);
		
		# Split 5'UTR & CDS on forward, 3'UTR & CDS on reverse
		if (
			$linedata->{exonStarts}->[$j] < $linedata->{cdsStart}
			and
			$linedata->{exonEnds}->[$j] >= $linedata->{cdsStart}
		) {
			# the start/stop codon is in this exon
			# we need to make the CDS out of a portion of this exon 
			$start = $linedata->{cdsStart};
			$stop  = $linedata->{exonEnds}->[$j];
		}
		
		# CDS only
		elsif (
			$linedata->{exonStarts}->[$j] >= $linedata->{cdsStart}
			and
			$linedata->{exonEnds}->[$j] <= $linedata->{cdsEnd}
		) {
			# entire exon is CDS
			$start = $linedata->{exonStarts}->[$j];
			$stop  = $linedata->{exonEnds}->[$j];
		}
	
		# Split 3'UTR & CDS on forward, 5'UTR & CDS on reverse
		elsif (
			$linedata->{exonStarts}->[$j] <= $linedata->{cdsEnd}
			and
			$linedata->{exonEnds}->[$j] > $linedata->{cdsEnd}
		) {
			# the stop/start codon is in this exon
			# we need to make the UTR out of a portion of this exon 
			$start = $linedata->{exonStarts}->[$j];
			$stop  = $linedata->{cdsEnd};
		}
	
		else {
			# UTR exon
			next;
		}
			
		# build the CDS object
		my $cds = Bio::SeqFeature::Lite->new(
			-seq_id        => $transcript->seq_id,
			-source        => $transcript->source,
			-start         => $start,
			-end           => $stop,
			-strand        => $transcript->strand,
			# -phase         => $linedata->{exonFrames}->[$j],
			-phase         => $phase,
			-primary_tag   => 'CDS',
			-primary_id    => $transcript->primary_id . ".cds$i", 
			-display_name  => $transcript->display_name . ".cds$i",
		);
		# the id and name still use $i for labeling to ensure numbering from 0
		
		# store this utr seqfeature in a temporary array
		push @cdss, $cds;
		
		# reset the phase for the next CDS
			# phase + (3 - (length % 3)), readjust to 0..2 if necessary
			# adapted from Barry Moore's gtf2gff3.pl script
		$phase = $phase + (3 - ( $cds->length % 3) );
		$phase -=3 if $phase > 2;
	}
	
	# associate found UTRs with the transcript
	foreach my $cds (@cdss) {
		$transcript->add_SeqFeature($cds);
	}
}



sub add_codons {
	
	my ($transcript, $linedata) = @_;
	
	# generate the start and stop codons
	my ($start_codon, $stop_codon);
	if ($transcript->strand == 1) {
		# forward strand
		
		# start codon
		$start_codon = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'start_codon',
				-start         => $linedata->{cdsStart},
				-end           => $linedata->{cdsStart} + 2,
				-strand        => 1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.start_codon',
				-display_name  => $transcript->display_name . '.start_codon',
		);
		
		# stop codon
		$stop_codon = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'stop_codon',
				-start         => $linedata->{cdsEnd} - 2,
				-end           => $linedata->{cdsEnd},
				-strand        => 1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.stop_codon',
				-display_name  => $transcript->display_name . '.stop_codon',
		);
	}
	
	else {
		# reverse strand
		
		# stop codon
		$stop_codon = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'stop_codon',
				-start         => $linedata->{cdsStart},
				-end           => $linedata->{cdsStart} + 2,
				-strand        => -1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.stop_codon',
				-display_name  => $transcript->display_name . '.stop_codon',
		);
		
		# start codon
		$start_codon = Bio::SeqFeature::Lite->new(
				-seq_id        => $transcript->seq_id,
				-source        => $transcript->source,
				-primary_tag   => 'start_codon',
				-start         => $linedata->{cdsEnd} - 2,
				-end           => $linedata->{cdsEnd},
				-strand        => -1,
				-phase         => 0,
				-primary_id    => $transcript->primary_id . '.start_codon',
				-display_name  => $transcript->display_name . '.start_codon',
		);
	}
	
	# associate with transcript
	$transcript->add_SeqFeature($start_codon);
	$transcript->add_SeqFeature($stop_codon);
}



sub print_current_gene_list {
	my ($gene2seqf, $gff_fh) = @_;
	
	# we need to sort the genes in genomic order before writing the GFF
	my %pos2seqf;
	print "  Sorting ", format_with_commas( scalar(keys %{ $gene2seqf }) ), 
		" top features....\n";
	foreach my $g (keys %{ $gene2seqf }) {
		
		# each value is an array of gene/transcripts
		foreach my $t ( @{ $gene2seqf->{$g} } ) {
		
			# get coordinates
			my $chr   = $t->seq_id;
			my $start = $t->start;
			
			# make sure start positions are unique, just in case
			while (exists $pos2seqf{$chr}{$start}) {
				$start++;
			}
			
			# store the seqfeature
			$pos2seqf{$chr}{$start} = $t;
		}
	}
	
	# print in genomic order
	print "  Writing features to GFF....\n";
	foreach my $chr (sort {$a cmp $b} keys %pos2seqf) {
		# sort by chromosome first
		# just simple ASCIIbetical sort
		
		foreach my $start (sort {$a <=> $b} keys %{ $pos2seqf{$chr} }) {
			# next sort by increasing start position
			
			# set gff version
			$pos2seqf{$chr}{$start}->version(3); 
			
			# print the seqfeature recursively
			$gff_fh->print( $pos2seqf{$chr}{$start}->gff_string(1) . "\n");
				# the gff_string method is undocumented in the POD, but is a 
				# valid method. Passing 1 should force a recursive action to 
				# print both parent and children.
		}
		
		# print directive to close out all previous genes
		$gff_fh->print("###\n"); 
	}
}



sub print_chromosomes {
	
	my $out_fh = shift;
	
	# open the chromosome file
	my $chromo_fh = open_to_read_fh($chromof) or die 
		"unable to open specified chromosome file '$chromof'!\n";
	
	# convert the chromosomes into GFF features
	while (my $line = $chromo_fh->getline) {
		next if ($line =~ /^#/);
		chomp $line;
		my ($chr, $end, $path) = split /\t/, $line;
		unless (defined $chr and $end =~ m/^\d+$/) {
			die " format of chromsome doesn't seem right! Are you sure?\n";
		}
		
		# generate seqfeature
		my $chrom = Bio::SeqFeature::Lite->new(
			-seq_id        => $chr,
			-source        => 'UCSC', # using a generic source here
			-primary_tag   => $chr =~ m/scaffold/i ? 'scaffold' : 'chromosome',
			-start         => 1,
			-end           => $end,
			-primary_id    => $chr,
			-display_name  => $chr,
		);
		
		# print the gff
		$chrom->version(3);
		$out_fh->print( $chrom->gff_string . "\n" );
	}
	
	# finished
	$out_fh->print( "###\n" );
	$chromo_fh->close;
}




__END__

=head1 NAME ucsc_table2gff3.pl



=head1 SYNOPSIS

   ucsc_table2gff3.pl --ftp <text> --db <text>
   
   ucsc_table2gff3.pl [--options] --table <filename>
  
  Options:
  --ftp [refgene|ensgene|xenorefgene|known|all]
  --db <text>
  --host <text>
  --table <filename>
  --status <filename>
  --sum <filename>
  --ensname <filename>
  --enssrc <filename>
  --kgxref <filename>
  --chromo <filename>
  --source <text>
  --(no)chr
  --(no)gene
  --(no)cds
  --(no)utr
  --(no)codon
  --(no)gz
  --version
  --help

=head1 OPTIONS

The command line flags and descriptions:

=over 4

=item --ftp [refgene|ensgene|xenorefgene|known|all]

Request that the current indicated tables and supporting files be 
downloaded from UCSC via FTP. Four different tables may be downloaded, 
including refGene, ensGene, xenoRefGene mRNA gene prediction tables, and 
the UCSC known gene table (if available). Specify all to download all 
four tables. A comma delimited list may also be provided.

=item --db <text>

Specify the genome version database from which to download the requested 
table files. See L<http://genome.ucsc.edu/FAQ/FAQreleases.html> for a 
current list of available UCSC genomes. Examples included hg19, mm9, and 
danRer7.

=item --host <text>

Optionally provide the host FTP address for downloading the current 
gene table files. The default is 'hgdownload.cse.ucsc.edu'.

=item --table <filename>

Provide the name of a UCSC gene or gene prediction table. Tables known 
to work include the refGene, ensGene, xenoRefGene, and UCSC knownGene 
tables. The file may be gzipped. When converting multiple tables, use 
this option repeatedly for each table.

=item --status <filename>

Optionally provide the name of the refSeqStatus table file. This file 
provides additional information for the refSeq-based gene prediction 
tables, including refGene, xenoRefGene, and knownGene tables. The 
file may be gzipped.

=item --sum <filename>

Optionally provide the name of the refSeqSummary file. This file 
provides additional information for the refSeq-based gene prediction 
tables, including refGene, xenoRefGene, and knownGene tables. The 
file may be gzipped.

=item --ensname <filename>

Optionally provide the name of the ensemblToGeneName file. This file 
provides a key to translate the Ensembl unique gene identifier to the 
common gene name. The file may be gzipped.

=item --enssrc <filename>

Optionally provide the name of the ensemblSource file. This file 
provides a key to translate the Ensembl unique gene identifier to the 
type of transcript, provided by Ensembl as the source. The file may be 
gzipped.

=item --kgxref <filename>

Optionally provide the name of the kgXref file. This file 
provides additional information for the UCSC knownGene gene table.
The file may be gzipped.

=item --chromo <filename>

Optionally provide the name of the chromInfo text file. Chromosome 
and/or scaffold features will then be written at the beginning of the 
output GFF file (when processing a single table) or written as a 
separate file (when processing multiple tables). The file may be gzipped.

=item --source <text>

Optionally provide the text to be used as the GFF source. The default is 
automatically derived from the source table file name, if recognized, or 
'UCSC' if not recognized.

=item --(no)chr

When downloading the current gene tables from UCSC using the --ftp 
option, indicate whether (or not) to include the chromInfo table. 
The default is true. 

=item --(no)gene

Specify whether (or not) to assemble mRNA transcripts into genes. This 
will create the canonical gene->mRNA->(exon,CDS) heirarchical structure. 
Otherwise, mRNA transcripts are kept independent. The gene name, when 
available, are always associated with transcripts through the Alias tag. 
The default is true.

=item --no(cds)

Specify whether (or not) to include CDS features in the output GFF file. 
The default is true.

=item --(no)utr

Specify whether (or not) to include three_prime_utr and five_prime_utr 
features in the transcript heirarchy. If not defined, the GFF interpreter 
must infer the UTRs from the CDS and exon features. The default is true.

=item --(no)codon

Specify whether (or not) to include start_codon and stop_codon features 
in the transcript heirarchy. The default is false.

=item --(no)gz

Specify whether the output file should be compressed with gzip.

=item --version

Print the version number.

=item --help

Display the POD documentation

=back

=head1 DESCRIPTION

This program will convert a UCSC gene or gene prediction table file into a
GFF3 format file. It will build canonical gene->transcript->[exon, CDS,
UTR] heirarchical structures. It will attempt to identify non-coding genes
as to type using the gene name as inference. Various additional
informational attributes may also be included with the gene and transcript
features, which are derived from supporting table files.

Four table files are supported. Gene prediction tables, including refGene, 
xenoRefGene, and ensGene, are supported. The UCSC knownGene gene table, if 
available, is also supported. Supporting tables include refSeqStatus, 
refSeqSummary, ensemblToGeneName, ensemblSource, and kgXref. 

The latest table files may be automatically downloaded using FTP from 
UCSC or other host. Since these files are periodically updated, this may 
be the best option. Alternatively, individual files may be specified 
through command line options. Files may be obtained manually through FTP, 
HTTP, or the UCSC Table Browser.

If provided, chromosome and/or scaffold features may also be written to a 
GFF file. If only one table is being converted, then the chromosome features 
are prepended to the GFF file; otherwise, a separate chromosome GFF file is 
written.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Howard Hughes Medical Institute
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  
