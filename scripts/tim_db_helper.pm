package tim_db_helper;

use strict;
require Exporter;
use Carp qw(carp cluck croak confess);
use File::Spec;
use Bio::DB::SeqFeature::Store;
use Statistics::Lite qw(
	sum
	mean
	median
	min
	max
	range
	stddevp
);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_data_helper qw(
	generate_tim_data_structure 
	parse_list
);
use tim_db_helper::config;
our $VERSION = '1.9.0';

# check for wiggle support
our $WIGGLE_OK = 0;
eval {
	require tim_db_helper::wiggle;
	tim_db_helper::wiggle->import;
};
unless ($@) {
	$WIGGLE_OK = 1;
}; 
$@ = undef;

# check for BigWig support
our $BIGWIG_OK = 0;
eval { 
	require tim_db_helper::bigwig;
	tim_db_helper::bigwig->import;
};
unless ($@) {
	$BIGWIG_OK = 1;
}; 
$@ = undef;

# check for BigBed support
our $BIGBED_OK = 0;
eval { 
	require tim_db_helper::bigbed;
	tim_db_helper::bigbed->import;
};
unless ($@) {
	$BIGBED_OK = 1;
}; 
$@ = undef;

# check for Bam support
our $BAM_OK = 0;
eval { 
	require tim_db_helper::bam;
	tim_db_helper::bam->import;
};
unless ($@) {
	$BAM_OK = 1;
}; 
$@ = undef;



# define reusable variables
our $TAG_EXCEPTIONS; # for repeated use with validate_included_feature()
our %total_read_number; # for rpm calculations

# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw();
our @EXPORT_OK = qw(
	open_db_connection
	get_dataset_list 
	validate_dataset_list 
	process_and_verify_dataset 
	check_dataset_for_rpm_support 
	get_new_feature_list 
	get_new_genome_list 
	validate_included_feature 
	get_chromo_region_score 
	get_region_dataset_hash 
	get_chromosome_list 
);


# The true statement
1; 

=head1 NAME

tim_db_helper

=head1 DESCRIPTION

These are helper subroutines to work with microarray and/or next generation
sequencing data stored in a Bioperl SeqFeature Store database. The included
subroutines allow to connect to the database and retrieve data relative to
various genomic features within the database, including genes, transcripts,
transcription start sites, genomic bins, etc. The data may collected and
summarized in a variety of methods, including mean, median, enumeration,
etc.

The data may be stored in the database in a variety of mechanisms. The
simplest is to store the data directly into the database as the score value
of the GFF features stored in the database.   Alternatively, the data may
be stored in a binary wig files and referenced in the attribute tag of the
data's feature. Two types of wig files are supported: a scaled binary file
(.wib file) supported by the module C<Bio::Graphics::Wiggle>, and a binary
BigWig file (.bw file) supported by the module C<Bio::DB::BigWig>. The
BigWig file format is much preferred as it maintains spatial resolution of
the original data and does not lose precision by scaling to 8-bit values,
unlike the .wib file format.

While these functions may appear to be simply a rehashing of the methods
and functions in Bio::DB::SeqFeature::Store, they either provide a simpler
function to often used database methodologies or are designed to work
intimately with the tim data format file and data structures (see
C<tim_file_helper.pm>). One key advantage to these functions is the ability
to work with datasets that are stranded (transcriptome data, for example).

A note on the database. While the Bio::DB::SeqFeature::Store database is
designed to work with GFF3 - formatted data, these functions make some
assumptions that we are working with non-standard GFF3 data. Specifically,
the GFF3 feature's method (third column) is supposed to be a standard
Sequence Ontology term. I'm instead using this column to identify different
datasets. For example, technically one should use a method of
'microarray_oligo' and use the feature name as the identifier of the
dataset, such as 'RSC_ypd_244k'. Instead, these functions assume the method
is being used as the unique identifier, so 'RSC_ypd_244k' is set as both
the method and the feature name in the source GFF3 files. While it would be
possible to use the display_name entirely to select unique datasets, it is
more convenient and accessible to identify through the method.

Historically, this module was initially written to use Bio::DB::GFF for 
database usage. It has since been re-written to use Bio::DB::SeqFeature::Store.

Complete usage and examples for the functions are provided below.

=head1 USAGE

Call the module at the beginning of your perl script and include the module 
names to export. 

  
  use tim_db_helper qw(
	  get_new_feature_list 
	  get_feature_dataset 
  );
  

This will export the indicated subroutine names into the current namespace. 
Their usage is detailed below. The configuration object may also be imported 
into the program's namespace as C<$TIM_CONFIG> to allow access to the local 
database configuration.


=over

=cut






################################################################################
################           General subroutines             #####################
################################################################################


### Open a connection to the SeqFeature Store MySQL database

=item open_db_connection

This module will open a connection to a BioPerl style database.
It returns an object that represents the connection. Several 
different types of databases are supported.

=over 4

=item Bio::DB::SeqFeature::Store database

These may be represented by a relational database (e.g. MySQL database), 
a SQLite database file (file.sqlite or file.db), or a single GFF3 file 
(file.gff) that can be loaded into an in-memory database. In-memory databases 
should only be used with small files as they demand a lot of memory.

Parameters for connecting to a relational database are stored in the BioToolBox 
configuration file, C<biotoolbox.cfg>. These include database adaptors, 
user name, password, etc. Information regarding the configuration file may 
be found within the file itself. 

=item Bio::DB::Sam database 

A self-contained database represented by a sorted, indexed Bam file 
(file.bam). See http://samtools.sourceforge.net for more details. Files 
may be either local or remote (prefixed with http:// or ftp://).

=item Bio::DB::BigWig database

A self-contained database of scores represented by a BigWig (file.bw). See
http://genome.ucsc.edu/goldenPath/help/bigWig.html for more information.
Files may be either local or remote (prefixed with http:// or ftp://).

=item Bio::DB::BigWigSet database

A local or remote directory of one or more BigWig files that can treated 
collectively as a database. A special text file may be included in the 
directory to assign metadata to each BigWig file, including attributes such 
as type, source, display name, strand, etc. See L<Bio::DB::BigWigSet> for 
more information on the formatting of the metadata file.

=item Bio::DB::BigBed database

A self-contained database of regions represented by a BigBed (file.bb). See
http://genome.ucsc.edu/goldenPath/help/bigBed.html for more information.
Files may be either local or remote (prefixed with http:// or ftp://).

=back

Pass the name of a relational database or the path of the database file to 
the subroutine. The opened database object is returned. If it fails, then 
an error message should be generated and nothing is returned.

Example:

	my $db_name = 'cerevisiae';
	my $db = open_db_connection($db_name);
	
	my $file = 'file.bam';
	my $db = open_db_connection($file);


=cut

sub open_db_connection {
	my $database = shift;
	unless ($database) {
		cluck 'no database name passed!';
		return;
	}
	
	# first check if it is a database reference
	my $db_ref = ref $database;
	if ($db_ref =~ /^Bio::DB/) {
		# the provided database is already an open database object
		# nothing to open, return as is
		
		# determine the name if possible
		my $db_name;
		if ($db_ref =~ /^Bio::DB::SeqFeature::Store/) {
			# a SeqFeature database, using any DBI adapter
			$db_name = $database->{'dbh'}->{'name'}; 
				# dig through the object internals to identify the original 
				# name of the database
				# this should be relatively well documented through DBI
				# but could break in the future since it's not official API
		}
		elsif ($db_ref eq 'Bio::DB::Sam') {
			# a Bam database
			$db_name = $database->{'bam_path'};
		}
		else {
			# determining the database name from other sources is
			# either not possible or not easy, so won't bother unless
			# there is a really really good need
			$db_name = q(); # undefined
		}
		
		# return as appropriate either both object and name or just object
		return wantarray ? ($database, $db_name) : $database;
	}
	
	# determine type of database to connect to
	my $db;
	my $error;
	
	
	### Attempt to open the database
	# we go through a series of checks to determine if it is remote, local, 
	# an indexed big data file, SQLite file, etc
	# when all else fails, try to open a SQL connection
	
	# check if it is a remote file
	if ($database =~ /^http|ftp/i) {
		
		# a remote Bam database
		if ($database =~ /\.bam$/i) {
			# open using BigWig adaptor
			if ($BAM_OK) {
				$db = open_bam_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote Bam file" .
						" '$database'! $!\n";
				}
			}
			else {
				$error = " Bam database cannot be loaded because\n" . 
					" Bio::DB::Sam is not installed\n";
			}
		}
		
		# a remote BigBed database
		elsif ($database =~ /\.bb$/i) {
			# open using BigBed adaptor
			if ($BIGBED_OK) {
				$db = open_bigbed_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote BigBed file" .
						" '$database'! $!\n";
				}
			}
			else {
				$error = " BigBed database cannot be loaded because\n" . 
					" Bio::DB::BigBed is not installed\n";
			}
		}
		
		# a remote BigWig database
		elsif ($database =~ /\.bw$/i) {
			# open using BigWig adaptor
			if ($BIGWIG_OK) {
				$db = open_bigwig_db($database);
				unless ($db) {
					$error = " ERROR: could not open remote BigWig file" .
						" '$database'! $!\n";
				}
			}
			else {
				$error = " BigWig database cannot be loaded because\n" . 
					" Bio::DB::BigWig is not installed\n";
			}
		}
		
		# a presumed remote directory, presumably of bigwig files
		else {
			# open using BigWigSet adaptor
			if ($BIGWIG_OK) {
				$db = open_bigwigset_db($database);
				unless ($db) {
					$error = " ERROR: could not open presumed remote " .
						"BigWigSet directory '$database'! $!\n";
				}
			}
			else {
				$error = " Presumed BigWigSet database cannot be loaded because\n" . 
					" Bio::DB::BigWigSet is not installed\n";
			}
		}
	
	}
	
	# a directory, presumably of bigwig files
	elsif (-d $database) {
		# open using BigWigSet adaptor
		if ($BIGWIG_OK) {
			$db = open_bigwigset_db($database);
			unless ($db) {
				$error = " ERROR: could not open local BigWigSet " . 
					"directory '$database'! $!\n";
			}
		}
		else {
			$error = " Presumed BigWigSet database cannot be loaded because\n" . 
				" Bio::DB::BigWigSet is not installed\n";
		}
	}
	
	# check for a known file type
	elsif ($database =~ /gff3|bw|bb|bam|db|sqlite/i) {
		
		# first check that it exists
		if (-e $database and -r _) {
		
			# a single gff3 file that we can load into memory
			if ($database =~ /\.gff3?(?:\.gz)?$/i) {
				# open gff3 file using a memory adaptor
				print " Loading file into memory database...\n";
				eval {
					$db = Bio::DB::SeqFeature::Store->new(
						-adaptor => 'memory',
						-gff     => $database,
					);
				};
				unless ($db) {
					$error = " ERROR: could not load file '$database' into memory!\n";
				}
			}
			
			# a SQLite database
			elsif ($database =~ /\.(?:sqlite|db)$/i) {
				# open using SQLite adaptor
				eval {
					$db = Bio::DB::SeqFeature::Store->new(
						-adaptor  => 'DBI::SQLite',
						-dsn      => $database,
					);
				};
				unless ($db) {
					$error = " ERROR: could not open SQLite file '$database'! $!\n";
				}
			}
			
			# a Bam database
			elsif ($database =~ /\.bam$/i) {
				# open using BigWig adaptor
				if ($BAM_OK) {
					$db = open_bam_db($database);
					unless ($db) {
						$error = " ERROR: could not open local Bam file" .
							" '$database'! $!\n";
					}
				}
				else {
					$error = " Bam database cannot be loaded because\n" . 
						" Bio::DB::Sam is not installed\n";
				}
			}
			
			# a BigBed database
			elsif ($database =~ /\.bb$/i) {
				# open using BigBed adaptor
				if ($BIGBED_OK) {
					$db = open_bigbed_db($database);
					unless ($db) {
						$error = " ERROR: could not open local BigBed file" .
							" '$database'! $!\n";
					}
				}
				else {
					$error = " BigBed database cannot be loaded because\n" . 
						" Bio::DB::BigBed is not installed\n";
				}
			}
			
			# a BigWig database
			elsif ($database =~ /\.bw$/i) {
				# open using BigWig adaptor
				if ($BIGWIG_OK) {
					$db = open_bigwig_db($database);
					unless ($db) {
						$error = " ERROR: could not open local BigWig file" .
							" '$database'! $!\n";
					}
				}
				else {
					$error = " BigWig database cannot be loaded because\n" . 
						" Bio::DB::BigWig is not installed\n";
				}
			}
		}
		
		# file does not exist or can be read
		else {
			if (not -e _) {
				# file does not exist
				$error = " ERROR: file '$database' does not exist!\n";
			}
			else {
				# file must not be readable then
				$error = " ERROR: file '$database' can not be read!\n";
			}
		}
	}
	
	# unrecognized real file
	elsif (-e $database) {
		# file exists, I just don't recognize the extension
		$error = " File '$database' type is not recognized\n";
	}
	
	
	# otherwise assume name of a database
	unless ($db) {
		# open the connection using parameters from the configuration file
		# we'll try to use database specific parameters first, else use 
		# the db_default parameters
		my $adaptor = $TIM_CONFIG->param($database . '.adaptor') || 
			$TIM_CONFIG->param('default_db.adaptor');
		my $user = $TIM_CONFIG->param($database . '.user') || 
			$TIM_CONFIG->param('default_db.user');
		my $pass = $TIM_CONFIG->param($database . '.pass') ||
			$TIM_CONFIG->param('default_db.pass') || undef;
		
		# check for empty passwords
		# config::simple passes an empty array when nothing was defined
		if (ref $pass eq 'ARRAY' and scalar @$pass == 0) {$pass = undef}
		
		# set up the dsn
		# it can be specifically defined
		my $dsn = $TIM_CONFIG->param($database . '.dsn') || undef;
		unless (defined $dsn) {
			# or dsn can be generated with the dsn_prefix
			$dsn = $TIM_CONFIG->param($database . '.dsn_prefix') || 
				$TIM_CONFIG->param('default_db.dsn_prefix');
			$dsn .= $database;
		}
		
		# establish the database connection
		eval {
			# to prevent annoying error messages from B:DB:SF:S
			local $SIG{__WARN__} = sub {}; 
			
			# attempt a connection
			$db = Bio::DB::SeqFeature::Store->new(
				-adaptor => $adaptor,
				-dsn     => $dsn,
				-user    => $user,
				-pass    => $pass,
			);
		};
		
		unless ($db) {
			$error .= " ERROR: unknown $adaptor database '$database'\n";
		}
	}
	
	# conditional return
	if ($db) {
		# return as appropriate either both object and name or just object
		return wantarray ? ($db, $database) : $db;
	} 
	else {
		$error .= " no database could be found or connected!\n";
		warn $error;
		return;
	}
}



### Retrieve a list of the microrarray data sets from the db

=item get_dataset_list

This subroutine will retrieve a list of the available features stored in the 
database and returns a hash of the feature's GFF types, represented as 
"type:source", corresponding to the third and second GFF columns, respectively.
The hash is keyed with an incrementing number, and the value is the GFF type 
of the dataset. A hash is returned rather than a list to help facilitate 
presenting and having the user select an item from the list. The list of
available features are sorted asciibetically before numbered.

Pass either the name of the database or an established database object. 
Supported databases include both Bio::DB::SeqFeature::Store and 
Bio::DB::BigWigSet databases. 

By default, the list of feature types are filtered by the source. Features 
whose source are listed in the C<source_exclude> array of the 
C<biotoolbox.cfg> file are excluded from the final hash. These usually 
include sources from official genomic authorities, such as 'SGD', 'GeneDB', 
'UCSC', 'Ensembl', etc. In this way, only special features (e.g. microarray 
datasets) are included in the list. Filtering is not performed with 
Bio::DB::BigWigSet databases (it is generally not needed).

To include all features without filtering, pass a second true argument 
(1, 'all', etc.).

Example:

	my $db_name = 'cerevisiae';
	my %microarray_dataset = get_dataset_list($db_name);
	foreach (sort {$a <=> $b} keys %microarray_dataset) {
		# print number in ascending order, dataset name
		print $_, $microarray_dataset{$_};
	}
	
	my %all_features = get_dataset_list($db_name, 'all');


=cut

sub get_dataset_list {
	
	my $database = shift;
	my $use_all_features = shift;
	
	# Open a db connection 
	my ($db, $db_name) = open_db_connection($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# get sources to skip
		# usually these are features from an official genome authority
	my %source2skip;
	unless ($use_all_features) {
		foreach ($TIM_CONFIG->param($db_name . '.source_exclude') ) {
			# database specific source exclusions
			$source2skip{$_} = 1;
		}
		unless (keys %source2skip) {
			# no database specific exclusions, we'll read default then
			foreach ($TIM_CONFIG->param('default_db.source_exclude') ) {
				$source2skip{$_} = 1;
			}
		}
	}
		
	# process the database types, according to the type of database
	my %dataset;
	
	# a SeqFeature database
	my $db_ref = ref $db;
	if ($db_ref =~ m/Bio::DB::SeqFeature::Store/) {
		my $i = 1;
		foreach my $type (
			map $_->[1],
			sort {$a->[0] cmp $b->[0]} 
			map [$_->method, $_],
			$db->types
		) {
			# sort the types asciibetically by method
			
			my $source = $type->source;
			
			# add the type to the list
			# check and skip unwanted sources
			unless (exists $source2skip{$source}) {
				# keep if it's not on the unwanted list
				$dataset{$i} = $type;
				$i++;
			}
		}
	}
	
	# a BigWigSet database
	elsif ($db_ref eq 'Bio::DB::BigWigSet') {
		
		# get the metadata
		my $metadata = $db->metadata;
		
		# since a BigWigSet database has very few types, and they are all
		# data sources and not features, there is no need to filter the list
				
		# collect
		my %types;
		foreach my $file (keys %{ $metadata }) {
			# get the type for each file
			my ($primary, $source, $type, $name);
			
			# get the appropriate tags
			foreach my $attribute (keys %{ $metadata->{$file} } ) {
				if ($attribute =~ m/^primary_tag|method$/i) {
					$primary = $metadata->{$file}{$attribute};
				}
				elsif ($attribute =~ m/^source/i) {
					$source = $metadata->{$file}{$attribute};
				}
				elsif ($attribute =~ m/^type/i) {
					$type = $metadata->{$file}{$attribute};
				}
				elsif ($attribute =~ m/^display_name/i) {
					$name = $metadata->{$file}{$attribute};
				}
			}
				
			# assemble the type
				# as per the rest of biotoolbox convention, we prefer to use 
				# the primary_tag or type first, but then use things like
				# display_name
			if ($primary and $source) {
				$types{ "$primary:$source" } += 1;
			}
			elsif ($type) {
				$types{ $type } += 1;
			}
			elsif ($primary) {
				$types{ $primary } += 1;
			}
			elsif ($name) {
				$types{ $name } += 1;
			}
		}
		
		# put the types into the final dataset hash
		my $i = 1;
		foreach my $type (sort {$a cmp $b} keys %types) {
			$dataset{$i} = $type;
			$i++;
		}
	}
	
	# some other database
	else {
		carp " no dataset lists for database type $db_ref!\n";
	}
	return %dataset;
}


### Validate a list of microarray data sets against those in a db

=item validate_dataset_list

This subroutine will validate that a list of microarray data set names exist 
within the given database. This is to help with, for example, checking that 
a dataset name written on a commandline is spelled correctly and actually 
exists within the given database. This is why the above subroutine, 
get_dataset_list(), is so helpful as it avoids having to validate existance and
spelling. 

The subroutine is passed an array. The first element of the array must be 
either the name of the database or an established database object reference. 
The subsequent elements are the names of the datasets to be verified.

The subroutine returns a scalar string consisting of the names of the bad
dataset names (to be passed on to the user). The list is separated by a 
comma and space ", ". 

This subroutine could do with a much better interface and needs a re-write and
re-work of the API.

Example:
	
	my @dataset_names = qw(
		H3_ChIP_44k
		H3k4me3_ChIP_44k
		H3k7ac_ChIP_44k   
	); 
	# 3rd dataset should be H3k9ac_ChIP_44k
	my $db_name = 'cerevisiae';
	my $bad_dataset = validate_dataset_list(
		$db_name,
		@dataset_names
	);
	print $bad_dataset; # prints 'H3k7ac_ChIP_44k'

=cut


sub validate_dataset_list {
	my $database = shift;
	
	# verify passed data
	unless ($database) {
		cluck "no database passed!\n";
		return;
	}
	unless (scalar @_ > 0) { 
		carp "no datasets to validate!\n";
		return;
	}
	
	# collect a list of available datasets
	my %dataset_list = get_dataset_list($database, 1);
		# the 1 forces get_dataset_list to collect all feature types
		# without filtering
	return unless %dataset_list; # no need to continue if we don't have a list
	
	# transform dataset list hash into something more useable
	my %dataset_checklist;
	foreach my $item (values %dataset_list) {
		# store the gff type into our dataset checklist
		$dataset_checklist{$item} = 1; # the value is unimportant
		
		# break the gff type into type:source
		if ($item =~ /^([\w\-]+):.+/) {
			# store the actual type or method, ignore source
			$dataset_checklist{$1} = 1;
		}
	}
	unless (%dataset_checklist) {
		carp "database has no features to validate against!\n";
		return;
	}
	
	
	# now go through the list of datasets to check the name
	my @baddatasets; # an array of the names that fail validation
	foreach my $dataset (@_) {
		# we may have combined datasets indicated by a &
		if ($dataset =~ /&/) {
			foreach (split /&/, $dataset) {
				unless (exists $dataset_checklist{$_}) {
					push @baddatasets, $_;
				}
			}
		} else { # only one dataset
			unless (exists $dataset_checklist{$dataset}) {
				push @baddatasets, $dataset;
			}
		}
	}
	
	# return the name of bad datasets
	return join(", ", @baddatasets);
}



### Process and verify a dataset

=item process_and_verify_dataset()

This subroutine will process a dataset list. It will verify that the 
dataset exists, either in the presented database, or if a local file that 
the file exists and is readable. For file-based datasets, it will prepend 
the 'file:' prefix that is necessary for the get dataset or score 
methods in tim_db_helper.

If no dataset names are passed, then an interactive list will be 
presented to the user for selection. The list will include features 
present in the database for the user to select. One or more features 
may be selected. If the single dataset option is set to true, then only 
one feature is accepted. The user response is validated before 
returning the list. 

To use this subroutine, pass an anonymous array with the following keys 
and values. Not every key is required.

  db       => The name of the database or a reference to an 
              established BioPerl database object. Typically, a 
              Bio::DB::SeqFeature::Store database is used.
  dataset  => Pass either a single dataset name as a scalar or an 
              anonymous array reference of a list of dataset names. 
              These may have been provided as a command line option and 
              need to be verified. If nothing is passed, then a list of 
              possible datasets will be presented to the user to be 
              chosen.
  prompt   => Provide a phrase to be prompted to the user to help in 
              selecting datasets from the list. If none is provided, a 
              generic prompt will be used.
  single   => A Boolean value (1 or 0) indicating whether only a single 
              dataset is allowed when selecting datasets from a 
              presented list. If true, only one dataset choice is 
              accepted. If false, one or more dataset choices are 
              allowed.

The subroutine will return a list of the accepted datasets. It will print 
bad dataset names to standard error.

=cut

sub process_and_verify_dataset {
	
	# Retrieve passed values
	my $arg_ref = shift; # the passed argument values as a hash reference
	
	# Check for single option
	my $single = $arg_ref->{'single'} || 0;
	
	# Collect the datasets
	my @datasets;
	if (exists $arg_ref->{'dataset'} and defined $arg_ref->{'dataset'}) {
		
		# check if it's an anonymous array of datasets
		if (ref $arg_ref->{'dataset'} eq 'ARRAY') {
			@datasets = @{ $arg_ref->{'dataset'} };
		}
		else {
			push @datasets, $arg_ref->{'dataset'};
		}
	}
	
	
	# Open database object
	my $db = open_db_connection( $arg_ref->{'db'} ) or 
		confess "no database name or connection!!\n";
	
	
	# Initialize main output arrays
	my @good_datasets;
	my @bad_datasets;
	
	# Check provided datasets
	if (@datasets) {
		
		# check for multiple comma-delimited datasets
		my @list_to_check;
		foreach my $item (@datasets) {
			if ($item =~ /,/) {
				# this one has a comma, therefore it has more than dataset
				push @list_to_check, split(/,/, $item);
			}
			else {
				# a singleton
				push @list_to_check, $item;
			}
		}
		
		# now verify the datasets
		foreach my $dataset (@list_to_check) {
			
			# check for a remote file
			if ($dataset =~ /^(?: http | ftp) .+ \. (?: bam | bw | bb) $/xi) {
				# a remote file
				# assume it is good, no verification here though
				# it will either work or won't work
				push @good_datasets, $dataset;
			}
			
			# a local file
			elsif ($dataset =~ /\.(?:bw|bb|bam)$/i) {
				# presume we have a local bigfile or aligment file
				
				# user may have requested two or more files to be merged
				# these should be combined with an ampersand
				# check each one 
				my @files;
				foreach my $file (split /\&/, $dataset) {
					if (-e $file) {
						# file exists
						push @files, "file:$file";
					}
					else {
						# file doesn't exist! can't use this set of files
						@files = ();
						last;
					}
				}
				if (@files) {
					push @good_datasets, join("&", @files);
				}
				else {
					push @bad_datasets, $dataset;
				}
			}
			
			# a feature type in a database
			else {
				# must be a database feature type
			
				# check for a database
				unless ($db) {
					cluck " dataset '$dataset' is a presumed database feature ",
						"but no database was passed!\n";
					return;
				}
				
				# validate the given dataset
				my $bad = validate_dataset_list($db, $dataset);
				if ($bad) {
					push @bad_datasets, $dataset;
				}
				else {
					push @good_datasets, $dataset;
				}
			}
		}
	}
	
	# User must select datasets
	else {
		# dataset not specified
		# present the dataset list to the user and get an answer
		
		# check for a database
		unless ($db) {
			cluck " no database provided to select datasets!\n";
			return;
		}
				
		# get the dataset list
		my %datasethash = get_dataset_list($db);
		
		# present the list
		print "\n These are the available data sets in the database:\n";
		foreach (sort {$a <=> $b} keys %datasethash) {
			# print out the list of microarray data sets
			print "  $_\t$datasethash{$_}\n"; 
		}
		
		# prompt the user
		if ($arg_ref->{'prompt'}) {
			# provided custom prompt
			print $arg_ref->{'prompt'};
		}
		else {
			# generic prompt
			print " Enter the number of the data set you would like to analyze  ";
		}
		
		# get answer from the user
		my $answer = <STDIN>;
		chomp $answer;
		my @answer_list = parse_list($answer);
		
		# take the first one if requested
		if ($single) {
			unless (scalar @answer_list == 1) {
				splice(@answer_list, 1);
			}
		}
		
		# verify the answer list
		foreach my $answer (@answer_list) {
			
			# check for merged datasets
			if ($answer =~ /&/) {
				# a merged dataset
				my @list = split /&/, $answer;
				my $check = 1;
				
				# check all are good
				foreach (@list) {
					unless (exists $datasethash{$_}) {
						$check = 0;
					}
				}
				
				# if all are good
				if ($check) {
					push @good_datasets, 
						join( "&", map { $datasethash{$_} } @list);
				}
				else {
					push @bad_datasets, $answer;
				}
			}
			
			else {
				# a single dataset
				# check if it is good
				
				if (exists $datasethash{$answer}) {
					push @good_datasets, $datasethash{$answer};
				} 
				else {
					push @bad_datasets, $datasethash{$answer};
				}
			}
		}
	}
	
	# Print bad results
	if (@bad_datasets) {
		print " The following datasets could not be verified:\n";
		foreach (@bad_datasets) {
			print "      $_\n";
		}
	}
	
	# Return good results
	if ($single) {
		return $good_datasets[0];
	}
	else {
		return @good_datasets;
	}
}





### Process and verify a dataset

=item check_dataset_for_rpm_support()

This subroutine will check a dataset for RPM, or Reads Per Million mapped, 
support. Only two types of database files support this, Bam files and 
BigBed files. If the dataset is either one of these, or the name of a 
database feature which points to one of these files, then it will 
calculate the total number of mapped alignments (Bam file) or features 
(BigBed file). It will return this total number. If the dataset does 
not support RPM (because it is not a Bam or BigBed file, for example), 
then it will return undefined.

Pass this subroutine one or two values. The first is the name of the 
dataset. Ideally it should be validated using process_and_verify_dataset() 
and have an appropriate prefix (file, http, or ftp). If it does not 
have a prefix, then it is assumed to be a database feature. The second 
passed feature is the name of a BioPerl database or an opened database 
object. This database will be checked for the indicated dataset, and 
the first returned feature checked for an attribute pointing to a 
supported file.

=cut

sub check_dataset_for_rpm_support {
	
	# get passed dataset and databases
	my $dataset = shift;
	my $database = shift;
	
	# check that we haven't done this already
	
	
	# Check the dataset to see if supports RPM or RPKM method
	# if so, then calculate the total number of reads
	# this uses the global variable $rpkm_read_sum
	my $rpm_read_sum;
	
	if (exists $total_read_number{$dataset}) {
		# this dataset has already been summed
		# no need to do it again
		
		$rpm_read_sum = $total_read_number{$dataset};
	}
	
	elsif ($dataset =~ /\.bam$/) {
		# a bam file dataset
		
		if ($BAM_OK) {
			# tim_db_helper::bam was loaded ok
			# sum the number of reads in the dataset
			$rpm_read_sum = sum_total_bam_alignments($dataset);
		}
		else {
			carp " Bam support is not available! " . 
				"Is Bio::DB::Sam installed?\n";
			return;
		}
	}
	
	elsif ($dataset =~ /\.bb$/) {
		# a bigbed file dataset
		
		if ($BIGBED_OK) {
			# tim_db_helper::bigbed was loaded ok
			# sum the number of features in the dataset
			$rpm_read_sum = sum_total_bigbed_features($dataset);
		}
		else {
			carp " BigBed support is not available! " . 
				"Is Bio::DB::BigBed installed?\n";
			return;
		}
	}
	
	elsif ($dataset !~ /^file|ftp|http/i) {
		# a database feature
		# this feature might point to a bam or bigbed file
		
		# Check the database
		my $db; # the database object to be used
		if (defined $database) {
			my $db_ref = ref $database;
			if ($db_ref =~ /^Bio::DB/) {
				$db = $database;
			}
			else {
				# the name of a database was passed, create a database connection
				$db = open_db_connection($database);
				unless ($db) {
					carp " No database to check feature '$dataset' for RPM support!\n";
					return;
				}
			}
		}
		else {
			carp " No database to check feature '$dataset' for RPM support!\n";
			return;
		}
		
		# get a sample of the features from the database
		my @features;
		if ($dataset =~ /&/) {
			# in case we have a combined datasets
			@features = $db->features(-type => [ split /&/, $dataset ]);
		}
		else {
			@features = $db->features(-type => $dataset);
		}
		unless (@features) {
			 # no feature found
			 # basically nothing to do
			 return;
		}
		
		# look for the database file in the attributes
		if ($features[0]->has_tag('bamfile')) {
			# specifying a bam file
			my ($bamfile) = $features[0]->get_tag_values('bamfile');
			
			if ($BAM_OK) {
				# tim_db_helper::bam was loaded ok
				# sum the number of reads in the dataset
				$rpm_read_sum = sum_total_bam_alignments($bamfile);
			}
			else {
				carp " Bam support is not available! " . 
					"Is Bio::DB::Sam installed?\n";
				return;
			}
		}
		
		elsif ($features[0]->has_tag('bigbedfile')) {
			# specifying a bigbed file
			my ($bedfile) = $features[0]->get_tag_values('bigbedfile');
			
			if ($BIGBED_OK) {
				# tim_db_helper::bigbed was loaded ok
				# sum the number of features in the dataset
				$rpm_read_sum = sum_total_bigbed_features($bedfile);
			}
			else {
				carp " BigBed support is not available! " . 
					"Is Bio::DB::BigBed installed?\n";
				return;
			}
		}
		
		else {
			# can't find an appropriate dataset
			# nothing to do
			return;
		}
	}
	
	else {
		# some other non-supported dataset
		return;
	}
	
	# return the sum value if we've made it this far
	$total_read_number{$dataset} = $rpm_read_sum;
	return $rpm_read_sum;
}







################################################################################
################           Feature subroutines             #####################
################################################################################


### Generate a new list of features


=item get_new_feature_list 

This subroutine will generate a new feature list collected from the database. 
Once the list of genomic features is generated, then data may be collected
for each item in the list. 

The subroutine will generate and return a data hash as described in 
tim_file_helper.pm. The data table will have two or three columns. The 
feature name and type:source are listed in columns one and two, respectively.
If the features have an Alias tag, then a third column is included with 
a comma delimited list of the feature aliases.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db       => The name of the database or a reference to an 
              established database object. 
  features => A scalar value containing a name representing the 
              type(s) of feature(s) to collect. This name will be 
              parsed into an actual list with the internal subroutine 
              _features_to_classes(). Refer to that documentation for 
              a list of appropriate features.
  Optional: 
  dubious  => A boolean value (1 or 0) indicating whether genes 
              flagged in the database as 'dubious' should be 
              included. The default is false (not kept).

The subroutine will return a reference to the data hash. It will print 
status messages to STDOUT. 

Example

	my $db_name = 'cerevisiae';
	my %data = get_new_feature_list( {
		'db'        => $db_name,
		'features'  => 'genes',
		'dubious'   => 0,
	} );


=cut


sub get_new_feature_list {

	# Retrieve passed values
	my $arg_ref = shift; # the passed argument values as a hash reference
	
	# Open a db connection 
	my ($db, $db_name) = open_db_connection($arg_ref->{'db'});
	unless ($db) {
		carp 'no database connected!';
		return;
	}

	# Verify a SeqFeature::Store database
	my $db_ref = ref $db;
	unless ($db_ref =~ /^Bio::DB::SeqFeature::Store/) {
		carp "Database type $db_ref doesn't support generating feature lists!\n";
		return;
	}
	
	
	# Translate the features into a list of classes
	my @classes = _features_to_classes($arg_ref->{'features'});
	unless (@classes) {
		carp "no or unknown features passed!";
		return;
	}
	
	# Generate data structures
	my $new_data = generate_tim_data_structure(
		$arg_ref->{'features'},
		'Name',
		'Type'
	);
	unless ($new_data) {
		cluck " cannot generate tim data structure!\n";
		return;
	}
	my $feature_table = $new_data->{'data_table'}; 
	
	# name of the database
	$new_data->{'db'} = $db_name; 
	
	# List of types
	if (scalar @classes > 1) {
		$new_data->{1}->{'include'} = join(",", @classes);
	}
	
	# Collect the genes from the database
	print "   Searching for " . join(", ", @classes) . "\n";
	my @featurelist; # an array of found feature objects in the database
	@featurelist = $db->features(
			-types => \@classes
	); 
	unless (@featurelist) {
		# there should be some features found in the database
		carp "no features found in database!";
		return;
	}
	print "   Found " . scalar @featurelist . " features in the database.\n";
	
	
	# Get the names of chromosomes to avoid
	my @excluded_chromosomes = 
		$TIM_CONFIG->param("$db_name\.chromosome_exclude");
	unless (@excluded_chromosomes) {
		@excluded_chromosomes = 
			$TIM_CONFIG->param('default_db.chromosome_exclude');
	}
	my %excluded_chr_lookup = map {$_ => 1} @excluded_chromosomes;
	
	
	# Check for aliases
	for (my $i = 0; $i < 50; $i++) {
		# we're checking the first 50 or so features looking for an Alias tag
		# checking that many because not all features may have the tag
		# we like to have Aliases, because it makes interpreting gene names
		# a little easier
		# or all of them if there are 
		last unless (defined $featurelist[$i]);
		
		if ($featurelist[$i]->has_tag('Alias')) {
			
			# add an Alias column to the data table
			push @{ $feature_table->[0] }, 'Aliases';
			$new_data->{2} = {
					'name'  => 'Aliases',
					'index' => 2,
			};
			$new_data->{'number_columns'} = 3;
			last;
		}
	}
	
	# Process the features
	FEATURE_COLLECTION_LIST:
	foreach my $feature (@featurelist) {
		
		# skip genes from excluded chromosomes
		if (exists $excluded_chr_lookup{ $feature->seq_id } ) {
			next FEATURE_COLLECTION_LIST;
		}
		
		# skip anything that matches the tag exceptions
		unless ( validate_included_feature($feature) ) {
			next FEATURE_COLLECTION_LIST;
		}
		
		
		# Record the feature information
		my @data = (
			$feature->display_name, 
			$feature->type 
		);
		
		# Add alias info if available
		if (exists $new_data->{2}) {
			# we only look for Alias info if we have a column for it
			if ($feature->has_tag('Alias')) {
				push @data, join(q( ), $feature->get_tag_values('Alias'));
			}
			else {
				push @data, '.'; # internal null value
			}
		}
		
		# Record information
		push @{$feature_table}, \@data;
		$new_data->{'last_row'} += 1;
	}
	
	
	# print result of search
	print "   Kept " . $new_data->{'last_row'} . " features.\n";
	
	# sort the table
	my @feature_table_sorted;
	my $header = shift @{$feature_table}; # move header
	@feature_table_sorted = sort { 
		# sort first by type, then by name
		( $a->[1] cmp $b->[1] ) and ( $a->[0] cmp $b->[0] )
	} @{$feature_table}; 
	unshift @feature_table_sorted, $header;
	
	# put the feature_table into the data hash
	$new_data->{'data_table'} = \@feature_table_sorted;
	
	# return the new data structure
	return $new_data;
}



### Generate a new list genomic windows

=item get_new_genome_list 

This subroutine will generate a new list of genomic windows. The genome
is split into intervals of a specified size that is moved along the 
genome in specified step sizes.

The subroutine will generate and return a data hash as described in 
tim_file_helper.pm. The data table will have 3 columns, including 
Chromosome, Start, and Stop.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db       => The name of the database or a reference to an 
              established database object. 
  Optional: 
  win      => A scalar value containing an integer representing the
              size of the window in basepairs. The default value 
              is defined in biotoolbox.cfg file.
  step     => A scalar value containing an integer representing the
              step size for advancing the window across the genome. 
              The default is the window size.

The subroutine will return a reference to the data hash. It will print 
status messages to STDOUT. 

Example

	my $db_name = 'cerevisiae';
	my $window_size = 500;
	my $step_size = 250;
	my %data = get_new_genome_list( {
		'db'        => $db_name,
		'win'       => $window_size,
		'step'      => $step_size,
	} );


=cut

sub get_new_genome_list {

	# Collect the passed arguments
	my $arg_ref = shift; 
	
	
	# Open a db connection 
	my ($db, $db_name) = open_db_connection($arg_ref->{'db'});
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	
	# Determine win and step sizes
	my ($win, $step);
	if ($arg_ref->{'win'}) {
		$win = $arg_ref->{'win'};
	}
	else {
		$win = 
			$TIM_CONFIG->param("$db_name\.window") ||
			$TIM_CONFIG->param('default_db.window');
		print "  Using default window size $win bp\n";
	}
	if ($arg_ref->{'step'}) {
		$step = $arg_ref->{'step'};
	}
	else {
		$step = $win;
	}
	
	
	# Generate data structures
	my $new_data = generate_tim_data_structure(
		'genome',
		'Chromosome',
		'Start',
		'Stop'
	);
	unless ($new_data) {
		cluck " cannot generate tim data structure!\n";
		return;
	}
	my $feature_table = $new_data->{'data_table'}; 
	
	# Begin loading basic metadata information
	$new_data->{'db'}      = $db_name; # the db name
	$new_data->{1}{'win'}  = $win; # under the Start metadata
	$new_data->{1}{'step'} = $step;
	
	
	
	# Collect the chromosomes
	# include option to exclude those listed in biotoolbox.cfg that
	# we don't want
	my @chromosomes = get_chromosome_list($db, 1);
	unless (@chromosomes) {
		carp " no sequence IDs were found in the database!\n";
		return;
	}
	
	
	# Collect the genomic windows
	print "   Generating $win bp windows in $step bp increments\n";
	foreach (@chromosomes) {
		
		# name and length as sub-array in each element
		my ($chr, $length) = @{$_};
		
		for (my $start = 1; $start <= $length; $start += $step) {
			# set the end point
			my $end = $start + $win - 1; 
			
			if ($end > $length) {
				# fix end to the length of the chromosome
				$end = $length;
			} 
			
			# add to the output list
			push @{$feature_table}, [ $chr, $start, $end,];
			$new_data->{'last_row'}++;
		}
	}
	print "   Kept " . $new_data->{'last_row'} . " windows.\n";
	
	# Return the data structure
	return $new_data;
}


=item validate_included_feature

This subroutine will validate a database feature to make sure it is 
useable. It will check feature attributes and compare them against 
a list of attributes and values to be avoided. The list of unwanted 
attributes and values is stored in the BioToolBox configuration file 
biotoolbox.cfg. 

Pass the subroutine a Bio::DB::SeqFeature::Store database feature. 
It will return true (1) if the feature passes validation and false 
(undefined) if it contains an excluded attribute and value.

=cut

sub validate_included_feature {
	
	# feature to check
	my $feature = shift;
	
	# get the list of feature exclusion tags
	unless (defined $TAG_EXCEPTIONS) {
		$TAG_EXCEPTIONS = $TIM_CONFIG->get_block('exclude_tags');
	}
	
	# Check the tag exceptions
	# we will check for the presence of the exception tag
	# if the feature tag value matches the exception tag value
	# then the feature should be excluded
	foreach my $key (keys %{ $TAG_EXCEPTIONS }) {
		if ($feature->has_tag($key)) {
			# first check that the feature has this tag
			my @tag_values = $feature->get_tag_values($key);
			if (ref $TAG_EXCEPTIONS->{$key} eq 'ARRAY') {
				# there's more than one exception value!
				# need to check them all
				foreach my $exception (@{ $TAG_EXCEPTIONS->{$key} }) {
					if (grep {$_ eq $exception} @tag_values) {
						# this feature should be excluded
						return;
					}
				}
			}
			else {
				# single tag exception value
				if (grep {$_ eq $TAG_EXCEPTIONS->{$key}} @tag_values) {
					# this feature should be excluded
					return;
				}
			}
		}
	}
	
	# if we've made it thus far, the feature is good
	return 1;
}






################################################################################
##################           Score subroutines             #####################
################################################################################



### Get a dataset score for a single region

=item get_chromo_region_score 

This subroutine will retrieve a dataset value for a single specified 
region in the genome. The region is specified with chromosomal coordinates:
chromosome name, start, and stop. It will collect all dataset values within the
window, combine them with the specified method, and return the single value.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db       => The name of the database or a reference to an 
              established database object. 
  dataset  => The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  method   => The method used to combine the dataset values found
              in the defined region. Acceptable values include 
              sum, mean, median, range, stddev, min, max, rpm, 
              and rpkm. See _get_segment_score() documentation 
              for more info.
  chromo   => The name of the chromosome (reference sequence)
  start    => The start position of the region on the chromosome
  stop     => The stop position of the region on the chromosome
  end      => Alias for stop
  
  Optional:
  strand   => The strand of the region (-1, 0, or 1) on the 
              chromosome. The default is 0, or unstranded.
  value    => Specify the type of value to collect. Acceptable 
              values include score, count, or length. The default 
              value type is score. 
  log      => Boolean value (1 or 0) indicating whether the dataset 
              values are in log2 space or not. If undefined, the 
              dataset name will be checked for the presence of the 
              phrase 'log2' and set accordingly. This argument is
              critical for accurately mathematically combining 
              dataset values in the region.
  stranded => Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, or all.
              Default is 'all'.
         	  
The subroutine will return the region score if successful.

Examples

	my $db = open_db_connection('cerevisiae');
	my $score = get_chromo_region_score( {
		'db'      => $db,
		'method'  => 'mean',
		'dataset' => $dataset,
		'chr'     => $chromo,
		'start'   => $startposition,
		'stop'    => $stopposition,
		'log'     => 1,
	} );
	


=cut

sub get_chromo_region_score {
	
	# retrieve passed values
	my $arg_ref = shift; 
	
	# check the data source
	unless ($arg_ref->{'dataset'}) {
		confess " no dataset requested!";
	}
	
	# Open a db connection 
	my $db = open_db_connection( $arg_ref->{'db'} ) or 
		confess "no database name or connection!!\n";
	
	# get coordinates
	my $chromo = $arg_ref->{'chromo'} || $arg_ref->{'seq'} || 
		$arg_ref->{'seq_id'} || undef;
	my $start = $arg_ref->{'start'} || undef;
	my $stop = $arg_ref->{'stop'} || $arg_ref->{'end'} || undef;
	my $strand = $arg_ref->{'strand'} || 0;
	unless ($chromo and $start and $stop) {
		cluck "one or more genomic region coordinates are missing!";
		return;
	};
	
	# define default values as necessary
	my $value_type = $arg_ref->{'value'} || 'score';
	my $stranded = $arg_ref->{'stranded'} || 'all';
	my $log = $arg_ref->{'log'} || undef;
	unless (defined $log) {
		# we need to know whether we are working with a log2 dataset or not
		# as it will adversely affect the math!
		if ($arg_ref->{'dataset'} =~ /log2/i) {
			# if we're smart we'll encode the log2 status in the dataset name
			# but chances are, we're not that smart
			$log = 1;
		} else {
			# otherwise assume it is non-log
			# unsafe, but what else to do? we'll put the onus on the user
			$log = 0;
		}
	}
	
	# get the scores for the region
	# pass to internal subroutine to combine dataset values
	return _get_segment_score(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$arg_ref->{'dataset'},
				$value_type,
				$arg_ref->{'method'}, 
				$stranded,  
				$log,
	);
}





### Retrieve hash of dataset values for a region


=item get_region_dataset_hash 

This subroutine will retrieve dataset values or feature attributes from
features located within a defined region and return them as a hash.
The (start) positions will be the hash keys, and the corresponding dataset 
values or attributes will be the hash values. The region is defined based on 
a genomic feature in the database. The region corresponding to the entire 
feature is selected by default. Alternatively, it may be adjusted by passing 
appropriate arguments.

Different dataset values may be collected. The default is to collect 
score values of the dataset features found within the region (e.g. 
microarray values). Alternatively, a count of found dataset features
may be returned, or the lengths of the found dataset features. When
lengths are used, the midpoint position of the feature is used in the
returned hash rather than the start position.

The returned hash is keyed by relative coordinates and their scores. For 
example, requesting a region from -200 to +200 of a feature (using the 
start and stop options, below) will return a hash whose keys are relative 
to the feature start position, i.e. the keys will >= -200 and <= 200. 
Absolute coordinates relative to the reference sequence or chromosome 
may be optionally returned instead.

The subroutine is passed a reference to an anonymous hash containing the 
arguments. The keys include

  Required:
  db       => The name of the annotation database or a reference to 
              an established database object. 
  dataset  => The name of the dataset in the database to be 
              collected. The name should correspond to a feature 
              type in the database, either as type or type:source. 
              The name should be verified using the 
              subroutine validate_dataset_list() prior to passing.
              Multiple datasets may be given, joined by '&', with no
              spaces. Alternatively, specify a data file name. 
              A local file should be prefixed with 'file:', while 
              a remote file should be prefixed with the transfer 
              protocol (ftp: or http:).
  name     => The name of the genomic feature.
  type     => The type of the genomic feature.
  Optional:
  ddb      => The name of the data-specific database or a reference 
              to an established database. Use when the data and 
              annotation are in separate databases.
  chromo   => The chromosome or sequence name (seq_id). This may be 
              used instead of name and type to specify a genomic 
              segment. This must be used with start and stop options, 
              and optionally strand options.
  start    => Indicate an integer value representing the start  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "stop".
  stop|end => Indicate an integer value representing the stop  
              position of the region relative to the feature start.
              Use a negative value to begin upstream of the feature.
              Must be combined with "start".
  extend   => Indicate an integer value representing the number of 
              bp the feature's region should be extended on both
              sides.
  position => Indicate the relative position of the feature from 
              which the "start" and "stop" positions are calculated.
              Three values are accepted: "5", which denotes the 
              5' end of the feature, "3" which denotes the 
              3' end of the feature, or "4" which denotes the 
              middle of the feature. This option is only used in 
              conjunction with "start" and "stop" options. The 
              default value is "5".
  strand   => For those features or regions that are NOT 
              inherently stranded (strand 0), artificially set the 
              strand. Three values are accepted: -1, 0, 1. This 
              will overwrite any pre-existing value (it will not, 
              however, migrate back to the database).
  stranded => Indicate whether the dataset values from a specific 
              strand relative to the feature should be collected. 
              Acceptable values include sense, antisense, or all.
              Default is all.
  value    => Indicate which attribute will be returned. Acceptable 
              values include "score", "count", or "length". The  
              default behavior will be to return the score values.
  avoid    => Boolean value to indicate that other features of the 
              same type should be avoided. This only works if name 
              and type was provided. Any positioned scores which 
              overlap the other feature(s) are not returned. The 
              default is false (return all values).
  absolute => Boolean value to indicate that absolute coordinates 
              should be returned, instead of transforming to 
              relative coordinates, which is the default.
          	  
The subroutine will return the hash if successful.

Example

	my $db = open_db_connection('cerevisiae');
	my %region_scores = get_region_dataset_hash( {
		'db'      => $db,
		'dataset' => $dataset,
		'name'    => $name,
		'type'    => $type,
	} );
	


=cut

sub get_region_dataset_hash {
	
	# retrieve passed values
	my $arg_ref = shift; 
	
	### Initialize parameters
	
	# Open a db connection 
	my $db = open_db_connection( $arg_ref->{'db'} ) or 
		confess "no database name or connection!!\n";
	
	# Open the data database if provided
	my $ddb;
	if (defined $arg_ref->{'ddb'}) {
		$ddb = open_db_connection( $arg_ref->{'ddb'} ) or
			confess "requested data database could not be opened!\n";
	}
	else {
		# use the same database for both
		$ddb = $db;
	}
	
	# check the data source
	unless ($arg_ref->{'dataset'}) {
		confess " no dataset requested!";
	}
	
	# confirm options and check we have what we need 
	my $name   = $arg_ref->{'name'}   || undef;
	my $type   = $arg_ref->{'type'}   || undef;
	my $chromo = $arg_ref->{'chromo'} || undef;
	my $start  = $arg_ref->{'start'}  || undef;
	my $stop   = $arg_ref->{'stop'}   || $arg_ref->{'end'} || undef;
	my $strand = $arg_ref->{'strand'} || undef;
	unless (
		(defined $name and defined $type) or 
		(defined $chromo and defined $start and defined $stop)
	) {
		cluck "the feature name and type or genomic coordinates are missing!";
		return;
	};
	
	# assign defaults
	my $stranded       = $arg_ref->{'stranded'} || 'all';
	my $value_type     = $arg_ref->{'value'}    || 'score';
	my $relative_pos   = $arg_ref->{'position'} || 5;
	
	
	# the final coordinates
	my $fref_pos; # to remember the feature reference position
	my $fchromo;
	my $fstart;
	my $fstop;
	my $fstrand;
	
	
	
	### Define the chromosomal region segment
	# we will use the primary database to establish the intitial feature
	# and determine the chromosome, start and stop
	
	
	# Extend a named database feature
	if (
		defined $name and 
		defined $type and 
		$arg_ref->{'extend'}
	) {
		# if an extension is specified to the feature region
		# first define the feature
		my @features = $db->features( 
				-name  => $name,
				-type  => $type,
		);
		if (scalar @features > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $feature = shift @features; 
		
		# record the feature reference position and strand
		if ($relative_pos == 5 and $feature->strand >= 0) {
			$fref_pos = $feature->start;
		}
		elsif ($relative_pos == 3 and $feature->strand >= 0) {
			$fref_pos = $feature->end;
		}
		elsif ($relative_pos == 5 and $feature->strand < 0) {
			$fref_pos = $feature->end;
		}
		elsif ($relative_pos == 3 and $feature->strand < 0) {
			$fref_pos = $feature->start;
		}
		elsif ($relative_pos == 4) {
			# strand doesn't matter here
			$fref_pos = $feature->start + int(($feature->length / 2) + 0.5);
		}
		
		# record final coordinates
		$fchromo = $feature->seq_id;
		$fstart  = $feature->start - $arg_ref->{'extend'};
		$fstop   = $feature->end + $arg_ref->{'extend'};
		$fstrand = defined $strand ? $strand : $feature->strand;
	} 
		
	# Specific start and stop coordinates of a named database feature
	elsif (
			defined $name and
			defined $type and 
			defined $start and
			defined $stop
	) {
		
		# first define the feature to get endpoints
		my @features = $db->features( 
				-name  => $name,
				-type => $type,
		);
		if (scalar @features > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $feature = shift @features; 
		
		# determine the cooridnates based on the identified feature
		if ($relative_pos == 5 and $feature->strand >= 0) {
			# feature is on forward, top, watson strand
			# set segment relative to the 5' end
			
			# record final coordinates
			$fref_pos  = $feature->start;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->start + $start;
			$fstop     = $feature->start + $stop;
			$fstrand   = defined $strand ? $strand : $feature->strand;
		}
		
		elsif ($relative_pos == 5 and $feature->strand < 0) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 5' end
			
			# record final coordinates
			$fref_pos = $feature->end;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->end - $stop;
			$fstop     = $feature->end - $start;
			$fstrand   = defined $strand ? $strand : $feature->strand;
		}
		
		elsif ($relative_pos == 3 and $feature->strand >= 0) {
			# feature is on forward, top, watson strand
			# set segment relative to the 3' end
			
			# record final coordinates
			$fref_pos = $feature->end;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->end + $start;
			$fstop     = $feature->end + $stop;
			$fstrand   = defined $strand ? $strand : $feature->strand;
		}
		
		elsif ($relative_pos == 3 and $feature->strand < 0) {
			# feature is on reverse, bottom, crick strand
			# set segment relative to the 3' end
			
			# record final coordinates
			$fref_pos = $feature->start;
			$fchromo   = $feature->seq_id;
			$fstart    = $feature->start - $stop;
			$fstop     = $feature->start - $start;
			$fstrand   = defined $strand ? $strand : $feature->strand;
		}
		
		elsif ($relative_pos == 4) {
			# feature can be on any strand
			# set segment relative to the feature middle
			
			# record final coordinates
			$fref_pos = $feature->start + int(($feature->length / 2) + 0.5);
			$fchromo   = $feature->seq_id;
			$fstart    = $fref_pos + $start;
			$fstop     = $fref_pos + $stop;
			$fstrand   = defined $strand ? $strand : $feature->strand;
		}
	}
	
	# an entire named database feature
	elsif (
		defined $name and 
		defined $type
	) {
		my @features = $db->features( 
				-name  => $name,
				-type => $type,
		);
		if (scalar @features > 1) {
			# there should only be one feature found
			# if more, there's redundant or duplicated data in the db
			# warn the user, this should be fixed
			warn " Found more than one feature of '$type => $name' in " .
				"the database!\n Using the first feature only!\n";
		}
		my $feature = shift @features; 
		
		# determine the strand
		$fstrand   = defined $strand ? $strand : $feature->strand;
		
		# record the feature reference position and strand
		if ($relative_pos == 5 and $fstrand >= 0) {
			$fref_pos = $feature->start;
		}
		elsif ($relative_pos == 3 and $fstrand >= 0) {
			$fref_pos = $feature->end;
		}
		elsif ($relative_pos == 5 and $fstrand < 0) {
			$fref_pos = $feature->end;
		}
		elsif ($relative_pos == 3 and $fstrand < 0) {
			$fref_pos = $feature->start;
		}
		elsif ($relative_pos == 4) {
			# strand doesn't matter here
			$fref_pos = $feature->start + int(($feature->length / 2) + 0.5);
		}
		
		# record final coordinates
		$fchromo   = $feature->seq_id;
		$fstart    = $feature->start;
		$fstop     = $feature->end;
	}
	
	# a genomic region
	elsif (
		defined $chromo and
		defined $start and 
		defined $stop
	) {
		# coordinates are easy
		$fchromo   = $chromo;
		$fstart    = $start;
		$fstop     = $stop;
		
		# determine the strand
		$fstrand   = defined $strand ? $strand : 0; # default is no strand
		
		# record the feature reference position and strand
		if ($relative_pos == 5 and $fstrand >= 0) {
			$fref_pos = $fstart;
		}
		elsif ($relative_pos == 3 and $fstrand >= 0) {
			$fref_pos = $fstop;
		}
		elsif ($relative_pos == 5 and $fstrand < 0) {
			$fref_pos = $fstop;
		}
		elsif ($relative_pos == 3 and $fstrand < 0) {
			$fref_pos = $fstart;
		}
		elsif ($relative_pos == 4) {
			# strand doesn't matter here
			$fref_pos = $fstart + int( ( ($fstop - $fstart + 1) / 2) + 0.5);
		}
	}
	
	# or else something is wrong
	else {
		confess " programming error! not enough information provided to" .
			" identify database feature!\n";
	}
	
	
	### Data collection
	my %datahash = _get_segment_score(
		$ddb, # using the data database here
		$fchromo,
		$fstart,
		$fstop,
		$fstrand, 
		$arg_ref->{'dataset'}, 
		$value_type,
		'indexed', # method
		$stranded, 
		0, # log value
	);
	
	
	### Check for conflicting features
	if (exists $arg_ref->{'avoid'} and $arg_ref->{'avoid'} and $type) {
		# we need to look for any potential overlapping features of the 
		# same type and remove those scores
		
		# get the overlapping features of the same type
		my @overlap_features = $db->features(
			-seq_id  => $fchromo,
			-start   => $fstart,
			-end     => $fstop,
			-type    => $type
		);
		if (@overlap_features) {
			# there are one or more feature of the same type in this 
			# region
			# one of them is likely the one we're working with
			# but not necessarily - user may be looking outside original feature
			# the others are not what we want and therefore need to be 
			# avoided
			foreach my $feat (@overlap_features) {
				
				# skip the one we want
				next if ($feat->display_name eq $name);
				
				# now eliminate those scores which overlap this feature
				foreach my $position (keys %datahash) {
					
					# delete the scored position if it overlaps with 
					# the offending feature
					if (
						$position >= $feat->start and
						$position <= $feat->end
					) {
						delete $datahash{$position};
					}
				}
			}
		}
		
	}
	
	
	
	### Convert the coordinates to relative positions
		# previous versions of this function that used Bio::DB::GFF returned 
		# the coordinates as relative positions, e.g. -200..200
		# to maintain this compatibility we will convert the coordinates to 
		# relative positions
		# most downstream applications of this function expect this, and it's
		# a little easier to work with. Just a little bit, though....
	if (exists $arg_ref->{'absolute'} and $arg_ref->{'absolute'}) {
		# do not convert to relative positions
		return %datahash;
	}
	else {
		my %relative_datahash;
		if ($fstrand >= 0) {
			# forward strand
			foreach my $position (keys %datahash) {
				# relative position is real position - reference
				$relative_datahash{ $position - $fref_pos } = $datahash{$position};
			}
		}
		elsif ($fstrand < 0) {
			# reverse strand
			foreach my $position (keys %datahash) {
				# the relative position is -(real position - reference)
				$relative_datahash{ $fref_pos - $position } = $datahash{$position};
			}
		}
		
		# return the collected dataset hash
		return %relative_datahash;
	}
}




=item get_chromosome_list

This subroutine will collect a list of chromosomes or reference sequences 
in a Bio::DB database and return the list along with their sizes in bp. 
Many BioPerl-based databases are supported, including 
Bio::DB::SeqFeature::Store, Bio::DB::Sam, Bio::DB::BigWig, 
Bio::DB::BigWigSet, and Bio::DB::BigBed, or any others that support the 
"seq_ids" method. See the L<open_db_connection> subroutine for more 
information.

Pass the subroutine either the name of the database, the path to the 
database file, or an opened database object.

Optionally pass a second value, a boolean argument to limit and exclude 
unwanted chromosomes as defined by the "chromosome_exclude" option in 
the BioToolBox configuration file, C<biotoolbox.cfg>. A true value limits 
chromosomes, and false includes all chromosomes. The default is to return 
all chromosomes. Sometimes some sequences are simply not wanted in 
analysis, like the mitochondrial chromosome or unmapped contigs.

The subroutine will return an array, with each element representing each 
chromosome or reference sequence in the database. Each element is an anonymous 
array of two elements, the chromosome name and length in bp.

Example
	
	my $db = open_db_connection('cerevisiae');
	# get all chromosomes in the database
	my @chromosomes = get_chromosome_list($db);
	foreach (@chromosomes) {
		my $name = $_->[0];
		my $length = $_->[1];
		print "chromosome $name is $length bp\n";
	}

=cut

sub get_chromosome_list {
	
	# options
	my $database = shift;
	my $limit = shift || 0;
	
	# Open a db connection 
	my ($db, $db_name) = open_db_connection($database);
	unless ($db) {
		carp 'no database connected!';
		return;
	}
	
	# Check for BigWigSet database
	# these need to be handled a little differently
	if (ref $db eq 'Bio::DB::BigWigSet') {
		# BigWigSet databases are the only databases that don't 
		# support the seq_ids method
		# instead we have to look at one of the bigwigs in the set
		my $bw_file = ($db->bigwigs)[0];
		$db = open_db_connection($bw_file);
	}
	
	# Get chromosome exclusion list
	# these are chromosomes that we do not want to include in the final 
	# list
	# they must be explicitly requested to ignore
	my %excluded_chr_lookup;
	if ($limit) {
		my @excluded_chromosomes = 
			$TIM_CONFIG->param("$db_name\.chromosome_exclude");
		unless (@excluded_chromosomes) {
			@excluded_chromosomes = 
				$TIM_CONFIG->param('default_db.chromosome_exclude');
		}
		%excluded_chr_lookup = map {$_ => 1} @excluded_chromosomes;
	}
		
	
	# Collect chromosome lengths
	# we want to maintain the original order of the chromosomes so we 
	# won't be putting it into a hash
	# instead an array of arrays
	my @chrom_lengths;
	
	# Database specific approaches to collecting chromosomes
	# I used to have one generic approach, but idiosyncrasies and potential 
	# bugs make me use different approaches for better consistency
	
	# Bigfile
	if (ref $db eq 'Bio::DB::BigWig' or ref $db eq 'Bio::DB::BigBed') {
		foreach my $chr ($db->seq_ids) {
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
			# get chromosome size
			my $length = $db->length($chr);
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# Bam
	elsif (ref $db eq 'Bio::DB::Sam') {
		for my $tid (0 .. $db->n_targets - 1) {
			# each chromosome is internally represented in the bam file as 
			# a numeric target identifier
			# we can easily convert this to an actual sequence name
			# we will force the conversion to go one chromosome at a time
			
			# sequence info
			my $chr    = $db->target_name($tid);
			my $length = $db->target_len($tid);
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}
	}
	
	# SeqFeature::Store or other Bioperl
	else {
		foreach my $chr ($db->seq_ids) {
			
			# check for excluded chromosomes
			if (exists $excluded_chr_lookup{$chr} ) {
				next;
			}
			
			# generate a segment representing the chromosome
			# due to fuzzy name matching, we may get more than one back
			my @segments = $db->segment($chr);
			# need to find the right one
			my $segment;
			while (@segments) {
				$segment = shift @segments;
				last if $segment->seq_id eq $chr;
			}
			
			# check segment
			unless ($segment) {
				carp " No genome segment for seq_id $chr!!?? skipping\n";
				next;
			}
			
			# get the chromosome length
			my $length = $segment->length;
			
			# store
			push @chrom_lengths, [ $chr, $length ];
		}	
	}
	
	
	# Return
	unless (@chrom_lengths) {
		carp " no chromosome sequences identified in database!\n";
		return;
	}
	return @chrom_lengths;
}




################################################################################
###############           Internal subroutines             #####################
################################################################################


### Internal subroutine to convert a feature category into a list of classes



=item _features_to_classes

This internal subroutine provides a conveniant look up and conversion of a 
single-word description of a category of features into a list of actual
feature classes in the database. For example, the word 'gene' may include
all ORFs, snRNAs, snoRNAs, and ncRNAs.

Pass the subroutine the feature category name as a scalar value. The 
actual list of feature types will be collected and returned as an array. 
Multiple values may be passed as a comma-delimited string (no spaces).

The aliases and feature lists are specified in the tim_db_helper 
configuration file, biotoolbox.cfg. Additional lists and aliases 
may be placed there. The lists are database specific, or they can be 
added to the default database.

If the passed category name does not match an alias in the config file, 
then it is assumed to be a feature in the database. No further validation 
will be done (if it's not valid, simply no features would be returned from 
the database). 

Also, feature types may be passed as the GFF's method:source, in which case 
they are assumed to be valid and not checked.

=cut

sub _features_to_classes {
	my $feature = shift;
	my @types;
		
	my $alias2types = $TIM_CONFIG->get_block('features');
	if (exists $alias2types->{$feature} ) {
		# looks like the feature is an alias for a list of features
		# defined in the config file
		if (ref $alias2types->{$feature} eq 'ARRAY') {
			# there's a list of them
			@types = @{ $alias2types->{$feature} };
		}
		else {
			# only one
			$types[0] = $alias2types->{$feature};
		}
	}
	else { 
		# We'll assume this is a specific feature in the database.
		# It may be provided as type:source or type.
		# We'll pass these on directly to the originating subroutine
		# more than one may be passed delimited by commas, but no spaces
		@types = split /,/, $feature;
	} 
	
	return @types;
}




### Internal subroutine to retrieve the scores from an established region object

=item _get_segment_score

This internal subroutine is used to collect the dataset scores for an 
established genomic region. It works with a variety of data sources. 

First, the data may be stored directly in the Bio::DB::SeqFeature::Store 
(using the original GFF score value). Second, the feature may reference 
a data file as an attribute (e.g., wigfile=/path/to/file.wib). Finally, 
the name(s) of a data file may be passed from which to collect the data. 
Supported data files include BigWig (.bw), BigBed (.bb), and Bam (.bam). 
A Bio::DB::BigWigSet database is also supported.

The subroutine is passed an array of ten specific values, all of which 
must be defined and presented in this order. These values include
  
  [0] The opened database object that may be used to generate the 
      the segment and contain the data to collect. If the dataset 
      request is from a big file (bigWig, bigBed, Bam), then this 
      database will likely not be used. Otherwise a 
      Bio::DB::SeqFeature::Store or Bio::DB::BigWigSet database 
      object should be passed.
  [1] The chromosome or seq_id of the segment to examine
  [2] The start position of the segment
  [3] The stop or end position of the segment
  [4] The strand of the region (or original feature) -1, 0, or 1.
  [5] The dataset name for filename. Multiple datasets may be included, 
      delimited with an ampersand (&). Multiple datasets are merged into 
      one, unless excluded by strand. Local data source files should be 
      prepended with 'file:', while remote data source files should be 
      prepended with the transfer protocol (http: or ftp:).
  [6] The data type to be collecting. In most cases, the score value 
      is used, but others may be collected. Accepted values include
         
         score
         count
         length
         
  [7] The method of combining all of the dataset values found in the 
      segment into a single value. Accepted values include
         
         sum
         mean
         median
         min
         max
         range (returns difference between max and min)
         stddev (returns the standard deviation of a population)
         indexed (returns hash of postion => score)
         rpm (returns reads per million mapped, only valid with 
              bam and bigbed databases)
         rpkm (same as rpm but normalized for length in kb)
         
  [8] The strandedness of acceptable data. Genomic segments 
      established from an inherently stranded database feature 
      (e.g. ORF) have a specific strandedness. If the dataset strand 
      does not match the specified criteria for strandedness, then it 
      is ignored. If the dataset does not have a specified strand, 
      then it is used regardless of the specified criteria for 
      strandedness. Currently, only transcription data is stranded. 
      Accepted values include
         
         sense       take only values on the same strand
         antisense   take only values on the opposite strand
         all         take all values
         
  [9] The log status of the dataset. Many microarray datasets are in 
      log2 format, and must be de-logged to perform accurate 
      mathematical calculations, such as mean or median. Supply a 
      boolean (0 or 1) value.
      
  
The subroutine returns a single numeric value if appropriate dataset values
are identified in the genomic segment. If not, then a non-value (.) is 
returned. Use of a period as a non-value avoids un-defined errors in 
some subsequent programmatic manipulations and provides a convenient human 
visual marker.

=cut

sub _get_segment_score {
	
	# get passed arguments
	my (
		$db,
		$chromo,
		$start,
		$stop,
		$strand, 
		$dataset, 
		$value_type,
		$method, 
		$strandedness, 
		$log
	) = @_;
	
	# define
	my @scores; # array of collected scores
	my %pos2data; # hash of position to scores
	my $dataset_type; # remember what type of database the data is from
	my $iterator; # seqfeature stream object for reiterating db features
	my $db_type = ref $db; # source of the originating db 
	
	my @datasetlist = split /[&,]/, $dataset; 
		# multiple datasets may be combined into a single search, for example
		# transcriptome data that is on f and r strands. These are given as
		# ampersand or comma delimited lists
	
	
	### Determine where we are going to get the data
		# first check whether the provided dataset(s) look like a data file
		# next check whether the database segment object came from a BigWigSet
		# finally assume it is a SeqFeature database object
		# then look for a wigfile, bigwigfile, or bamfile attribute
		# finally then just take the score directly from the database objects
	
	### Data source files provided
	if ($datasetlist[0] =~ /^file|http|ftp/) {
		
		# collect the data according to file type
		
		# BigWig Data file
		if ($datasetlist[0] =~ /\.bw$/i) {
			# file is in bigwig format
			# this uses the Bio::DB::BigWig adaptor
			
			# check that we have bigwig support
			if ($BIGWIG_OK) {
				# get the dataset scores using tim_db_helper::bigwig
				
				# the data collection depends on the method
				if ($value_type eq 'score' and 
					$method =~ /min|max|mean|sum|count/
				) {
					# we can use the low-level, super-speedy, summary method 
					# warn " using collect_bigwig_score() with file\n";
					return collect_bigwig_score(
						$chromo,
						$start,
						$stop,
						$method,
						@datasetlist
					);
				}
				
				elsif ($value_type eq 'count' and $method eq 'sum') {
					# we can use the low-level, super-speedy, summary method 
					# warn " using collect_bigwig_score() with file\n";
					return collect_bigwig_score(
						$chromo,
						$start,
						$stop,
						'count', # special method
						@datasetlist
					);
				}
				
				elsif ($method eq 'indexed') {
					# collect hash of position => scores
					# warn " using collect_bigwig_position_score() with file\n";
					return collect_bigwig_position_scores(
						$chromo,
						$start,
						$stop,
						@datasetlist
					);
				}
				
				else {
					# use the longer region collection method
					# warn " using collect_bigwig_scores() with file\n";
					@scores = collect_bigwig_scores(
						$chromo,
						$start,
						$stop,
						@datasetlist
					);
					$dataset_type = 'bw';
				}
			}
			else {
				croak " BigWig support is not enabled! " . 
					"Is Bio::DB::BigWig installed?\n";
			}
		}		
		
		# BigBed Data file
		elsif ($datasetlist[0] =~ /\.bb$/i) {
			# data is in bigbed format
			# this uses the Bio::DB::BigBed adaptor
			
			# check that we have bigbed support
			if ($BIGBED_OK) {
				# get the dataset scores using tim_db_helper::bigbed
				
				if ($method eq 'indexed') {
					# warn " using collect_bigbed_position_scores() with file\n";
					return collect_bigbed_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
				}
				
				else {
					# warn " using collect_bigbed_scores() with file\n";
					@scores = collect_bigbed_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
					$dataset_type = 'bb';
				}
			}
			else {
				croak " BigBed support is not enabled! " . 
					"Is Bio::DB::BigBed installed?\n";
			}
		}
		
		# BAM data file
		elsif ($datasetlist[0] =~ /\.bam$/i) {
			# data is in bam format
			# this uses the Bio::DB::Sam adaptor
			
			# check that we have Bam support
			if ($BAM_OK) {
				# get the dataset scores using tim_db_helper::bam
				
				if ($method eq 'indexed') {
					# warn " using collect_bam_position_scores() with file\n";
					return collect_bam_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
				}
				else {
					# warn " using collect_bam_scores() with file\n";
					@scores = collect_bam_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type, 
						@datasetlist
					);
					$dataset_type = 'bam';
				}
			}
			else {
				croak " Bam support is not enabled! " . 
					"Is Bio::DB::Sam installed?\n";
			}
		}
		
		# Unsupported Data file
		else {
			confess " Unsupported file type for file '$datasetlist[0]!\n";
		}
		
	}
	
	
	### BigWigSet database
	elsif ($db_type =~ m/^Bio::DB::BigWigSet/) {
		# calling features from a BigWigSet database object
		
		# we may be able to take advantage of a special low-level 
		# super-speedy interface based on the BigWigSet summary feature
		
		# the data collection depends on the method
		if ($value_type eq 'score' and 
			$method =~ /min|max|mean|sum|count/
		) {
			# we can use the low-level, super-speedy, summary method 
			# warn " using collect_bigwigset_score()\n";
			return collect_bigwigset_score(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				$method,
				@datasetlist
			);
		}
		elsif ($value_type eq 'count' and $method eq 'sum') {
			# we can use the low-level, super-speedy, summary method 
			# warn " using collect_bigwigset_score()\n";
			return collect_bigwigset_score(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				'count', # special method
				@datasetlist
			);
		}
		
		elsif ($value_type eq 'score' and $method eq 'indexed') {
			# want positioned score data
			# warn " using collect_bigwigset_position_score()\n";
			return collect_bigwigset_position_scores(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				@datasetlist
			);
		}
		
		else {
			# simply collect a list of the scores
			# warn " using collect_bigwigset_scores()\n";
			@scores = collect_bigwigset_scores(
				$db,
				$chromo,
				$start,
				$stop,
				$strand, 
				$strandedness, 
				@datasetlist
			);
			$dataset_type = 'bw';
		}
	}
		
	
	### SeqFeature database
	elsif ($db_type =~ m/^Bio::DB::SeqFeature/) {
		# a SeqFeature database
		# normal collection
		$iterator = $db->get_seq_stream(
			-seq_id      => $chromo,
			-start       => $start,
			-end         => $stop,
			-primary_tag => [@datasetlist],
		);
	}
	
	
	### Some other database?
	else {
		# some other Bio::DB database????
		# until I code in every single possibility
		# let's just try a basic features method using whatever the 
		# default type is and hope for the best
		$iterator = $db->get_seq_stream(
			-seq_id      => $chromo,
			-start       => $start,
			-end         => $stop,
		);
	}
		
	
	
	### Process database SeqFeature objects
	if ($iterator) {
		# We have a seqfeature object stream
		# First check whether we're dealing with a datafile pointed to by
		# an attribute tag 
		# Failing that, assume it's the seqfeature objects themselves we want
		
		# collect the first feature
		my $feature = $iterator->next_seq;
		return _return_null($method) unless $feature;
		
		# deal with features that might not be from the chromosome we want
		# sometimes chromosome matching is sloppy and we get something else
		while ($feature->seq_id ne $chromo) {
			$feature = $iterator->next_seq || undef;
		}
		return _return_null($method) unless $feature;
		
		## Wig Data
		if ( $feature->has_tag('wigfile') ) {
			# data is in wig format, or at least the first datapoint is
			
			# determine the type of wigfile
			my ($wigfile) = $feature->get_tag_values('wigfile');
			
			## Bio::Graphics wib file
			if ($wigfile =~ /\.wib$/) {
				# data is in old-style binary wiggle format
				# based on the Bio::Graphics::Wiggle adaptor
				
				# get the full list of features to pass off to the 
				# helper subroutine
				my @features;
				push @features, $feature;
				while (my $f = $iterator->next_seq) {
					push @features, $f;
				}
				
				# check that we have wiggle support
				if ($WIGGLE_OK) {
					# get the dataset scores using tim_db_helper::wiggle
					
					if ($method eq 'indexed') {
						# warn " using collect_wig_position_scores() from tag\n";
						return collect_wig_position_scores(
							$start,
							$stop,
							$strand, 
							$strandedness, 
							@features
						);
					}
					else {
						# warn " using collect_wig_scores() from tag\n";
						@scores = collect_wig_scores(
							$start,
							$stop,
							$strand, 
							$strandedness, 
							@features
						);
						$dataset_type = 'wig';
					}
				}
				else {
					croak " Wiggle support is not enabled! " . 
						"Is Bio::Graphics::Wiggle installed?\n";
				}
			}
			
			## BigWig file
			elsif ($wigfile =~ /\.bw$/) {
				# data is in bigwig format
				# this uses the Bio::DB::BigWig adaptor
				
				# collect the wigfile paths
				# also check strand while we're at it
				my @wigfiles;
				while ($feature) {
					
					# check if we can take this feature
					if (
						$strandedness eq 'all' # stranded data not requested
						or $feature->strand == 0 # unstranded data
						or ( 
							# sense data
							$strand == $feature->strand 
							and $strandedness eq 'sense'
						) 
						or (
							# antisense data
							$strand != $feature->strand  
							and $strandedness eq 'antisense'
						)
					) {
						# we can take this file, it passes the strand test
						my ($file) = $feature->get_tag_values('wigfile');
						push @wigfiles, "file:$file";
					}
					
					# prepare for next
					$feature = $iterator->next_seq || undef;
				}
				
				# if no wigfiles are found, return empty handed
				# should only happen if the strands don't match
				return _return_null($method) unless (@wigfiles);
				
				# check that we have bigwig support
				if ($BIGWIG_OK) {
					
					# the data collection depends on the method
					if ($value_type eq 'score' and 
						$method =~ /min|max|mean|sum|count/
					) {
						# we can use the low-level, super-speedy, summary method 
						# warn " using collect_bigwig_score() from tag\n";
						return collect_bigwig_score(
							$chromo,
							$start,
							$stop,
							$method,
							@wigfiles
						);
					}
					
					elsif ($value_type eq 'count' and $method eq 'sum') {
						# we can use the low-level, super-speedy, summary method 
						# warn " using collect_bigwig_score() from tag\n";
						return collect_bigwig_score(
							$chromo,
							$start,
							$stop,
							'count', # special method
							@wigfiles
						);
					}
					
					elsif ($method eq 'indexed') {
						# warn " using collect_bigwig_position_scores() from tag\n";
						return collect_bigwig_position_scores(
							$chromo,
							$start,
							$stop,
							@wigfiles
						);
					}
					
					else {
						# use the longer region collection method
						# warn " using collect_bigwig_scores() from tag\n";
						@scores = collect_bigwig_scores(
							$chromo,
							$start,
							$stop,
							@wigfiles
						);
						$dataset_type = 'bw';
					}
				}
				else {
					croak " BigWig support is not enabled! " . 
						"Is Bio::DB::BigWig installed?\n";
				}
			}
			else {
				croak " Unrecognized wigfile attribute '$wigfile'!" . 
					" Unable to continue!\n";
			}
		}
		
		
		## BigWig Data
		elsif ( $feature->has_tag('bigwigfile') ) {
			# data is in bigwig format
			# this uses the Bio::DB::BigWig adaptor
			
			# collect the wigfile paths
			# also check strand while we're at it
			my @wigfiles;
			while ($feature) {
				
				# check if we can take this feature
				if (
					$strandedness eq 'all' # stranded data not requested
					or $feature->strand == 0 # unstranded data
					or ( 
						# sense data
						$strand == $feature->strand 
						and $strandedness eq 'sense'
					) 
					or (
						# antisense data
						$strand != $feature->strand  
						and $strandedness eq 'antisense'
					)
				) {
					# we can take this file, it passes the strand test
					my ($file) = $feature->get_tag_values('bigwigfile');
					push @wigfiles, "file:$file";
				}
				
				# prepare for next
				$feature = $iterator->next_seq || undef;
			}
			
			# if no wigfiles are found
			# should only happen if the strands don't match
			return _return_null($method) unless (@wigfiles);
			
			# check that we have bigwig support
			if ($BIGWIG_OK) {
				# get the dataset scores using tim_db_helper::bigwig
				
				# the data collection depends on the method
				if ($value_type eq 'score' and 
					$method =~ /min|max|mean|sum|count/
				) {
					# we can use the low-level, super-speedy, summary method 
					return collect_bigwig_score(
						$chromo,
						$start,
						$stop,
						$method,
						@wigfiles
					);
				}
				
				elsif ($value_type eq 'count' and $method eq 'sum') {
					# we can use the low-level, super-speedy, summary method 
					return collect_bigwig_score(
						$chromo,
						$start,
						$stop,
						'count', # special method
						@wigfiles
					);
				}
				
				elsif ($method eq 'indexed') {
					return collect_bigwig_position_scores(
						$chromo,
						$start,
						$stop,
						@wigfiles
					);
				}
				
				else {
					# use the longer region collection method
					@scores = collect_bigwig_scores(
						$chromo,
						$start,
						$stop,
						@wigfiles
					);
					$dataset_type = 'bw';
				}
			}
			else {
				croak " BigWig support is not enabled! " . 
					"Is Bio::DB::BigWig installed?\n";
			}
		}
		
		
		## BigBed Data
		elsif ( $feature->has_tag('bigbedfile') ) {
			# data is in bigbed format
			# this uses the Bio::DB::BigBed adaptor
			
			# collect the bedfile paths
			my @bedfiles;
			while ($feature) {
				my ($file) = $feature->get_tag_values('bigbedfile');
				push @bedfiles, "file:$file";
				
				# prepare for next
				$feature = $iterator->next_seq || undef;
			}
			
			# check that we have bigbed support
			if ($BIGBED_OK) {
				# get the dataset scores using tim_db_helper::bigbed
				
				if ($method eq 'indexed') {
					# warn " using collect_bigbed_position_scores() from tag\n";
					return collect_bigbed_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type,
						@bedfiles
					);
				}
				else {
					# warn " using collect_bigbed_scores() from tag\n";
					@scores = collect_bigbed_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type,
						@bedfiles
					);
					$dataset_type = 'bb';
				}
			}
			else {
				croak " BigBed support is not enabled! " . 
					"Is Bio::DB::BigBed installed?\n";
			}
		}
		
		# Bam Data
		elsif ( $feature->has_tag('bamfile') ) {
			# data is in bam format
			# this uses the Bio::DB::Sam adaptor
			
			# collect the bamfile paths
			my @bamfiles;
			while ($feature) {
				my ($file) = $feature->get_tag_values('bamfile');
				push @bamfiles, "file:$file";
				
				# prepare for next
				$feature = $iterator->next_seq || undef;
			}
			
			# check that we have bam support
			if ($BAM_OK) {
				# get the dataset scores using tim_db_helper::bigbed
				
				if ($method eq 'indexed') {
					# warn " using collect_bam_position_scores() from tag\n";
					return collect_bam_position_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type,
						@bamfiles
					);
				}
				else {
					# warn " using collect_bam_scores() from tag\n";
					@scores = collect_bam_scores(
						$chromo,
						$start,
						$stop,
						$strand, 
						$strandedness, 
						$value_type,
						@bamfiles
					);
					$dataset_type = 'bam';
				}
			}
			else {
				croak " Bam support is not enabled! " . 
					"Is Bio::DB::Sam installed?\n";
			}
		}
		
		
		## Database Data
		else {
			# Working with data stored directly in the database
			# this is more straight forward in collection
			
			# Walk through the datapoints
			# warn " using database\n";
			while ($feature) {
			
				# Check which data to take based on strand
				if (
					$strandedness eq 'all' # all data is requested
					or $strand == 0 # region is unstranded
					or $feature->strand == 0 # unstranded data
					or ( 
						# sense data
						$strand == $feature->strand 
						and $strandedness eq 'sense'
					) 
					or (
						# antisense data
						$strand != $feature->strand  
						and $strandedness eq 'antisense'
					)
				) {
					# we have acceptable data to collect
				
					# data is in the database
					# much easier to collect
					
					# store data in either indexed hash or score array
					if ($method eq 'indexed') {
					
						# determine position to record
						my $position;
						if ($feature->start == $feature->end) {
							# just one position recorded
							$position = $feature->start;
						}
						else {
							# calculate the midpoint
							$position = int( 
								($feature->start + $feature->end) / 2
							);
						}
						
						# store the appropriate value
						if ($value_type eq 'score') {
							push @{ $pos2data{$position} }, $feature->score;
						}
						elsif ($value_type eq 'count') {
							$pos2data{$position} += 1;
						}
						elsif ($value_type eq 'length') {
							push @{ $pos2data{$position} }, $feature->length;
						}
					}
					
					else {
						# just store the score in the array
						
						# store the appropriate value
						if ($value_type eq 'score') {
							push @scores, $feature->score;
						}
						elsif ($value_type eq 'count') {
							push @scores, 1;
						}
						elsif ($value_type eq 'length') {
							push @scores, $feature->length;
						}
					}
				}
				
				# prepare for next
				$feature = $iterator->next_seq || undef;
			}
			
			# post-process the collected position->score values 
			# combine multiple values recorded at the same position
			if (
				$method eq 'indexed' and 
				($value_type eq 'score' or $value_type eq 'length')
			) {
				# each 'value' is an array of one or more scores or lengths 
				# from the datapoints collected above
				# the mean value is the best we can do right now for 
				# combining the data
				# really would prefer something else
				# we don't have a true method to utilize
				foreach my $position (keys %pos2data) {
					$pos2data{$position} = mean( @{$pos2data{$position}} );
				}
			}
			
			$dataset_type = 'db';
		} 
	
	} # end database collection
	
	
	
	### Determine region score from collected scores
	# We have collected the positioned scores
	# Now return the appropriate values
	
	# indexed scores
	if ($method eq 'indexed') {
		# requested indexed position scores
		# we will simply return the data hash
		# regardless whether there are scores or not
		return %pos2data;
	}
	
	# single region score
	else {
		# requested a single score for this region
		# we need to combine the data
		my $region_score;
		
		# first deal with log2 values if necessary
		if ($log) {
			@scores = map {2 ** $_} @scores;
		}
		
		# check that we have scores
		return _return_null($method) unless (@scores);
		
		# determine the region score according to method
		# we are using subroutines from Statistics::Lite
		if ($method eq 'median') {
			# take the median value
			$region_score = median(@scores);
		}
		elsif ($method eq 'mean') {
			# or take the mean value
			$region_score = mean(@scores);
		} 
		elsif ($method eq 'range') {
			# or take the range value
			# this is 'min-max'
			$region_score = range(@scores);
		}
		elsif ($method eq 'stddev') {
			# or take the standard deviation value
			# we are using the standard deviation of the population, 
			# since these are the only scores we are considering
			$region_score = stddevp(@scores);
		}
		elsif ($method eq 'min') {
			# or take the minimum value
			$region_score = min(@scores);
		}
		elsif ($method eq 'max') {
			# or take the maximum value
			$region_score = max(@scores);
		}
		elsif ($method eq 'count') {
			# count the number of values
			$region_score = scalar(@scores);
		}
		elsif ($method eq 'sum') {
			# sum the number of values
			$region_score = sum(@scores);
		}
		elsif ($method =~ /rpk?m/) {
			# convert to reads per million mapped
			# this is only supported by bam and bigbed db, checked above
			
			# total the number of reads if necessary
			unless (exists $total_read_number{$dataset} ) {
				
				# check the type of database
				if ($dataset_type eq 'bam') {
					# a bam database
					
					$total_read_number{$dataset} = 
						sum_total_bam_alignments($dataset);
					print "\n [total alignments: ", 
						format_with_commas( $total_read_number{$dataset} ), 
						"]\n";
				}
				
				elsif ($dataset_type eq 'bb') {
					# bigBed database
					
					$total_read_number{$dataset} = 
						sum_total_bigbed_features($dataset);
					print "\n [total features: ", 
						format_with_commas( $total_read_number{$dataset} ), 
						"]\n";
				}
			}	
			
			# calculate the region score according to the method
			if ($method eq 'rpkm') {
				$region_score = 
					( sum(@scores) * 10^9 ) / 
					( ($stop - $start + 1) * $total_read_number{$dataset} );
			}
			elsif ($method eq 'rpm') {
				$region_score = 
					( sum(@scores) * 10^6 ) / $total_read_number{$dataset};
			}
			else {
				# this dataset doesn't support rpm methods
				# use the sum method instead
				$region_score = sum(@scores);
			}
		}
		else {
			# somehow bad method snuck past our checks
			confess " unrecognized method '$method'!";
		}
	
		# convert back to log2 if necessary
		if ($log) { 
			$region_score = log($region_score) / log(2);
		}
		
		# finished
		return $region_score;
	}
}


=item _return_null

Internal method for returning a 0 or internal null '.' character based 
on the method being used.

=cut

sub _return_null {
	my $method = shift;
	
	# return based on the method
	if ($method eq 'sum') { 
		return 0;
	}
	elsif ($method eq 'count') { 
		return 0;
	}
	else {
		# internal null value
		return '.';
	}
}




=back

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

