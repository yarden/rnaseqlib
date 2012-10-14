package tim_db_helper::bigbed;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::BigBed;
our $VERSION = '1.9.0';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigbed_scores
	collect_bigbed_position_scores
	open_bigbed_db
	sum_total_bigbed_features
);

# Hashes of opened file objects
our %OPENED_BEDFILES; # opened bigbed file objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????

# Hash of Bigfile chromosomes
our %BIGBED_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BIGBED_CHROMOS{bigfile}{chromos}

# The true statement
1; 



### Modules ###



### Collect BigBed scores only
sub collect_bigbed_scores {
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect BigBed scores!\n";
	}
	my ($chromo, $start, $stop, $strand, $stranded, $method, @bed_features) = @_;
		# method can be score, count, or length
	
	# initialize the score array
	# this will record score, count, or lengths per the method
	my @scores;
	
	# look at each bedfile
	# usually there is only one, but for stranded data there may be 
	# two bedfiles (+ and -), so we'll check each wig file for strand info
	foreach my $bedfile (@bed_features) {
	
		# open the bedfile
		my $bb = open_bigbed_db($bedfile) or 
			croak "Unable to open bigBed file '$bedfile'! $!\n";
			
		# first check that the chromosome is present
		unless (exists $BIGBED_CHROMOS{$bedfile}{$chromo}) {
			next;
		}
		
		# collect the features overlapping the region
		my $bb_stream = $bb->features(
			-seq_id   => $chromo, 
			-start    => $start, 
			-end      => $stop,
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$stranded eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$strand == $bed->strand 
					and $stranded eq 'sense'
				) 
				or (
					# antisense data
					$strand != $bed->strand  
					and $stranded eq 'antisense'
				)
			) {
				# we have acceptable data to collect
			
				# store the appropriate datapoint
				if ($method eq 'score') {
					push @scores, $bed->score;
				}
				elsif ($method eq 'count') {
					push @scores, 1;
				}
				elsif ($method eq 'length') {
					push @scores, $bed->length;
				}
			}
		}
	}

	# return collected data
	return @scores;
}




### Collect positioned BigBed scores
sub collect_bigbed_position_scores {
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect BigBed position scores!\n";
	}
	my ($chromo, $start, $stop, $strand, $stranded, $method, @bed_features) = @_;
		# method can be score, count, or length
	
	# set up hash, either position => count or position => [scores]
	my %bed_data;
	
	# look at each bedfile
	# usually there is only one, but for stranded data there may be 
	# two bedfiles (+ and -), so we'll check each wig file for strand info
	foreach my $bedfile (@bed_features) {
	
		# check for opened bedfile
		my $bb = open_bigbed_db($bedfile) or 
			croak "Unable to open bigBed file '$bedfile'! $!\n";
			
		# first check that the chromosome is present
		unless (exists $BIGBED_CHROMOS{$bedfile}{$chromo}) {
			next;
		}
		
		# collect the features overlapping the region
		my $bb_stream = $bb->features(
			-seq_id   => $chromo, 
			-start    => $start, 
			-end      => $stop,
			-iterator => 1,
		);
		
		# process each feature
		while (my $bed = $bb_stream->next_seq) {
			
			# First check whether the strand is acceptable
			if (
				$stranded eq 'all' # all data is requested
				or $bed->strand == 0 # unstranded data
				or ( 
					# sense data
					$strand == $bed->strand 
					and $stranded eq 'sense'
				) 
				or (
					# antisense data
					$strand != $bed->strand  
					and $stranded eq 'antisense'
				)
			) {
				# we have acceptable data to collect
			
				# determine position to record
				my $position;
				if ($bed->start == $bed->end) {
					# just one position recorded
					$position = $bed->start;
				}
				else {
					# calculate the midpoint
					$position = int( 
						( ($bed->start + $bed->end) / 2) + 0.5
					);
				}
				
				# check the position
				next unless (
					# want to avoid those whose midpoint are not technically 
					# within the region of interest
					$position >= $start and $position <= $stop
				);
				
				# store the appropriate datapoint
				# for score and length, we're putting these into an array
				if ($method eq 'score') {
					# perform addition to force the score to be a scalar value
					push @{ $bed_data{$position} }, $bed->score + 0;
				}
				elsif ($method eq 'count') {
					$bed_data{$position} += 1;
				}
				elsif ($method eq 'length') {
					# I hope that length is supported, but not sure
					# may have to calculate myself
					push @{ $bed_data{$position} }, $bed->length;
				}
			}
		}
	}

	# combine multiple datapoints at the same position
	if ($method eq 'score' or $method eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %bed_data) {
			$bed_data{$position} = mean( @{$bed_data{$position}} );
		}
	}
	
	# return collected data
	return %bed_data;
}



### Open a bigBed database connection
sub open_bigbed_db {
	
	my $bedfile = shift;
	
	# check whether the file has been opened or not
	if (exists $OPENED_BEDFILES{$bedfile} ) {
		# this file is already opened, use it
		return $OPENED_BEDFILES{$bedfile};
	}
	
	else {
		# this file has not been opened yet, open it
		my $path = $bedfile;
		$path =~ s/^file://; # clean up file prefix if present
		my $bb;
		eval {
			$bb = Bio::DB::BigBed->new($path);
		};
		return unless $bb;
		
		# store the opened object for later use
		$OPENED_BEDFILES{$bedfile} = $bb;
		
		# collect the chromosomes for this bigBed file
		%{ $BIGBED_CHROMOS{$bedfile} } = map { $_ => 1 } $bb->seq_ids;
		
		return $bb;
	}
}



### Sum the total number of features in the bigBed file
sub sum_total_bigbed_features {
	
	# Passed arguments;
	my $bb_file = shift;
	unless ($bb_file) {
		carp " no BigBed file or BigBed db object passed!\n";
		return;
	}
	
	
	# Open BigBed file if necessary
	my $bed;
	my $bb_ref = ref $bb_file;
	if ($bb_ref =~ /Bio::DB::BigBed/) {
		# we have an opened bigbed db object
		$bed = $bb_file;
	}
	else {
		# we have a name of a sam file
		$bed = open_bigbed_db($bb_file);
		return unless ($bed);
	}
	
	# Count the number of alignments
	my $total_read_number = 0;
	
	# loop through the chromosomes
	my @chroms = $bed->seq_ids;
	foreach my $chrom (@chroms) {
		
		# use an iterator to speed things up
		my $iterator = $bed->features(
			-seq_id    => $chrom,
			-iterator  => 1,
		);
		
		# count the number of bed features we go through
		while (my $f = $iterator->next_seq) {
			# we don't actually do anything with them
			$total_read_number++;
		}
	}
	
	return $total_read_number;
}



__END__

=head1 NAME

tim_db_helper::bigbed

=head1 DESCRIPTION

This module supports the use of bigBed file in the biotoolbox scripts.
It is used to collect the dataset scores from a binary 
bigBed file (.bb). The file may be local or remote.

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the bed feature, then the data for that bed 
feature is collected.  

For loading bigbed files into a Bio::DB database, see the biotoolbox perl 
script 'big_filegff3.pl'.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened bigbed file object is stored in a global 
hash in case it is needed again.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::BigBed to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::bigbed;

It will automatically export the name of the subroutines. 

=over

=item collect_bigbed_scores

This subroutine will collect only the data values from a binary bigbed file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed seven or more arguments in the following order:
    
    1) The chromosome or seq_id
    2) The start position of the segment to collect 
    3) The stop or end position of the segment to collect 
    4) The strand of the original feature (or region), -1, 0, or 1.
    5) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       or "all". Only those scores which match the indicated 
       strandedness are collected.
    6) The method or type of data collected. 
       Acceptable values include 'score' (returns the bed feature 
       score), 'count' (returns the number of bed features found), or 
       'length' (returns the length of the bed features found).  
    7) The paths, either local or remote, to one or more BigBed files.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bigbed_position_scores

This subroutine will collect the score values from a binary bigBed file 
for the specified database region keyed by position. 

The subroutine is passed the same arguments as collect_bigbed_scores().

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for score or length data methods) or sum 
(for count methods) is returned.

=item open_bigbed_db()

This subroutine will open a BigBed database connection. Pass either the 
local path to a bigBed file (.bb extension) or the URL of a remote bigBed 
file. It will return the opened database object.

=item sum_total_bigbed_features()

This subroutine will sum the total number of bed features present in a 
BigBed file. This may be useful, for example, in calculating fragments 
(reads) per million mapped values when the bigbed file represents 
sequence alignments.

Pass either the name of a bigBed file (.bb), either local or remote, or an 
opened BigBed database object. A scalar value of the total number of features 
is returned.

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



