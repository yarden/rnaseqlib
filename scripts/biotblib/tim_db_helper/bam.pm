package tim_db_helper::bam;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(mean);
use Bio::DB::Sam;
our $VERSION = '1.8.4';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bam_scores
	collect_bam_position_scores
	open_bam_db
	sum_total_bam_alignments
);

# Hashes of opened file objects
our %OPENED_BAMFILES; # opened bam file objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????

# Hash of Bigfile chromosomes
our %BAM_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BAM_CHROMOS{bigfile}{chromos}


# The true statement
1; 



### Modules ###



### Collect Bam scores only
sub collect_bam_scores {
	
	# set the do_index boolean to false
	# return the scores
	return _collect_bam_data(0, @_);
}




### Collect positioned Bam scores
sub collect_bam_position_scores {
	
	# collect the raw data
	# set the do_index boolean to true
	my %bam_data = _collect_bam_data(1, @_);
	
	# grab the method from the passed arguments
	my $method = $_[3];
	
	# combine multiple datapoints at the same position
	if ($method eq 'length') {
		# each value is an array of one or more datapoints
		# we will take the simple mean
		foreach my $position (keys %bam_data) {
			$bam_data{$position} = mean( @{$bam_data{$position}} );
		}
	}
	
	# return collected data
	return %bam_data;
}




### Actual collection of scores
sub _collect_bam_data {
	
	# pass the required information
	unless (scalar @_ >= 8) {
		confess " At least eight arguments must be passed to collect Bam data!\n";
	}
	my (
		$do_index, 
		$chromo,
		$start,
		$stop,
		$strand, 
		$stranded, 
		$value_type, 
		@bam_features
	) = @_;
		# method can be score, count, or length
	
	# initialize score structures
	# which one is used depends on the $do_index boolean variable
	my %pos2data; # either position => count or position => [scores]
	my @scores; # just scores
	
	# look at each bamfile
	# usually there is only one, but there may be more than one
	foreach my $bamfile (@bam_features) {
	
		## Open the Bam File
		my $bam = open_bam_db($bamfile);
			
		# first check that the chromosome is present
		unless (exists $BAM_CHROMOS{$bamfile}{$chromo}) {
			next;
		}
		
		
		## Set the code to filter alignments based on strand 
		my $filter;
		if ($stranded eq 'sense' and $strand == 1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == 1 ? 1 : 0;
			};
		}
		elsif ($stranded eq 'sense' and $strand == -1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == -1 ? 1 : 0;
			};
		}
		elsif ($stranded eq 'antisense' and $strand == 1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == -1 ? 1 : 0;
			};
		}
		elsif ($stranded eq 'antisense' and $strand == -1) {
			$filter = sub {
				my $a = shift;
				return $a->strand == 1 ? 1 : 0;
			};
		}
		else {
			# no strand requested, take all
			$filter = sub {
				return 1;
			};
		}
		
		
		## Collect the data according to the requested value type
		# we will either use simple coverage (score method) or
		# process the actual alignments (count or length)
		
		## Coverage
		if ($value_type eq 'score') {
			# collecting scores, or in this case, basepair coverage of 
			# alignments over the requested region
			
			my $coverage;
			if ($stranded eq 'sense' or $stranded eq 'antisense') {
				# Cannot currently collect stranded data with the coverage 
				# method. I will keep the filter in here anyway to future-proof 
				# in case Lincoln ever adds this support (don't hold your 
				# breath!)
				
				($coverage) = $bam->features(
					-type     => 'coverage',
					-seq_id   => $chromo,
					-start    => $start,
					-end      => $stop,
					-filter   => $filter,
				);
			}
			else {
				# no stranded data wanted
				($coverage) = $bam->features(
					-type     => 'coverage',
					-seq_id   => $chromo,
					-start    => $start,
					-end      => $stop,
				);
			}
			
			# convert the coverage data
			# by default, this should return the coverage at 1 bp resolution
			if ($coverage) {
				@scores = $coverage->coverage;
				
				# check whether we need to index the scores
				if ($do_index) {
					for (my $i = $start; $i <= $stop; $i++) {
						# move the scores into the position score hash
						$pos2data{$i} += shift @scores;
					}
				}
			}
		}
		
		
		## Alignments
		else {
			# either collecting counts or length
			# working with actual alignments
			
			# get the alignments in a stream
			my $iterator;
			if ($stranded eq 'sense' or $stranded eq 'antisense') {
				# include the filter subroutine
				$iterator = $bam->features(
					-type     => 'match',
					-seq_id   => $chromo,
					-start    => $start,
					-end      => $stop,
					-filter   => $filter,
					-iterator => 1,
				);
			}
			else {
				# no stranded data wanted
				$iterator = $bam->features(
					-type     => 'match',
					-seq_id   => $chromo,
					-start    => $start,
					-end      => $stop,
					-iterator => 1,
				);
			}
			
			# process the alignments
				# want to avoid those whose midpoint are not technically 
				# within the region of interest
			if ($do_index) {
				# record position and score
				
				if ($value_type eq 'count') {
					while (my $a = $iterator->next_seq) {
						
						# enumerate at the alignment's midpoint
						my $position = 
							int( ( ($a->start + $a->end) / 2) + 0.5);
						
						# check midpoint is within the requested region
						if (
							$position >= $start and 
							$position <= $stop
						) {
							$pos2data{$position} += 1;
						}
					}
				}
				elsif ($value_type eq 'length') {
					while (my $a = $iterator->next_seq) {
						
						# record length at the alignment's midpoint
						my $position = 
							int( ( ($a->start + $a->end) / 2) + 0.5);
						
						# check midpoint is within the requested region
						if (
							$position >= $start and 
							$position <= $stop
						) {
							push @{ $pos2data{$position} }, $a->length;
						}
					}
				}
			}
			
			else {
				# record scores only, no position
				
				if ($value_type eq 'count') {
					while (my $a = $iterator->next_seq) {
						# enumerate 
						push @scores, 1;
					}
				}
				elsif ($value_type eq 'length') {
					while (my $a = $iterator->next_seq) {
						# record length
						push @scores, $a->length;
					}
				}
			
			}
		
		}
	}

	
	## Return collected data
	if ($do_index) {
		return %pos2data;
	}
	else {
		return @scores;
	}
}


### Open a bigWig database connection
sub open_bam_db {
	
	my $bamfile = shift;
	
	# check if we have seen this bam file before
	if (exists $OPENED_BAMFILES{$bamfile} ) {
		# this file is already opened, use it
		return $OPENED_BAMFILES{$bamfile};
	}
	
	else {
		# this file has not been opened yet, open it
		
		# check the path
		my $path = $bamfile;
		$path =~ s/^file://; # strip the file prefix if present
		
		# open the bam database object
		my $bam;
		eval {
			$bam = Bio::DB::Sam->new(
					-bam         => $path,
					-autoindex   => 1,
			);
		};
		return unless $bam;
		
		# store the opened object for later use
		$OPENED_BAMFILES{$bamfile} = $bam;
			
		# collect the chromosomes for this bam
		%{ $BAM_CHROMOS{$bamfile} } = map { $_ => 1 } $bam->seq_ids;
		
		# done
		return $bam;
	}
}



### Determine total number of alignments in a bam file
sub sum_total_bam_alignments {
	
	# Passed arguments;
	my $sam_file = shift;
	my $min_mapq = shift || 0; # by default we take all alignments
	my $paired   = shift || 0; # by default we assume all alignments are single-end
	unless ($sam_file) {
		carp " no Bam file or bam db object passed!\n";
		return;
	}
	
	
	# Open Bam file if necessary
	my $sam;
	my $sam_ref = ref $sam_file;
	if ($sam_ref =~ /Bio::DB::Sam/) {
		# we have an opened sam db object
		$sam = $sam_file;
	}
	else {
		# we have a name of a sam file
		$sam = open_bam_db($sam_file);
		return unless ($sam);
	}
	
	# Count the number of alignments
	my $total_read_number = 0;
	
	# loop through the chromosomes
	for my $tid (0 .. $sam->n_targets - 1) {
		# each chromosome is internally represented in the bam file as 
		# a numeric target identifier
		# we can easily convert this to an actual sequence name
		# we will force the conversion to go one chromosome at a time
		
		# sequence name
		my $seq_id = $sam->target_name($tid);
		
		# process the reads according to single or paired-end
		# paired end alignments
		if ($paired) {
			$sam->fetch($seq_id, 
				sub {
					my $a = shift;
					
					# check paired alignment
					return if $a->unmapped;
					return unless $a->proper_pair;
					return if $a->qual < $min_mapq;
					
					# we're only counting forward reads of a pair 
					return if $a->strand != 1; 
					
					# count this fragment
					$total_read_number++;
				}
			);
		}
		
		# single end alignments
		else {
			$sam->fetch($seq_id, 
				sub {
					my $a = shift;
					
					# check paired alignment
					return if $a->unmapped;
					return if $a->qual < $min_mapq;
					
					# count this fragment
					$total_read_number++;
				}
			);
		}
	}
	
	# done
	return $total_read_number;
}




__END__




=head1 NAME

tim_db_helper::bam

=head1 DESCRIPTION

This module is used to collect the dataset scores from a binary 
bam file (.bam) of alignments. Bam files may be local or remote, 
and are usually prefixed with 'file:', 'http://', of 'ftp://'.

Collected data values may be restricted to strand by specifying the desired 
strandedness, 
depending on the method of data collection. Collecting scores, or basepair 
coverage of alignments over the region of interest, does not currently support 
stranded data collection (as of this writing). However, enumerating 
alignments (count method) and collecting alignment lengths do support 
stranded data collection.

Currently, paired-end bam files are treated as single-end files. There are 
some limitations regarding working with paired-end alignments that don't 
work well (search, strand, length, etc). If paired-end alignments are to 
be analyzed, they should be processed into another format (BigWig or BigBed). 
See the biotoolbox scripts 'bam2gff_bed.pl' or 'bam2wig.pl' for solutions.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened bam file object is stored in a global 
hash in case it is needed again.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::Sam to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::bam;

It will automatically export the name of the subroutines. 

=over

=item collect_bam_scores

This subroutine will collect only the data values from a binary bam file 
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
    6) The type of data collected. 
       Acceptable values include 'score' (returns the basepair coverage
       of alignments over the region of interest), 'count' (returns the 
       number of alignments found at each base position in the region, 
       recorded at the alignment's midpoint), or 'length' (returns the 
       mean lengths of the alignments found at each base position in 
       the region, recorded at the alignment's midpoint). 
    7) The paths, either local or remote, to one or more Bam files.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bam_position_scores

This subroutine will collect the score values from a binary bam file 
for the specified database region keyed by position. 

The subroutine is passed the same arguments as collect_bam_scores().

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. The feature midpoint is used 
as the key position. When multiple features are found at the same 
position, a simple mean (for length data methods) or sum 
(for count methods) is returned.

=item open_bam_db()

This subroutine will open a Bam database connection. Pass either the 
local path to a Bam file (.bam extension) or the URL of a remote Bam 
file. It will return the opened database object.

=item sum_total_bam_alignments()

This subroutine will sum the total number of properly mapped alignments 
in a bam file. Pass the subroutine one to three arguments. 
    
    1) The name of the Bam file which should be counted. Alternatively,  
       an opened Bio::DB::Sam object may also be given. Required.
    2) Optionally pass the minimum mapping quality of the reads to be 
       counted. The default is 0, or all alignments are counted.
    3) Optionally pass a boolean value (1 or 0) indicating whether 
       the Bam file represents paired-end alignments. Only proper 
       alignment pairs are counted. The default is to treat all 
       alignments as single-end.
       
The subroutine will return the number of alignments.

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

