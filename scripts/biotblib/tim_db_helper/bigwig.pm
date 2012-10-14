package tim_db_helper::bigwig;

# modules
require Exporter;
use strict;
use Carp;
use Statistics::Lite qw(min max mean);
use Bio::DB::BigWig qw(binMean binStdev);
use Bio::DB::BigWigSet;
our $VERSION = '1.9.0';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_bigwig_score
	collect_bigwig_scores
	collect_bigwig_position_scores
	collect_bigwigset_score
	collect_bigwigset_scores
	collect_bigwigset_position_scores
	open_bigwig_db
	open_bigwigset_db
);


# Hashes of opened file objects
our %OPENED_BIGFILES; # opened bigwig file objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????

# Hash of Bigfile chromosomes
our %BIGWIG_CHROMOS;
	# sometimes user may request a chromosome that's not in the bigfile
	# that could lead to an exception
	# we will record the chromosomes list in this hash
	# $BIGWIG_CHROMOS{bigfile}{chromos}

# The true statement
1; 



### Modules ###

### Collect single BigWig score
sub collect_bigwig_score {
	
	# Pass the required information
	unless (scalar @_ >= 5) {
		confess " At least five arguments must be passed to collect BigWig score!\n";
	}
	my ($chromo, $start, $stop, $method, @wig_files) = @_;
	
	# Confirm the method 
	unless ($method =~ /^min|max|mean|count|sum|stddev$/) {
		confess " can not use method $method!\n";
	}
	
	# Collecting summary features
	# we will collect a summary object for each requested wig feature  
	my @summaries;
	
	# Walk through each requested feature
	# There is likely only one
	foreach my $wig (@wig_files) {
		
		# open a new db object
		my $bw = open_bigwig_db($wig) or 
			croak " Unable to open bigWig file '$wig'! $!\n";
		
		# first check that the chromosome is present
		if (exists $BIGWIG_CHROMOS{$wig}{$chromo}) {
			my (@features) = $bw->features(
				-seq_id    => $chromo,
				-start     => $start,
				-end       => $stop,
				-type      => 'summary',
			);
			
			# check
			confess(' BigWig database returned ' . scalar(@features) . 
				'summary features! Expected only 1!\n') 
				if scalar(@features) > 1;
			
			# keep the summary
			push @summaries, $features[0]->score;
		}
	}
	
	# now process the summary features
	return _process_summaries($method, @summaries);
}





### Collect multiple BigWig scores
sub collect_bigwig_scores {
	
	# pass the required information
	unless (scalar @_ >= 4) {
		confess " At least four arguments must be passed to collect BigWig scores!\n";
	}
	my ($chromo, $start, $stop, @wig_files) = @_;
	
	# initialize 
	my @scores; 
	
	# Walk through each requested feature
	# There is likely only one
	foreach my $wig (@wig_files) {
	
		# Open the BigWig file
		my $bw = open_bigwig_db($wig) or 
			croak " Unable to open bigWig file '$wig'! $!\n";;
		
		# first check that the chromosome is present
		unless (exists $BIGWIG_CHROMOS{$wig}{$chromo}) {
			next;
		}
		
		# initialize a feature stream for this segment
		my $iterator = $bw->features(
			-seq_id     => $chromo,
			-start      => $start,
			-end        => $stop,
			-type       => 'region',
			-iterator   => 1,
		);
		
		# collect the scores
		while (my $f = $iterator->next_seq) {
			push @scores, $f->score;
		}
	} 
	
	# return collected data
	return @scores;
}




### Collect positioned BigWig scores
sub collect_bigwig_position_scores {
	
	# pass the required information
	unless (scalar @_ >= 4) {
		confess " At least four arguments must be passed to collect BigWig position scores!\n";
	}
	my ($chromo, $start, $stop, @wig_files) = @_;
	
	# initialize 
	my %pos2data; # hash of position => score
	my %duplicates; # hash of duplicate positions, position => number
	
	# Walk through each requested feature
	# There is likely only one
	foreach my $wig (@wig_files) {
	
		# Open the BigWig file
		my $bw = open_bigwig_db($wig) or 
			croak " Unable to open bigWig file '$wig'! $!\n";;
		
		# first check that the chromosome is present
		unless (exists $BIGWIG_CHROMOS{$wig}{$chromo}) {
			next;
		}
		
		# initialize a feature stream for this segment
		my $iterator = $bw->features(
			-seq_id     => $chromo,
			-start      => $start,
			-end        => $stop,
			-type       => 'region',
			-iterator   => 1,
		);
		
		# collect the scores
		while (my $f = $iterator->next_seq) {
			# process the feature
			_process_position_score_feature($f, \%pos2data, \%duplicates);
		}
	} 
	
	# check for duplicate positions
	if (%duplicates) {
		_remove_duplicate_positions(\%pos2data, \%duplicates);
	}
	
	# return collected data
	return %pos2data;
}



sub collect_bigwigset_score {
	
	# Pass the required information
	unless (scalar @_ >= 8) {
		confess " At least eight arguments must be passed to collect BigWigSet score!\n";
	}
	my ($db, $chromo, $start, $stop, $strand, $stranded, $method, @types) = @_;
	
	# Confirm the method 
	unless ($method =~ /^min|max|mean|count|sum|stddev$/) {
		confess " can not use method $method!\n";
	}
	
	# Confirm the chromosome
	my $first_path = ($db->bigwigs)[0];
	unless ( exists $BIGWIG_CHROMOS{$first_path}{$chromo} ) {
		# chromosome is not present, at least in the first bigwig in the set
		# return null
		return $method =~ /sum|count/ ? 0 : '.';
	}
	
	# Reset which feature_type to collect from the database
	# this is normally set to region when we opened the bigwigset db
	# but it may be changed by a previous method
	# we now want summary feature_type
	$db->feature_type('summary');
	
	# Collecting summary features
	# we will collect a summary object for each requested wig feature  
	my @summaries;
	
	# Work through all the requested feature types
	# the BigWigSet feature request doesn't work well with multiple features
	# so we'll do them one at a time
	foreach my $type (@types) {
	
		# Collect the summary features
		$type =~ s/\:.+$//; # strip the source if present
			# features method only works with primary_tag, not full 
			# primary_tag:source type
		my @features = $db->features(
			-seq_id   => $chromo,
			-start    => $start,
			-end      => $stop,
			-type     => $type,
		);
		# if the type doesn't work, then try display_name instead
		unless (@features) {
			@features = $db->features(
				-seq_id   => $chromo,
				-start    => $start,
				-end      => $stop,
				-name     => $type,
			);
		}
			# since we're collecting summary features, we will only get one 
			# per bigwig file that matches the request
			# no need for a seqfeature stream
		
		# Determine which features to take based on strandedness
		
		# Stranded features
		if (
			$strand != 0 and
			($stranded eq 'sense' or $stranded eq 'antisense')
		) {
			# we will have to check the strand for each object
			# feature objects we collect don't have the standard strand set
			# instead, we will have to get the attribute tag named strand
			
			# check each feature
			foreach my $f (@features) {
				
				# get the feature strand
				my $fstrand = 0; # default
				if ($f->has_tag('strand') ) {
					($fstrand) = $f->get_tag_values('strand');
				}
				
				# collect summary if strand is appropriate
				if (
					$fstrand == 0 or
					(
						# sense data
						$strand == $fstrand 
						and $stranded eq 'sense'
					) 
					or (
						# antisense data
						$strand != $fstrand  
						and $stranded eq 'antisense'
					)
				) {
					# we have acceptable data to collect
					push @summaries, $f->score;
				}
			}
		}
		
		# Non-stranded features
		else {
			# take all the features found
			
			# keep all the summaries
			foreach my $f (@features) {
				push @summaries, $f->score;
			}
		}
	}
	
	# now process the summary features
	return _process_summaries($method, @summaries);
}




sub collect_bigwigset_scores {
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect BigWigSet scores!\n";
	}
	my ($db, $chromo, $start, $stop, $strand, $stranded, @types) = @_;
	
	# Confirm the chromosome
	my $first_path = ($db->bigwigs)[0];
	unless ( exists $BIGWIG_CHROMOS{$first_path}{$chromo} ) {
		# chromosome is not present, at least in the first bigwig in the set
		# return nothing
		return;
	}
	
	# Reset which feature_type to collect from the database
	# this is normally set to region when we opened the bigwigset db
	# but it may be changed by a previous method
	# we now want region feature_type
	$db->feature_type('region');
	
	# initialize collection array
	my @scores; 
	
	# Go through each feature type requested
	foreach my $type (@types) {
		
		# Collect the feature scores
		$type =~ s/\:.+$//; # strip the source if present
			# features method only works with primary_tag, not full 
			# primary_tag:source type
			# since the default feature_type for the bigwigset database is 
			# region, we will get lots of features returned, one for each 
			# datapoint
			# use an iterator to process them
		my $iterator = $db->get_seq_stream(
			-seq_id   => $chromo,
			-start    => $start,
			-end      => $stop,
			-type     => $type,
		);
		
		# check that we have a feature stream
		my $feature = $iterator->next_seq;
		unless ($feature) {
			# uh oh, no feature! perhaps we didn't correctly identify the 
			# correct bigwig in the set
			# try again using display_name instead
			$iterator = $db->get_seq_stream(
				-seq_id   => $chromo,
				-start    => $start,
				-end      => $stop,
				-name     => $type,
			);
			$feature = $iterator->next_seq;
		}
		
		# Determine which features to take based on strandedness
		
		# Stranded features
		if (
			$strand != 0 and
			($stranded eq 'sense' or $stranded eq 'antisense')
		) {
			# we will have to check the strand for each object
			# feature objects we collect don't have the standard strand set
			# instead, we will have to get the attribute tag named strand
			
			# check each feature
			while ($feature) {
				
				# get the feature strand
				my $fstrand = 0; # default
				if ($feature->has_tag('strand') ) {
					($fstrand) = $feature->get_tag_values('strand');
				}
				
				# collect score if strand is appropriate
				if (
					$fstrand == 0 or
					(
						# sense data
						$fstrand == $strand 
						and $stranded eq 'sense'
					) 
					or (
						# antisense data
						$fstrand != $strand  
						and $stranded eq 'antisense'
					)
				) {
					# we have acceptable data to collect
					push @scores, $feature->score;
				}
				
				# prepare for the next feature
				$feature = $iterator->next_seq;
			}
		}
		
		# Non-stranded features
		else {
			# take all the features found
			while ($feature) {
				push @scores, $feature->score;
				$feature = $iterator->next_seq;
			}
		}
	}
	
	# Finished
	return @scores;
}



sub collect_bigwigset_position_scores {
	
	# pass the required information
	unless (scalar @_ >= 7) {
		confess " At least seven arguments must be passed to collect BigWigSet position scores!\n";
	}
	my ($db, $chromo, $start, $stop, $strand, $stranded, @types) = @_;
	
	# Confirm the chromosome
	my $first_path = ($db->bigwigs)[0];
	unless ( exists $BIGWIG_CHROMOS{$first_path}{$chromo} ) {
		# chromosome is not present, at least in the first bigwig in the set
		# return nothing
		return;
	}
	
	# Reset which feature_type to collect from the database
	# this is normally set to region when we opened the bigwigset db
	# but it may be changed by a previous method
	# we now want region feature_type
	$db->feature_type('region');
	
	# initialize collection hash, position => score
	my %pos2data; 
	my %duplicates;
	
	# Go through each feature type requested
	foreach my $type (@types) {
		
		# Collect the feature scores
		$type =~ s/\:.+$//; # strip the source if present
			# features method only works with primary_tag, not full 
			# primary_tag:source type
			# since the default feature_type for the bigwigset database is 
			# region, we will get lots of features returned, one for each 
			# datapoint
			# use an iterator to process them
		my $iterator = $db->get_seq_stream(
			-seq_id   => $chromo,
			-start    => $start,
			-end      => $stop,
			-type     => $type,
		);
		
		# check that we have a feature stream
		my $feature = $iterator->next_seq;
		unless ($feature) {
			# uh oh, no feature! perhaps we didn't correctly identify the 
			# correct bigwig in the set
			# try again using display_name instead
			$iterator = $db->get_seq_stream(
				-seq_id   => $chromo,
				-start    => $start,
				-end      => $stop,
				-name     => $type,
			);
			$feature = $iterator->next_seq;
		}
		
		# Determine which features to take based on strandedness
		
		# Stranded features
		if (
			$strand != 0 and
			($stranded eq 'sense' or $stranded eq 'antisense')
		) {
			# we will have to check the strand for each object
			# feature objects we collect don't have the standard strand set
			# instead, we will have to get the attribute tag named strand
			
			# Check each feature
			
			# Stranded features
			while ($feature) {
				
				# get the feature strand
				my $fstrand = 0; # default
				if ($feature->has_tag('strand') ) {
					($fstrand) = $feature->get_tag_values('strand');
				}
				
				# collect score if strand is appropriate
				if (
					$strand == 0 or
					(
						# sense data
						$fstrand == $strand 
						and $stranded eq 'sense'
					) 
					or (
						# antisense data
						$fstrand != $strand  
						and $stranded eq 'antisense'
					)
				) {
					# acceptable data point
					# process
					_process_position_score_feature(
						$feature, \%pos2data, \%duplicates);
				}
				
				# prepare next
				$feature = $iterator->next_seq;
			}
		}
		
		# Non-stranded features
		else {
			# take all the features found
			while ($feature) {
				# process
				_process_position_score_feature(
					$feature, \%pos2data, \%duplicates);
				$feature = $iterator->next_seq;
			}
		}
	}
	
	# Remove duplicates
	if (%duplicates) {
		_remove_duplicate_positions(\%pos2data, \%duplicates);
	}
	
	
	# Finished
	return %pos2data;
}




### Open a bigWig database connection
sub open_bigwig_db {
	
	my $wigfile = shift;
	my $bw;
	
	# check whether the file has been opened or not
	if (exists $OPENED_BIGFILES{$wigfile} ) {
		# this file is already opened, use it
		$bw = $OPENED_BIGFILES{$wigfile};
	}
	
	else {
		# this file has not been opened yet, open it
		my $path = $wigfile;
		$path =~ s/^file://; # clean up file prefix if present
		eval {
			$bw = Bio::DB::BigWig->new( -bigwig => $path);
		};
		return unless $bw;
		
		# store the opened object for later use
		$OPENED_BIGFILES{$wigfile} = $bw;
		
		# collect the chromosomes for this bigwig
		%{ $BIGWIG_CHROMOS{$wigfile} } = map { $_ => 1 } $bw->seq_ids;
	}
	
	return $bw;
}




### Open a bigWigSet database connection
sub open_bigwigset_db {
	
	my $directory = shift;
	
	# check for trailing slash
	# this seems to interfere with generating the list of files, leading 
	# to duplicates: both raw files as well as contents from metadata
	$directory =~ s/\/$//; # strip trailing slash
	
	# open the database connection 
	# we're using the region feature type because that's what the rest of 
	# tim_db_helper modules expect and work with
	my $bws;
	eval {
		$bws = Bio::DB::BigWigSet->new(
					-dir            => $directory,
					-feature_type   => 'region',
		);
	};
	return unless $bws;
	
	# check that we haven't just opened a new empty bigwigset object
	my @paths = $bws->bigwigs;
	
	if (@paths) {
		# we have bigwig files, must be a valid bigwigset directory
		
		# collect the chromosomes from the first bigwig
		# we will assume all of the bigwigs have the same chromosomes!
		my $bw = $bws->get_bigwig($paths[0]);
		
		# collect the chromosomes for this bigwig
		%{ $BIGWIG_CHROMOS{$paths[0]} } = map { $_ => 1 } $bw->seq_ids;
		
		return $bws;
	}
	else {
		# no valid bigWig files, not valid
		return;
	}
}




### Internal subroutine for processing summary features
# for personal use only
sub _process_summaries {
	my ($method, @summaries) = @_;
	
	## Process the collected summaries
	my $value;
	
	# No summaries
	if (scalar @summaries == 0) {
		# nothing was found!
		
		# return empty handed
		if ($method =~ /sum|count/) {
			$value = 0;
		}
		else {
			# internal null
			$value = '.';
		}
	}
	
	# One summary
	elsif (scalar @summaries == 1) {
		# great! only one summary returned
		
		# return based on the method
		if ($method eq 'mean') {
			$value = binMean( $summaries[0] );
		}
		elsif ($method eq 'sum') {
			$value = $summaries[0]->{sumData};
		}
		elsif ($method eq 'min') {
			$value = $summaries[0]->{minVal};
		}
		elsif ($method eq 'max') {
			$value = $summaries[0]->{maxVal};
		}
		elsif ($method eq 'count') {
			$value = $summaries[0]->{validCount};
		}
		elsif ($method eq 'stddev') {
			$value = binStdev( $summaries[0] );
		}
		else {
			confess " unknown method $method!\n";
		}
	}
	
	# Multiple summaries
	else {
		# more than one summary
		# this will take a little more work
		
		# return based on the method
		if ($method eq 'mean') {
			my $sum = 0;
			my $count = 0;
			foreach my $s (@summaries) {
				$sum += $s->{sumData};
				$count += $s->{validCount};
			}
			$value = $count ? $sum/$count : '.';
		}
		elsif ($method eq 'sum') {
			my $sum = 0;
			foreach my $s (@summaries) {
				$sum += $s->{sumData};
			}
			$value = $sum;
		}
		elsif ($method eq 'min') {
			$value = min( map { $_->{minVal} } @summaries);
		}
		elsif ($method eq 'max') {
			$value = max( map { $_->{maxVal} } @summaries);
		}
		elsif ($method eq 'count') {
			my $count = 0;
			foreach my $s (@summaries) {
				$count += $s->{validCount};
			}
			$value = $count;
		}
		elsif ($method eq 'stddev') {
			carp " can not determine correct stddev value with " . 
				scalar(@summaries) . " bigwig summaries!\n";
			# take the first one only as something to report
			$value = binStdev( $summaries[0] );
		}
		else {
			confess " unknown method $method!\n";
		}
	}
	
	# Done
	unless (defined $value) {
		# one last check to avoid undefined values
		$value = $method =~ /sum|count/ ? 0 : '.';
	}
	return $value;
}



### Internal subroutine for processing position and score from features
# for personal use only
sub _process_position_score_feature {
	
	# collect the feature and hashes
	my ($f, $pos2data, $duplicates) = @_;
	
	# process across the length of this feature
		# for most wig features this is almost certainly a single point
		# but just in case this is spanned data, we will record the 
		# value at every position
	for (my $pos = $f->start; $pos <= $f->end; $pos++) {
		
		# check for duplicate positions
		# duplicates should not normally exist for a single wig file
		# but we may be working with multiple wigs that are being combined
		if (exists $pos2data->{$pos} ) {
			
			if (exists $duplicates->{$pos} ) {
				# we have lots of duplicates at this position!
				
				# append an incrementing number at the end
				$duplicates->{$pos} += 1; # increment first
				my $new = $pos . '.' . $duplicates->{$pos};
				$pos2data->{$new} = $f->score;
			}
			else {
				# first time duplicate
				
				# record this one
				my $new = $pos . '.1';
				$pos2data->{$new} = $f->score;
				
				# remember
				$duplicates->{$pos} = 1;
			}
		}
		else {
			$pos2data->{$pos} = $f->score;
		}
	}
}


### Internal subroutine for removing duplicate positions from pos2data hash
# for personal use only
sub _remove_duplicate_positions {
	
	# collect the feature and hashes
	my ($pos2data, $duplicates) = @_;
	
	# remove the duplicates
		# we will combine all of them with a simple mean, what else to do?
	foreach my $pos (keys %{ $duplicates } ) {
		my $num = $duplicates->{$pos};
		
		# collect all the values
		my @values;
		push @values, $pos2data->{$pos}; # initial value
		for my $i (1..$num) {
			push @values, $pos2data->{ "$pos\.$i" };
			delete $pos2data->{ "$pos\.$i" };
		}
		$pos2data->{$pos} = mean(@values);
	}
}



__END__

=head1 NAME

tim_db_helper::bigwig

=head1 DESCRIPTION

This module supports the use of bigwig file in the biotoolbox scripts.
It is used to collect the dataset scores from a binary 
bigWig file (.bw), or from a directory of bigWig files, known as a 
BigWigSet. BigWig files may be local or remote.

In either case, the file is read using the Bio::DB::BigWig module, and 
the values extracted from the region of interest. 

Stranded data collection is not supported with bigWig files. However, 
since the BigWigSet database supports metadata attributes for each 
included bigWig, it has the potential for collecting stranded data. To 
do so, each bigWig metadata must include the strand attribute. 

For loading bigwig files into a Bio::DB database, see the biotoolbox perl 
script 'big_file2gff3.pl'. This will prepare either a GFF3 file for loading 
into a Bio::DB::SeqFeature::Store database, or a Bio::DB::BigWigSet 
database.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened bigwig file object is stored in a global 
hash in case it is needed again.

When a single score is requested for a region, then a special low-level 
statistical method is employed to significantly reduce data collection 
times. Up to a ten fold improvement or better has been observed over the 
simple point-by-point collection, depending on the size of the region 
requested.

=head1 USAGE

The module requires Lincoln Stein's Bio::DB::BigWig to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::bigwig;

It will automatically export the name of the subroutines. 

=over

=item collect_bigwig_score

This subroutine will collect a single value from a binary bigWig file. 
It uses the low-level summary method to collect the statistical 
information and is therefore significantly faster than the other 
methods, which rely upon parsing individual data points across the 
region.

The subroutine is passed five or more arguments in the following order:
    
    1) The chromosome or seq_id
    2) The start position of the segment to collect 
    3) The stop or end position of the segment to collect 
    4) The method of collecting the data. Acceptable values include 
       mean, min, max, sum, count, and stddev. 
    5) One or more paths to bigWig files from which to collect the data

The object will return either a valid score. When nothing is found, it 
will return 0 for methods sum and score, or a null '.' value.

=item collect_bigwigset_score

Similar to collect_bigwig_score() but using a BigWigSet database of 
BigWig files. Unlike individual BigWig files, BigWigSet features support 
stranded data collection if a strand attribute is defined in the metadata 
file. 

The subroutine is passed eight or more arguments
    
    1) The opened BigWigSet database object
    2) The chromosome or seq_id
    3) The start position of the segment to collect from
    4) The stop or end position of the segment to collect from
    5) The strand of the segment to collect from
    6) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       or "all". Only those scores which match the indicated 
       strandedness are collected.
    7) The method of collecting the data. Acceptable values include 
       mean, min, max, sum, count, and stddev. 
    8) One or more database feature types for the data 

=item collect_bigwig_scores

This subroutine will collect only the score values from a binary BigWig file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed four or more arguments in the following order:
    
    1) The chromosome or seq_id
    2) The start position of the segment to collect 
    3) The stop or end position of the segment to collect 
    4) One or more paths to bigWig files from which to collect the data

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_bigwigset_scores

Similar to collect_bigwig_scores() but using a BigWigSet database of 
BigWig files. Unlike individual BigWig files, BigWigSet features support 
stranded data collection if a strand attribute is defined in the metadata 
file. 

The subroutine is passed seven or more arguments
    
    1) The opened BigWigSet database object
    2) The chromosome or seq_id
    3) The start position of the segment to collect from
    4) The stop or end position of the segment to collect from
    5) The strand of the segment to collect from
    6) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       or "all". Only those scores which match the indicated 
       strandedness are collected.
    7) One or more database feature types for the data 

=item collect_bigwig_position_scores

This subroutine will collect the score values from a binary BigWig file 
for the specified database region keyed by position. 

The subroutine is passed four or more arguments in the following order:
    
    1) The chromosome or seq_id
    2) The start position of the segment to collect 
    3) The stop or end position of the segment to collect 
    4) One or more paths to bigWig files from which to collect the data

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. Note that only one value is 
returned per position, regardless of the number of dataset features 
passed. Usually this isn't a problem as only one dataset is examined at a 
time.

=item collect_bigwigset_position_score

Similar to collect_bigwig_position_scores() but using a BigWigSet database 
of BigWig files. Unlike individual BigWig files, BigWigSet features support 
stranded data collection if a strand attribute is defined in the metadata 
file. 

The subroutine is passed seven or more arguments
    
    1) The opened BigWigSet database object
    2) The chromosome or seq_id
    3) The start position of the segment to collect from
    4) The stop or end position of the segment to collect from
    5) The strand of the segment to collect from
    6) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       or "all". Only those scores which match the indicated 
       strandedness are collected.
    7) One or more database feature types for the data 

=item open_bigwig_db()

This subroutine will open a BigWig database connection. Pass either the 
local path to a bigWig file (.bw extension) or the URL of a remote bigWig 
file. It will return the opened database object.

=item open_bigwigset_db()

This subroutine will open a BigWigSet database connection using a directory 
of BigWig files and one metadata index file, as described in 
Bio::DB::BigWigSet. Essentially, this treats a directory of BigWig files as 
a single database with each BigWig file representing a different feature 
with unique attributes (type, source, strand, etc). 

Pass the subroutine a scalar value representing the local path to the 
directory. It presumes a feature_type of 'region', as expected by the other 
tim_db_helper subroutines and modules. It will return the opened database 
object.


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



