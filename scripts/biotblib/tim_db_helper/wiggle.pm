package tim_db_helper::wiggle;

# modules
require Exporter;
use strict;
use Carp;
use Bio::Graphics::Wiggle;
our $VERSION = '1.7.0';


# Exported names
our @ISA = qw(Exporter);
our @EXPORT = qw(
	collect_wig_scores
	collect_wig_position_scores
);


# Hashes of opened file objects
our %OPENED_WIGFILES; # opened wigfile objects
	# in empirical testing, this doesn't really seem to speed things up
	# like I thought it would
	# oh well, keep it anyway????


# The true statement
1; 



### Modules ###

sub collect_wig_scores {
	
	# we will actually call collect_wig_position_scores()
	# but only return the values
	
	my %wig_data = collect_wig_position_scores(@_);
	
	# return the values
	return values %wig_data;
}



sub collect_wig_position_scores {
	
	# pass the required information
	my ($start, $stop, $strand, $stranded, @wig_features) = @_;
	
	# set up hash, position => score
	my %wig_data;
	
	# look at each wigfile
	# usually there is only one, but for stranded data there may be 
	# two wigfiles (+ and -), so we'll check each wig file for strand info
	foreach my $feature (@wig_features) {
	
		# Check which data to take based on strand
		if (
			$stranded eq 'all' # stranded data not requested
			or $feature->strand == 0 # unstranded data
			or ( 
				# sense data
				$strand == $feature->strand 
				and $stranded eq 'sense'
			) 
			or (
				# antisense data
				$strand != $feature->strand  
				and $stranded eq 'antisense'
			)
		) {
			# we have acceptable data to collect
			
			# collect from wigfile if present
			if ($feature->has_tag('wigfile') ) {
				
				# get wigfile name
				my @wigfiles = $feature->get_tag_values('wigfile');
				my $wigfile = shift @wigfiles;
				confess " no wigfile passed!\n" unless $wigfile;
				
				# check for opened wigfile
				my $wig;
				if (exists $OPENED_WIGFILES{$wigfile} ) {
					# this file is already opened, use it
					$wig = $OPENED_WIGFILES{$wigfile};
				}
				else {
					# this file has not been opened yet, open it
					unless (-e $wigfile) {
						confess " Binary wiggle file '$wigfile' does not exist!\n";
					}
					$wig = Bio::Graphics::Wiggle->new($wigfile,0);
					unless ($wig) {
						confess " unable to open data wigfile '$wigfile'";
					}
					
					# store the opened object for later use
					$OPENED_WIGFILES{$wigfile} = $wig;
				}
				
				# adjust as necessary to avoid wig errors
				if ($start < $wig->start) {
					# adjust the start position
					$start = $wig->start;
				}
				elsif ($start > $wig->end) {
					# nothing we can do here, no values
					return;
				}
				if ($stop > $wig->end) {
					# adjust the end position
					$stop = $wig->end;
				}
				elsif ($stop < $wig->start) {
					# nothing we can do here, no values
					return;
				}
				
				# collect the wig values
				my $scores_ref = $wig->values($start => $stop);
				
				# re-associate position with the scores
				my $step = $wig->step;
				my $pos = $start;
				foreach my $s (@{ $scores_ref }) {
					#print Dumper($s);
					if (defined $s) {
						# the binary wig file (.wib) is usually set up with 
						# a step of 1 bp, even if the original wig file was not
						# this can result in lots of undefined values at the 
						# positions where there was no original data
						# hence the defined check here
						# store a real value in the hash keyed under the position
						$wig_data{$pos} = $s;
					}
					
					# adjust position by the step size, 
					# regardless whether defined or not
					$pos += $step;
				}
			}
		}
	}	
	
	# return the wig data hash
	return %wig_data;
}




__END__


=head1 NAME

tim_db_helper::wiggle

=head1 DESCRIPTION

This module is used to collect the dataset scores from a binary 
wig file (.wib) that is referenced in the database. Typically, a single 
feature representing the dataset is present across each chromosome. The 
feature should contain an attribute ('wigfile') that references the 
location of the binary file representing the dataset scores. The file is 
read using the Bio::Graphics::Wiggle module, and the values extracted from the 
region of interest. 

Scores may be restricted to strand by specifying the desired strandedness. 
For example, to collect transcription data over a gene, pass the strandedness 
value 'sense'. If the strand of the region database object (representing the 
gene) matches the strand of the wig file data feature, then the data is 
collected.

For loading wig files into a Bio::DB database, see the perl script 
'wiggle2gff3.pl' included with the Bio::Graphics distribution, as well as 
Bio::Graphics::Wiggle::Loader.

To speed up the program and avoid repetitive opening and 
closing of the files, the opened wig file object is stored in a global 
hash in case it is needed again.
 
=head1 USAGE

The module requires Lincoln Stein's Bio::Graphics to be installed. 

Load the module at the beginning of your program.

	use tim_db_helper::wiggle;

It will automatically export the name of the subroutines. 

=over

=item collect_wig_scores

This subroutine will collect only the score values from a binary wig file 
for the specified database region. The positional information of the 
scores is not retained, and the values are best further processed through 
some statistical method (mean, median, etc.).

The subroutine is passed three or more arguments in the following order:
    
    1) The start position of the segment to collect from
    2) The stop or end position of the segment to collect from
    3) The strand of the original feature (or region), -1, 0, or 1.
    4) A scalar value representing the desired strandedness of the data 
       to be collected. Acceptable values include "sense", "antisense", 
       "none" or "no". Only those scores which match the indicated 
       strandedness are collected.
    5) One or more database feature objects that contain the reference 
       to the wib file. They should contain the attribute 'wigfile'.

The subroutine returns an array of the defined dataset values found within 
the region of interest. 

=item collect_wig_position_scores

This subroutine will collect the score values from a binary wig file 
for the specified database region keyed by position. 

The subroutine is passed the same arguments as collect_wig_scores().

The subroutine returns a hash of the defined dataset values found within 
the region of interest keyed by position. Note that only one value is 
returned per position, regardless of the number of dataset features 
passed.

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




