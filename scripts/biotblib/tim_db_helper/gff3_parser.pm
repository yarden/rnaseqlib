package tim_db_helper::gff3_parser;

use strict;
require Exporter;
use Carp qw(carp cluck croak confess);
use Bio::SeqFeature::Lite;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_file_helper qw(
	open_to_read_fh
);
my $VERSION = '1.8.1';

sub new {
	my $class = shift;
	my $self = {};
	bless $self, $class;
	
	# check for a file name
	if (@_) {
		$self->parse_file(shift);
	}
	
	# done
	return $self;
}

sub parse_file {
	my $self = shift;
	
	# check file
	my $file = shift;
	unless ($file) {
		cluck("no file name passed!\n");
		return;
	}
	unless ($file =~ /\.gff3?(?:\.gz)?$/) {
		cluck("file doesn't look like a GFF file!\n");
		return;
	}
	
	# open file
	my $fh = open_to_read_fh($file) or 
		croak("unable to open file '$file'!\n");
	
	$self->{'fh'} = $fh;
}


sub top_features {
	my $self = shift;
	
	my $in_fh = $self->{fh};
	
	# initialize
	my @top_features;
	my @orphan_features;
	my %loaded; # hash of loaded features
	
	# Collect the top features
	# we will continue to read the file until we either reach the file end
	# or we reach a close features directive (###)
	
	# Each line will be processed into a SeqFeature object, and then checked 
	# for parentage. Child objects with a Parent tag will be appropriately 
	# associated with its parent, or put into an orphanage. Any orphans 
	# left will be checked a final time for parents; if a parent can't be 
	# found, it will be lost. Features without a parent are assumed to be 
	# top-level features.
	
	while (my $line = $in_fh->getline) {
		
		chomp $line;
		
		# process comment and pragma lines
		if ($line =~ /^##gff-version ([\d\.]+)$/) {
			# version pragma, check it
			if ($1 != 3) {
				croak " Input GFF version is $1 not 3!\n" ;
			}
			next;
		}
		elsif ($line =~ /^###/) {
			# close features directive line, 
			# all subfeatures have been loaded for this chromosome/scaffold
			last;
		}
		elsif ($line =~ /^##FASTA$/) {
			# FASTA pragma
			# go no further
			last;
		}
		elsif ($line =~ /^#/) {
			# some other unrecognized pragma or a comment line
			# skip
			next;
		}
		elsif ($line =~ /^>/) {
			# fasta header line, skip
			next;
		}
		elsif ($line =~ /^[agctn]+$/i) {
			# fasta sequence, skip
			next;
		}
		
		
		# generate the SeqFeature object for this GFF line
		my $feature = $self->from_gff3_string($line);
		
		### Process the feature
		
		# check the ID
		my $id = $feature->primary_id;
		if ($id) {
			# this ID should be unique in the GFF file
			# complain if it isn't
			if (exists $loaded{$id}) {
				carp " Feature ID '$id' occurs more than once in file! skipping\n";
				next;
			} 
			else {
				$loaded{$id} = $feature;
			}
		}
		# if the feature didn't have an ID, we'll just assume it is
		# a child of another feature, otherwise it may get lost
		
		# look for parents and children
		if ($feature->has_tag('Parent')) {
			# must be a child
			my ($parent) = $feature->get_tag_values('Parent');
			if (exists $loaded{$parent}) {
				# we've seen this id
				# associate the child with the parent
				$loaded{$parent}->add_SeqFeature($feature);
			}
			else {
				# can't find the parent, maybe not loaded yet?
				# put 'em in the orphanage
				push @orphan_features, $feature;
			}
		}
		else {
			# must be a parent
			push @top_features, $feature;
		}
	}
	# Finished loading the GFF lines (so far....)
	
	# check for orphans
	if (@orphan_features) {
		
		my $lost_count = 0;
		while (my $feature = shift @orphan_features) {
			
			# find the parent
			my ($parent) = $feature->get_tag_values('Parent');
			if (exists $loaded{$parent}) {
				# we've seen this id
				# associate the child with the parent
				$loaded{$parent}->add_SeqFeature($feature);
			}
			else {
				# can't find the parent
				# forget about them, just another statistic
				$lost_count++;
			}
		}
		
		# report
		if ($lost_count) {
			carp " $lost_count features could not be associated with" . 
				" reported parents!\n";
		}
	}
	
	# finished (for now)
	return @top_features;
}



sub from_gff3_string {
	
	my $self = shift;
	
	# get the string
	my $string = shift;
	unless ($string) {
		cluck("must pass a string!\n");
		return;
	}
	chomp $string;
	
	# parse the string
	my (
		$seq_id, 
		$source, 
		$primary_tag, 
		$start, 
		$end, 
		$score, 
		$strand, 
		$phase, 
		$group_tags
	) = split(/\t/, $string);
	
	# generate the basic SeqFeature
	my $feature = Bio::SeqFeature::Lite->new(
		-seq_id         => $seq_id,
		-start          => $start,
		-end            => $end,
	);
	
	# add more attributes
	if ($primary_tag ne '.') {
		$feature->primary_tag($primary_tag);
	}
	if ($source ne '.') {
		$feature->source($source);
	}
	if ($score ne '.') {
		$feature->score($score);
	}
	if ($phase =~ /^[012]$/) {
		$feature->phase($phase);
	}
	
	# add strand
	if ($strand eq '+') {
		$feature->strand(1);
	}
	elsif ($strand eq '-') {
		$feature->strand(-1);
	}
	else {
		$feature->strand(0);
	}
	
	# process groups
	my @groups = split(/\s*;\s*/, $group_tags);
	foreach my $group (@groups) {
		
		my ($tag, $value) = split /=/, $group;
		$tag = $self->unescape($tag);
		my @values = map { $self->unescape($_) } split(/,/, $value);
		
		# determine the appropriate attribute based on tag name
		if ($tag eq 'Name') {
			$feature->display_name($values[0]);
		}
		elsif ($tag eq 'ID') {
			$feature->primary_id($values[0]);
		}
		else {
			foreach (@values) {
				$feature->add_tag_value($tag, $_);
			}
		}
	}
	
	# finished
	return $feature;
}


sub unescape {
  # Borrowed unashamedly from bioperl Bio::Tools::GFF
  # which in turn was borrowed from Bio::DB::GFF
  
  my $self = shift;
  
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}



__END__

=head1 NAME

tim_db_helper::gff3_parser

=head1 DESCRIPTION

This module parses a GFF3 file into SeqFeature objects. Children features 
are associated with parents as sub SeqFeature objects, assuming the Parent 
tag is included and correctly identifies the unique ID tag of the parent. 
Any feature without a Parent tag is assumed to be a parent. Children 
features referencing a parent feature that has not been loaded may be 
lost. 

Embedded Fasta sequences are ignored, as are most comment and pragma lines.

Close directives (###) in the GFF3 file are highly encouraged to limit 
parsing, otherwise the entire file will be slurped into memory. Refer to 
the GFF3 definition at http://www.sequenceontology.org for more details.

The SeqFeature objects that are returned are Bio::SeqFeature::Lite objects. 
Refer to that documentation for more information.

=head1 SYNOPSIS

  use tim_db_helper::gff3_parser;
  my $filename = 'file.gff3';
  
  my $parser = tim_db_helper::gff3_parser->new($filename) or 
  	die "unable to open gff file!\n";
  
  while (my @top_features = $parser->top_features() ) {
  	while (@top_features) {
  		my $feature = shift @top_features;
  		# each $feature is a Bio::SeqFeature::Lite object
  		my @children = $feature->get_SeqFeatures();
  	}
  }

=head1 METHODS

=over

=item new()

=item new($file)

Initialize a new gff3_parser object.

Optionally pass the name of the GFF3 file, and it will be automatically 
opened by calling parse_file().

=item parse_file($file)

Pass the name of a GFF3 file to be parsed. The file must have a .gff or 
.gff3 extension, and may optionally be gzipped (.gz extension). 

=item top_features()

This method will return an array of the top (parent) features defined in 
the GFF3 file. The file will be progressively parsed from the beginning 
until either a close features pragma (###) or the end of the file is 
reached. Features containing a Parent attribute are associated with the 
corresponding feature, if it was loaded. 

When close pragmas are present in the file, call this method repeatedly 
to finish parsing the remainder of the file.

=item from_gff3_string($string)

This method will parse a GFF3 formatted string or line of text and  
return a Bio::SeqFeature::Lite object.

=item unescape($text)

This method will unescape special characters in a text string. Certain 
characters, including ";" and "=", are reserved for GFF3 formatting and 
are not allowed, thus requiring them to be escaped.

=back

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

This package is free software; you can redistribute it and/or modify
it under the terms of the GPL (either version 1, or at your option,
any later version) or the Artistic License 2.0.  

