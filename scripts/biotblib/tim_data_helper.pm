package tim_data_helper;

### modules
require Exporter;
use strict;
use Carp;


### Variables
# Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
	generate_tim_data_structure
	verify_data_structure
	index_data_table
	find_column_index
	parse_list
	format_with_commas
);
our $VERSION = '1.8.5';



### The True Statement
1; 



#################   The Subroutines   ###################

### Generate a new empty tim data structure
sub generate_tim_data_structure {
	
	# Collect the feature
	my $feature = shift;
	
	# Collect the array of dataset headers
	my @datasets = @_;
	
	# Initialize the hash structure
	my %tim_data = (
		'program'        => $0,
		'feature'        => $feature,
		'db'             => q(),
		'gff'            => 0,
		'bed'            => 0,
		'number_columns' => 0,
		'last_row'       => 0,
		'headers'        => 1,
		'other'          => [],
		'data_table'     => [],
	);
	
	# Add the datasets
	my $index = 0;
	foreach my $dataset (@datasets) {
		
		# the metadata structure
		$tim_data{$index} = {
			'name'        => $dataset,
			'index'       => $index,
		};
		
		# the table header
		$tim_data{'data_table'}->[0][$index] = $dataset;
		
		# increment number columns
		$tim_data{'number_columns'} += 1;
		
		$index++;
	}
	
	# Finished
	return \%tim_data;
}





### Subroutine to check for missing required values in the data hash
sub verify_data_structure {
	my $datahash_ref = shift;
	
	# check for data table
	unless (
		defined $datahash_ref->{'data_table'} and 
		ref $datahash_ref->{'data_table'} eq 'ARRAY'
	) {
		carp " No data table in passed data structure!";
		return;
	}
	
	# check for column index
	if (defined $datahash_ref->{'number_columns'}) {
		my $number = scalar @{ $datahash_ref->{'data_table'}->[0] };
		if ($datahash_ref->{'number_columns'} != $number) {
			warn " number of data columns [$number] doesn't match " . 
				"metadata value [" . $datahash_ref->{'number_columns'} . "]!\n";
			# fix it for them
			$datahash_ref->{'number_columns'} = $number;
		}
	}
	else {
		$datahash_ref->{'number_columns'} = 
			scalar @{ $datahash_ref->{'data_table'}->[0] };
	}
	
	# check metadata
	for (my $i = 0; $i < $datahash_ref->{'number_columns'}; $i++) {
		unless (
			$datahash_ref->{$i}{'name'} eq 
			$datahash_ref->{'data_table'}->[0][$i]
		) {
			carp " incorrect or missing metadata!\n  Column header names don't" 
				. " match metadata name values for index $i\n" . 
				"  compare '" . $datahash_ref->{$i}{'name'} . "' with '" .
				$datahash_ref->{'data_table'}->[0][$i] . "'\n";
			return;
		}
	}
	
	# check for last row index
	if (defined $datahash_ref->{'last_row'}) {
		my $number = scalar( @{ $datahash_ref->{'data_table'} } ) - 1;
		if ($datahash_ref->{'last_row'} != $number) {
			warn " data table last_row index [$number] doesn't match " . 
				"metadata value [" . $datahash_ref->{'last_row'} . "]!\n";
			# fix it for them
			$datahash_ref->{'last_row'} = $number;
		}
	}
	else {
		# define it for them
		$datahash_ref->{'last_row'} = 
			scalar( @{ $datahash_ref->{'data_table'} } ) - 1;
	}
	
	# check for proper gff structure
	if ($datahash_ref->{'gff'}) {
		# if any of these checks fail, we will reset the gff version to 
		# the default of 0, or no gff
		my $gff_check = 1; # start with assumption it is true
		
		# check number of columns
		if ($datahash_ref->{'number_columns'} != 9) {
			$gff_check = 0;
		}
		
		# check column indices
		if (
			exists $datahash_ref->{0} and
			$datahash_ref->{0}{'name'} !~ 
			m/^chr|chromo|seq|refseq|ref_seq|seq|seq_id/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{1} and
			$datahash_ref->{1}{'name'} !~ m/^source/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{2} and
			$datahash_ref->{2}{'name'} !~ m/^type|method/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{3} and
			$datahash_ref->{3}{'name'} !~ m/^start/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{4} and
			$datahash_ref->{4}{'name'} !~ m/^stop|end/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{5} and
			$datahash_ref->{5}{'name'} !~ m/^score|value/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{6} and
			$datahash_ref->{6}{'name'} !~ m/^strand/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{7} and
			$datahash_ref->{7}{'name'} !~ m/^phase/i
		) {
			$gff_check = 0;
		}
		if (
			exists $datahash_ref->{8} and
			$datahash_ref->{8}{'name'} !~ m/^group|attribute/i
		) {
			$gff_check = 0;
		}
		
		# update gff value as necessary
		if ($gff_check == 0) {
			# reset metadata
			$datahash_ref->{'gff'} = 0;
			$datahash_ref->{'headers'} = 1;
			
			# remove the AUTO key from the metadata
			for (my $i = 0; $i < $datahash_ref->{'number_columns'}; $i++) {
				if (exists $datahash_ref->{$i}{'AUTO'}) {
					delete $datahash_ref->{$i}{'AUTO'};
				}
			}
		}
	}
	
	# check for proper BED structure
	if ($datahash_ref->{'bed'}) {
		# if any of these checks fail, we will reset the bed flag to 0
		# to make it not a bed file format
		my $bed_check = 1; # start with assumption it is correct
		
		# check number of columns
		if (
			$datahash_ref->{'number_columns'} < 3 and 
			$datahash_ref->{'number_columns'} > 12 
		) {
			$bed_check = 0;
		}
		
		# check column index names
		# we're assuming these names are derived from tim_file_helper or 
		# something like it
		if (
			exists $datahash_ref->{0} and
			$datahash_ref->{0}{'name'} !~ 
			m/^chr|chromo|seq|refseq|ref_seq|seq|seq_id/i
		) {
			$bed_check = 0;
		}
		if (
			exists $datahash_ref->{1} and
			$datahash_ref->{1}{'name'} !~ m/^start/i
		) {
			$bed_check = 0;
		}
		if (
			exists $datahash_ref->{2} and
			$datahash_ref->{2}{'name'} !~ m/^stop|end/i
		) {
			$bed_check = 0;
		}
		
		# the remaining columns are tricky, as they may or may not be 
		# named as I expect, especially if it was generated de novo
		# so only check these if the original file extension was bed
		if (
			exists $datahash_ref->{'extension'} and 
			$datahash_ref->{'extension'} =~ /bed|bdg/i
		) {
			if (
				exists $datahash_ref->{3} and
				$datahash_ref->{3}{'name'} !~ m/^name|id|score/i
				# for bed this should be name or ID
				# for bedgraph this should be score
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{4} and
				$datahash_ref->{4}{'name'} !~ m/^score|value/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{5} and
				$datahash_ref->{5}{'name'} !~ m/^strand/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{6} and
				$datahash_ref->{6}{'name'} !~ m/^thickstart/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{7} and
				$datahash_ref->{7}{'name'} !~ m/^thickend/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{8} and
				$datahash_ref->{8}{'name'} !~ m/^itemrgb/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{9} and
				$datahash_ref->{9}{'name'} !~ m/^blockcount/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{10} and
				$datahash_ref->{10}{'name'} !~ m/^blocksizes/i
			) {
				$bed_check = 0;
			}
			if (
				exists $datahash_ref->{11} and
				$datahash_ref->{11}{'name'} !~ m/^blockstarts/i
			) {
				$bed_check = 0;
			}
		}
		
		# reset the BED tag value as appropriate
		if ($bed_check) {
			$datahash_ref->{'bed'} = $datahash_ref->{'number_columns'};
		}
		else {
			# reset metadata
			$datahash_ref->{'bed'} = 0;
			$datahash_ref->{'headers'} = 1;
			
			# remove the AUTO key from the metadata
			for (my $i = 0; $i < $datahash_ref->{'number_columns'}; $i++) {
				if (exists $datahash_ref->{$i}{'AUTO'}) {
					delete $datahash_ref->{$i}{'AUTO'};
				}
			}
		}
	}
	
	# check proper SGR file structure
	if (
		$datahash_ref->{'extension'} =~ /sgr/i or
		$datahash_ref->{'filename'} =~ /sgr/i
	) {
		# there is no sgr field in the data structure
		# so we're just checking for the extension
		# we will change the extension as necessary if it doesn't conform
		if (
			$datahash_ref->{'number_columns'} != 3 or
			$datahash_ref->{0}{'name'} !~ /^chr|seq|ref/i or
			$datahash_ref->{1}{'name'} !~ /^start|position/i
		) {
			# doesn't smell like a SGR file
			# change the extension so the write subroutine won't think it is
			# make it a text file
			$datahash_ref->{'extension'} =~ s/sgr/txt/i;
			$datahash_ref->{'filename'}  =~ s/sgr/txt/i;
			$datahash_ref->{'headers'} = 1;
			
			# remove the AUTO key from the metadata
			for (my $i = 0; $i < $datahash_ref->{'number_columns'}; $i++) {
				if (exists $datahash_ref->{$i}{'AUTO'}) {
					delete $datahash_ref->{$i}{'AUTO'};
				}
			}
		}
	}
	
	return 1;
}





#### Index a data table
sub index_data_table {
	
	# get the arguements
	my ($data_ref, $increment) = @_;
	
	# check data structure
	unless (defined $data_ref) {
		carp " No data structure passed!";
		return;
	}
	unless ( verify_data_structure($data_ref) ) {
		return;
	}
	if (exists $data_ref->{'index'}) {
		warn " data structure is already indexed!\n";
		return 1;
	}
	
	# check column indices
	my $chr_index = find_column_index($data_ref, '^chr|seq|refseq');
	my $start_index = find_column_index($data_ref, '^start');
	unless (defined $chr_index and $start_index) {
		carp " unable to find chromosome and start dataset indices!\n";
		return;
	}
	
	# define increment value
	unless (defined $increment) {
		# calculate default value
		if (exists $data_ref->{$start_index}{'win'}) {
			# in genome datasets, window size metadata is stored with the 
			# start position
			# increment is window size x 20
			# seems like a reasonable compromise between index size and efficiency
			$increment = $data_ref->{$start_index}{'win'} * 20;
		}
		else {
			# use some random made-up default value that could be totally 
			# inappropriate, maybe we should carp a warning instead
			$increment = 100;
		}
	}
	$data_ref->{'index_increment'} = $increment;
	
	# generate index
	my $table_ref = $data_ref->{'data_table'};
	my %index;
	for (my $row = 1; $row <= $data_ref->{'last_row'}; $row++) {
		
		# the index will consist of a complex hash structure
		# the first key will be the chromosome name
		# the first value will be the second key, and is the integer of 
		# the start position divided by the increment
		# the second value will be the row index number 
		
		# calculate the index value
		my $index_value = int( $table_ref->[$row][$start_index] / $increment );
		
		# check and insert the index value
		unless (exists $index{ $table_ref->[$row][$chr_index] }{ $index_value} ) {
			# insert the current row, which should be the first occurence
			$index{ $table_ref->[$row][$chr_index] }{ $index_value } = $row;
		}
	}
	
	# associate the index hash
	$data_ref->{'index'} = \%index;
	
	# success
	return 1;
}





### Subroutine to find a column
sub find_column_index {
	my ($data_ref, $name) = @_;
	
	# the $name variable will be used as a regex in identifying the name
	# fix it so that it will possible accept a # character at the beginning
	# without a following space, in case the first column has a # prefix
	$name =~ s/ \A (\^?) (.+) \Z /$1#?$2/x;
	
	# walk through each column index
	my $index;
	for (my $i = 0; $i < $data_ref->{'number_columns'}; $i++) {
		# check the names of each column
		if ($data_ref->{$i}{'name'} =~ /$name/i) {
			$index = $i;
			last;
		}
	}
	return $index;
}




### Parse string into list
sub parse_list {
	# this subroutine will parse a string into an array
	# it is designed for a string of numbers delimited by commas
	# a range of numbers may be specified using a dash
	# hence 1,2,5-7 would become an array of 1,2,5,6,7
	
	my $string = shift;
	return unless defined $string;
	if ($string =~ /[^\d,\-\s\&]/) {
		carp " the string contains characters that can't be parsed\n";
		return;
	}
	$string =~ s/\s+//g; 
	my @list;
	foreach (split /,/, $string) {
		# check for a range
		if (/\-/) { 
			my ($start, $stop) = split /\-/;
			# add each item in the range to the list
			for (my $i = $start; $i <= $stop; $i++) {
				push @list, $i;
			}
			next;
		} 
		else {
			# either an ordinary number or an "&"ed list of numbers 
			push @list, $_;
		}
	}
	return @list;
}



### Format a number into readable comma-delimited by thousands number
sub format_with_commas {
	# for formatting a number with commas
	my $number = shift;
	if ($number =~ /[^\d,\-\.]/) {
		carp " the string contains characters that can't be parsed\n";
		return $number;
	}
	
	# check for decimals
	my ($integers, $decimals);
	if ($number =~ /^\-?(\d+)\.(\d+)$/) {
		$integers = $1;
		$decimals = $2;
	}
	else {
		$integers = $number;
	}
	
	# format
	my @digits = split //, $integers;
	my @formatted;
	while (@digits) {
		if (@digits > 3) {
			unshift @formatted, pop @digits;
			unshift @formatted, pop @digits;
			unshift @formatted, pop @digits;
			unshift @formatted, ',';
		}
		else {
			while (@digits) {
				unshift @formatted, pop @digits;
			}
		}
	}
	
	# finished
	return join("", @formatted) . $decimals;
}




__END__

=head1 NAME

tim_data_helper

=head1 DESCRIPTION

These are general subroutines for working with data, and specifically the 
tim data structure. It also provides a catchall location for common 
subroutines that don't fit in either tim_file_helper or tim_db_helper.

=head1 TIM DATA STRUCTURE

The tim data structure is a complex data structure that is commonly used 
throughout the biotoolbox scripts, thus simplifying data input/output and 
manipulation. The primary structure is a hash with numerous keys. The actual 
data table is represented as an array of arrays.  Metadata for the columns 
(datasets) are stored as hashes. 

This whole data structure is intended to eventually become the structure for 
a blessed object and the basis for a class of object-oriented methods that 
work on the structure. One of these days I'll find time to implement that and 
rewrite all of the biotoolbox scripts.... Maybe for that mythical 2.0 release.

The description of the primary keys in the tim data structure are described 
here.

=over 4

=item program

This includes the scalar value from the Program header
line and represents the name of the program that generated the
data file.

=item db

This includes the scalar value from the Database header
line and is the name of the database from which the file data
was generated.

=item feature

This includes the scalar value from the Feature header
line and describes the type of features in the data file.

=item gff

This includes a scalar value of the source GFF file version, obtained 
from either the GFF file pragma or the file extension. 
The default value is 0 (not a GFF file). As such, it may be treated 
as a boolean value.

=item bed

If the source file is a BED file, then this tag value is set to the 
number of columns in the original BED file, an integer of 3 to 12. 
The default value is 0 (not a BED file). As such, it may be treated 
as a boolean value.

=item number_columns

This includes an integer representing the total 
number of columns in the data table. It is automatically calculated
from the data table and updated each time a column is added.

=item last_row

This represents the integer for the index number (0-based) of the last
row in the data table. It is calculated automatically from the data 
table.

=item other

This key points to an anonymous array of additional, unrecognized 
header lines in the parsed file. For example, metadata from older 
file formats or general comments not suitable for other locations. 
The entire line is added to the array, and is rewritten before the 
column metadata is written. The line ending character is automatically 
stripped when it is added to this array upon file loading, and 
automatically added when writing out to a text file.

=item filename

The original path and filename of the file opened and parsed. (Just in 
case you forgot ;) Joking aside, missing extensions may have been added 
to the filename by the different functions upon opening (a convenience for 
users) in the case that they weren't initially provided. The actual 
complete name will be found here.

=item basename

The base name of the original file name, minus the extension(s). 
Useful when needing to assign a new file name based on the current 
file name.

=item extension

The known extension(s) of the original file name. Known extensions 
currently include '.txt, .gff, .gff3, .bed, .sgr' as well as 
their gzip equivalents.

=item path

The parent directories of the original file. The full filename can be 
regenerated by concatenating the path, basename, and extension.

=item headers

A boolean flag (1 or 0) to indicate whether headers are present or not. 
Some file formats, e.g. BED, GFF, etc., do not explicitly have column 
headers; the headers flag should be set to false in this case. Standard 
tim data formatted text files should be set to true.

=item <column_index_number>

Each column will have a metadata index. Usually this is read from 
the column's metadata line. The key will be the index number (0-based) 
of the column. The value will be an anonymous hash consisting of 
the column metadata. For metadata header lines from a parsed file, these
will be the key=value pairs listed in the line. There should always 
be two mandatory keys, 'name' and 'index'. 

=item data_table

This key will point to an anonymous array of arrays, representing the
tab-delimited data table in the file. The primary array 
will be row, representing each feature. The secondary array will be 
the column, representing the descriptive and data elements for each 
feature. Any value can be looked up by 
$data_structure_ref->{'data_table'}->[$row][$column]. The first row
should always contain the column (dataset) names, regardless whether 
the original data file had dataset names (e.g. GFF or BED files).

=item index

This is an optional index hash to provide faster lookup to specific 
data table rows for genomic bin features. The index is generated 
using the index_data_table() function. The index is comprised of 
another hash data structure, where the first key represents the 
chromosome name, and the second key represents an index value. 
The index value is the integer (or whole number rounding down) of 
the start value divided by the index_increment value. For example, 
with a genomic bin feature at chr1:10691..10700 and an index_increment
value of 100, the index value would be {chr1}{106}. The value of 
that key would be the index number of that row, or more specifically, 
the row index for the first occurence of that index_value (which 
would've been genomic bin feature chr1:10601..10610). Hence, the 
index will get you relatively close to your desired genomic 
position within the data_table, but you will still need to step 
through the features (rows) starting at the indexed position 
until you find the row you want. That should save you a little bit 
of time, at least. The index is not stored upon writing to a 
standard tim data text file.

=item index_increment

This is a single number representing the increment value to calculate 
the index value for the index. It is generated along with the index 
by the index_data_table() function. The index_increment value is not 
stored upon writing to a standard tim data text file.

=back

=head1 USAGE

Call the module at the beginning of your perl script. Include the name(s) 
of the subroutines to import.
  
  use tim_data_helper qw(generate_tim_data_structure);
  

The specific usage for each subroutine is detailed below.


=over

=item generate_tim_data_structure()

As the name implies, this generates a new empty data structure as described 
above. Populating the data table and metadata is the responsibility of the 
end user.

Pass the module an array. The first element should be the name of the features 
in the data table. This is an arbitrary, but required, value. The remainder 
of the array should be the name(s) of the columns (datasets). A rudimentary 
metadata hash for each dataset is generated (consisting only of name and 
index). The name is also entered into the first row of the data table (row 0, 
the header row).

It will return the reference to the tim_data_structure.

Example
	
	my $main_data_ref = generate_tim_data_structure(qw(
		genomic_intevals
		Chromo
		Start
		Stop
	));

=item verify_data_structure()

This subroutine verifies the data structure. It checks items such as the
presence of the data table array, the number of columns in the data table
and metadata, the metadata index of the last row, the presence of basic
metadata, and verification of dataset names for each column. For data 
structures with the GFF or BED tags set to true, it will verify the 
format, including column number and column names; if a check fails, it 
will reset the GFF or BED key to false. It will automatically correct 
some simple errors, and complain about others.

Pass the data structure reference. It will return 1 if successfully 
verified, or false if not.

=item index_data_table()

This function creates an index hash for genomic bin features in the 
data table. Rather than stepping through an entire data table of 
genomic coordinates looking for a specific chromosome and start 
feature (or data row), an index may be generated to speed up the 
search, such that only a tiny portion of the data_table needs to be 
stepped through to identify the correct feature.

This function generates two additional keys in the tim data structure 
described above, L</index> and L</index_increment>. Please refer to 
those items in L<TIM DATA STRUCTURE> for their description and 
usage.

Pass this subroutine one or two arguments. The first is the reference 
to the data structure. The optional second argument is an integer 
value to be used as the index_increment value. This value determines 
the size and efficiency of the index; small values generate a larger 
but more efficient index, while large values do the opposite. A 
balance should be struck between memory consumption and speed. The 
default value is 20 x the feature window size (determined from the 
metadata). Therefore, finding the specific genomic coordinate 
feature should take no more than 20 steps from the indexed position. 
If successful, the subroutine returns a true value.

Example

	my $main_data_ref = load_tim_data_file($filename);
	index_data_table($main_data_ref) or 
		die " unable to index data table!\n";
	...
	my $chr = 'chr9';
	my $start = 123456;
	my $index_value = 
		int( $start / $main_data_ref->{index_increment} ); 
	my $starting_row = $main_data_ref->{index}{$chr}{$index_value};
	for (
		my $row = $starting_row;
		$row <= $main_data_ref->{last_row};
		$row++
	) {
		if (
			$main_data_ref->{data_table}->[$row][0] eq $chr and
			$main_data_ref->{data_table}->[$row][1] <= $start and
			$main_data_ref->{data_table}->[$row][2] >= $start
		) {
			# do something
			# you could stop here, but what if you had overlapping
			# genomic bins for some odd reason?
		} elsif (
			$main_data_ref->{data_table}->[$row][0] ne $chr
		) {
			# no longer on same chromosome, stop the loop
			last;
		} elsif (
			$main_data_ref->{data_table}->[$row][1] > $start
		) {
			# moved beyond the window, stop the loop
			last;
		}
	}
		
=item find_column_index()

This subroutine helps to find the index number of a dataset or column given 
only the name. This is useful if the file contents are not in a standard 
order, for example a typical tim data text file instead of a GFF or BED file.

Pass the subroutine two arguments: 1) The reference to the data structure, and 
2) a scalar text string that represents the name. The string will be used in 
regular expression pattern, so Perl REGEX notation may be used. The search 
is performed with the case insensitive flag. The index position of the first 
match is returned.

Example

	my $main_data_ref = load_tim_data_file($filename);
	my $chromo_index = find_column_index($main_data_ref, "^chr|seq");
	
=item parse_list()

This subroutine parses a scalar value into a list of values. The scalar is 
a text string of numbers (usually column or dataset indices) delimited by 
commas and/or including a range. For example, a string "1,2,5-7" would become 
an array of [1,2,5,6,7].

Pass the module the scalar string.

It will return the array of numbers.

Example
	
	my $index_request = '1,2,5-7';
	my @indices = parse_list($index_request); # returns [1,2,5,6,7]


=item format_with_commas()

This subroutine process a large number (e.g. 4327908475) into a human-friendly 
version with commas delimiting the thousands (4,327,908,475).

Pass the module a scalar string with a number value.

It will return a scalar value containing the formatted number.

Example
	
	my $count = '4327908475';
	print " The final count was " . format_with_commas($count) . "\n";



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



