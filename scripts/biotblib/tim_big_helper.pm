package tim_big_helper;

### modules
require Exporter;
use strict;
use Carp qw(carp cluck);
use File::Temp;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use tim_db_helper qw(get_chromosome_list);
use tim_db_helper::config;



### Export
our @ISA = qw(Exporter);
our @EXPORT = qw(
);
our @EXPORT_OK = qw(
	wig_to_bigwig_conversion
	bed_to_bigbed_conversion
	generate_chromosome_file
);


our $VERSION = '1.9.0';

1;


### Wig to BigWig file conversion
sub wig_to_bigwig_conversion {
	
	# Collect passed arguments
	my $argument_ref = shift;
	unless ($argument_ref) {
		cluck "no arguments passed!";
		return;
	}
	
	# wigfile
	my $wigfile = $argument_ref->{'wig'} || q();
	unless ($wigfile) {
		cluck "no wig file passed!";
		return;
	}
	
	
	# Identify bigwig conversion utility
	my $bw_app_path = $argument_ref->{'bwapppath'} || q();
	unless ($bw_app_path) {
		# check for an entry in the configuration file
		$bw_app_path = $TIM_CONFIG->param('applications.wigToBigWig') || 
			undef;
	}
	unless ($bw_app_path) {
		# try checking the system path as a final resort
		$bw_app_path = `which wigToBigWig`;
		chomp $bw_app_path;
	}
	unless ($bw_app_path) {
		carp " Utility 'wigToBigWig' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}
	
	
	# Generate list of chromosome sizes if necessary
	my $chromo_file = $argument_ref->{'chromo'} || q();
	unless ($chromo_file) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		my $database = $argument_ref->{'db'} || q();
		unless ($database) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$chromo_file = generate_chromosome_file($database);
		unless ($chromo_file) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	
	# Generate the bw file name
	# we can substitute one of three possible names for bw
	my $bw_file = $wigfile;
	$bw_file =~ s/\.(?:bed|bdg|bedgraph|wig)$/.bw/;
	
	
	# Generate the bigwig file 
	print " converting $wigfile to bigWig....\n";
	if ($bw_app_path =~ /wigToBigWig$/) {
		# include the -clip option in case there are any positions 
		# out of bounds of the chromosome
		# it will just warn instead of fail
		system $bw_app_path, '-clip', $wigfile, $chromo_file, $bw_file;
	}
	elsif ($bw_app_path =~ /bedGraphToBigWig$/) {
		# this doesn't have the -clip option, too bad
		system $bw_app_path, $wigfile, $chromo_file, $bw_file;
	}
	
	# check the result
	if (-e $bw_file and -s $bw_file) {
		# conversion successful
		if ($chromo_file =~ /^chr_sizes\w{5}/) {
			# we no longer need our temp chromosome file
			unlink $chromo_file;
		}
		return $bw_file;
	}
	else {
		carp " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bw_file) {
			# 0-byte file was created
			unlink $bw_file;
		}
		if ($chromo_file =~ /^chr_sizes\w{5}/) {
			# leave the temp chromosome file as a courtesy
			carp " Leaving temporary chromosome file '$chromo_file'\n";
		}
		return;
	}
}


### Bed to BigBed file conversion
sub bed_to_bigbed_conversion {
	
	# Collect passed arguments
	my $argument_ref = shift;
	unless ($argument_ref) {
		carp "no arguments passed!";
		return;
	}
	
	# bedfile
	my $bedfile = $argument_ref->{'bed'} || undef;
	unless ($bedfile) {
		carp "no bed file passed!";
		return;
	}
	
	
	# identify bigbed conversion utility
	my $bb_app_path = $argument_ref->{'bbapppath'} || undef;
	unless ($bb_app_path) {
		# check for an entry in the configuration file
		$bb_app_path = $TIM_CONFIG->param('applications.bedToBigBed') || 
			undef;
	}
	unless ($bb_app_path) {
		# try checking the system path as a final resort
		$bb_app_path = `which bedToBigBed`;
		chomp $bb_app_path;
	}
	unless ($bb_app_path) {
		carp " Utility 'bedToBigBed' not specified and can not be found!" . 
			" Conversion failed!\n";
		return;
	}
	
	
	# Generate list of chromosome sizes if necessary
	my $chromo_file = $argument_ref->{'chromo'} || undef;
	unless ($chromo_file) {
		# a pre-generated list of chromosome sizes was not provided
		# need to generate one from the database
		my $database = $argument_ref->{'db'} || q();
		unless ($database) {
			carp " No requisite database or chromosome info file provided!" .
				" Conversion failed\n";
			return;
		}
		$chromo_file = generate_chromosome_file($database);
		unless ($chromo_file) {
			carp " Cannot generate chromosome info file! Conversion failed\n";
			return;
		}
	}
	
	
	# Generate the bb file name
	my $bb_file = $bedfile;
	$bb_file =~ s/\.bed$/.bb/;
	
	
	# Generate the bigBed file using Jim Kent's utility
	print " converting $bedfile to BigBed....\n";
	system $bb_app_path, $bedfile, $chromo_file, $bb_file;
	
	
	# Check the result
	if (-e $bb_file and -s $bb_file) {
		# conversion successful
		if ($chromo_file =~ /^chr_sizes\w{5}/) {
			# we no longer need our temp chromosome file
			unlink $chromo_file;
		}
		return $bb_file;
	}
	else {
		print " Conversion failed. You should try manually and watch for errors\n";
		if (-e $bb_file) {
			# 0-byte file was created
			unlink $bb_file;
		}
		if ($chromo_file =~ /^chr_sizes\w{5}/) {
			# leave the temp chromosome file as a courtesy
			print " Leaving temporary chromosome file '$chromo_file'\n";
		}
		return;
	}
}



sub generate_chromosome_file {
	
	my $database = shift;
	print " generating chromosome file....\n";
	
	# generate chromosome lengths file
	my @chromosomes = get_chromosome_list($database);
	unless (@chromosomes) {
		carp " no chromosome sequences identified in database!\n";
		return;
	}
	
	# prepare temp file
	my $chr_fh = new File::Temp(
		'UNLINK'   => 0,
		'TEMPLATE' => 'chr_sizesXXXXX',
	);
	my $chromo_file = $chr_fh->filename;

	# write out
	foreach my $chr (@chromosomes) {
		# chromosome name and size
		$chr_fh->print( $chr->[0] . "\t" . $chr->[1] . "\n");
	}
	$chr_fh->close;
	
	return $chromo_file;
}




__END__

=head1 NAME

tim_big_helper

=head1 DESCRIPTION

This module helps in the conversion of wig and bed files to bigWig and 
bigBed files, respectively. It uses external applications to 
accomplish this, taking care of generating a chromosome file from a 
database if necessary. 

Two exported subroutines are available for wig and bed conversions. 

=head1 USAGE

Load the module at the beginning of your program and include the name or 
names of the subroutines to export. None are automatically exported.

	use tim_big_helper qw(wig_to_bigwig_conversion);


=over

=item wig_to_bigwig_conversion()

This subroutine will convert a wig file to a bigWig file. See the UCSC 
documentation regarding wig (http://genome.ucsc.edu/goldenPath/help/wiggle.html)
and bigWig (http://genome.ucsc.edu/goldenPath/help/bigWig.html) file formats. 
It uses Jim Kent's wigToBigWig utility to perform the conversion. 
This must be present on the system for the conversion to succeed. 

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "<chromosome name> 
<size in bases>". This file may be specified, or automatically generated if 
given a Bio::DB database name (preferred to ensure genome version 
compatibility).

The function returns the name of the bigWig file, which will be the 
input wig file basename with the BigWig ".bw". Note that the it does 
not check for success of writing the bigwig file. Check STDERR for errors 
in bigwig file generation.

Pass the function an anonymous hash of arguments, including the following:

  Required:
  wig         => The name of the wig source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bwapppath   => Provide the full path to Jim Kent's wigToBigWig 
                 utility. This parameter may instead be defined in the 
                 configuration file "biotoolbox.cfg". 

Example

	my $wig_file = 'example_wig';
	my $bw_file = wig_to_bigwig_conversion( {
			'wig'   => $wig_file,
			'db'    => $database,
	} );
	if (-e $bw_file) {
		print " success! wrote bigwig file $bw_file\n";
		unlink $wig_file; # no longer necessary
	}
	else {
		print " failure! see STDERR for errors\n";
	};

=item bed_to_bigbed_conversion

This subroutine will convert a bed file to a bigBed file. See the UCSC 
documentation regarding bed (http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED)
and bigBed (http://genome.ucsc.edu/goldenPath/help/bigBed.html) file formats. 
It uses Jim Kent's bedToBigBed utility to perform the conversion. This 
must be present on the system for the conversion to succeed. 

The conversion requires a list of chromosome name and sizes in a simple text 
file, where each line is comprised of two columns, "<chromosome name> 
<size in bases>". This file may be specified, or automatically generated if 
given a Bio::DB database name (preferred to ensure genome version 
compatibility).

The function returns the name of the bigBed file, which will be the 
input bed file basename with the extension ".bb". Note that the it does 
not check for success of writing the bigbed file. Check STDERR for errors 
in bigbed file generation.

Pass the function an anonymous hash of arguments, including the following:

  Required:
  bed         => The name of the bed source file. 
  db          => Provide an opened database object from which to generate 
                 the chromosome sizes information.
  Optional: 
  chromo      => The name of the chromosome sizes text file, described 
                 above, as an alternative to providing the database name.
  bbapppath   => Provide the full path to Jim Kent's bedToBigBed  
                 utility. This parameter may instead be defined in the 
                 configuration file "biotoolbox.cfg". 

Example

	my $wig_file = 'example_wig';
	my $bw_file = wig_to_bigwig_conversion( {
			'wig'   => $wig_file,
			'db'    => $database,
	} );
	if (-e $bw_file) {
		print " success! wrote bigwig file $bw_file\n";
		unlink $wig_file; # no longer necessary
	}
	else {
		print " failure! see STDERR for errors\n";
	};

=item generate_chromosome_file

This subroutine will generate a chromosome sizes files appropriate for 
the big file conversion utilities from an available database. It is a 
two column text file, the first column is the chromosome name, and the 
second column is the length in bp. The file is written in the 
current directory with a name of "chr_sizesXXXXX", where X are random 
characters as defined by File::Temp. 

The chromosome names and lengths are obtained from a Bio::DB 
database using the L<tim_db_helper::get_chromosome_list()> 
subroutine.

Pass the subroutine a database name, path to a supported database file, 
or opened Bio::DB object.

The file will be written, closed, and the filename returned.

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


