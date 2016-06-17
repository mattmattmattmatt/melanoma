#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"input=s",
		"output_prefix=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{input} || !$OPT{output_prefix});

	   
=pod

=head1 SYNOPSIS

split_by_chr.pl -input input_file -output_prefix output_prefix [options]

Required flags: -input -output_prefix

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

split_by_chr.pl -> Splits an input file into by chromosome files

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Dan Andrews

=head1 EXAMPLE

split_by_chr.pl -input input_file -output_prefix output_stub

=cut

my $input_file = $OPT{input};
my $output_file_stub = $OPT{output_prefix};

# Work through input reads and write to chromosome specific files

my %CHR_FILEHANDLES;

open (my $INPUT, $input_file) 
    or die "Unable to open input file [$input_file]";

while (<$INPUT>){
	next if /^#/;
    my @cols = split /\t/;

    my $chr = $cols[0];
    $chr =~ s/chr//;

    if (! defined $CHR_FILEHANDLES{$chr}) {

		my $output_file = $output_file_stub ."\." . $chr;
	
		open($CHR_FILEHANDLES{$chr}, ">$output_file")
		    or die "Unable to open filehandle [$output_file]";
    }

    print {$CHR_FILEHANDLES{$chr}} $_;
}

close($INPUT);

# Close chromosome-specific filehandles
foreach my $chr (keys %CHR_FILEHANDLES) {
    close($CHR_FILEHANDLES{$chr})
}


