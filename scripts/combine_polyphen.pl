#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use modules::Utils;
use Bio::EnsEMBL::Registry;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
			"nomap_in=s",
			"nomap_out=s",	   		
			"map_in=s",
			"map_out=s",
	   		"debug=s",
	   		);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{nomap_in} || !$OPT{map_in});



=pod

=head1 SYNOPSIS

combine_polyphen.pl -map_in map_input_file -nomap_in nomap_input_file -map_out map_output_file(default=nomap_in.new) -nomap_out nomap_output_file(default=map_in.new) -debug [options]

Required flags: -map_in -nomap_in

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

combine_polyphen.pl -> clean up the polyphen overlap files

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

/home/matt/work/pipeline/conf/parsers/combine_polyphen.pl -map_in NCBIM37.uniprot_mapping.1 -nomap_in NCBIM37.uniprot_nomapping.1

=cut

my $map_infile = $OPT{map_in};
my $nomap_infile = $OPT{nomap_in};
my $map_outfile = defined $OPT{map_out}?$OPT{map_out}:$map_infile.'.new';
my $nomap_outfile = defined $OPT{nomap_out}?$OPT{nomap_out}:$nomap_infile.'.new';




my %data;

open (MAPPED, $map_infile) || modules::Exception->throw("Cannot open file $map_infile");

while (<MAPPED>) {
	chomp;
	my @fields = split(" ");
	push @{$data{$fields[0]}{$fields[1]}{mapped}}, $fields[3]; 
}

open (UNMAPPED, $nomap_infile) || modules::Exception->throw("Cannot open file $nomap_infile");

while (<UNMAPPED>) {
	chomp;
	my @fields = split(" ");
	push @{$data{$fields[0]}{$fields[1]}{unmapped}}, $fields[3];
}

open(MAPOUT,">$map_outfile") || modules::Exception->throw("Cannot open output file $map_outfile");
open(UNMAPOUT,">$nomap_outfile") || modules::Exception->throw("Cannot open output file $nomap_outfile");

for my $chr (keys %data) {
	COORD:
	for my $coord (sort {$a<=>$b} keys %{$data{$chr}}) {
		if (exists $data{$chr}{$coord}{mapped}) {
			for my $entry (@{$data{$chr}{$coord}{mapped}}) {
				if ($entry =~ /CCDS\d+$/) {
					print MAPOUT join(" ",
								$chr,
								$coord,
								$coord,
								$entry
								) ."\n";
					next COORD;
				}
			}
			#Otherwise just print the first entry
			my ($first_entry) = @{$data{$chr}{$coord}{mapped}};
			print MAPOUT join(" ",
								$chr,
								$coord,
								$coord,
								$first_entry
								) . "\n";
			 
		} elsif (exists $data{$chr}{$coord}{unmapped}) {
			my ($first_entry) = @{$data{$chr}{$coord}{unmapped}};
			print UNMAPOUT join(" ",
								$chr,
								$coord,
								$coord,
								$first_entry
								) . "\n";		
		} else {
			modules::Exception->throw("ERROR: Shouldn't be here");
		}
	}
}




