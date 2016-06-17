#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "input=s",
	   "debug=s",
	   "cosmic_out=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{input});



=pod

=head1 SYNOPSIS

parse_cosmic.pl -input input_file -cosmic_out (default=GRCh37.cosmic) [options]

Required flags: -input

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_cosmic.pl -> Script to parse cosmic db dump

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./parse_cosmic.pl -input CosmicCompleteExport_v62_291112.tsv

=cut

my $cosmic_file = $OPT{input};
my $cosmic_output = defined $OPT{mel_out}?$OPT{mel_out}:'GRCh37.cosmic';

my $snv_out = $cosmic_output . '.overlap.snv';
my $indel_out = $cosmic_output . '.overlap.indel';
my $gene_out = $cosmic_output . '.gene';


if ( !-e $cosmic_file ) {
	modules::Exception->throw("File $cosmic_file doesn't exist");	
}

my %cosmic_data = ();
my %genes = ();

open(FILE,"$cosmic_file") || modules::Exception->throw("Can't open file $cosmic_file\n");

while (<FILE>) {
	chomp;
	next if /^Gene name/;
	my @fields = split("\t");
	
	my ($genename,$ens,undef,undef,undef,undef,$prim,$sub,$hist,undef,undef,undef,$base_str,$aa_str,$mut_type,undef,undef,undef,$coord_str,$strand,$somatic_str,undef,undef,$tumour_type) = @fields;
	next unless $coord_str =~ /:/;
	
	my ($chr,$start_coord,$end_coord) = $coord_str =~ /([0-9XY]+):(\d+)-(\d+)/;
	if ($chr == 23) {
		$chr = 'X';
	} elsif ($chr == 24) {
		$chr = 'Y';
	}
	my $mutant_type;
	my $mut_group;
	my $hist_group;
	my $base_change;
	my $aa_change;
	
	my $somatic_confirmed = $somatic_str =~ /omatic/?'somatic':'unknown';
	
	#print join("\t", $genename, $ens, $prim,$sub,$hist, $base_str,$aa_str,$mut_type,$coord_str,$strand,$somatic_str,$tumour_type) . "\n";
	my $cosmic_str;
	
	#SNV handling
	if ($mut_type =~ /Substitution/) {
		$mutant_type = 'SNV';
		$mut_group = 'snv';
		
		
		
		if ($base_str =~ /([ACGT]+\>[ACGTYRBMK]+)$/) {
			$base_change = $1;
			$base_change =~ s/>/->/;
		} else {
			$base_change = $base_str;
			#modules::Exception->warning("ERROR: Problem with bp $base_str\n $_\n");			
			
		}
		
		if ($aa_str =~ /p\.([A-Z\*])\d+\>?([A-Z\*])$/) {
			my $first_aa = $1;
			my $second_aa = $2;
			if ($first_aa eq $second_aa) {
				$aa_change = 'SYN';		
			} else {
				$aa_change = $first_aa . '->' . $second_aa;
			} 
		} elsif (/coding silent/) {
			$aa_change = 'SYN';
		} else {
			$aa_change = $aa_str;
			#modules::Exception->warning("ERROR: Problem with aa $aa_str\n $_\n");			
		}
		
	} elsif ($mut_type =~ /Insertion/) {
		$mutant_type = 'INS';
		$mut_group = 'indel';
		
		if ($base_str =~ /(ins[0-9]+)/) {
			$base_change = $1;
		} elsif ($base_str =~ /ins([ACTGN]+)/) {
			$base_change = $1;
		} elsif ($base_str =~ /([ACGT]+>[ACGT]+)$/) {
			$base_change = $1;
		} elsif ($base_str =~ />([ACGT]+)$/) {
			$base_change = $1;
		} else {
			$base_change = $base_str;
			#modules::Exception->warning("ERROR: Problem with bp $base_str\n $_\n");	
		}
		
		if ($_ =~ /Frameshift/) {
			$aa_change = 'frameshift';
		} elsif ($aa_str =~ /ins([A-Z]+)/) {
			$aa_change = $1;
		} elsif ($aa_str =~ />([A-Z]+)/) {
			$aa_change = $1;
		}elsif ($aa_str =~ /(ins\d+)/) {
			$aa_change = $1;
		} elsif ($aa_str =~ /p\.([A-Z*]+).*([A-Z]+)/) {
			$aa_change = $1.'..'.$2;
		}else {
			#modules::Exception->warning("ERROR: Problem with aa $aa_str\n $_\n");
			$aa_change = $aa_str;
		}
		
		$end_coord-- if $end_coord>$start_coord;
		
	} elsif ($mut_type =~ /Deletion/) {
		$mutant_type = 'DEL';
		$mut_group = 'indel';
		
		if ($base_str =~ /(del[0-9]+)/) {
			$base_change = $1;
		} elsif ($base_str =~ /del([ACGT]+)/) {
			$base_change = $1;
		} elsif ($base_str =~ /([ACGT]+>[ACGT]+)$/) {
			$base_change = $1;
		} elsif ($base_str =~ />([ACGT]+)$/) {
			$base_change = $1;
		} else {
			$base_change = $base_str;
			#modules::Exception->warning("ERROR: Problem with bp $base_str\n");
		}
		
		if ($_ =~ /Frameshift/) {
			$aa_change = 'frameshift';
		} elsif ($aa_str =~ /del([A-Z*]+)/) {
			$aa_change = $1;
		} elsif ($aa_str =~ /p\.([A-Z*]+)/ && (length($base_change) == 3|| $base_change eq 'del3')) {
			$aa_change = $1;			
		} elsif ($aa_str =~ /p\.([A-Z*]+).*([A-Z]+)/) {
			$aa_change = $1.'..'.$2;
		} else {
			$aa_change = $aa_str;
			#modules::Exception->warning("ERROR: Problem with aa $aa_str\n");
		}
	} else {
		#Don't bother with strange mutation types (unknown, complex, etc)
		next;
	}
	
	
	$genename =~ s/_ENST\d+//;
	$genename =~ s/_CCDS//;
	my @gene_names = split(",",$genename);
	for my $name (@gene_names) {
		$genes{uc($name)}++;
	}
	$cosmic_str = join('^^^', $mutant_type, $base_change, $strand, $aa_change, $genename , $hist, $tumour_type, $somatic_confirmed);
	#print "$chr $start_coord $end_coord $cosmic_str\n" if $start_coord != $end_coord;
	
	#Want to group entries with different strands or tumour types as the same entries
	$cosmic_data{$mut_group}{$chr}{$start_coord}{$end_coord}{$cosmic_str}++;
}

open(SNV,">$snv_out") || modules::Exception->throw("ERROR: Can't open file $snv_out for output");

for my $chr (sort keys %{$cosmic_data{snv}}) {
	for my $start_coord (sort {$a<=>$b} keys %{$cosmic_data{snv}{$chr}}) {
		for my $end_coord (sort {$a<=>$b} keys %{$cosmic_data{snv}{$chr}{$start_coord}}) {
			my $final_cosmic_str =  &Collapse(keys %{$cosmic_data{snv}{$chr}{$start_coord}{$end_coord}});
			print SNV join("\t",$chr,$start_coord,$end_coord,$final_cosmic_str) . "\n";
		}
	}
}

close SNV;

open(INDEL,">$indel_out") || modules::Exception->throw("ERROR: Can't open file $indel_out for output");

for my $chr (sort keys %{$cosmic_data{indel}}) {
	for my $start_coord (sort {$a<=>$b} keys %{$cosmic_data{indel}{$chr}}) {
		for my $end_coord (sort {$a<=>$b} keys %{$cosmic_data{indel}{$chr}{$start_coord}}) {
			my $final_cosmic_str = &Collapse(keys %{$cosmic_data{indel}{$chr}{$start_coord}{$end_coord}});
			print INDEL join("\t",$chr,$start_coord,$end_coord,$final_cosmic_str) . "\n";
		}
	}
}

close INDEL;

open(GENE,">$gene_out") || modules::Exception->throw("Can't open file $gene_out\n");
print GENE join("\n",sort keys %genes) . "\n";
close GENE;


#Combines the entries with same coordinates and variant type into single line
sub Collapse {
	my @str_entries = @_;
	
	my %combined_fields;
	
	for my $str_entry (@str_entries) {
		my @fields = split('\^\^\^',$str_entry);
		
		for (my $a = 0; $a < @fields; $a++) {
			$combined_fields{$a}{$fields[$a]}++;	
		}
		
	}

	my $combined_str;
	
	for my $count (sort {$a<=>$b} keys %combined_fields) {
		$combined_str .= "^^^" . join(",",keys %{$combined_fields{$count}});
	}
	
	$combined_str =~ s/^\^\^\^//;	
	return $combined_str;
	
}
	
