#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Set::IntSpan;
use modules::Utils;
use modules::Exception;
use modules::SystemCall;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "exon=s",
	   "output=s",
	   "exon_pileup=s",
	   "mindepthcutoff=i",
	   "full",
	   "chr=s",
	   "debug"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{exon} || !$OPT{exon_pileup});



=pod

=head1 SYNOPSIS

exon_coverage.pl -chr chr -exon exonfile -exon_pileup exon_pileup_file -mindepthcutoff min_depth_to_be_considered_covered -full print_all_entries(default=not_full_cover_entries) -output output_file [options]

Required flags: -exon -exon_pileup

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

exon_coverage.pl -> Generate coverage for the exon from a run

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./exon_coverage.pl -mindepthcutoff 4 -exon /home/matt/work/pipeline/conf/mouse/NCBIM37/exon/040411/NCBIM37.exon.1 -exon_pileup /home/matt/work/analysis/pipeline/runs/IGL00054_11/IGL00054_exon_pileup_out_11 -output /home/matt/work/analysis/pipeline/runs/IGL00054_15/overlap/IGL00054_15_run_coverage_match -chr 1

=cut

my $output;
if (defined $OPT{output}) {
	$output = $OPT{output};
	my $basedir = dirname($output);
} else {
	$output = "summary.txt";
}

if ($OPT{chr} eq 1) {
	system("rm $output") if -e $output;
}

if ($OPT{chr}) {
	open(FILE,">>$output") || modules::Exception->throw("Can't open file to write $output\n");
} else {
	open(FILE,">$output") || modules::Exception->throw("Can't open file to write $output\n");
}

my $mindepth = defined $OPT{mindepthcutoff}?$OPT{mindepthcutoff}:1;

#Bin size for exon lookup index
my $binsize = 10000;

my $exon = $OPT{exon};
if (!-e $exon) {
	modules::Exception->throw("ERROR: Cannot open exon file $OPT{exon}");
}

my $pileup = $OPT{exon_pileup};
if (!-e $pileup) {
	modules::Exception->throw("ERROR: Cannot open pileup file $OPT{exon_pileup}");
}

#Store the final stats
my %exon_stats;

#Store total bases not covered
my $missing_bases = 0;

#Store total bases covered by exons
my $total_bases = 0;

my %lines = ();

#Get the depth stats
&Depth_Stats();
my $summary;
my $missing_percent = sprintf("%.2f",$missing_bases/$total_bases*100);

if ($OPT{chr}) {
	$summary = "\nUNCOVERED BASES for chr$OPT{chr}: $missing_bases/$total_bases ($missing_percent%)\n\n";
} else {
	$summary = "\nUNCOVERED BASES for genome: $missing_bases/$total_bases ($missing_percent%)\n\n";
}

print FILE join("\n", sort {my @a = split("\t",$a); my @b = split("\t",$b); $a[1] <=> $b[1]} keys %lines);

print FILE $summary;
print STDERR $summary;


#Generate the depth statistics
sub Depth_Stats {
	open(EXON,"$exon") || die "Can't open file $exon\n";

	while (<EXON>) {
		my @fields = split;
		my $coord_index_start = int($fields[1]/$binsize);
		my $coord_index_end = int($fields[2]/$binsize);
		while ($coord_index_start <= $coord_index_end) {
			#Add initial entries to our data structure so we get entries for every exon
			push @{$exon_stats{$fields[0]}{$coord_index_start}{"$fields[1]-$fields[2]"}{exons}},$fields[3];
			$coord_index_start++;
		}
	}
	
	open(PILEUP,$pileup) || die "Can't open pileup file $pileup";
	my %exon_depths = ();
	
	while (<PILEUP>) {
		if ($OPT{chr}) {
			next unless /^$OPT{chr}\s+/;
		}
		
		my @fields = split;
		my %exon_info = &Coord_to_Exons($fields[0],$fields[1]);
		for my $range ( keys %exon_info ) {
			#Account for exons that share the same exons
			my $exon_str = join(",",@{$exon_info{$range}});
			#store the depth for each individual base
		    $exon_depths{$exon_str}{$range}{$fields[1]} = $fields[3];
		}
	}
	
	
	
	my $exon_count = keys %exon_depths;
	
	my $tmp_count = 0;
	
	for my $chr ( sort keys %exon_stats ) {
	    for my $index ( sort {$a<=>$b} keys %{$exon_stats{$chr}} ) {
	        for my $exon_range ( sort {my ($astart) = $a =~ /^(\d+)/; my ($bstart) = $b =~ /^(\d+)/; $astart <=> $bstart } keys %{$exon_stats{$chr}{$index}} ) {
	        	my @exon_depths = ();
	        	my $exon_str = join(",",@{$exon_stats{$chr}{$index}{$exon_range}{exons}});
	        	my $exon_cover_set;
	        	my $cover_str;
				my $median;
				my $min;
				my $max;
	        	my ($start,$end) = $exon_range =~ /(\d+)-(\d+)/;
	        	my $exon_length = $end-$start;
	        	$total_bases += $exon_length;
	        	
	        	#check there is some coverage
	        	if (exists $exon_depths{$exon_str}{$exon_range}) {
	        		my $exon_set = Set::IntSpan->new($exon_range);
	        		my $exon_cover_set = Set::IntSpan->new();
	        		my $exon_cover_flag = 0;
		        	for my $exon_coord (keys %{$exon_depths{$exon_str}{$exon_range}}) {
		        		push @exon_depths, $exon_depths{$exon_str}{$exon_range}{$exon_coord};
		        		if ($exon_depths{$exon_str}{$exon_range}{$exon_coord} >= $mindepth) {
							#Build up the Set object to get the ranges later on
			        		$exon_cover_set = $exon_cover_set->union(Set::IntSpan->new("$exon_coord-$exon_coord"));
			        		$exon_cover_flag = 1;
			        	}
		        	}
		        	#Calculate the median
		       		$median = modules::Utils->median(\@exon_depths);
    		   		$min = modules::Utils->min(\@exon_depths);
       				$max = modules::Utils->max(\@exon_depths);
       				
       				#Classify the coverage
		        	if ($exon_set->equal($exon_cover_set)) {
						$cover_str = "FULL_COVER";
					} elsif ($exon_cover_flag) {
						#Otherwise get the range that isn't covered
						my $uncovered_set = $exon_set - $exon_cover_set;
						$cover_str = $uncovered_set->run_list();
						$missing_bases += $uncovered_set->cardinality();
					} else {
						$cover_str = "ONLY_LOW_COVER";
						$missing_bases += $exon_length;
					}
					
					
	        	} else {
	        		$median = "N/A";
					$min = "N/A";
					$max = "N/A";
					$cover_str = "NO_COVER";
					$missing_bases += $exon_length;
	        	}    
	        	
	        	my $line = join("\t",$chr, $start, $end, $exon_str, $median, $min, $max, $cover_str) if $cover_str ne 'FULL_COVER';
       			$lines{$line}++;
	        }
	    }
	}
}


#Lookup exon name based on coord
sub Coord_to_Exons {
	my ($chr,$pileup_coord) = @_;
	my %exons = ();
	
	my $group_index = int($pileup_coord/$binsize);
	for my $range (sort { my ($astart) = $a =~ /^(\d+)/; my ($bstart) = $b =~ /^(\d+)/; $astart <=> $bstart } keys %{$exon_stats{$chr}{$group_index}})	{
		my ($start,$end) = $range =~ /(\d+)-(\d+)/;
		if ($start > $pileup_coord) {
			#Break the loop as we're out of range; speeds things up
			last;
		}
		if ($pileup_coord >= $start && $pileup_coord <= $end) {
			for my $exon ( @{$exon_stats{$chr}{$group_index}{$range}{exons}} ) {
			    push @{$exons{"$start-$end"}},$exon;
			}
		}
	}
	return %exons;
	
}




