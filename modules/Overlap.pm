package modules::Overlap;

use strict;
use modules::Exception;
use Data::Dumper;

=pod

=head1 SYNOPSIS

overlapcoordinates_fast.pl will read the 2 input files (ref, coord) and check which elements in ref overlap elements in coord (-fail will report the regions that have no match)

Required flags: -ref -coord

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

    -fail  reports only non-matching elements

    -all   reports everything, whether matched or not

    -full  reports only matching that are 100% overlapping

    -full_rev reports only elements in -coord, that are fully overlapped by elements in -ref

    -count Will only report the number of elements in coord that overlap with each element in ref

    -binsize Will use bins of this size to split up our lookup table (default == 1000)

    -us Will replace all spaces in the strings in the coord descriptions with '_'

    -individual_lines Instead of joining all overlaps on one line, they get split into multiple lines
		
    -fiftyp Makes matches only if 50% of the ref region is overlapped by a coord region
    
    -sort_coord Sort the overlapping entries for each line with multiple overlapping coord entries by start coord
    
    -sort_name Sort the overlapping entries for each line with multiple overlapping coord entries by name
    
    -ignore Ignore coord entries containing this text or comma delimited string of text entries
    
    -output outfile name
    
    -append append to existing file
    
    -silent only_generate_overlap_object
    
    -exact only_report_exact_coordinate_matches
    
    -max_distance max_distance_for_start_and_end_coord_to_be_away
    
=head1 NAME

Overlap.pm

=head1 DESCRIPTION

August, 2007

This script will read each element from both the supplied ref and coord and report which, if any elements are overlapping between them

The input files should have the following format

chr start end (anything after these fields will be kept as the descriptive string to be reported upon match)

Usually, ref will be a list of genes, etc, while coord is a list of clone, for example

=head1 AUTHOR

Matthew Field

=cut

#Class the finds overlap between coordinate sets
sub new
{
	my ( $class, $logfile ) = @_;
	my $self = bless {}, $class;
	return $self;
}

#Function that does the actual overlap
sub overlap
{
	my $self = shift;
	my %args = @_;
	
	
	#print Dumper \%args;
	if ( !defined $args{-ref} || !defined $args{-coord} )
	{
		modules::Exception->throw(
"Need to define -ref (reference file to overlap with coordinate set) and -coord (file containing coordinate set)"
		);
	}
	
	
	my $delim      = '^^^';
	my $ref_file   = $args{-ref};
	my $coord_file = $args{-coord};
	if ( !-e $ref_file )
	{
		modules::Exception->throw("Problem with ref file $ref_file");
	}

	if ( !-e $coord_file )
	{
		modules::Exception->throw("Problem with coord file $coord_file");
	}

	#If there are entries to ignore matches
	my %coord_ignore = ();
	if ( $args{-ignore} )
	{
		%coord_ignore = map { $_ => 1 } split( ",", $args{-ignore} );
	}

	#Given that we are generally checking single coordinates keep as 1 by default; when we start looking at larger ranges we want to up this value to 1000
	my $binsize = defined $args{-binsize} ? $args{-binsize} : 1000;
	my $underscore = defined $args{-us} ? 1 : 0;

	#Load coord -  This is usually an annotation file (eg dbsnp file)
	open( COORD, $coord_file ) || die "Can't open file $coord_file \n";
	my %coord_lists;
	my %coord_locs;
	my %overlaps_return_struct;
	my $cnt           = 1;
	my $overlap_count = 0;

	my %dup_coords = ();

	while (<COORD>)
	{
		#ignore blank lines
		next unless /\S/;
		next if /OVERLAP COUNT/;
		chomp;
		my ( $chr, $start, $end, $flagthis );
		if ( $_ =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $3;
			$flagthis = $4;
		}
		elsif ( $_ =~ /^(\S+)\s+(\d+)\s+(\d+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $3;
			$flagthis = $cnt++;
		}
		elsif ( $_ =~ /^(\S+)\s+(\d+)\s+(\S+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $start;
			$flagthis = $3;
		}
		elsif ( $_ =~ /^(\S+)\s+(\d+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $start;
			$flagthis = $cnt++;
		}
		else
		{
			modules::Exception->throw(
						"ERROR: Problem with format of -coord file $coord_file at line $_" );
		}

		if ( $end < $start )
		{
			modules::Exception->throw(
						  "ERROR: End coord less than start coord for line $_");
		}

		if ( $args{-chr} )
		{
			next unless $chr eq $args{-chr};
		}

		if ($underscore)
		{
			$flagthis =~ s/ /_/g;
		}

		#If we're seen the key before note it's a dup entry
		if ( exists $coord_locs{$flagthis} )
		{
			$dup_coords{$flagthis} = 1;
		}

		#save the info for this element in a hash
		$coord_locs{$flagthis}{$chr}{set}{"$start-$end"} = 1;

		#Add this element to a hash table to speed up the lookups later on
		my $start_ind = int( $start / $binsize );
		my $stop_ind  = int( $end / $binsize );
		foreach ( $start_ind - 1 .. $stop_ind + 1 )
		{
			$coord_lists{$chr}{$_}{$flagthis} = 1;
		}
		

	}
	close(COORD);

	my @overlap_entries = ();

	#my $count = 0;
	##Load the stuff in ref - this is the input list (eg filter_basic snps)
	open( REF, $ref_file ) || die "Can't open file $ref_file \n";
	while (<REF>)
	{
		next unless /\S/;
		next if /OVERLAP COUNT/;
		chomp;
		my ( $chr, $start, $end, $flagthis );
		if ( $_ =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $3;
			$flagthis = $4;
		}
		elsif ( $_ =~ /^(\S+)\s+(\d+)\s+(\d+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $3;
			$flagthis = '';
		}
		elsif ( $_ =~ /^(\S+)\s+(\d+)\s+(\S+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $start;
			$flagthis = $3;
		}
		elsif ( $_ =~ /^(\S+)\s+(\d+)/ )
		{
			$chr      = $1;
			$start    = $2;
			$end      = $start;
			$flagthis = '';
		}
		else
		{
			modules::Exception->throw(
				"ERROR: Problem with format of -ref file $ref_file at line $_");
		}

		if ( $end < $start )
		{
			modules::Exception->throw(
						  "ERROR: End coord less than start coord for line $_");
		}

		if ( $args{-chr} )
		{
			next unless $chr eq $args{-chr};
		}
		my %overlaps =
		  $self->_getOverlapList( $chr, $start, $end, $binsize, \%coord_lists,
								  \%coord_locs, \%args );

		#print Dumper \%overlaps;

		my @ovlps = ();
		if ( exists( $args{-sort_coord} ) )
		{
			@ovlps = sort
			{
				my ($astart) = $coord_locs{$a}{$chr}{set} =~ /(\d+)/; my ($bstart) = $coord_locs{$b}{$chr}{set} =~ /(\d+)/; $astart <=> $bstart
			} keys %overlaps;
		}
		elsif ( exists ($args{-sort_name} ) ) {
			@ovlps = sort
			{
				$a cmp $b
			} keys %overlaps;
		}
		
		else
		{
			@ovlps = keys %overlaps;
		}

		my $outstr = '';

		foreach my $entry (@ovlps)
		{
			# We split the coord values if we're ignoring entries; this allows 1 to 1 matching
			#Different handling for different formats of coord strings:
			# 's1' -> s1 is treated as a single string entry
			# 's1,s2,,,sN' -> s1 to sN are treated as individual entries
			#Here we need to chop up the entry and see if any novel ones remain 
			if ( keys %coord_ignore > 1) {
				my @entries = split(",", $entry);
			    my $proceed = 0;
			    my @new_entries = ();
				for my $local_entry ( @entries ) {
					#Add entries not in ignore list
				    if (! exists $coord_ignore{$local_entry}) {
				    	push @new_entries, $local_entry;
				    }
				}
				#Next if no entries remain after removing the ignore list entries
				next unless @new_entries;			
				$entry = join(",",@new_entries);
			} elsif (keys %coord_ignore == 1) {
				#Here we have a single ignore entry so just skip if it exists
				if ( exists $coord_ignore{$entry} ) {
					next;
				}
			} 
			
			#Flag the coord entry as a dup
			if ( exists $dup_coords{$entry} )
			{
				$entry .= "-DUP";
			}

			if ( exists( $args{-individual_lines} ) )
			{

				#Here we report every single overlap coord entry
				push @overlap_entries, "$chr\t$start\t$end\t$flagthis\t$entry";
				push @{ $overlaps_return_struct{PASS}{$chr}{$start}{$end}
					  {$flagthis} }, $entry;
				$overlap_count++;
			}
			else
			{

				#Build up the output string with a delimiter
				if ( $outstr ne '' )
				{
					$outstr = $outstr . "$delim$entry";
				}
				else
				{
					$outstr = $entry;
				}
				push @{ $overlaps_return_struct{PASS}{$chr}{$start}{$end}
					  {$flagthis} }, $entry;
			}
		}

		if ( !defined $args{-individual_lines} && !exists( $args{-fail} ) )
		{
			if ( length($outstr) > 0 )
			{
				push @overlap_entries, "$chr\t$start\t$end\t$flagthis\t$outstr";
				$overlap_count++;
			}
			else
			{
				if ( exists( $args{-all} ) )
				{
					my $non_mark =
					  defined $args{-all} ? $args{-all} : "NO_MATCH";
					push @overlap_entries,
					  "$chr\t$start\t$end\t$flagthis\t$non_mark";
					push @{ $overlaps_return_struct{FAIL}{$chr}{$start}{$end}
						  {$flagthis} }, $non_mark;
				}
			}
		}

		#only print the unmatched stuff if asked for
		if ( exists( $args{-fail} ) && $outstr eq '' )
		{
			push @overlap_entries, "$chr\t$start\t$end\t$flagthis\tNO_MATCH";
			push
			  @{ $overlaps_return_struct{FAIL}{$chr}{$start}{$end}{$flagthis} },
			  "NO_MATCH";
			$overlap_count++;
		}

	}


	if (@overlap_entries)
	{
		if ( $args{-output} )
		{
			open( OUT, ">$args{-output}" )
			  || die "Can't open file to write $args{-output}\n";
			print OUT join( "\n", @overlap_entries ), "\n";
			close(OUT);
		}
		elsif ( $args{-append} )
		{
			open( OUT, ">>$args{-append}" )
			  || die "Can't open file to append $args{-append}\n";
			print OUT join( "\n", @overlap_entries ), "\n";
			close(OUT);
		}
		elsif (!$args{-silent}) {
			print join( "\n", @overlap_entries ), "\n";
		} 
		
	}

	#Remove the pass entries if we only want fails
	if ( exists $args{-fail} )
	{
		delete $overlaps_return_struct{'PASS'};
	}
	return ( \%overlaps_return_struct, $overlap_count );
}

#Generate hash containing overlapping coord entries for the input ref entry
sub _getOverlapList
{
	my $self = shift;
	my ( $chr, $start, $end, $binsize, $reflist, $reflocs, $args ) = @_;
	my %return_hash;

	#find all the ref entries that overlap with the given region
	my $strt = int( $start / $binsize );
	my $stp  = int( $end / $binsize );
	#print "INDEX $binsize $start $end $strt $stp\n";
	foreach ( $strt .. $stp )
	{
		foreach my $entry ( keys %{ $reflist->{$chr}{$_} } )
		{
			ENTRY: foreach my $coords ( keys %{ $reflocs->{$entry}{$chr}{set} } )
			{
				my ( $ref_start, $ref_end ) = split( '-', $coords );

				#Flag sets whether physical overlap required
				my $no_other_args = 1; 
				if (exists $args->{-full} || exists $args->{-full_rev} || exists $args->{-exact} || exists $args->{-fiftyp}) {
					$no_other_args = 0;
				}
				
				
				#Check if it's standard overlap
				if ( $ref_start <= $end && $ref_end >= $start )
				{
					#Just looking for any overlap
					$return_hash{$entry} = 1;
					
					#Check if there are any special overlap requirements and remove match if conditions fail
					if ( exists $args->{-full} )
					{
						#check the "ref" (list being looped) is fully overlapped by "coord" (the reference hash)
						if ( $ref_start > $start || $ref_end < $end )
						{
							delete $return_hash{$entry};
						}
					}

					if ( exists $args->{-full_rev} )
					{

						#check the "coord" (the reference hash) is fully overlapped by "ref" (list being looped)
						if ( $start > $ref_start || $end < $ref_end )
						{
							delete $return_hash{$entry};
						}
					}
					
					if ( exists $args->{-exact} )
					{
						#check the "coord" (the reference hash) is an exact coordinate match with "ref" (list being looped)
						if ( $start != $ref_start || $end != $ref_end )
						{
							delete $return_hash{$entry};
						}
					}
					
					if ( exists $args->{-fiftyp} )
					{    
						
						my $overlap_bases = 0;
						my $ref_size = $ref_end - $ref_start;
						my $coord_size = $end - $start;

						if ($start<=$ref_start && $end>=$ref_end) {
							$overlap_bases = $ref_size;
						} elsif ($start>=$ref_start && $end<=$ref_end) {
							$overlap_bases = $coord_size;
						} elsif ($start<=$ref_start && $end<=$ref_end) {
							$overlap_bases = $end-$ref_start;
						} elsif ($start>=$ref_start && $end>=$ref_end) {
							$overlap_bases = $ref_end-$start;
						}
						
						if ($overlap_bases < (0.5 * $coord_size) ) {
							delete $return_hash{$entry};
						}
						
						#need at least half of the region to overlap
#						my $my_size  = $end - $start;
#						my $ref_size = $ref_end - $ref_start;
#						my $left     = $end - $ref_start;
#						my $right    = $ref_end - $start;
#						my $mx       = ( $left > $right ) ? $left : $right;
#						$mx = ( $mx < $ref_size ) ? $mx : $ref_size;
#
#				   		print "Some Overlap: R $right L $left $mx $my_size ". $mx/$my_size." \n";
#						if ( $mx < ( 0.5 * $my_size ) )
#						{
#							#print "Overlap\n";
#							delete $return_hash{$entry};
#						}
					}
					
					#This overlap uses distance
					if ( exists $args->{-max_distance} )
					{
						my $max_distance = $args->{-max_distance};
						#check both the start and end are not farther away than max_distance
						if(abs($start-$ref_start) > $max_distance || abs($end-$ref_end) > $max_distance) {
							delete $return_hash{$entry};
						}
					}
					
					#Don't bother with the rest of the entry coordinates once we've found a match
					if (exists $return_hash{$entry}) {
						last ENTRY;
					}
				} 
				#This overlap type doesn't require a physical overlap always so need to check here as well if no other arguments passed that require overlap
				elsif ( exists $args->{-max_distance} && $no_other_args)
				{
					my $max_distance = $args->{-max_distance};
					#check both the start and end are not farther away than max_distance
					if(abs($start-$ref_start) < $max_distance && abs($end-$ref_end) < $max_distance) {
						$return_hash{$entry} = 1;
					}
				}
				

#				if ( exists $args->{-full} )
#				{
#
#					#check the "ref" (list being looped) is fully overlapped by "coord" (the reference hash)
#					if ( $ref_start <= $start && $ref_end >= $end )
#					{
#						$return_hash{$entry} = 1;
#					}
#				}
#				elsif ( exists $args->{-exact} )
#				{
#
#					#check the "coord" (the reference hash) is an exact coordinate match with "ref" (list being looped)
#					if ( $start == $ref_start && $end == $ref_end )
#					{
#						$return_hash{$entry} = 1;
#					}
#				}
#				elsif ( exists $args->{-full_rev} )
#				{
#
#					#check the "coord" (the reference hash) is fully overlapped by "ref" (list being looped)
#					if ( $start <= $ref_start && $end >= $ref_end )
#					{
#						$return_hash{$entry} = 1;
#					}
#				}
#				
#				else
#				{
#					if ( exists $args->{-fiftyp} )
#					{    #need at least half of the region to overlap
#						if ( $ref_start <= $end && $ref_end >= $start )
#						{
#							my $my_size  = $end - $start;
#							my $ref_size = $ref_end - $ref_start;
#							my $left     = $end - $ref_start;
#							my $right    = $ref_end - $start;
#							my $mx       = ( $left > $right ) ? $left : $right;
#							$mx = ( $mx < $ref_size ) ? $mx : $ref_size;
#
#					   #print "Some Overlap: $mx $my_size ". $mx/$my_size." \n";
#							if ( $mx > ( 0.5 * $my_size ) )
#							{
#
#								#print "Overlap\n";
#								$return_hash{$entry} = 1;
#							}
#						}
#					}
#					elsif ( $ref_start <= $end && $ref_end >= $start )
#					{
#						#Just looking for any overlap
#						$return_hash{$entry} = 1;
#						#Don't bother with the rest of the entry coordinates once we've found a match
#						last ENTRY;
#					}
#				}
			}
		}
	}

	return %return_hash;
}

1;
