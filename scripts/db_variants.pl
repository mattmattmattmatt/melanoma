#! /usr/bin/perl -w

use strict;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::BulkInsert;
use modules::Exception;
use modules::Pipeline;
use modules::Utils;
use modules::QualityEncoding;
use modules::SystemCall;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(
			\%OPT,                   
			"help|h",
			"man|m",                 
			"runid=i",
			"snv_pileup=s",          
			"indel_pileup=s",
			"overlap_outfile_snv=s", 
			"overlap_outfile_indel=s",
			"infile_snv=s",          
			"infile_indel=s",
			"chr=s",
			"sample_name=s",
			"writeDB=i",             
			"coord=i",
			"debug",
);
pod2usage( -verbose => 2 ) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{sample_name} || !$OPT{snv_pileup} || !$OPT{indel_pileup} || !$OPT{overlap_outfile_snv} || !$OPT{overlap_outfile_indel} || !$OPT{infile_snv} || !$OPT{infile_indel} );

=pod

=head1 SYNOPSIS

db_variants.pl -runid runid -snv_pileup snv_pileup_file  -indel_pileup indel_pileup_file  -writeDB(default=1) -overlap_outfile_snv output_overlap_snv_file -infile_snv input_overlap_file -debug output_commands [options]

Required flags: -runid -snv_pileup -indel_pileup -snv_normal_pileup -indel_normal_pileup -overlap_outfile_snv -overlap_outfile_indel -infile_snv -infile_indel 

=head1 OPTIONS

	-help  brief help message
	
	-man   full documentation

=head1 NAME

db_variants.pl -> Add all variants to the db

=head1 DESCRIPTION

Feb 15, 2011

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./db_variants.pl

=cut

my $delim             = '^^^';
my $runid             = $OPT{runid};
my $write_to_db       = defined $OPT{writeDB} ? $OPT{writeDB} : 1;
my $pileup_file_snv   = $OPT{snv_pileup};
my $pileup_file_indel = $OPT{indel_pileup};
my $sample_name = $OPT{sample_name};

my $tumour_flag = modules::Pipeline->get_tumour_flag(-sample_name=>$sample_name);
my $pileup_normal_file_snv;
my $pileup_normal_file_indel;

if ($tumour_flag) {
	($pileup_normal_file_snv = $pileup_file_snv) =~ s/pileup/normal.pileup/;
	($pileup_normal_file_indel = $pileup_file_indel) =~ s/pileup/normal.pileup/;
	if ( !-e $pileup_normal_file_snv ) {
		modules::Exception->throw("File $pileup_normal_file_snv doesn't exist");	
	}
	if ( !-e $pileup_normal_file_indel ) {
		modules::Exception->throw("File $pileup_normal_file_indel doesn't exist");	
	}
	
}

my $output_snv = $OPT{overlap_outfile_snv};
my $input_snv  = $OPT{infile_snv};
if ( !-e $input_snv ) {
	modules::Exception->throw("File $input_snv doesn't exist");
}

my $output_indel = $OPT{overlap_outfile_indel};
my $input_indel  = $OPT{infile_indel};
if ( !-e $input_indel ) {
	modules::Exception->throw("File $input_indel doesn't exist");
}

#Get the quality stores from a previous step
my %snv_info = ();
my $run;
my $filter_snv;
my $filter_deletion;
my $filter_insertion;

if ($write_to_db) {
	($run) = modules::Adaptors::Run->search( 'id' => $runid );
	($filter_snv) = modules::Adaptors::Filter->search( 'name' => 'filter_snv' );
	($filter_insertion) =
	  modules::Adaptors::Filter->search( 'name' => 'filter_insertion' );
	($filter_deletion) =
	  modules::Adaptors::Filter->search( 'name' => 'filter_deletion' );

	unless (    defined $run
			 && defined $filter_snv
			 && defined $filter_deletion
			 && defined $filter_insertion )
	{
		modules::Exception->throw(
					 "Unable to retrieve run and filter objects from database");
	}
}

open( SNV, "$input_snv" ) || modules::Exception->throw("Can't open file $input_snv\n");
while (<SNV>) {
	chomp;
	next unless /^[0-9XYM]/;
	my @fields = split("\t");
	
	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $fields[0];
	}
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $fields[1];
	}
	
	my ( $type, $base_change, $quality, $clr, $snv_class, $depth ) = split( '\^\^\^', $fields[3] );
	$quality =~ s/Q//;
	$depth   =~ s/D//;
	
	$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{clr} = $clr;
	$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{quality} = $quality;
	$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{quality_reads} = $depth;
	$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{class} = $snv_class;
	
}

my $sys = modules::SystemCall->new();

#If it's the first chromosome remove the file as we're appending; should happen with db cleanup but if run fails in this step it wouldn't work
if ($OPT{chr} eq 1 || !defined $OPT{chr}) {
	$sys->run("rm $output_snv") if (-e $output_snv);
	$sys->run("rm $output_indel") if (-e $output_indel);
}

#Append since we're running on per chr basis
open( SNVOUT, ">>$output_snv" ) || modules::Exception->throw("Can't open file $output_snv\n");
open( INDELOUT, ">>$output_indel" ) || modules::Exception->throw("Can't open file $output_indel\n");

#Sometimes we may not have snvs and indels for a particular chr
if ( !-e $pileup_file_snv && !-e $pileup_file_indel ) {
	modules::Exception->warning("SKIP STEP db_variants: input pileup file $pileup_file_snv and $pileup_file_indel  doesn't exist\n");
	exit;
}


my $count         = 0;
my $print_counter = 0;

open( my $PF_SNV, $pileup_file_snv ) || modules::Exception->throw("Unable to open file [$pileup_file_snv]");

#my @snv_filter_inserts = ();

my @snv_rows;
my @snvfilter_rows;
my %snv_dbid;
my %seen_chr;

my $homozygous_min_snv_frequency   = 0.8;
my $heterozygous_min_snv_frequency = 0.3;

while (<$PF_SNV>) {

	print STDERR "Processing base " . $print_counter . "\n" if ( $print_counter % 10000000 == 0 && $print_counter != 0 ); 
	$print_counter++;
	chomp;
	my @cols = split /\t/;
	
	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $cols[0];
	}
	
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $cols[1];
	}
	
	#Skip if it failed the filter_vcf step
	next unless exists $snv_info{ $cols[0] }{ $cols[1] } ;	

	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $cols[1];
	}
	my $ref_base = uc( $cols[2] );

	if ( $ref_base eq 'N' || $ref_base eq '*' ) {    # skip this base if it is an N or an indel
		    #	print STDERR "-- skip N base\n";
		next;
	}

	my $orig_snv_string = $cols[4];

	$cols[4] =~ s/[\.\,]/$ref_base/g;
	$cols[4] =~ s/[\~\$]//g;
	$cols[4] =~ s/\^.//g;

	while ( $cols[4] =~ /([\+\-])(\d+)/ )
	{
		my $sign   = "\\" . $1;
		my $length = $2;

		$cols[4] =~ s/$sign\d+.{$length}//;
	}

	my $read_depth = $cols[3];

	$cols[4] =~ s/(.)/$1\$delim/g;
	my @bases = split /\$delim/, $cols[4];
	my %bases;

	# call a consensus base for this snv
	foreach my $base (@bases)
	{
		next if $base eq '*';    #Skip deletions annotated this way
		$bases{ uc($base) }++;
	}
	my @sorted_bases = sort { $bases{$b} <=> $bases{$a} } keys %bases;

	my $ref_base_freq =
	  defined $bases{$ref_base} ? sprintf( "%.3f",$bases{ uc($ref_base) } / $read_depth) : '0.0';

	my @variant_bases;
	foreach my $sorted_base (@sorted_bases) {
		if ( defined $ref_base && ( $sorted_base eq '.' || uc($sorted_base) eq $ref_base ) ) {
			next;
		} else {
			push @variant_bases, uc($sorted_base);
		}
	}

	if ( scalar @variant_bases == 0 ) {    # Not a variant
		next;
	}

	my $commonest_variant_base      = $variant_bases[0];
	my $next_commonest_variant_base = $variant_bases[1] if scalar @variant_bases > 1;

	my $commonest_variant_base_freq = sprintf( "%.3f", $bases{$commonest_variant_base} / $read_depth );
	my $next_commonest_variant_base_freq = sprintf( "%.3f", $bases{$next_commonest_variant_base} / $read_depth ) if defined $next_commonest_variant_base;

	my $commonest_base = $sorted_bases[0];
	my $var_base       = uc($commonest_variant_base);

	$cols[4] = $orig_snv_string;

	$count++;

	$cols[5] =~ s/(.)/$1\t/g;
	my @chars = split /\t/, $cols[5];

	#    pop @chars;

	my @qualities;
	for ( my $i = 0 ; $i < scalar @chars ; $i++ ) {
		$qualities[$i] = ord( $chars[$i] ) - 33;
	}

	$cols[5] = join( "\:", @qualities );

	#     Homozygous snvs

	my $zygosity;

	if ( defined $commonest_variant_base_freq && $commonest_variant_base_freq >= $homozygous_min_snv_frequency ) {
		$zygosity = 'hom';
	} elsif ((defined $ref_base_freq
			&& defined $commonest_variant_base_freq
			&& ( $ref_base_freq + $commonest_variant_base_freq ) >=
			$homozygous_min_snv_frequency
			&& $ref_base_freq >= $heterozygous_min_snv_frequency
			&& $commonest_variant_base_freq >= $heterozygous_min_snv_frequency )
		  || (defined $next_commonest_variant_base_freq
			  && defined $commonest_variant_base_freq
			  && ( $commonest_variant_base_freq +
				   $next_commonest_variant_base_freq ) >=
			  $homozygous_min_snv_frequency
			  && $next_commonest_variant_base_freq >=
			  $heterozygous_min_snv_frequency
			  && $commonest_variant_base_freq >= $heterozygous_min_snv_frequency)){
		$zygosity = 'het';
	} else {
		$zygosity = 'other';
	}

	#Get the quality score
	my $base_change_key = $cols[2] . '->' . $var_base;
	if ( !exists $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key} ) {
		#Here mpileup calls a different varbase then our pileup module so go with the mpileup change
		for my $mpileup_base_change (keys %{ $snv_info{ $cols[0] }{ $cols[1] } } ) {
			($var_base) = $mpileup_base_change =~ /(\S)$/;
		}
		$base_change_key = $cols[2] . '->' . $var_base;
	}

	if (! exists $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{quality}) {
		modules::Exception->throw("ERROR: $cols[0] $cols[1] $base_change_key\n");
	}

	my $quality_score = $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{quality};
	my $clr_str = $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{clr};
	my $snv_class = $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{class};

	### Write to database
	my $median = modules::Utils->median( \@qualities );

	#Different depending on whether tumour or normal
	my $normal_base_str;
	my $tumour_base_str;
	
	#Slightly different for normal and tumour samples
	if ($tumour_flag) {
		#Get the normal base str
		$tumour_base_str = $cols[4];
		my @lines = split("\n",`grep -w $cols[1] $pileup_normal_file_snv`);
		
		if (@lines) {
			my @fields = split("\t",$lines[0]);
			$normal_base_str = $fields[4];
		} 
		
	} else {		
		$tumour_base_str = 'N/A';
		$normal_base_str = $cols[4];
	
	}
	
	my  $snv_str = join( "\t",
						$cols[0],
						$cols[1],
						$cols[1],
						'SNV'. $delim . $cols[2] . '->'. $var_base
						  . $delim
						  . $snv_class
						  . $delim . 'Q'
						  . $quality_score 
						  . $delim 
						  . $clr_str 
						  . $delim . 'D' 
						  . $read_depth 
						  . $delim 
						  . $zygosity
						  . $delim
						  . $bases{$commonest_variant_base} . '/'
						  . $read_depth
						  . $delim . 'MQ' 
						  . $median 
						  . $delim . 'TBS:'
						  . $tumour_base_str
						  . $delim . 'NBS:'
						  . $normal_base_str
						  . $delim
						  . $cols[5]
						  );

	print SNVOUT $snv_str . "\n";


	#Make sure there's no single quotes as it will break the mysql inserts
	$normal_base_str =~ s/'//g;
	$normal_base_str =~ s/"//g;
	$tumour_base_str =~ s/'//g;
	$tumour_base_str =~ s/"//g;

	
	
	
	if ($write_to_db) {
		
		#default
		my $clr_score = 0;
		
		#if it's a tumour sample
		if ($clr_str =~ /(\d+)/) {
			$clr_score = $1;
		}
		
		my %snv = (
					'run_id'               => $runid,
					'chr'                  => $cols[0],
					'coord'                => $cols[1],
					'ref_base'             => $cols[2],
					'var_base'             => $var_base,
					'ref_base_freq'        => $ref_base_freq,
					'var_base_freq'        => $commonest_variant_base_freq,
					'median_quality_score' => $median,
					'snv_score'            => $quality_score,
					'clr_score'            => $clr_score,
					'tumour_base_string'          => $tumour_base_str,					
					'normal_base_string'          => $normal_base_str,
					'read_depth'           => $read_depth,
					'snv_class'			   => $snv_class
		);

		$seen_chr{ $cols[0] }++;
		push @snv_rows, \%snv;

		#my $inserted_snv = modules::Adaptors::snv->insert(\%snv);

		my %snv_filter = (
			'filter_id'   => $filter_snv->id,
			'filtermatch' => 1,
			'filterpass'  => 1,
			'chr'         => $snv{chr},         # These will need to be removed
			'coord'       => $snv{coord},       # later on; use for lookup
		);

		push @snvfilter_rows, \%snv_filter;

	} else {
		print "$snv_str";
	}
}

#Insert the snv_filters in batch mode

if ( $write_to_db && @snv_rows ) {
	# Insert snv rows
	modules::Adaptors::BulkInsert->insert( -table_name => 'snvs',
										   -data       => \@snv_rows );

	# Get snv database ids, now that they are in the database
	foreach my $seen_chr ( keys %seen_chr ) {
		my $snv_iterator = modules::Adaptors::SNV->search( run_id => $runid, chr => $seen_chr );
		while ( my $snv_from_db = $snv_iterator->next ) {
			$snv_dbid{ $snv_from_db->chr }->{ $snv_from_db->coord } =
			  $snv_from_db->id;
		}
	}

	# Set snv_filter snv_ids using our snv_dbid lookup
	for ( my $i = 0 ; $i < scalar @snvfilter_rows ; $i++ ) {
		my $chr   = $snvfilter_rows[$i]->{chr};
		my $coord = $snvfilter_rows[$i]->{coord};

		if ( defined $snv_dbid{$chr} && defined $snv_dbid{$chr}->{$coord} ) {
			$snvfilter_rows[$i]->{snv_id} = $snv_dbid{$chr}->{$coord};
			delete $snvfilter_rows[$i]->{chr};
			delete $snvfilter_rows[$i]->{coord};
		} else {
			modules::Exception->throw("snvfilter_row [chr:$chr,coord:$coord] points to a non-existant snv.  Something bad has happened." );
		}
	}

	# Insert snv_filter rows
	modules::Adaptors::BulkInsert->insert( -table_name => 'snvs_filters',
										   -data       => \@snvfilter_rows );
}

close($PF_SNV);
close(SNVOUT);

#Now indels...

#First parse the pileup to get the base info; do this first as we don't always see support from mpileup output
my %indel_info = ();
open( my $PF_INDEL, $pileup_file_indel ) || modules::Exception->throw("Unable to open file [$pileup_file_indel]");
while(<$PF_INDEL>) {
	print STDERR "Processing base " . $print_counter . "\n" if ( $print_counter % 10000000 == 0 && $print_counter != 0 ); 
	$print_counter++;
	chomp;

	my ($chr,$start,$ref,$total_depth,$base_str,$qual_str) = split /\t/;
	
	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $chr;
	}
	
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $start;
	}
	
	#Account for start of indel being one more than vcf coordinate
	my $indel_start = $start + 1;
	$indel_info{$chr}{$indel_start}{base_str} = $base_str;
	$indel_info{$chr}{$indel_start}{qual_str} = $qual_str;
	$indel_info{$chr}{$indel_start}{total_reads} = $total_depth;
}




my @insertion_rows        = ();
my @deletion_rows = ();
my @indel_filter_rows = ();
my %seen_chr_indel          = ();
my %indel_ids         = ();
my $indel_flag = 0;


open( INDEL, "$input_indel" ) || modules::Exception->throw("Can't open file $input_indel\n");

while (<INDEL>) {
	chomp;
	next unless /^[0-9XYM]/;
	my ($chr,$start,$end,$rest) = split("\t");
	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $chr;
	}
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $start;
	}
	
	$seen_chr_indel{$chr}++;
	
	my ( $type, $bases, $quality, $clr_str, $indel_class, $depth ) = split( '\^\^\^', $rest );
	$quality =~ s/Q//;
	$depth =~ s/D//;
	
	my $zygosity = 'N/A'; 
	my $mpileup_str = 'N/A';
	my $quality_str = 'N/A';
	my $ref_freq = 0;
	my $var_freq = 0;
	my $var_count = 0;
	my $median = 0;
	my $total_reads = 0;

	if (exists $indel_info{$chr}{$start}) {
		$mpileup_str = $indel_info{$chr}{$start}{base_str};
		$quality_str = $indel_info{$chr}{$start}{qual_str};
		
		my $mpileup_pattern;
		my $length = length($bases);
		if ($type eq 'DEL') {
			$mpileup_pattern = '-'.$length.$bases;
		} else {
			$mpileup_pattern = '+'.$length.$bases;
		}
		
		my $pattern = quotemeta($mpileup_pattern);
		
		#Count up the indel occurrences
		while ($mpileup_str =~ /$pattern/gi) { 
			$var_count++
		}

		my $ref_count = 0;
		for my $char (split("",$mpileup_str)) {
			$ref_count++ if $char eq ',' || $char eq '.';
		}
		
		#refs are listed as .-2TT or ,-2TT so subtract these
		$ref_count = $ref_count - $var_count;
		
		$total_reads = $indel_info{$chr}{$start}{total_reads};
		$ref_freq = $ref_count == 0?0.0:sprintf( "%.3f",$ref_count/$total_reads);
		$var_freq = $var_count == 0?0.0:sprintf( "%.3f",$var_count/$total_reads);

		if ($var_freq > $homozygous_min_snv_frequency) {
			$zygosity = 'hom';
		} elsif ($var_freq > $heterozygous_min_snv_frequency) {
			$zygosity = 'het';
		}  else {
			$zygosity = 'other';
		}

		my @chars = split("",$indel_info{$chr}{$start}{qual_str});

	#    pop @chars;

		my @qualities;
		for ( my $i = 0 ; $i < scalar @chars ; $i++ ) {
			$qualities[$i] = ord( $chars[$i] ) - 33;
		}
		$median = modules::Utils->median( \@qualities );
		$quality_str = join( "\:", @qualities );
		
		
	} else {
		print STDERR "No pileup info for $chr $start $rest\n";
	}

	#Different depending on whether tumour or normal
	my $normal_base_str;
	my $tumour_base_str;
	
	#Slightly different for normal and tumour samples
	if ($tumour_flag) {
		#Get the normal base str
		$tumour_base_str = $mpileup_str;
		my $pileup_start = $start -1 ;
		my @lines = split("\n",`grep -w $pileup_start $pileup_normal_file_indel`);
		
		if (@lines) {
			my @fields = split("\t",$lines[0]);
			$normal_base_str = $fields[4];
		} 
		
	} else {		
		$tumour_base_str = 'N/A';
		$normal_base_str = $mpileup_str;
	
	}
	
	
	my $indel_str = join( "\t",
								$chr,
								$start,
								$end,
								$rest . 
								$delim . $zygosity . 
								$delim . $var_count .'/' . $total_reads .
								$delim . 'MQ' . $median . 
								$delim . 'TBS:' . $tumour_base_str .
								$delim . 'NBS' . $normal_base_str .
								$delim . $quality_str
								) . "\n";
	
	print INDELOUT "$indel_str";
	
	if ($write_to_db) {
		$normal_base_str =~ s/'//g;
		$normal_base_str =~ s/"//g;
		$tumour_base_str =~ s/'//g;
		$tumour_base_str =~ s/"//g;
		
		#default
		
		my $clr_score = 0;
		if ($clr_str =~ /(\d+)/) {
			$clr_score = $1;
		}	
		
		my %indel = (
				  'run_id'         => $runid,
				  'chr'            => $chr,
				  'start_coord'    => $start,
				  'end_coord'      => $end,
				  'var_score'      => $quality,
				  'clr_score'      => $clr_score,
				  'read_depth'     => $depth,
				  'ref_base_freq'  => $ref_freq,
				  'var_base_freq'  => $var_freq,
				  'normal_base_string' => $normal_base_str,
				  'tumour_base_string' => $tumour_base_str,
				  'median_quality_score' => $median,
				  'var_type'       => $type,   
				  'var_caller'     => 'samtools_indel',
				  'variant_class'  => $indel_class
				);

		my %indel_filter = (
						'filtermatch' => 1,
						'filterpass'  => 1,
						'chr' => $chr,    # These will need to be removed
						'start_coord' => $start,
						'end_coord'   => $end
						);
		
		if ($type eq 'DEL') {
			$indel_filter{filter_id} = $filter_deletion->id;
			push @deletion_rows, \%indel;
		} elsif ($type eq 'INS') {
			$indel{'inserted_bases'} = $bases;
			push @insertion_rows, \%indel;
			$indel_filter{'filter_id'} = $filter_insertion->id;
		} else {
			modules::Exception->throw("ERROR: Indel $chr $start has var type $type");
		}
		
		push @indel_filter_rows, \%indel_filter;
		$indel_flag = 1;
	} else {
		print "$indel_str";
	}
	
	
}


#If there are entries to write
if ( $indel_flag && $write_to_db )
{

	#print Dumper \@indel_rows;

	if (@insertion_rows) {
		# Insert indel rows
		modules::Adaptors::BulkInsert->insert( -table_name => 'variants',
										   -data       => \@insertion_rows );
	}
	
	if (@deletion_rows) {
		# Insert indel rows
		modules::Adaptors::BulkInsert->insert( -table_name => 'variants',
										   -data       => \@deletion_rows );
	}
	
	# Get indel database ids, now that they are in the database
	foreach my $seen_chr ( keys %seen_chr_indel )
	{
		my $indel_iterator =
		  modules::Adaptors::Variant->search( run_id => $runid,
											  chr    => $seen_chr );
		while ( my $indel_from_db = $indel_iterator->next )
		{
			$indel_ids{ $indel_from_db->chr }{ $indel_from_db->start_coord }
			  { $indel_from_db->end_coord } = $indel_from_db->id;
		}
	}

	# Set variant_filter variant_ids using our variant_db lookup
	for ( my $i = 0 ; $i < scalar @indel_filter_rows ; $i++ )
	{
		my $chr         = $indel_filter_rows[$i]->{chr};
		my $start_coord = $indel_filter_rows[$i]->{start_coord};
		my $end_coord   = $indel_filter_rows[$i]->{end_coord};

		if (    defined $indel_ids{$chr}
			 && defined $indel_ids{$chr}{$start_coord}
			 && defined $indel_ids{$chr}{$start_coord}{$end_coord} )
		{
			$indel_filter_rows[$i]->{variant_id} =
			  $indel_ids{$chr}{$start_coord}{$end_coord};

			#Remove the lookup keys
			delete $indel_filter_rows[$i]->{chr};
			delete $indel_filter_rows[$i]->{start_coord};
			delete $indel_filter_rows[$i]->{end_coord};
		}
		else
		{
			modules::Exception->throw(
				   "indel_filter_row [chr:$chr,coord:$start_coord] points to a "
					 . "non-existant indel.  Something bad has happened." );
		}
	}

	# Insert indel_filter rows
	modules::Adaptors::BulkInsert->insert( -table_name => 'variants_filters',
										   -data       => \@indel_filter_rows );
}

close INDELOUT;
close $PF_INDEL;
