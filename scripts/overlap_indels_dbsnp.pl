#! /usr/bin/perl -w

use strict;
use modules::Adaptors::Variant;
use modules::Adaptors::Filter;
use modules::Adaptors::Variant_Filter;
use modules::Overlap;
use modules::Exception;
use modules::Adaptors::BulkInsert;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "writeDB=i",
		   "coord_file=s",
		   "overlap_outfile=s",
		   "ref_file=s",
		   "indel_type=s",
		   "runid=i",
		   "chr=s",
		   "match=s"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{coord_file} || !$OPT{overlap_outfile} || !$OPT{ref_file});

	   
=pod

=head1 SYNOPSIS

overlap_indels_dbsnp.pl -match determine_if_match_means_pass_or_fail(options=PASS/FAIL;default=FAIL) -indel_type comma_delim_list_of_Var_types(default=INS,DEL)-runid runid -ref_filtername db_filtername -overlap_outfile output_overlap_file -ref_file input_overlap_file -coord_file formatted_snp_file -chr chrname -writeDB write_results_to_DB(default=0) [options]

Required flags: -coord_filtername -coord_file -runid -overlap_outfile

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

overlap_indels_dbsnp.pl -> Script to run a overlaps for indels for dbsnp

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE



=cut

my $coord_filter = defined $OPT{coord_filtername}?$OPT{coord_filtername}:'filter_dbsnp_indel';

#Dbsnp snv flag requires us to insert the allele freq into the database if the info is there
my $dbsnp_indel_flag = 0;
if ($coord_filter =~ /filter_dbsnp_indel/) {
	$dbsnp_indel_flag = 1;
}

my $coord_file = $OPT{coord_file};
my $runid = $OPT{runid};
my $chromosome = defined $OPT{chr}?$OPT{chr}:'';
my $match = defined $OPT{match}?$OPT{match}:'FAIL';
my $filterpass;
if ($match ne 'PASS' && $match ne 'FAIL') {
	modules::Exception->throw("Match parameter must be set to PASS or FAIL (default=FAIL)");
} else {
	$filterpass = $match eq 'FAIL'?0:1;
}

my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;

my $indel_type = defined $OPT{indel_type}?$OPT{indel_type}:'INS,DEL';

my $overlap_args;

my $overlap = modules::Overlap->new();

#Need to create tmp files for overlap mo
my $outfile_match = $OPT{overlap_outfile};
#Now call for nonmatching; need to change some of the arguments
(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;

my $ref_file =  $OPT{ref_file};

#print "REF $ref_file OUT $outfile_match\n";
my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
#Build up the arguments for the overlap call
my %args = (
			-ref=>$ref_file,
			-coord=>$coord_file,
			-silent=>1, #Don't need the file
			-exact=>1 #Only allow exact coordinate matches
			);

if ($OPT{chr}) {
	$args{-chr} = $chromosome;
	#If step died in the middle remove the appended file
	if ($OPT{chr} eq 1) {
		
		system("rm $outfile_match") if -e $outfile_match;
		system("rm $outfile_nomatch") if -e $outfile_nomatch;
	}
} 



#Call for matching
($overlap_match,$match_count) = $overlap->overlap(%args);

$args{-fail} = 1;
($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);


my ($filter) = modules::Adaptors::Filter->search('name' => $coord_filter);

unless ( defined $filter ) {
    modules::Exception->throw("Unable to retrieve experiment and filter objects from database");
}

# Fetch INDELS from database

my %search_params = ('run_id' => $runid);

if ($OPT{chr}){
	$search_params{'chr'} = $chromosome;	
}

my @db_indels;
my @indel_types = split(",",$indel_type);

#Get only the relevant variant types for overlap
for my $indel_type ( @indel_types ) {
	$search_params{'var_type'} = $indel_type;
    my @indel_type = modules::Adaptors::Variant->search(\%search_params);
    push (@db_indels,@indel_type);
	delete $search_params{'var_type'};
}


# Sort INDELs

my @sorted_indels = sort {$a->start_coord <=> $b->start_coord} @db_indels;

my @indel_filter_inserts;

#Here we overwrite the default match/nomatch files as coordinate matching isn't enough with indel
open(MATCH,">>$outfile_match") || modules::Exception->throw("Can't open file to write $outfile_match\n");
open(NOMATCH,">>$outfile_nomatch") || modules::Exception->throw("Can't open file to write $outfile_nomatch\n");

foreach my $newindel  (@sorted_indels) {
	my %indel_filter = (
    					variant_id    => $newindel->id,
		      			filter_id => $filter->id,
		      		  );
		
	if (exists $overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}) {
		
		#If the indel is potentially matched set as non matched by default and change once we find the exact match
		$indel_filter{'filtermatch'} = 0;
		$indel_filter{'filterpass'} = $filterpass == 0?1:0;
		
		my $dbsnp_match = 0;
		
		#Check if allele freq is available for snvs
	    for my $indel_entry ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}} ) {
	    	next if $dbsnp_match; #As we're only matching by coordinates quit if we've found a match 
    		my ($indel_type,$indel_bases,$qual) = split('\^\^\^',$indel_entry);
    		#If there are multiple entries		
			my @matches = ();
			my %matches = ();
			my %rs = ();
			my $rs;

	    	for my $match_entry (@{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}{$indel_entry}}) {
	    		#Need to confirm that indel matches wrt all indel type and bases
	    		my ($match_type,$local_rs,$match_bases) = split('\^\^\^',$match_entry);
	    		next unless $match_type eq $indel_type;
	    		
    			#For dels we need to ensure length is the same; we do this by checking the bases match
    			my ($just_match_bases) = $match_bases =~ /^([ATCG]+)/;
    			if ($match_type eq 'DEL') {    				
    				next unless $just_match_bases eq $indel_bases;
    			}
    			$indel_filter{'filtermatch'} = 1;
				$indel_filter{'filterpass'} = $filterpass == 1?1:0;
				#Set the rs as matched even if there are no allele frequencies reported
				$rs = $local_rs;
				
				#If there's allele frequency info
		    	if ($match_bases =~ /\(([0-9\.]+):([0-9\.]+)\)/) {
		    		if ($match_type eq 'INS' && $just_match_bases eq $indel_bases) {
		    			#make sure the inserted bases match
						$matches{$1} = $2;
						$rs{$1} = $local_rs;
		    		} elsif ($match_type eq 'DEL') {
		    			#deletion we already know bases match
		    			$matches{$1} = $2;
		    			$rs{$1} = $local_rs;
		    		} 
		    	} 
		    	
				push @matches, $match_entry;
				$dbsnp_match = 1;
	    	}
	    	
		    if (keys %matches) {
		    	#Get the highest ref_allele value for multiple matches; we use allele freq in determining pass
		    	my ($match_ref_allele) = reverse sort {$a<=>$b} keys %matches; 
		    	$indel_filter{'attribute'} = 'rs='.$rs{$match_ref_allele} . ';ref_allele_freq=' . $match_ref_allele . ';var_allele_freq=' . $matches{$match_ref_allele};
		    } elsif ($dbsnp_match) {
		    	$indel_filter{'attribute'} = 'rs='.$rs;
		    } 
	    	
	    	if (@matches) {
		    	my $match_entries = join(",",@matches);
				print MATCH join("\t",
									$newindel->chr,
									$newindel->start_coord,
									$newindel->end_coord,
									$indel_entry,
									$match_entries
									) . "\n";
	    	}
	    	
	    }
		
		
		#Never found an allele and indel type match so write to nomatch file
		if (!$dbsnp_match) {
			my ($indel_entry) = keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}};
			print NOMATCH join("\t",
								$newindel->chr,
								$newindel->start_coord,
								$newindel->end_coord,
								$indel_entry,
								'NO_MATCH'
								) . "\n";
		} 
		
	} elsif (exists $overlap_nomatch->{FAIL}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}) {
		my ($indel_entry) = keys %{$overlap_nomatch->{FAIL}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}};
		print NOMATCH join("\t",
								$newindel->chr,
								$newindel->start_coord,
								$newindel->end_coord,
								$indel_entry,
								'NO_MATCH'
								) ."\n";
		#If the indel didn't match
		$indel_filter{'filtermatch'} = 0;
		$indel_filter{'filterpass'} = $filterpass == 0?1:0;
	} else {
		my $snp_str = $newindel->chr." ".$newindel->start_coord;
		print  modules::Exception->throw("ERROR: Can't find snp $snp_str in the fail or pass list\n");
	}
		
		
	if (!exists $indel_filter{attribute}) {
		$indel_filter{attribute} = 'N/A';
	}
	
	#Only insert snp_filter entries that pass or are dbsnp entries
	push @indel_filter_inserts, \%indel_filter;
	
}

close MATCH;
close NOMATCH;

modules::Adaptors::BulkInsert->insert(-table_name=>'variants_filters', -data=>\@indel_filter_inserts) if $write_to_db && @indel_filter_inserts;
print STDERR "Total match count: $match_count\n";
print STDERR "Total nomatch count: $nomatch_count\n";

