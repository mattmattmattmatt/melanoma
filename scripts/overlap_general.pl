#! /usr/bin/perl -w

use strict;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;;
use modules::Adaptors::Filter;
use modules::Adaptors::Variant_Filter;
use modules::Overlap;
use modules::Exception;
use modules::Adaptors::BulkInsert;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Env qw($ENSEMBL_REGISTRY);
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "writeDB=i",
		   "coord_file_snv=s",
		   "coord_file_indel=s",
		   "overlap_outfile_snv=s",
		   "overlap_outfile_indel=s",
		   "ref_file_snv=s",
		   "ref_file_indel=s",
		   "coord_filtername=s",
		   "runid=i",
		   "chr=s",
		   "indel_type=s",
		   "match=s",
		   "snv",
		   "exact",
		   "indel",
		   "add_attribute"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{coord_filtername} || (!$OPT{coord_file_snv} && !$OPT{coord_file_indel}) || (!$OPT{overlap_outfile_snv} && !$OPT{overlap_outfile_indel}) || (!$OPT{ref_file_snv} && !$OPT{ref_file_indel}) );

	   
=pod

=head1 SYNOPSIS

overlap_general.pl -snv overlap_snvs -indel overlap_indels -add_attribute add_matching_value_to_snvs_filters -exact only_match_exact_coord_matches(for_deletions) -indel_type indel_enum_to_search(default=INS,DEL) -match determine_if_match_means_pass_or_fail(options=PASS/FAIL;default=FAIL) -runid runid -overlap_outfile_snv snv_output_overlap_file -overlap_outfile_indel indel_output_overlap_file -ref_file_snv snv_input_overlap_file -ref_file_indel indel_input_overlap_file -coord_file_snv formatted_snv_file -coord_file_indel formatted_indel_file -chr chrname -writeDB write_results_to_DB(default=0)[options]

Required flags: -runid -coord_filtername -runid (-overlap_outfile_snv || -overlap_outfile_indel) (-ref_file_snv || -ref_file_indel) (-coord_file_snv || -coord_file_indel)

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

overlap_general.pl -> Script to run a generic overlap

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./overlap_general.pl 

=cut

if (!$OPT{snv} && !$OPT{indel}) {
	modules::Exception->throw("ERROR: Need to overlap either -snv or -indel");
}

my $coord_filter = $OPT{coord_filtername};

if ($OPT{indel} && $coord_filter =~ /filter_dbsnp/) {
	modules::Exception->throw("ERROR: Use overlap_indels_dbsnp.pl to overlap indels and dbsnp");
}

#Do all the common operations first
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

my ($filter) = modules::Adaptors::Filter->search('name' => $coord_filter);
	
unless ( defined $filter ) {
	modules::Exception->throw("Unable to retrieve experiment and filter objects from database");
}

my $filter_exon = 0;
if ($coord_filter eq 'filter_exon') {
	$filter_exon = 1;
}
	
my $filter_splice = 0;
if ($coord_filter eq 'filter_splicesite') {
	$filter_splice = 1;
}	

#Now get the peptide info from ensembl to see if we need to combine any variants
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($ENSEMBL_REGISTRY);

# get the adaptors
my $gene_adaptor = $registry->get_adaptor('human', "Core", "Gene");

if ($OPT{snv}) {

	if (!defined $OPT{coord_file_snv} || !defined $OPT{ref_file_snv} || !defined $OPT{overlap_outfile_snv}) {
		modules::Exception->throw("ERROR: For snv need to define -coord_file_snv, -ref_file_snv, and -overlap_outfile_snv");
	}

	my $filter_dbsnp = 0;
	if ($coord_filter =~ /filter_dbsnp/) {
		$filter_dbsnp = 1;
	}
	
	my $coord_file = $OPT{coord_file_snv};
	
	
	#print "IGNORE STRING: $ignore_string\n";
	
	my $overlap_args;
	
	my $overlap = modules::Overlap->new();
	my $outfile_match = $OPT{overlap_outfile_snv};
	my $ref_file =  $OPT{ref_file_snv};
	
	#print "REF $ref_file OUT $outfile_match\n";
	my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
	#Build up the arguments for the overlap call
	my %args = (
				-ref=>$ref_file,
				-coord=>$coord_file,
				);
	

	if ($OPT{exact}) {
		$args{-exact} = 1;
	}
	
	if ($OPT{chr}) {
		$args{-chr} = $chromosome;
		$args{-append} = $outfile_match;
	} else {
		$args{-output} = $outfile_match;
	}
	
	#Call for matching
	($overlap_match,$match_count) = $overlap->overlap(%args);
	
	
	#Now call for nonmatching; need to change some of the arguments
	(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;
	
	if ($OPT{chr}) {
		$args{-append} = $outfile_nomatch;
	} else {
		$args{-output} = $outfile_nomatch;
	}
	
	$args{-fail} = 1;
	($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);	
	
	# Fetch SNPs from database
	
	my %search_params = ('run_id' => $runid);
	
	if ($OPT{chr}){
		$search_params{'chr'} = $chromosome;	
	}
	
	my @snvs = modules::Adaptors::SNV->search(\%search_params);
		
	# Sort snvs
	
	my @sorted_snvs = sort {$a->coord <=> $b->coord} @snvs;
	
	my @snv_filter_inserts;
	
	foreach my $newsnv  (@sorted_snvs) {
		my %snv_filter = (
	    					snv_id    => $newsnv->id,
			      			filter_id => $filter->id,
			      		  );
		
		#Make sure to set the attribute if it's dbsnp as this is where we record allele frequency
		if ($filter_dbsnp) {
				$snv_filter{'attribute'} = 'N/A';
		}

		if ($filter_exon) {
			$snv_filter{'attribute'} = 'N/A';
		}
		
		if (exists $overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}) {
			#If the snv is matched
			$snv_filter{'filtermatch'} = 1;
			#Set the pass based on the filter context
			$snv_filter{'filterpass'} = $filterpass == 1?1:0;
			
			#Check if allele freq is available for snvs
			if ($filter_dbsnp) {		
				for my $end ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}} ) {
				    	my ($ref_base,$var_base) = $basic_entry =~ /([ACGT])->([ACGT])/;
				    	if (!defined $ref_base) {
				    		print STDERR "ERROR: $basic_entry\n";
				    	}
				    	
				    	#Want to report the highest ref allele freq with multiple entries
				    	my %matches = ();
				    	my %rs = ();
				    	my $rs;
				    	my $allele_freq = my $rs_found = 0;
				    	for my $match_entry (@{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}{$basic_entry}}) {
					    	if ($match_entry =~ /(rs\d+).*([ACGT])\(([0-9\.]+)\)->([ACGT])\(([0-9\.]+)\)/) {
					    		if ($ref_base eq $2 && $var_base eq $4) {
									$matches{$3} = $5;
									$rs{$3} = $1;
									$allele_freq = 1;
					    		}
					    	} elsif ($match_entry =~ /(rs\d+)/) {
					    		#If there are no allele freqs we still want the rs number
					    		$rs = $1;
					    		$rs_found = 1;
					    	} else {
					    		modules::Exception->throw("ERROR: SNV match has no rs number $newsnv->chr $newsnv->coord $match_entry");
					    	}
				    	}
					    if ($allele_freq) {
					    	#Get the highest ref_allele value for multiple matches; we use allele freq in determining pass
					    	my ($match_ref_allele) = reverse sort {$a<=>$b} keys %matches; 
					    	$snv_filter{'attribute'} = 'rs='.$rs{$match_ref_allele} . ';ref_allele_freq=' . $match_ref_allele . ';var_allele_freq=' . $matches{$match_ref_allele};
					    } elsif ($rs_found) {
					    	$snv_filter{'attribute'} = 'rs='.$rs;
					    }
				    }
				}
			} elsif ($filter_exon) {
				for my $end ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}} ) {
				    	for my $match_entry (@{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}{$basic_entry}}) {
					    	if ($match_entry =~ /(ENSG\d+)/) {
					    		my ($ens_gene_name) = $1;
								my $gene = $gene_adaptor->fetch_by_stable_id($ens_gene_name);
								next unless $gene;
								my $canonical_transcript = $gene->canonical_transcript();
								my $transcript_length = 0;
								my @exon_objs = @{$canonical_transcript->get_all_Exons()};
							    foreach my $exon (@exon_objs) {
							    	if ($exon->coding_region_start($canonical_transcript)) {
							    		$transcript_length += abs( $exon->coding_region_end($canonical_transcript) - $exon->coding_region_start($canonical_transcript)) + 1;
							    	}
							    }
							    
								my $strand = $canonical_transcript->strand();
								my $transcriptMapper = $canonical_transcript->get_TranscriptMapper();
								my @pepCoords = $transcriptMapper->genomic2pep($newsnv->coord,$newsnv->coord, $strand);
								if (@pepCoords > 1) {
									modules::Exception->throw("ERROR: SNV at $newsnv->chr $newsnv->coord has more than one peptide coordinate");
								}
								my $aa_length = $transcript_length/3 - 1; # Account for stop codon
								my $aa_pos = $pepCoords[0]->start; 
								$snv_filter{'attribute'} = 'aa_pos='.$aa_pos.';aa_len='.$aa_length;
					    	} else {
					    		modules::Exception->throw("ERROR: Exon entry at $newsnv->chr $newsnv->coord doesn't contain ENSEMBL gene id");
					    	}
				    	}
				    }
				}
			} elsif ($filter_splice) {
				my $splice_distance;
				for my $end ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}} ) {
				    	for my $match_entry (@{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}{$basic_entry}}) {
							($splice_distance) = $match_entry =~ /(\d+)$/;
				    	}
				    }
				}
				$snv_filter{'attribute'} = 'splice_dist='.$splice_distance; 
				
			} elsif ($OPT{add_attribute}) {
				#Here we want to add the matching field as an attibute
				for my $end ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}} ) {
				    	my @matches = @{$overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}{$basic_entry}};
						if (@matches > 1) {			
							my $snv_str = $newsnv->chr." ".$newsnv->coord;
							modules::Exception->throw("ERROR: Shouldn't have more than 1 match for $snv_str\n");
						} else {
							my $match = $overlap_match->{PASS}{$newsnv->chr}{$newsnv->coord}{$end}{$basic_entry}->[0];
							$match =~ s/-DUP//;
							$snv_filter{'attribute'} = $match;
						}
				    }
				}
			}
		} elsif (exists $overlap_nomatch->{FAIL}{$newsnv->chr}{$newsnv->coord}) {
			#If the snv didn't match
			$snv_filter{'filtermatch'} = 0;
			$snv_filter{'filterpass'} = $filterpass == 0?1:0;
		} else {
			my $snv_str = $newsnv->chr." ".$newsnv->coord;
			print  modules::Exception->throw("ERROR: Can't find snv $snv_str in the fail or pass list\n");
		}
			
		#Only insert snv_filter entries that pass or are dbsnp entries
		if ($snv_filter{'filterpass'} == 1 || $filter_dbsnp) {
			#print Dumper \%snv_filter;
			push @snv_filter_inserts, \%snv_filter;
		}	
	}
	
	modules::Adaptors::BulkInsert->insert(-table_name=>'snvs_filters', -data=>\@snv_filter_inserts) if $write_to_db && @snv_filter_inserts;
	print STDERR "Total match count: $match_count\n";
	print STDERR "Total nomatch count: $nomatch_count\n";

}

if ($OPT{indel}) {	
	if (!defined $OPT{coord_file_indel} || !defined $OPT{ref_file_indel} || !defined $OPT{overlap_outfile_indel}) {
		modules::Exception->throw("ERROR: For indel need to define -coord_file_indel, -ref_file_indel, and -overlap_outfile_indel");
	}
	
	my $coord_file = $OPT{coord_file_indel};
	my $indel_type = defined $OPT{indel_type}?$OPT{indel_type}:'INS,DEL';
	
	my $overlap_args;
	
	my $overlap = modules::Overlap->new();
	my $outfile_match = $OPT{overlap_outfile_indel};
	my $ref_file =  $OPT{ref_file_indel};
	
	#print "REF $ref_file OUT $outfile_match\n";
	my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
	#Build up the arguments for the overlap call
	my %args = (
				-ref=>$ref_file,
				-coord=>$coord_file,
				);
	
	#Used for deletion matching in filter_common
	if ($OPT{exact}) {
		$args{-exact} = 1;
	}
	
	if ($OPT{chr}) {
		$args{-chr} = $chromosome;
		$args{-append} = $outfile_match
	} else {
		$args{-output} = $outfile_match;
	}
	
	#Call for matching
	($overlap_match,$match_count) = $overlap->overlap(%args);
	
	#Now call for nonmatching; need to change some of the arguments
	(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;
	
	if ($OPT{chr}) {
		$args{-append} = $outfile_nomatch;
	} else {
		$args{-output} = $outfile_nomatch;
	}
	
	$args{-fail} = 1;
	($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);
		
	# Fetch INDELS from database
	
	my %search_params = ('run_id' => $runid);
	
	if ($OPT{chr}){
		$search_params{'chr'} = $chromosome;	
	}
	
	my @db_indels;
	my @indel_types = split(",",$indel_type);
	
	#Get only the relevant variant types for overlap
	for my $type ( @indel_types ) {
		$search_params{'var_type'} = $type;
	    my @indel_type = modules::Adaptors::Variant->search(\%search_params);
	    push (@db_indels,@indel_type);
		delete $search_params{'var_type'};
	}
	
	# Sort INDELs
	
	my @sorted_indels = sort {$a->start_coord <=> $b->start_coord} @db_indels;
	
	my @indel_filter_inserts;
	
	foreach my $newindel  (@sorted_indels) {
		my %indel_filter = (
	    					variant_id    => $newindel->id,
			      			filter_id => $filter->id,
			      		  );
			
		if (exists $overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}) {
			#With other inputs we can't check the allele so go by start and end coord only
			$indel_filter{'filtermatch'} = 1;
			$indel_filter{'filterpass'} = $filterpass == 1?1:0;
			if ($filter_exon) {
				for my $end ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}} ) {
				    	for my $match_entry (@{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}{$basic_entry}}) {
					    	if ($match_entry =~ /(ENSG\d+)/) {
					    		my ($ens_gene_name) = $1;
								my $gene = $gene_adaptor->fetch_by_stable_id($ens_gene_name);
								next unless $gene;
								my $canonical_transcript = $gene->canonical_transcript();
								my $transcript_length = 0;
								my @exon_objs = @{$canonical_transcript->get_all_Exons()};
							    foreach my $exon (@exon_objs) {
							    	if ($exon->coding_region_start($canonical_transcript)) {
							    		$transcript_length += abs( $exon->coding_region_end($canonical_transcript) - $exon->coding_region_start($canonical_transcript)) + 1;
							    	}
							    }
							    
								my $strand = $canonical_transcript->strand();
								my $transcriptMapper = $canonical_transcript->get_TranscriptMapper();
								my @pepCoords = $transcriptMapper->genomic2pep($newindel->start_coord,$newindel->start_coord, $strand);
								if (@pepCoords > 1) {
									modules::Exception->throw("ERROR: Indel at $newindel->chr $newindel->start_coord has more than one peptide coordinate");
								}
								my $aa_length = $transcript_length/3 - 1; # Account for stop codon
								my $aa_pos = $pepCoords[0]->start; 
								$indel_filter{'attribute'} = 'aa_pos='.$aa_pos.';aa_len='.$aa_length;
					    	} else {
					    		modules::Exception->throw("ERROR: Exon entry at $newindel->chr $newindel->start_coord doesn't contain ENSEMBL gene id");
					    	}
				    	}
				    }
				}
			} elsif ($filter_splice) {
				my $splice_distance;
				for my $end ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}} ) {
				    	for my $match_entry (@{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}{$basic_entry}}) {
							($splice_distance) = $match_entry =~ /(\d+)$/;
				    	}
				    	
				    }
				}
				$indel_filter{'attribute'} = 'splice_dist='.$splice_distance; 
				
			} elsif ($OPT{add_attribute}) {
				#Here we want to add the matching field as an attibute
				for my $end ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}} ) {
				    	my @matches = @{$overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}{$basic_entry}};
						if (@matches > 1) {			
							my $indel_str = $newindel->chr." ".$newindel->start_coord;
							modules::Exception->throw("ERROR: Shouldn't have more than 1 match for $indel_str\n");
						} else {
							$indel_filter{'attribute'} = $overlap_match->{PASS}{$newindel->chr}{$newindel->start_coord}{$end}{$basic_entry}->[0];
						}
				    }
				}
			}
		} elsif (exists $overlap_nomatch->{FAIL}{$newindel->chr}{$newindel->start_coord}{$newindel->end_coord}) {
			#If the indel didn't match
			$indel_filter{'filtermatch'} = 0;
			$indel_filter{'filterpass'} = $filterpass == 0?1:0;
		} else {
			#This should never happen as the indel should either match or not match
			if ($coord_filter =~ /filter_common/) {
				next;
			}
			my $indel_str = $newindel->chr." ".$newindel->start_coord;
			print  modules::Exception->throw("ERROR: Can't find indel $indel_str in the fail or pass list\n");
		}
			
		#Only insert variant_filter entries that pass or are dbsnp entries
		if ($indel_filter{'filterpass'} == 1) {
			push @indel_filter_inserts, \%indel_filter;
		}	
	}
	
	modules::Adaptors::BulkInsert->insert(-table_name=>'variants_filters', -data=>\@indel_filter_inserts) if $write_to_db && @indel_filter_inserts;
	print STDERR "Total match count: $match_count\n";
	print STDERR "Total nomatch count: $nomatch_count\n";
}