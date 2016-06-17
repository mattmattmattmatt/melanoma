#! /usr/bin/perl -w

use strict;
use modules::VEP;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::BulkInsert;
use modules::SystemCall;
use modules::Pipeline;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Env qw($ENSEMBL_REGISTRY);
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "overlap_outfile_snv=s",
		   "ref_file_snv=s",
		   "exonns_filter_name=s",
		   "tmpdir=s",
		   "runid=i",
		   "writeDB=i",
		   "chr=s",
		   "ref=s"
	    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{overlap_outfile_snv} || !$OPT{ref_file_snv} || !$OPT{runid} || !$OPT{tmpdir});


=pod

=head1 SYNOPSIS

filter_exon_ns.pl -ref_file_snv <overlap_infile_snv> -overlap_outfile_snv <overlap_outfile_snv> -chr chr_to_run_on -exonns_filter_name <filter_name for non-synonymous variants, eg. filter_exon_ns> -writeDB 1|0 -ref reference_genome(default=GRCh37)

Required flags: -runid -overlap_outfile_snv -ref_file_snv -tmpdir

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_exon_ns.pl -> Script to drive variant_effect_predictor.pl and parse details of exonic and non-synonymous variants for snvs and indels

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_exon_ns.pl

=cut

# Put command line options into the right places
my $pipe_conf = modules::Pipeline->get_pipe_conf();
my $clus_conf = modules::Pipeline->get_cluster_conf();


my $exonns_file = $OPT{overlap_outfile_snv};
#Need exon file to get ensembl codon information
(my $exon_file = $exonns_file) =~ s/_ns//;

my $ref_file = $OPT{ref_file_snv};
my $vep_executable = $pipe_conf->read('binaries','variant_predictor','binary');
my $vep_db_dir = $clus_conf->read('svn','conf_dir') . '/vep_index';
my $working_dir = $OPT{tmpdir};
my $runid = $OPT{runid};
my $exonns_filter_name = defined $OPT{exonns_filter_name}?$OPT{exonns_filter_name}:'filter_exon_ns';
my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;
my $ref = defined $OPT{ref}?$OPT{ref}:"GRCh37";
my $organism;

if ($ref =~ /GRCh/) {
	$organism = 'human';
} else {
	$organism = 'mouse';
}


my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	print "ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n";
	exit;
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $vep_input_file = $working_dir . '/vep.input';

open(my $EXONNS_FILE, $ref_file)
    or modules::Exception->throw("Unable to open input file [$ref_file]");
open(my $EXON_FILE, $exon_file)
    or modules::Exception->throw("Unable to open input file [$exon_file]");
        
    
open(my $VEP_INPUT_FILE, ">$vep_input_file")
    or modules::Exception->throw("Unable to open file for writing annovar input format [$vep_input_file]");

#Connect to database

my ($run) = modules::Adaptors::Run->search('id' => $runid);
my ($exonns_filter) = modules::Adaptors::Filter->search('name' => $exonns_filter_name);

unless (defined $run
	&& defined $exonns_filter){
    modules::Exception->throw("Did not get all necessary things from database [run_id:$runid][exonns:$exonns_filter_name]");
}

my %vep_input_info = ();

while (<$EXONNS_FILE>){
    chomp;

	#db_variants input file
	#1	55776	55776	SNV^^^C->T^^^Q202^^^TUM81^^^D34^^^het^^^16/34^^^MQ40^^^TBS:,TTt,,,.,T,,.,,tT,.T.Ttt.TtTT.tT,,NBS:.,.,.,.,.,..,,,,,....,,,.,.^^^41:40:41:41:40:41:41:41:41:41:34:34:41:41:40:38:34:39:41:41:41:38:32:37:41:41:36:39:39:37:28:26:34:35

    my ($chr, $start, $end, $attributes) = split /\t/;
	my (undef,$base_change) = split('\^\^\^',$attributes);
	my ($ref_base,$var_base) = $base_change =~ /(.*)->(.*)/;


    if (!defined $chr || ! defined $start || ! defined $end || ! defined $attributes){
		modules::Exception->warning("Did not parse input line successfully [$_]");
    }

    unless ($var_base =~ /[ATGC]/) {
		modules::Exception->warning("Did not parse attributes successfully [$attributes]");
    }

	my $base_str = $ref_base . '/' . $var_base;
	$vep_input_info{$chr}{$start} = join("\t", $chr, $start, $end, $base_str, '+' );
    #print $VEP_INPUT_FILE join("\t", $chr, $start, $end, $base_str, '+' ) . "\n";
}

close($EXONNS_FILE);

#Now get the peptide info from ensembl to see if we need to combine any variants
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($ENSEMBL_REGISTRY);

# get the adaptors
my $gene_adaptor = $registry->get_adaptor($organism, "Core", "Gene");

my %pep_to_coord = ();
my %coord_to_pep = ();
my %aa_length = ();

while (<$EXON_FILE>) {
	chomp;
	my @fields = split("\t");
	next unless /ENS/;
	my ($ens_gene_name) = $fields[4] =~ /(ENSG\d+)/;
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
	$aa_length{$ens_gene_name} = $transcript_length/3 - 1; # Account for stop codon
	my $strand = $canonical_transcript->strand();
	my $transcriptMapper = $canonical_transcript->get_TranscriptMapper();
	my @pepCoords = $transcriptMapper->genomic2pep($fields[1], $fields[1], $strand);
	foreach my $pepinfo (@pepCoords){
		push @{$pep_to_coord{$ens_gene_name . ':' . $pepinfo->start}}, "$fields[0]:$fields[1]";
		$coord_to_pep{$fields[0]}{$fields[1]} = $ens_gene_name . ':' . $pepinfo->start;
	}
	
}

close($EXON_FILE);
my %reported_coords = ();

for my $chr ( keys %vep_input_info ) {
    for my $coord ( keys %{$vep_input_info{$chr}} ) {
        next if exists $reported_coords{"$chr:$coord"};
        if (exists $coord_to_pep{$chr} && exists $coord_to_pep{$chr}{$coord}) {
                if (exists $pep_to_coord{$coord_to_pep{$chr}{$coord}} && @{$pep_to_coord{$coord_to_pep{$chr}{$coord}}} > 1) {
                        #Here we combine the entries from the same codon
                        my $full_ref = my $full_var;
                        my $max_coord =  0;
                        my $min_coord = 100000000000;

                        my @fields = split("\t",$vep_input_info{$chr}{$coord});

                        for my $vep_input ( @{$pep_to_coord{$coord_to_pep{$chr}{$coord}}} ) {
                                my ($vep_chr,$vep_coord) = split(':',$vep_input);
                                if ($vep_coord < $min_coord) {
                                        $min_coord = $vep_coord;
                                }
                                if ($vep_coord > $max_coord) {
                                        $max_coord = $vep_coord;
                                }

                                my @vep_fields = split("\t",$vep_input_info{$vep_chr}{$vep_coord});
                                my ($vep_ref,$vep_var) = split('/',$vep_fields[3]);
                                $full_ref .= $vep_ref;
                                $full_var .= $vep_var;
                                $reported_coords{"$vep_chr:$vep_coord"}++;
                        }

                        my $full_base = $full_ref . '/' . $full_var;

                        print $VEP_INPUT_FILE join("\t",$chr,$min_coord,$max_coord,$full_base,$fields[4]) . "\n";
                        
                        next;
                }
        }
        print $VEP_INPUT_FILE $vep_input_info{$chr}{$coord} . "\n";
        }
}

close($VEP_INPUT_FILE);


# Run annovar and parse result

my $vep 
    = modules::VEP->new('input_file'      => $vep_input_file,
				'executable_path' => $vep_executable,
				'db_dir'          => $vep_db_dir,
				'working_dir'     => $working_dir,
				'organism'		  => $organism
				);

$vep->run;

my $results = $vep->parse_result;

# Get a list of chromosomes included in the results

my %chrs;
foreach my $result (@$results){
    $chrs{$result->[1]}++;
}



# Work through results by chromosome. Write exon_ns matches to database and intermediate files.
open(my $EXONNS, ">$exonns_file") 
    or modules::Exception->throw("Unable to open file for writing");


my @exonns_inserts = ();

foreach my $chr (keys %chrs){

    my %chr_snvs;

    my $snv_iterator
	= modules::Adaptors::SNV->search(
										  'run_id' => $runid,
					    				  'chr' => $chr
					    				  );

    while (my $snv = $snv_iterator->next){
		$chr_snvs{$snv->coord} = $snv->id;
    }


    foreach my $result (@$results){
		next unless $result->[1] eq $chr;
	
		my ($snv_type,
		    $snv_chr,
		    $snv_start,
		    $snv_end,
		    $snv_gene,
		    $snv_exon,
		    $snv_ref_aa,
		    $snv_var_aa,
		    $poly_predict,
		    $poly_score,
		    $sift_predict,
		    $sift_score) = @$result;
	
		
		my $info_field = join('^^^',$snv_gene, $snv_exon, $poly_predict, $poly_score, $sift_predict, $sift_score);
		
		my $combined = $snv_start == $snv_end?0:1;
		my $snv_count = $snv_start;
		
		if ($snv_type =~ 'missense_variant' || $snv_type =~ 'stop_gained'){ # have a nonsynonymous change, gain a Stop or loose a Stop
		
		
#			if (!defined $coord_to_pep{$snv_chr}{$snv_count}) {
#				modules::Exception->throw("ERROR: No aa_pos for $snv_chr $snv_count");
#			}
		
			while ($snv_count <= $snv_end) {
			
				my $attribute;
				if ($combined) {
					$attribute = 'gene=' . $snv_gene 
						                    . ';exon=' . $snv_exon 
						                    . ';aa_change=' . $snv_ref_aa . '->' . $snv_var_aa
						                    . ';poly_pred='. $poly_predict
						                    . ';poly_score='. $poly_score
						                    . ';sift_pred='. $sift_predict
						                    . ';sift_score='.$sift_score
						                    . ';combined='.$snv_start.'-'.$snv_end;
					print $EXONNS 
								join("\t", 
				     			$snv_chr, 
				     			$snv_count, 
				     			$snv_count,
				     			'') 
								. join('^^^', 
				       					$info_field,
				       					$snv_ref_aa . '->' . $snv_var_aa, 'combined') . "\n";
				} else {
					$attribute = 'gene=' . $snv_gene 
						                    . ';exon=' . $snv_exon 
						                    . ';aa_change=' . $snv_ref_aa . '->' . $snv_var_aa
						                    . ';poly_pred='. $poly_predict
						                    . ';poly_score='. $poly_score
						                    . ';sift_pred='. $sift_predict
						                    . ';sift_score='.$sift_score;
						                    
					print $EXONNS 
								join("\t", 
				     			$snv_chr, 
				     			$snv_count, 
				     			$snv_count,
				     			'') 
								. join('^^^', 
				       					$info_field,
				       					$snv_ref_aa . '->' . $snv_var_aa) . "\n";
				}
			    
			    my %snv_filter_exonns = (
			    					'filter_id' => $exonns_filter->id,
						     		'snv_id'    => $chr_snvs{$snv_start},
						     		'filtermatch' => 1,
					       	 		'filterpass' => 1,
						     		'attribute' => $attribute
						                    );
		
				push @exonns_inserts, \%snv_filter_exonns;
		    	$snv_count++;
			   
			}
		}
    }
}

modules::Adaptors::BulkInsert->insert(-table_name=>'snvs_filters', -data=>\@exonns_inserts) if $write_to_db && @exonns_inserts;


close($EXONNS);


