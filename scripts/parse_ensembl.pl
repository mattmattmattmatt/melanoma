#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;
use modules::Utils;
use modules::Pipeline;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Env qw($ENSEMBL_REGISTRY_LOAD);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	 	"ref=s",
	   	"debug=s"
	   	);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

parse_ensembl.pl -ref reference_genome(default=GRCh37) [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_ensembl.pl -> Script to generate all the gene files required for the pipeline

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./parse_ensemlb.pl

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $pipe_conf = modules::Pipeline->get_pipe_conf();
my $ref = defined $OPT{ref}?$OPT{ref}:"GRCh37";
my $splice_length = $pipe_conf->read('cutoffs','splice_length_cutoff');

my $ens_link;
my $organism;
if ($ref =~ /GRCh/) {
	$organism = "human";
	$ens_link = "http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g="; 
} else {
	$organism = "mouse";
	$ens_link = "http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=";
}

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($ENSEMBL_REGISTRY_LOAD);

# get the adaptors
my $gene_adaptor = $registry->get_adaptor($organism, "Core", "Gene");
my $slice_adaptor = $registry->get_adaptor($organism, 'Core', 'Slice');
my @chr_slices = @{ $slice_adaptor->fetch_all('chromosome') };
my $go_adaptor = $registry->get_adaptor("Multi","Ontology","GOTerm");

#Stores the GO term mappings to prevent continually looking up the same db value
my %go_mapping = ();

my $conf_dir = "$svndir/conf/$organism/$ref";

#The output files

#File for filter_exon (per chr files)
my $exon_file = $ref . '.exon.overlap';

#Single exon file
my $exon_all_file = $ref . '.exon.overlap.all';

#File for filter_splicesite (per chr files)
my $splice_file = $ref . '.splice.overlap';

#File for exon coords (used for coverage stats)
my $coord_file = $ref . '.exon.coord';
open(COORD,">$coord_file") || modules::modules->throw("Can't open file to write $coord_file\n");

#gene merging file
my $all_gene_info = $ref. ".gene.all";
open(ALLINFO,">$all_gene_info") || modules::modules->throw("Can't open file to write $all_gene_info\n");

#Store all the data for sorting and writing to files
my %gene_data = ();



#URLs for linking
my $cosmic_url = "http://cancer.sanger.ac.uk/cosmic/gene/overview?ln=";
my $omim_link = 'http://omim.org/entry/';

#We only want to report canonical transcripts so keep track of the ones we see
my %transcripts = ();

#map gene names
my %gene_name_mapper = ();

#Use the latest COSMIC
my $cosmic_file;
if ($OPT{cosmic_file}) {
	$cosmic_file = $OPT{cosmic_file};
} else {
	opendir(COSMIC,"$conf_dir/cosmic/") || modules::Exception->throw("ERROR: Cannot open dbsnp directory $conf_dir/cosmic/");
	my @cosmic_files = grep {/^\d/} readdir COSMIC;
	closedir COSMIC;
	my ($cosmic_tmp) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @cosmic_files);
	$cosmic_file = $conf_dir . '/cosmic/'. $cosmic_tmp . '/'.$ref .'.cosmic.gene';
}
if ( ! -e $cosmic_file ) {
	modules::Exception->throw("File $cosmic_file doesn't exist");
}

#Use the latest Vogelstein
my $vogel_file;
if ($OPT{vogel_file}) {
	$vogel_file = $OPT{vogel_file};
} else {
	opendir(VOGEL,"$conf_dir/vogelstein/") || modules::Exception->throw("ERROR: Cannot open dbsnp directory $conf_dir/vogelstein/");
	my @vogel_files = grep {/^\d/} readdir VOGEL;
	closedir VOGEL;
	my ($vogel_tmp) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @vogel_files);
	$vogel_file .= $conf_dir . '/vogelstein/'. $vogel_tmp .'/'.$ref .'.vogelstein.gene';
}
if ( ! -e $vogel_file ) {
	modules::Exception->throw("File $vogel_file doesn't exist");
}

#COSMIC entries are joined by gene name; these non coordinate entries are only reported if there is no coordinate match
open(COSMIC,"$cosmic_file") || modules::Exception->throw("Can't open file $cosmic_file\n");

while (<COSMIC>) {
	chomp;
	next unless /\w/;
	my ($cosmic_gene) = $_;
	$gene_name_mapper{uc($cosmic_gene)}{cosmic} = $cosmic_url . $cosmic_gene;
}


#VOGELSTEIN entries are joined by gene name
open(VOGEL,"$vogel_file") || modules::Exception->throw("Can't open file $vogel_file\n");

while (<VOGEL>) {
	chomp;
	next unless /\w/;
	my ($vogel_gene,$info) = split(" ");
	$gene_name_mapper{uc($vogel_gene)}{vogel} = $info;
}


# traverse chromosomes
for my $chr_slice (@chr_slices) {
	
	#Get the genes for that chromosome
	my @genes = @{$gene_adaptor->fetch_all_by_Slice($chr_slice)};
	my $chr = $chr_slice->seq_region_name();
	next unless $chr =~ /[0-9XYM]/;
	if ($OPT{chr}) {
		next unless $OPT{chr} eq $chr;
	}
	my $count = 0;
	
	
	print "Processing " . scalar(@genes) . " gene IDs for chr $chr...\n";

	#traverse genes
    for my $gene (@genes) {
	    # let user know count
	    local $| = 1;
	    print "[$count/" . scalar(@genes) . "]\r";
	    $count++;
	
		#ENSG id
		my $gene_id = $gene->stable_id();
		
		if ($OPT{debug}) {
			next unless $OPT{debug} eq $gene_id;
		}

	
	    # get canonical transcript
	    my $canonical_transcript = $gene->canonical_transcript();
	
		#only protein coding biotypes	
		my $biotype = $canonical_transcript->biotype();
		#next unless $biotype eq 'protein_coding';
		
		#Stores the sequence info
		#my $sequence;
		my $coding = 0;
	    
		my %uniprot_ids = ();
	
		my @exon_objs = @{$canonical_transcript->get_all_Exons()};
		my $exon_count = 1;
	    foreach my $exon (@exon_objs) {
	    	my $coding_start = $exon->coding_region_start($canonical_transcript);
			next unless defined $coding_start; #Skip non coding exons
			$coding = 1;
			my $coding_end = $exon->coding_region_end($canonical_transcript);
			my $full_ensembl_name = $gene_id . '_exon'.$exon_count;
			push @{$gene_data{$chr}{$coding_start}{$coding_end}},$full_ensembl_name;
			
			
			#get uniprot id from this API call
			my @exon_sf = @{ $exon->get_all_supporting_features() };
			foreach my $sf (@exon_sf) {
				if ($sf->analysis()->display_label() =~ /uniprot/i){
		    		$uniprot_ids{$sf->hseqname()}++;
				}
			}
			
			$exon_count++;
				
	    }
	    
	    next unless $coding;
		
		
		#gene description
		my $desc = $gene->description();

	
		#external db links
		my @dblinks = @{$gene->get_all_DBLinks()};
	
		#hugo, go, ccds, and omim values
		my %hugo = ();
		my %go_terms = ();
		my %ccds = ();
		my %omims = ();
	
	
		#Get OMIM, CCDS, GO TERMS, and HGNC 
	  	foreach my $dbe (@dblinks){
	  		if ($dbe->dbname eq 'HGNC') {
	  			$hugo{$dbe->display_id}++;
	  		} elsif ($dbe->dbname eq 'MGI') {
	  			$hugo{$dbe->display_id}++;
	  		} elsif ($dbe->dbname eq 'CCDS') {
	  			$ccds{$dbe->display_id}++;
	  		} elsif ($dbe->dbname eq 'MIM_GENE') {
	  			my ($omim_id) = $dbe->display_id =~ /\[\*(\d+)/;
	  			$omims{"http://omim.org/entry/$omim_id"}++;
	  		} elsif ($dbe->dbname eq 'GO') {
	  			if (exists $go_mapping{$dbe->display_id}) {
	  				#use the lookup if we've already seen code
	  				$go_terms{$go_mapping{$dbe->display_id}}++;
	  			} else {
					#use the db if never seen code before
					my $go_obj = $go_adaptor->fetch_by_accession($dbe->display_id);
	  				$go_mapping{$dbe->display_id} = $go_obj->name;
	  				$go_terms{$go_obj->name}++;
	  			}
	  		}
	  		
	    	#print "\t".$dbe->dbname."\t".$dbe->display_id."\n";
	  	}
	  	
	  	#join by comma if multiple values
	  	my $ccds = keys %ccds?join(",",keys %ccds):"NO_CCDS";
		my $omim = keys %omims?join(",",keys %omims):"NO_OMIM";
		my $go_string = keys %go_terms?join("; ",keys %go_terms):"NO_GO";
	  	my $hugo;
	  	
	  	#Use the external name if no HGNC name
	  	if (!keys %hugo) {
	  		$hugo = $gene->external_name();
	  	} else {
	  		$hugo = join(",",keys %hugo);
	  	}
	  	
	  	my $hugo_lookup = uc($hugo);
	  	
	  	#Get refseq id from this API call
	  	my @ct_sf = @{ $canonical_transcript->get_all_supporting_features() };


		my %refseq_ids;

		foreach my $ct_sf (@ct_sf) {

	    	if ($ct_sf->analysis()->display_label() =~ /$organism cDNAs/i){
				$refseq_ids{$ct_sf->hseqname()}++
	    	}
		}
	  	
	  	my $refseq = keys %refseq_ids?join(",",keys %refseq_ids):"NO_REFSEQ";
	    my $uniprot_name = keys %uniprot_ids?join(",",keys %uniprot_ids):"NO_UNIPROT";
		my $vogelstein = exists $gene_name_mapper{$hugo_lookup}{vogel}?$gene_name_mapper{$hugo_lookup}{vogel}:"NO_VOGEL";
		my $cosmic = exists $gene_name_mapper{$hugo_lookup}{cosmic}?$gene_name_mapper{$hugo_lookup}{cosmic}:"NO_COSMIC";
	
		#my $tag = ">$gene_id";
		#print FASTA "$tag\n$sequence\n";
		my $ens_entry = $ens_link . $gene_id;
	    print ALLINFO join("\t",
	    					$ens_entry,
							$hugo,
							$cosmic,
							$vogelstein,
							$uniprot_name,
							$ccds,
							$refseq,
							$desc,
							$omim,
							$go_string
	    					) . "\n";

	}
}


open(EXONALL,">$exon_all_file") ||modules::Exception->throw("Can't open file to write $exon_all_file\n");

for my $chr (sort keys %gene_data) {
	my %exon_coord;
	my $exon_chr_file = $exon_file . '.' . $chr;
	my $splice_chr_file = $splice_file . '.' . $chr;
	open(EXON,">$exon_chr_file") || modules::Exception->throw("Can't open file to write $exon_chr_file\n");
	open(SPLICE,">$splice_chr_file") || modules::Exception->throw("Can't open file to write $splice_chr_file\n");
	
	for my $start (sort {$a<=>$b} keys %{$gene_data{$chr}}) {		
		
		for my $end (sort {$a<=>$b} keys %{$gene_data{$chr}{$start}}) {
			#Holds the overlap exon info
			
			
			my $exon_name = join(",",@{$gene_data{$chr}{$start}{$end}});
			(my $splice_name1 = $exon_name) =~ s/(exon\d+)/$1_splice1/g;
			(my $splice_name2 = $exon_name) =~ s/(exon\d+)/$1_splice2/g;
				
			
		
			    
		    print EXONALL join(" ",
		    				$chr,
		    				$start,
		    				$end,
		    				$exon_name
		    				). "\n";
		    
		    print EXON join(" ",
		    				$chr,
		    				$start,
		    				$end,
		    				$exon_name
		    				). "\n";
		    
		    my $distance = 10;
		    for my $splice_coord ($start-$splice_length..$start-1) {
			    print SPLICE join(" ",
			    				  $chr,
			    				  $splice_coord,
			    				  $splice_coord,
			    				  $splice_name1.'_'.$distance
			    				) . "\n";
		    	
		    	$distance--;
		    }
		    
		    for my $splice_coord ($end+1..$end+$splice_length) {
				print SPLICE join(" ",
			    				  $chr,
			    				  $splice_coord,
			    				  $splice_coord,
			    				  $splice_name2.'_'.$distance
			    				) . "\n";				
		    	$distance++;
		    }
		    				
		    				
		    
		    my $start_count = $start;
		    #Avoid duplicate coords; use hash and print at end
		   	while ($start_count <= $end) {
	    		$exon_coord{$start_count}++;
	    		$start_count++;
			}
			    
			
		}
	}
	
	for my $coord (sort {$a<=>$b} keys %exon_coord) {
		print COORD "$chr $coord\n";
	}
	
	close SPLICE;
	close EXON;
	
}


close EXONALL;
close COORD;
close ALLINFO;
#close FASTA;
