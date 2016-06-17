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
	   		"biomart_file=s",
	   		"ucsc_ccds_file=s",
	   		"ucsc_fasta_file=s",
			"uniprot_file=s",
			"uniprot_out=s",
			"ccds_out=s",
	   		"debug=s",
	   		"organism=s"
	   		);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{biomart_file} || !$OPT{ucsc_ccds_file} || !$OPT{uniprot_file} || !$OPT{ucsc_fasta_file});



=pod

=head1 SYNOPSIS

create_polyphen.pl -organism organism(default=mouse) -uniprot_out uniprot_output_mapping_file(default=uniprot.out) -ccds_out ccds_output_mapping_file(default=ccds.out) -ucsc_fasta_file fast_ccds_file -ucsc_ccds_file ncbi_ccds_file -uniprot_file uniprot_to_ens_mapping_file -biomart_file biomart_mapping_file(see README.conf for format) -debug [options]

Required flags: -biomart_file -ucsc_ccds_file -uniprot_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

create_polyphen.pl -> Script to generate a polyphen aa mapping chart

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

/home/matt/work/pipeline/conf/parsers/create_polyphen.pl -biomart_file ~/work/pipeline/conf/mouse/NCBIM37/gene/160611/NCBIM37.BiomartMapping.txt -ucsc_ccds_file ~/work/pipeline/conf/mouse/NCBIM37/ccds/040411/NCBIM37.ccdsGene.txt -uniprot_file ~/software/ngs/suites/polyphen/uniprot/mouse.seq -ucsc_fasta_file ~/work/pipeline/conf/mouse/NCBIM37/ccds/040411/NCBIM37.ccdsGene.fa -uniprot_out NCBIM37.uniprot_mapping.txt -ccds_out NCBIM37.uniprot_nomapping.txt


=cut

my $biomart_file = $OPT{biomart_file};
if ( ! -e $biomart_file ) {
	modules::Exception->throw("File $biomart_file doesn't exist");
}

my $ucsc_ccds_file = $OPT{ucsc_ccds_file};
if ( ! -e $ucsc_ccds_file ) {
	modules::Exception->throw("File $ucsc_ccds_file doesn't exist");
}

my $uniprot_file = $OPT{uniprot_file};
if ( ! -e $uniprot_file ) {
	modules::Exception->throw("File $uniprot_file doesn't exist");
}

my $organism = defined $OPT{organism}?$OPT{organism}:'mouse';
my $ucsc_fasta_file = $OPT{ucsc_fasta_file};
if ( ! -e $ucsc_fasta_file ) {
	modules::Exception->throw("File $ucsc_fasta_file doesn't exist");
}

my $ccds_out = defined $OPT{ccds_out}?$OPT{ccds_out}:"ccds.out";
my $uniprot_out = defined $OPT{uniprot_out}?$OPT{uniprot_out}:"uniprot.out";

open(NOMAP,">$ccds_out") || modules::Exception->throw("Can't open file to write $ccds_out\n");
open(MAP,">$uniprot_out") || modules::Exception->throw("Can't open file to write $uniprot_out\n");


my %gene_mapper = ();
my $registry = 'Bio::EnsEMBL::Registry';

print STDERR "Loading ensembl registry (this could be slow)\n";

my $reg_file = $ENV{SVNDIR} . '/conf/ensembl_registry.conf';

$registry->load_all($reg_file);
my $slice_adaptor = $registry->get_adaptor($organism, 'Core', 'Slice');

my $gene_adaptor = $registry->get_adaptor($organism, 'Core', 'Gene');


open(FILE,"$ucsc_fasta_file") || modules::Exception->throw("Can't open file $ucsc_fasta_file\n");

my $ccds_fasta;

#Keep track of ccds reported for bookkeeping at end
my %ccds_mapped = ();
my %all_ccds = ();

while (<FILE>) {
	chomp;
	if (/>(\S+)/) {
		$ccds_fasta = $1;

		#Don't use versions		
		if ($ccds_fasta =~ /\.[0-9]/) {
			$ccds_fasta =~ s/\.[0-9]//;
		}		
	} else {
		my $seq = $_;
		my $aa = modules::Utils->translate($seq);
		$gene_mapper{$ccds_fasta}{full_aa} = $aa;
		$aa =~ s/_//g;
		
		$gene_mapper{$ccds_fasta}{aa} = $aa;
		$gene_mapper{$ccds_fasta}{nt} = $seq;
	}
}

open(UCSC,"$ucsc_ccds_file") || die "Can't open file $ucsc_ccds_file\n";

#Get the exon coordinates from the ccds file
while (<UCSC>) {
	my @fields = split("\t");
	my @start_coords = split(",",$fields[9]);
	my @end_coords = split(",",$fields[10]);
	
	if (@start_coords != @end_coords) {
		modules::Exception->throw("ERROR: Incorrect number of exons");
	}
	
	my $strand = $fields[3];
	
	my @coord_str;
	(my $chr = $fields[2]) =~ s/chr//;
	for ( my $count = 0 ; $count < @start_coords ; $count++ ) {
		#Account for ucsc 0-based coords
	    my $buffered_start = $start_coords[$count] + 1;
	    my $buffered_end = $end_coords[$count];
	    push @coord_str, $buffered_start.'-'.$buffered_end;
	}
	
	my $coord_str = join(",",@coord_str);
	my $ccds_id = $fields[1];
	
	#Don't use version numbers
	if ($ccds_id =~ /\.[0-9]/) {
		$ccds_id =~ s/\.[0-9]//;
	}
	
	$all_ccds{$ccds_id} = 1;
	$gene_mapper{$ccds_id}{strand} = $strand;
	$gene_mapper{$ccds_id}{coord_str} = $coord_str;
	$gene_mapper{$ccds_id}{chr} = $chr;
}


#Map the two gene sets
my %uniprot_to_ccds = ();
my %ccds_to_uniprot = ();
my %uniprot_to_ensembl = ();

open(BIOMART,$biomart_file) || modules::Exception->throw("Can't open file $biomart_file");
#Get the biomart to ccds mapping info
while (<BIOMART>) {
	next unless /CCDS/;
	chomp;
	my ($ccds,$ucsc,$uniprot,$ens_gene,$genename,$genedesc,$go) =  split("\t");
	next if $uniprot eq 'NO_UNIPROT';
	my @uniprots = split(",",$uniprot);
	$ccds_to_uniprot{$ccds} = \@uniprots;
	for my $uniprot ( @uniprots ) {
		push @{$uniprot_to_ccds{$uniprot}}, $ccds;
		my @ensembls = split(",",$ens_gene);
		for my $ensembl ( @ensembls ) {
			push @{$uniprot_to_ensembl{$uniprot}}, $ensembl;    
		}    
	}
}


#Get the protein sequence from the database polyphen uses
open(UNIPROT,"$uniprot_file") || die "Can't open file $uniprot_file\n";
my $sp_accession;


#Keep track of sp info
my %sp_info = ();

my $uniprot_count = 0;
while (<UNIPROT>) {
	chomp;
	if (/>(\S+)/) {
		my @fields = split('\|',$1);
		$sp_accession = $fields[1];
	} else {
		if (exists $uniprot_to_ccds{$sp_accession}) {
			$sp_info{$sp_accession} = $_;
		}
	}
}

my $total = keys %gene_mapper;
print "ANALYSIS $total entries...\n\n";
my $count = 0;

CCDS:
for my $ccds (reverse sort { my ($anum) = $a =~ /(\d+)/; my ($bnum) = $b =~ /(\d+)/; $anum <=> $bnum } keys %gene_mapper) {
	
	if ($OPT{debug}) {
		my $found = 0;
		my @ccds = split(",",$OPT{debug});
		for my $ccds_debug ( @ccds ) {
			if ($ccds eq $ccds_debug) {
				$found = 1;
			} 	    
		}
		next unless $found;
	}
	
	if (exists $ccds_to_uniprot{$ccds}) {
		#Check the sp entries first
		for my $sp_accession (@{$ccds_to_uniprot{$ccds}}) {
			
			if (!defined $sp_info{$sp_accession}) {
				next;
			} 
			
			my $aa_uniprot_seq = $sp_info{$sp_accession};
			my $aa_uniprot_length = length($aa_uniprot_seq);
			
			#Check if the SP entry matches
			if ($gene_mapper{$ccds}{aa} eq $aa_uniprot_seq) {
    			&Report($gene_mapper{$ccds}{coord_str},$gene_mapper{$ccds}{strand},$gene_mapper{$ccds}{chr},$aa_uniprot_seq,$ccds,$sp_accession);
	    		print "PASS_CCDS: SP:$sp_accession:$aa_uniprot_length CCDS:$ccds $count\n";
	    		next CCDS;
	    	}
		}
		
		
		
		#Check the ens matches next
		for my $sp_accession (@{$ccds_to_uniprot{$ccds}}) {
			
			if (!defined $sp_info{$sp_accession}) {
				next;
			} 
			
			my $aa_uniprot_seq = $sp_info{$sp_accession};
			my $aa_uniprot_length = length($aa_uniprot_seq);
	    	#See if we can find an ensembl match if SP fails
	    	my ($match_ensembl,$ens_coord,$ens_strand,$ens_chr,$ens_name) = &Check_Ensembl($sp_accession,$aa_uniprot_seq);
	    		    	
			if ($match_ensembl) {
    			&Report($ens_coord,$ens_strand,$ens_chr,$aa_uniprot_seq,$ens_name,$sp_accession);	
			   	print "PASS_ENS: SP:$sp_accession:$aa_uniprot_length ENS:$ens_name CCDS:$ccds $count\n";
			   	next CCDS;
			}
	    	
		}
	    	
		#Here we found no match so it unmapped
		my $sp = join(",",@{$ccds_to_uniprot{$ccds}});
	    &Report($gene_mapper{$ccds}{coord_str},$gene_mapper{$ccds}{strand},$gene_mapper{$ccds}{chr},$gene_mapper{$ccds}{aa},$ccds. '_NO_MAP',$sp);
		print "FAIL_ALL: SP:$sp CCDS:$ccds $count\n";
	} else {
		#Here there are no SP entries for the ccds
		print "FAIL_ALL: SP:NO_SP CCDS:$ccds $count\n";
    	&Report($gene_mapper{$ccds}{coord_str},$gene_mapper{$ccds}{strand},$gene_mapper{$ccds}{chr},$gene_mapper{$ccds}{full_aa},$ccds. '_NO_MAP','NO_SP');
	}
	
}

#Report the matching entries; only report perfect matches
sub Report {
	my ($coord_str,$strand,$chr,$aa_uniprot,$mapped_entry,$sp) = @_;
	my @aa = split("",$aa_uniprot);	
	my $aa_length = @aa;
	my $aa_number = 1;
	my $aa_index_count = 0;
	my @aa_uniprot_coord = ();

	my %numbers_in_set = &Get_Numbers($coord_str);
	my @sorted_numbers = ();
	
	#Product strand affects 
	if ($strand eq '+') {
		@sorted_numbers  = sort {$a<=>$b} keys %numbers_in_set;
	} else {
		@sorted_numbers = reverse sort {$a<=>$b} keys %numbers_in_set;
	}

	#my $numbers = @sorted_numbers;
	#print "$coord_str $numbers\n";
	#print Dumper \@sorted_numbers;

	my $local_count = 0;
	for my $number ( @sorted_numbers ) {
	    if ($local_count % 3 == 0 && $local_count != 0) {
	    	#Always report things on the positive strand b/c snps are on the positive strand
	    	my @sorted_coord = sort {$a<=>$b} @aa_uniprot_coord;
			my $bp_seq;
			#First get the sequence
	    	for my $coord (@sorted_coord) {
	    		$bp_seq .= $slice_adaptor->fetch_by_region('chromosome',$chr, $coord, $coord)->seq();
	    	}
	    	my $position = 0;
			for my $coord (@sorted_coord) {
				my @bases = split("",$bp_seq);
				my $refbase = $bases[$position];
				if ($mapped_entry =~ /NO_MAP/) {
					#Here we weren't able to map the SP entry so we use the ccds aa sequence
					
					#Don't print huge proteins like TTN
					my $aa_print = $aa_uniprot;
					if ($aa_length > 2000) {
						$aa_print = "SEQ_TOO_LONG(".$aa_length." aa)";
					}
						
			
					print NOMAP "$chr $coord $coord $sp:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry:$aa_print\n";
					#print "$chr $coord $coord $sp_accession:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry:$aa_uniprot\n";
				} else {
		    		print MAP "$chr $coord $coord $sp:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry\n";
		    		#print "$chr $coord $coord $sp_accession:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry\n";
				}
		    	$position++;
			}
	    	@aa_uniprot_coord = ();
	    	$aa_number++;
	    	$aa_index_count++;
	    }
	    push @aa_uniprot_coord, $number;
	    $local_count++;
	}	
	$count++;
}

#check whether there is a corresponding ensembl entry when ccds mapping fails
sub Check_Ensembl {
	my ($sp_accession,$aa_uniprot) = @_;
	my $match_ensembl = 0;
    my @ensembl_coords = ();
    my $ensembl_coord_str;
    
    #here we try to use ensembl to resolve the difference in ccds/uniprot transcripts
	for my $ensembl_id ( @{$uniprot_to_ensembl{$sp_accession}} ) {
		#This entry crashes the script
		my $gene = $gene_adaptor->fetch_by_stable_id($ensembl_id);
		if (!$gene) {
			print "ERROR: Skip checking ensembl for $ensembl_id\n";
			next;
		}
	    my $slice = $slice_adaptor->fetch_by_gene_stable_id($ensembl_id);
		foreach my $transcript (@{$slice->get_all_Transcripts}){
			my $strand = $transcript->strand();
			my $name = $transcript->stable_id();
			my $chr = $transcript->seq_region_name();
			my $bp_ensembl_length = 0;
			my $exon_count = 0;
			
			my $aa_ensembl;
			if ($transcript->translate()) {
				$aa_ensembl = $transcript->translate()->seq();
			} else {
				next;
			}

			#If we've found a match get the coordinates
			if ($aa_ensembl eq $aa_uniprot) {
				foreach my $exon (@{$transcript->get_all_Exons}) {
					#my $estring = feature2string($exon);
					my $coding_start = $exon->coding_region_start($transcript);
					my $coding_end = $exon->coding_region_end($transcript);
					my $genome_start = $exon->seq_region_start();
					my $genome_end = $exon->seq_region_end();
					
					#Skip if not coding
	            	if (!defined $coding_start) {
	            		next;
	            	} 
				
	            	my $length = $coding_end - $coding_start + 1;
					my $genome_length = $genome_end - $genome_start + 1;
					if ($length != $genome_length) {			
						#Adjust the coordinates for UTRs
						if ($strand == 1) {
							if ($exon_count == 0) {
								$genome_start = $genome_end - $length + 1;
							} else {
								$genome_end = $genome_start + $length - 1;
							}
						} else {
							if ($exon_count == 0) {
								$genome_end = $genome_start + $length - 1;
							} else {
								$genome_start = $genome_end - $length + 1;
							}
						}
						$exon_count++;
					}
					push @ensembl_coords, $genome_start. '-' . $genome_end;		
				}
				$ensembl_coord_str = join(",",@ensembl_coords);
				my $return_strand =  $strand == 1?'+':'-';
				return (1,$ensembl_coord_str,$return_strand,$chr,$name);
			}
		}
	}
	return (0,'','','','');
}
	

#Get the numbers in a range of coordinates
sub Get_Numbers {
	my ($coord_str) = @_;
	my @num_groups = split(",",$coord_str);
	my %numbers;
	for my $num_group ( @num_groups ) {
	    my ($start,$end) = $num_group =~ /(\d+)-(\d+)/;
	    while ($start <= $end) {
	    	$numbers{$start}++;
	    	$start++;
	    }
	}
	return %numbers;
	
	
}


