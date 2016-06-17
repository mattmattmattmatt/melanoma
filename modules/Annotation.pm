package modules::Annotation;

use strict;
use Data::Dumper;
use modules::ConfigXML;

sub new {
	my ($class, @args) = @_;
	
	my $self = bless {}, $class;

    my %args = @args;

    my @required_args = (-annotation_file, 
    					-exon_coord_file, 
    					-splice_size, 
    					-report_xml, 
    					-sample_type, 
    					-mutant_type);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
   	$self->{binsize} =  1000000;   	
	$self->{splice_size} = $args{-splice_size};

	#Parse the xml to get the annotation columns
    my $config = modules::ConfigXML->new($args{-report_xml});
    my $sample_type = $args{-sample_type};
    my $mutant_type = $args{-mutant_type};
    if (!$config->exists($mutant_type,'sample_types',$sample_type)) {
    	modules::Exception->throw("ERROR: Cannot get annotation columns for $mutant_type $sample_type");
    }
    
    my @annotation_headers = ();
	if ($config->exists($mutant_type,'sample_types',$sample_type,'annotations')) {
		@annotation_headers = split(",",$config->read($mutant_type,'sample_types',$sample_type,'annotations'));
	} else {
		@annotation_headers = split(",",$config->read($mutant_type,'common','annotations'));
	}
    
	$self->{annotation_headers} = \@annotation_headers;

	$self->_create_map($args{-annotation_file}, $args{-exon_coord_file}, \@annotation_headers);
   	
    return $self;
}

sub _create_map {;
    my $self = shift;
    my ($anno_file,$exon_coord_file) = @_;
    my @annotation_headers = @{$self->{annotation_headers}};
    my %exon_map = ();
    my %coord_map = ();
    open(ANNO,"$anno_file") || die "Can't open file $anno_file\n";
    while (<ANNO>) {
    	chomp;
    	#TODO: Fix to match latest gene columns
    	my @fields = split("\t");
    	my $exon_link = $fields[0];
    	my ($exon) = $fields[0] =~ /(ENS.*)$/;
    	
    	
    	#print "E $exon\n";
    	
    	#The +1 account for the ensembl entry
		if (@fields != @annotation_headers+1) {
			modules::Exception->throw("ERROR: Column numbers don't match for annotation_file and annotation_headers");
		}    	
    	
    	$exon_map{$exon}{ensembl} = $exon_link;
    	
    	my $count = 1;
    	for my $annotation_header ( @annotation_headers ) {
    	    $exon_map{$exon}{$annotation_header} = $fields[$count];
    	    $count++;
    	}
    }	    
		    
    open(EXON,"$exon_coord_file") || die "Can't open file $exon_coord_file\n";
    
    my %exon_ranges = ();
    
	#Join the individual exon ranges
    while(<EXON>) {
    	my ($chr,$start,$end,$exon_str) = split;
    	
    	my @exons = split(",",$exon_str);
    	
    	for my $exon_full ( @exons ) {
	    	(my $exon_short = $exon_full) =~ s/_exon\d+//;
	    	
	    	
	    	if (exists $exon_ranges{$exon_short}) {
	    		if ($start < $exon_ranges{$exon_short}{min}) {
		    		$exon_ranges{$exon_short}{min} = $start;
	    		}
	    		if ($end > $exon_ranges{$exon_short}{max}) {
		    		$exon_ranges{$exon_short}{max} = $end;
	    		}
	    		
	    	} else {
	    		$exon_ranges{$exon_short}{min} = $start;
	    		$exon_ranges{$exon_short}{max} = $end;
	    	}
	    	$exon_ranges{$exon_short}{chr} = $chr;
    	    
    	}
    	
    }

    my $splice_size = $self->{splice_size};
      
    for my $exon (keys %exon_ranges) {
    	#next unless $exon eq 'ENSG00000100170';
		my $start = $exon_ranges{$exon}{min};
		my $end = $exon_ranges{$exon}{max};
		my $chr = $exon_ranges{$exon}{chr};
    	my $coord_str = $chr . ':'. $start . '-' . $end;
    	$exon_map{$exon}{coord} = $coord_str;

    	#Add the splice site buffer so we can get to the original gene
    	my $buffered_start = $start - $splice_size;
    	$buffered_start = 0 if $start < 0;
    	my $buffered_end = $end + $splice_size;
    	
    	my $start_bin = int($buffered_start / $self->{binsize});
    	my $end_bin = int($buffered_end / $self->{binsize});
		my $bin_count = $start_bin;
    	while ( $bin_count <= $end_bin ) {
	    	push @{$coord_map{$chr}{$bin_count}{"$buffered_start-$buffered_end"}}, $exon;
    		$bin_count++;	    
    	}
    }
     
    $self->{coord_map} = \%coord_map;
    $self->{exon_map} = \%exon_map;
}

sub get_gene_from_coord {
	my $self = shift;
	my %args = @_;

    my @required_args = ('-chr','-coord');

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		} 
    }
    
    my $chr = $args{-chr};
    my $coord = $args{-coord};
	my $coord_bin = int($coord / $self->{binsize});   
	my @annotation_headers = @{$self->{annotation_headers}};
	my $gene_name_col = $annotation_headers[0];
	
    my %gene_info = ();
	for my $range ( sort { my ($astart) = $a =~ /^(\d+)/; my ($bstart) = $b =~ /^(\d+)/; $astart <=> $bstart } keys %{$self->{coord_map}{$chr}{$coord_bin}} ) {
		my ($exon_start,$exon_end) = $range =~ /(\d+)-(\d+)/;
		if ($exon_start > $coord) {
			#Here we're out of range so move on
			last;
		}
		
		#Here we match the coordinate
		if ($exon_start <= $coord && $exon_end >= $coord) {
			my $exon_array = $self->{coord_map}{$chr}{$coord_bin}{$range};
			for my $exon ( @{$exon_array} ) {
				if (!defined $self->{exon_map}{$exon}{$gene_name_col}) {
					modules::Exception->throw("ERROR: exon $exon does not map correctly");
				} 
				
				#Not an annotation_header so have to add manually				
				$gene_info{ensembl}{$self->{exon_map}{$exon}{ensembl}}++;
								
				for my $annotation_header ( @annotation_headers ) {
				    $gene_info{$annotation_header}{$self->{exon_map}{$exon}{$annotation_header}}++;
				}
				
				#print Dumper \%gene_info;
				
#				$gene_info{gene_name}{$self->{exon_map}{$exon}{genename}}++;
#				$gene_info{gene_desc}{$self->{exon_map}{$exon}{genedesc}}++;
#				$gene_info{uniprot_name}{$self->{exon_map}{$exon}{uniprot}}++;
#				$gene_info{refseq}{$self->{exon_map}{$exon}{refseq}}++;
#				$gene_info{ucsc}{$self->{exon_map}{$exon}{ucsc}}++;
#				$gene_info{go}{$self->{exon_map}{$exon}{go}}++;
#				$gene_info{ensembl}{$exon}++;
#				$gene_info{omim}{$self->{exon_map}{$exon}{omim}}++;
				
			}
		}
	}
	return \%gene_info;
	
}


return 1;
