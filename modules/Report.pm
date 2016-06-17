package modules::Report;

use strict;
use modules::FilteredSNV;
use modules::Adaptors::SNV;
use modules::Adaptors::SNV_Row;
use modules::Adaptors::Variant_Row;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::BulkInsert;
use modules::FilteredVariant;
use modules::PileupLine;
use modules::Exception;
use modules::Utils;
use modules::ConfigXML;
use Data::Dumper;
use File::Basename;
use Bio::EnsEMBL::Registry;
use Env qw($ENSEMBL_REGISTRY);

my $EXON_NS_FILTER = 'filter_exon_ns';
my $EXON_FILTER = 'filter_exon';
my $DBSNP_FILTER = 'filter_dbsnp_snv';
my $DBSNP_FILTER_INDEL = 'filter_dbsnp_indel';
my $COSMIC_FILTER = 'filter_cosmic';
my $SPLICE_FILTER = 'filter_splicesite';

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

    my @required_args = (
			             -run, 
						 -gene_mapper,
						 -sample_type,
						 -gene_col_name
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
	
	#Run id
    $self->run($args{'-run'});
	#Gene name mapping
	$self->gene_mapper($args{'-gene_mapper'});
	#Sample group
	$self->sample_type($args{-sample_type});
	#Gene column name
	$self->gene_col_name($args{-gene_col_name});
	
	
	#Optional extra lists we've filtered out    
	my @filter_files = defined $args{-filter_files}?@{$args{-filter_files}}:();
	$self->filter_files(\@filter_files) if @filter_files;
	
	#Used to differentiate indel from snv reports when different behavior is required
	if (defined $args{-indel}) {
		$self->{variant_type} = 'indel';
	} else {
		$self->{variant_type} = 'snv';
		
		#These are snv specific
		if (defined $args{-polyphen_file}) {
			$self->_parse_polyphen($args{-polyphen_file});
		}	

		#This is for the polyphen entries that don't have scores
		if (defined $args{-polyphen_info_file}) {
			$self->_parse_polyphen_info($args{-polyphen_info_file})
		}
	}
	
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_all($ENSEMBL_REGISTRY);
	my $slice_adaptor = $registry->get_adaptor('human', 'core', 'Slice');
	$self->{adaptor} = $slice_adaptor;	
	
    return $self;
}

#Parse the polyphen file to add these scores to the columns -> these are the polyphen entries that have been run
sub _parse_polyphen {
	my $self = shift;
	my $polyphen_file = shift;
	my %polyphen = ();
	open(FILE,"$polyphen_file") || modules::Exception->throw("Can't open file $polyphen_file\n");
	while (<FILE>) {
		my ($chr,$coord,undef,$score,$prediction,$rest) = split("\t");
		my @rest_fields = split(':',$rest);
		my $nt_change = $rest_fields[2];
		$polyphen{$chr}{$coord}{$nt_change} = $prediction  . ',' . $score;
	}
	close FILE;
	$self->{polyphen} = \%polyphen; 
}

#Parse the polyphen info file -> this contains info to allow user to run polyphen themselves
sub _parse_polyphen_info {
	my $self = shift;
	my $polyphen_info_file = shift;
	my %polyphen_info = ();
	open(FILE,"$polyphen_info_file") || modules::Exception->throw("Can't open file $polyphen_info_file\n");
	while (<FILE>) {
		chomp;
		my ($chr,$coord,undef,$poly_info) = split("\t");
		$polyphen_info{$chr}{$coord} = $poly_info;
	}
	close FILE;
	$self->{polyphen_info} = \%polyphen_info; 
}

#Subroutine loads up all the filters used for input/reporting/pass,etc; these conditions are based on sample_type
sub load {
	my $self = shift;
	my %args = @_;

    my @required_args = (
			             -report_xml
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my %ignore_filters = ();
    if (defined $args{-ignore_filters}) {
    	%ignore_filters = %{$args{-ignore_filters}};
    }
    my $sample_type = $self->{sample_type};

    #Parse the xml
    my $config = modules::ConfigXML->new($args{-report_xml});

	#print Dumper $config;

	my $variant_type = $self->{variant_type};

	my $conf_header;
	if ($variant_type eq 'snv') {
		$conf_header = 'snv';
	} else {
		$conf_header = 'indel';
	}


	if (!$config->exists($conf_header,'sample_types',$sample_type)) {
		modules::Exception->throw("ERROR: no sample_type config entry for $sample_type");
	}
    

    #Load up the data struct
    
    #Input filter
    
    my @input_filters = ();
    if ($config->exists($conf_header,'sample_types',$sample_type,'input_filters')) {
    	@input_filters = split(",",$config->read($conf_header,'sample_types',$sample_type,'input_filters'));
    } else {
    	@input_filters = split(",",$config->read($conf_header,'common','input_filters'));
    }
    
    my @input_filters_db;
    for my $input_filter ( @input_filters ) {
    	if (!exists $ignore_filters{$input_filter}) {
	    	my ($filter_db) = modules::Adaptors::Filter->search(name => $input_filter);
	    	if (!defined $filter_db) {
	    		modules::Exception->throw("ERROR: filter $input_filter doesn't exist in db");
	    	}
	    	push @input_filters_db, $filter_db;
    	}
	}
    
    $self->{primary_filters} = \@input_filters_db;
    
    #Extra filters to report
    my @extra_filters = ();
    if ($config->exists($conf_header,'sample_types',$sample_type,'extra_columns')) {
    	@extra_filters = split(",",$config->read($conf_header,'sample_types',$sample_type,'extra_columns'));
    } else {
    	@extra_filters = split(",",$config->read($conf_header,'common','extra_columns'));
    }
    
    my @extra_filters_db;
    for my $extra_filter ( @extra_filters ) {
       	if (!exists $ignore_filters{$extra_filter}) {
			my ($filter_db) = modules::Adaptors::Filter->search(name => $extra_filter);
			if (!defined $filter_db) {
	    		modules::Exception->throw("ERROR: filter $extra_filter doesn't exist in db");
	    	}
		    push @extra_filters_db, $filter_db;
       	}
	}	
    $self->{filter_list_order} = \@extra_filters_db; 
    
    #Filter with snv_filter info
    my @filter_info = ();
    if ($config->exists($conf_header,'sample_types',$sample_type,'filter_info')) {
    	@filter_info = split(",",$config->read($conf_header,'sample_types',$sample_type,'filter_info'));
    } else {
    	@filter_info = split(",",$config->read($conf_header,'common','filter_info'));
    }
    
    my @filter_info_db;
    for my $filter_info ( @filter_info ) {
       	if (!exists $ignore_filters{$filter_info}) {
			my ($filter_db) = modules::Adaptors::Filter->search(name => $filter_info);
			if (!defined $filter_db) {
	    		modules::Exception->throw("ERROR: filter $filter_info doesn't exist in db");
	    	}
		    push @filter_info_db, $filter_db;
       	}
	}	
    $self->{filter_info} = \@filter_info_db; 
	
	
	
	#Report filter data structure
	
	
	my %report_set = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'summary_report')) {
		%report_set = map {$_ => 1} split(",",$config->read($conf_header,'sample_types',$sample_type,'summary_report'));
	} else {
		%report_set = map {$_ => 1} split(",",$config->read($conf_header,'common','summary_report'));	
	}
	for my $report_filter (keys %report_set) {
		if (exists $ignore_filters{$report_filter}) {
			delete $report_set{$report_filter};
		}
	}
	$self->{report_set} = \%report_set;    
	
	#Filter pass conditions 
	my $rules;
	if ($config->exists()) {
		$rules = $config->read($conf_header,'sample_types',$sample_type,'pass_conditions','rule');
	} else {
		$rules = $config->read($conf_header,'common','pass_conditions','rule');		
	}
	my @rules = ();
	#Get the rules	
	if (ref($rules) eq 'ARRAY'){ # Cope with steps being a single step or an array of steps
    	@rules = @$rules;
	} else {
    	@rules = ($rules);	
	}
	
	my %rule_set = ();
	for my $rule ( @rules ) {
	    my ($name,$condition) = split(':',$rule);
	    if (!exists $ignore_filters{$name}) {
	    	$rule_set{$name} = $condition;
	    }
	}
	
	$self->{filter_set} = \%rule_set;

	my @annotation_headers = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'annotations')) {
		@annotation_headers = split(",",$config->read($conf_header,'sample_types',$sample_type,'annotations'));
	} else {
		@annotation_headers = split(",",$config->read($conf_header,'common','annotations'));
	}

	my @run_headers = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'run_headers')) {
		@run_headers = split(",",$config->read($conf_header,'sample_types',$sample_type,'run_headers'));
	} else {
		@run_headers = split(",",$config->read($conf_header,'common','run_headers'));
	}

	my @extra_headers = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'extra_headers')) {
		@extra_headers = split(",",$config->read($conf_header,'sample_types',$sample_type,'extra_headers'));
	} else {
		@extra_headers = split(",",$config->read($conf_header,'common','extra_headers'));
	}
	
	
	$self->{extra_headers} = \@extra_headers;
	$self->{run_headers} = \@run_headers;
	$self->{annotation_headers} = \@annotation_headers;

	#Header lookup for columns
	$self->headers();
	
    
    
    
}

#Add the additional filter files
sub filter_files {
	my $self = shift;
	my ($filter_files) = @_;
	my %filter_files = ();
	#Handle any additional gene lists we want to filter
	my $count = 1;
	for my $filter_file ( @{$filter_files} ) {
	    if ( !-e $filter_file ) {
	    	modules::Exception->throw("File $filter_file doesn't exist");	
	    }
	    #Default filter name is the base of the filename
	    my $filter_name = basename($filter_file);
	
	    my %filter_list = ();
	    open(FILE,"$filter_file") || die "Can't open file $filter_file\n";
	    while (<FILE>) {
	    	chomp;
	    	my ($entry) = $_ =~ /^(\S+)/;
	    	next unless $entry;
	    	$filter_list{uc($entry)} = 1;
	    }
	    $filter_files{$filter_name} = \%filter_list;
		$count++;
	}
	$self->{filter_files} = \%filter_files;
	return $self->{filter_files};
}


#Add the genemapper
sub gene_mapper {
	 my ($self, $gene_mapper) = @_;

    if (defined $gene_mapper) {
		$self->{'gene_mapper'} = $gene_mapper;
    } elsif (! defined $self->{'gene_mapper'}) {
		modules::Exception->throw("genemapper not set");
    }

    return $self->{'gene_mapper'};
}

#Map all the headers
sub headers {
	my $self = shift;
	my %header_col_lookup;
	
	my @filter_names;
						
    foreach my $filter (@{$self->{filter_list_order}}) {
		push @filter_names, $filter->name;
    }

	my @annotation_headers = @{$self->{annotation_headers}};
	my @run_headers = @{$self->{run_headers}};
	my @extra_headers = @{$self->{extra_headers}};


	my @common_headers = (
						  @run_headers,			      
					      @filter_names,
					      @extra_headers,
					      @annotation_headers
						  );
	

	for (my $i = 0; $i < scalar @common_headers; $i++){
    	$header_col_lookup{$common_headers[$i]} = $i;
    	
	}
		
	$self->{headers} = \@common_headers;
	$self->{header_lookup} = \%header_col_lookup;
}

#Get the runid
sub run {
    my ($self, $run) = @_;

    if (defined $run) {
		$self->{'run'} = $run;
    } elsif (! defined $self->{'run'}) {
		modules::Exception->throw("run not set");
    }

    return $self->{'run'};
}


#Get the sample_typeid
sub sample_type {
    my ($self, $sample_type) = @_;

    if (defined $sample_type) {
		$self->{'sample_type'} = $sample_type;
    } elsif (! defined $self->{'sample_type'}) {
		modules::Exception->throw("sample_type not set");
    }

    return $self->{'sample_type'};
}

#Get the gene_col_name
sub gene_col_name {
	my ($self, $gene_col_name) = @_;

    if (defined $gene_col_name) {
		$self->{'gene_col_name'} = $gene_col_name;
    } elsif (! defined $self->{'gene_col_name'}) {
		modules::Exception->throw("gene_col_name not set");
    }

    return $self->{'gene_col_name'};
}

#Generate the file containing the polyphen info -> NOW DONE FROM ENSEMBL
#sub generate_polyphen_input {
#	my $self = shift;
#	my %args = @_;
#		
#	my %sorted_pass_rows = %{$self->{sorted_pass_rows}};
#	#Now print the appropriate fields
#	my $divider_count = @{$self->{headers}} - 1;
#
#	my $count = 0;
#
#	#Need this for polyphen	
#	my %header_col_lookup = %{$self->{header_lookup}};
#	my $uniprot_col_num = $header_col_lookup{'uniprot_name'};
#	my $ref_base_col_num = $header_col_lookup{'ref_base'};
#	my $var_base_col_num = $header_col_lookup{'var_base'};
#	my $aa_change_col_num = $header_col_lookup{'aa_change'};
#	my $dbsnp_col_num = $header_col_lookup{'filter_dbsnp_snv'};
#	my %polyphen_data = ();
#	
#
#	for my $type (keys %sorted_pass_rows) {
#		
#		for my $splice ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}} ) {
#		    for my $depth ( sort {$b<=>$a} keys %{$sorted_pass_rows{$type}{$splice}} ) {
#		    	for my $chr ( sort keys %{$sorted_pass_rows{$type}{$splice}{$depth}} ) {
#		    		for my $coord ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}} ) {
#		    			#We don't calculate polyphen for this type of entry -> NOW HANDLED WITH NEW RARE_ALLELE FILTER
#		    			#if ($sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$dbsnp_col_num] eq 'RARE_ALLELE') {
#		    			#	next;
#		    			#}
#		    			
#		    			$polyphen_data{$chr}{$coord}{ref} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$ref_base_col_num];
#		    			$polyphen_data{$chr}{$coord}{var} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$var_base_col_num];
#		    			$polyphen_data{$chr}{$coord}{uniprot} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$uniprot_col_num];
#						$polyphen_data{$chr}{$coord}{aachange} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$aa_change_col_num];
#						$count++;
#		    		}
#				}
#			}
#		}
#	}
#    return (\%polyphen_data,$count);
#}

#Print the final TSV that goes to the spreadsheet and the match and no match files
sub print_to_files {
    my ($self) = shift;
    my %args = @_;
    
    
    my @required_args = (
    					'-tsv_file',
    					'-pass_file'
			 			);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $output_filename = $args{-tsv_file};
    open(my $OUT, ">$output_filename")
	or modules::Exception->throw("Unable to open output file [$output_filename]");

	my $pass_file = $args{-pass_file};
	open(my $PASS, ">$pass_file")
    or modules::Exception->throw("Unable to open output file [$pass_file]");

	(my $fail_file = $pass_file) =~ s/match/nomatch/;

	open(my $FAIL, ">$fail_file") 
	or modules::Exception->throw("Unable to open output file [$fail_file]");


	# print the headers
    print $OUT join("\t",@{$self->{headers}}) . "\n";
	
	my %sorted_pass_rows = %{$self->{sorted_pass_rows}};
	my %rare_allele_rows = %{$self->{sorted_allele_rows}};
	
	
	#Now print the appropriate fields
	my $divider_count = @{$self->{headers}} - 1;

	my %header_col_lookup = %{$self->{header_lookup}};
	
	my $variant_type = $self->{variant_type};

	if ($variant_type eq 'snv') {

		for my $type (keys %sorted_pass_rows) {
			my $divider = "$type PASS\t" . "---\t" x $divider_count;
			chop $divider;
			print $OUT "$divider\n";
			
			for my $splice ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}} ) {
			    for my $depth ( sort {$b<=>$a} keys %{$sorted_pass_rows{$type}{$splice}} ) {
			    	for my $chr ( sort keys %{$sorted_pass_rows{$type}{$splice}{$depth}} ) {
			    		for my $coord ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}} ) {
			    			print $OUT join("\t",@{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}}) . "\n";
			    			#Make a useful key for overlap file  
			    			
			    			my $rest = join("^^^","PASS",
			    								$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele}] .'->'. $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele}],
			    								'ref_count:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele_count}],
		 	 									'var_count:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele_count}],
		 	 									'snv_score:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{snv_score}],
		 	 									'clr_score:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{clr_score}],
		 	 									'final_status:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{final_status}],
		 	 									'normal_alleles:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{normal_alleles}]	
			    								);
			    			$rest =~ s/ /_/g;
			    			
			    			print $PASS join("\t",
			    							$chr,
			    							$coord,
			    							$coord,
			    							$rest
			    							)."\n";
						}
					}
				}
			}
		}
		
		for my $type (keys %rare_allele_rows) {
			my $divider = "$type RARE_ALLELE(<2%)\t" . "---\t" x $divider_count;
			chop $divider;
			print $OUT "$divider\n" if exists $rare_allele_rows{$type};
			
			for my $splice ( sort {$a<=>$b} keys %{$rare_allele_rows{$type}} ) {
			    for my $depth ( sort {$b<=>$a} keys %{$rare_allele_rows{$type}{$splice}} ) {
			    	for my $chr ( sort keys %{$rare_allele_rows{$type}{$splice}{$depth}} ) {
			    		for my $coord ( sort {$a<=>$b} keys %{$rare_allele_rows{$type}{$splice}{$depth}{$chr}} ) {
			    			print $OUT join("\t",@{$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}}) . "\n";
			    			#Make a useful key for overlap file  
			    			my $rest = join("^^^","RARE_ALLELE",
			    								$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele}] .'->'. $rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele}],
			    								'ref_count:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele_count}],
		 	 									'var_count:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele_count}],
		 	 									'snv_score:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{snv_score}],
		 	 									'clr_score:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{clr_score}],
		 	 									'final_status:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{final_status}],
		 	 									'normal_alleles:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{normal_alleles}]	
			    								);
			    			$rest =~ s/ /_/g;
			    			
			    			print $PASS join("\t",
			    							$chr,
			    							$coord,
			    							$coord,
			    							$rest
			    							)."\n";
						}
					}
				}
			}
		}
		
		
		# Print a dividing line between 'above the line' and 'below the line'
		
		
		my $final_divider = "FAIL LINE\t" . "---\t" x $divider_count;
		chop($final_divider);
		print $OUT "$final_divider\n";
	    
	
		my @fail_rows = @{$self->{fail_rows}};
		for my $field_array ( @fail_rows) {
		    my $line = join("\t",@{$field_array});
		    print $OUT $line,"\n";
		 	my $rest = join("^^^",  'FAIL',
		 							$field_array->[$header_col_lookup{ref_allele}]."->".$field_array->[$header_col_lookup{var_allele}],
		 	 						'ref_count:'.$field_array->[$header_col_lookup{ref_allele_count}],
		 	 						'var_count:'.$field_array->[$header_col_lookup{var_allele_count}],
		 	 						'snv_score:'.$field_array->[$header_col_lookup{snv_score}],
		 	 						'clr_score:'.$field_array->[$header_col_lookup{clr_score}],
		 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}],
		 	 						'normal_alleles:'.$field_array->[$header_col_lookup{normal_alleles}]
		 	 						);
		 	$rest =~ s/ /_/g;
		 	print $FAIL join("\t",$field_array->[$header_col_lookup{chr}], $field_array->[$header_col_lookup{coord}], $field_array->[$header_col_lookup{coord}], $rest) ."\n";
		}
	
	} else {
		for my $type (keys %sorted_pass_rows) {
			my $divider = "$type PASS\t" . "---\t" x $divider_count;
			chop $divider;
			print $OUT "$divider\n";
			
			for my $splice ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}} ) {
			    for my $depth ( sort {$b<=>$a} keys %{$sorted_pass_rows{$type}{$splice}} ) {
			    	for my $chr ( sort keys %{$sorted_pass_rows{$type}{$splice}{$depth}} ) {
			    		for my $start_coord ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}} ) {
			    			for my $end_coord ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}} ) {
			    				for my $var_base ( sort keys %{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}} ) {
			    			
					    			print $OUT join("\t",@{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}}) . "\n";
					    			#Make a useful key for overlap file  
					    			my $rest = join("^^^","PASS",
					    									$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_type}],
															$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele}]."->".$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele}],
		 	 												'ref_count:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele_count}],
		 	 												'var_count:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele_count}],
		 	 												'snv_score:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_score}],
		 	 												'clr_score:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{clr_score}],
		 	 												'final_status:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{final_status}],
		 	 												'normal_alleles:'.$sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{normal_alleles}]
					    								);
					    			$rest =~ s/ /_/g;
					    			print $PASS join("\t",
					    							$chr,
					    							$start_coord,
					    							$end_coord,
					    							$rest
					    							)."\n";
			    				}
			    			}
						}
					}
				}
			}
		}
		
		for my $type (keys %rare_allele_rows) {
			my $divider = "$type RARE_ALLELE(<2%)\t" . "---\t" x $divider_count;
			chop $divider;
			print $OUT "$divider\n" if exists $rare_allele_rows{$type};
			
			for my $splice ( sort {$a<=>$b} keys %{$rare_allele_rows{$type}} ) {
			    for my $depth ( sort {$b<=>$a} keys %{$rare_allele_rows{$type}{$splice}} ) {
			    	for my $chr ( sort keys %{$rare_allele_rows{$type}{$splice}{$depth}} ) {
			    		for my $start_coord ( sort {$a<=>$b} keys %{$rare_allele_rows{$type}{$splice}{$depth}{$chr}} ) {
			    			for my $end_coord ( sort {$a<=>$b} keys %{$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}} ) {
			    				for my $var_base ( sort keys %{$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}} ) {
			    			
					    			print $OUT join("\t",@{$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}}) . "\n";
					    			#Make a useful key for overlap file  
					    			my $rest = join("^^^","RARE_ALLELE",
					    									$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_type}],
															$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele}]."->".$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele}],
		 	 												'ref_count:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele_count}],
		 	 												'var_count:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele_count}],
		 	 												'snv_score:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_score}],
		 	 												'clr_score:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{clr_score}],
		 	 												'final_status:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{final_status}],
		 	 												'normal_alleles:'.$rare_allele_rows{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{normal_alleles}]
					    							);
					    			$rest =~ s/ /_/g;
					    			print $PASS join("\t",
					    							$chr,
					    							$start_coord,
					    							$end_coord,
					    							$rest
					    							)."\n";
			    				}
			    			}
						}
					}
				}
			}
			
		}
		
		
		
		my $final_divider = "FAIL LINE\t" . "---\t" x $divider_count;
		chop($final_divider);
		print $OUT "$final_divider\n";
				
		my @fail_rows = @{$self->{fail_rows}};
		for my $field_array ( @fail_rows) {
			my $line = join("\t",@{$field_array});
		    print $OUT $line,"\n";
			my $rest = join("^^^",  'FAIL',
									$field_array->[$header_col_lookup{var_type}],
									$field_array->[$header_col_lookup{ref_allele}]."->".$field_array->[$header_col_lookup{var_allele}],
		 	 						'ref_count:'.$field_array->[$header_col_lookup{ref_allele_count}],
		 	 						'var_count:'.$field_array->[$header_col_lookup{var_allele_count}],
		 	 						'snv_score:'.$field_array->[$header_col_lookup{var_score}],
		 	 						'clr_score:'.$field_array->[$header_col_lookup{clr_score}],
		 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}],
		 	 						'normal_alleles:'.$field_array->[$header_col_lookup{normal_alleles}]
								);
			$rest =~ s/ /_/g;
			print $FAIL join("\t",$field_array->[$header_col_lookup{chr}], $field_array->[$header_col_lookup{start_coord}], $field_array->[$header_col_lookup{end_coord}], $rest) ."\n";
		}
		
#		my @fail_rows = @{$self->{fail_rows}};
#		for my $line ( sort {my @afields = split("\t",$a); my @bfields = split ("\t",$b); $afields[0] cmp $bfields[0] || $afields[1] <=> $bfields[1]} @fail_rows) {
#		    print $OUT $line,"\n";
#		    my @cols = split("\t",$line);
#		 	my $rest = join("^^^",  @cols[3,4,5,6,7,8,9]);
#		 	$rest =~ s/ /_/g;
#		 	print $FAIL join("\t",$cols[0], $cols[1], $cols[2], $rest) ."\n";
#		}
	}
	
	close($OUT);
    return 1;
}

#Generate the pass fail status
sub generate_pass_fail {
	my $self = shift;
	my %args = @_;
	my $debug = defined $args{-debug}?1:0;

	my $variant_type = $self->{variant_type};
	my @lines = ();
	
	if ($variant_type eq 'snv') {
	    # Get the snvs to report on
	    my %snv_ids = $self->_get_variants($args{-chr},$args{-start},$args{-end});
		#Get the formatted lines
		@lines = $self->_generate_snv_lines(-snv_ids=>\%snv_ids);
	} else {
		my %indel_ids = $self->_get_variants($args{-chr},$args{-start},$args{-end});
		@lines = $self->_generate_indel_lines(-indel_ids=>\%indel_ids);
	}

	$self->_sort_pass_fail(-lines=>\@lines);
}

#sort the lines into pass/fail and put into data structures
sub _sort_pass_fail {
	my $self = shift;
	my %args = @_;
	
	my @required_args = (
						'-lines'
			 			);

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    	
    my $variant_type = $self->{variant_type};
    
    my @pass_rows;
    my @rare_allele_rows;
    my %fail_rows;
	my @lines = @{$args{-lines}};
	my %filter_set = %{$self->{filter_set}};
	my %report_filters = %{$self->{report_set}};
	my %header_col_lookup = %{$self->{header_lookup}};
	my %filter_files = keys %{$self->{filter_files}}?%{$self->{filter_files}}:();
	
	#Keep track of genes filtered out for reporting			
	my %filtered_genes;			
	
	my $dbsnp_column_number;
	if ($variant_type eq 'snv') {
		$dbsnp_column_number = $header_col_lookup{filter_dbsnp_snv};
	} else {
		$dbsnp_column_number = $header_col_lookup{filter_dbsnp_indel};
	}
	
	my %filter_count = ();
	my $line_count = @lines;
	
	if ($variant_type eq 'snv') {
		$filter_count{filter_snv}{pass} = $line_count; 	
	} else {
		$filter_count{filter_indel}{pass} = $line_count;
	}
	
ROW:	
	for my $line (@lines) {


		my $dbsnp_fail = 0;
		chomp $line;
		my @cols = split ("\t",$line);




		my $gene;
		if (!defined $cols[$header_col_lookup{$self->{gene_col_name}}]) {
			#Cosmic overlaps don't need to overlap genes so don't report this warning
			#modules::Exception->warning("No gene for snv entry with line ($line)");
		} else {
			$gene = $cols[$header_col_lookup{$self->{gene_col_name}}];
		}
	
		#First get the filter pass counts
		for my $report_filter ( keys %report_filters ) {
		    #Keep a tally of the number of records passing each filter
			if (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] eq 'PASS') {
				$filter_count{$report_filter}{pass}++;
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] eq 'FAIL') {
				$filter_count{$report_filter}{fail}++;
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] == 1) {
				$filter_count{$report_filter}{pass}++;
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] == 0) {
				$filter_count{$report_filter}{fail}++;
		    }
		}
		
		if ($gene) {
			#Handle any additional file filters of gene names
			for my $file_filter ( keys %filter_files ) {
				if (exists $filter_files{$file_filter}->{$gene}) {
					$filtered_genes{$gene}++;
					$filter_count{$file_filter}{fail}++;
				} else {
					$filter_count{$file_filter}{pass}++;
				}
			}
		}	
		
		if ($cols[$dbsnp_column_number] eq 'PASS') {
			$cols[$dbsnp_column_number] = 'NOVEL';
		}
		
		
		#Entries with a cutoff
	    foreach my $filter (keys %filter_set) {
			
			#Handle the special cases first (currently only snvs)
			#These are all snv specific filters
			if ($filter_set{$filter} =~ /allele_freq=([0-9\.]+)/) {
				my $cutoff = $1;
				my $allele_freq_str = 	$cols[$header_col_lookup{dbsnp_allele_freq}];

				#special filtering for dbsnp human where rare alleles are still a pass
				if ($cols[$dbsnp_column_number] eq 'NOVEL') {
					#Skip if we're already passed or don't have the data
					next;
				} else {
					if ($allele_freq_str =~ /:R([0-9\.]+)/) {
						my ($dbsnp_allele_freq) = $allele_freq_str =~ /:R([0-9\.]+)/;
						if ($dbsnp_allele_freq > $cutoff) {
							my $local_line = join("\t",@cols);
							$local_line =~ s/OVERALL_PASS/filter_dbsnp_snv FAIL/;
							$fail_rows{$local_line} = 1;
							next ROW;
						} else {
							$cols[$dbsnp_column_number] = 'RARE_ALLELE';
							next;
				    	}
						
					} else {
						#matching cases with no allele info (reads)
						my $local_line = join("\t",@cols);
						$local_line =~ s/OVERALL_PASS/filter_dbsnp_snv FAIL/;
						$fail_rows{$local_line} = 1;
						next ROW;
					}
								
				}
			} elsif ($filter eq 'base_change') {
				#filter out base changes we don't want
				my @changes = split(",",$filter_set{$filter});
				my $current_change = $cols[$header_col_lookup{ref_allele}] . '->' . $cols[$header_col_lookup{var_allele}];
				my $matched = 0;
				
				for my $change (@changes) {
					if ($change =~ /^\!/) {
						$change =~ s/^\!//;
						#If it matches a change we're filtering out
						if ($current_change eq $change) {
							my $local_line = join("\t",@cols);
							$local_line =~ s/OVERALL_PASS/base_change FAIL/;
							$fail_rows{$local_line} = 1;
							next ROW;
						}
					} else {
						modules::Exception->throw("ERROR: Can only filter out base changes");
					}
				}
			}  elsif (!exists $header_col_lookup{$filter}) {
				#Catchall; shouldn't get here as all the special cases should have been handled already
				modules::Exception->throw("ERROR: Problem with pass/fail analysis for filter $filter");
			} elsif ($cols[$header_col_lookup{$filter}] =~ /FAIL/) {
				#If it's a pass/fail field (eg filter_dbsnp_snv)
				my $local_line = join("\t",@cols);
				$local_line =~ s/OVERALL_PASS/$filter FAIL/;
				$fail_rows{$local_line} = 1;
				next ROW;
			} elsif ($cols[$header_col_lookup{$filter}] =~ /^\d/ && $cols[$header_col_lookup{$filter}] < $filter_set{$filter}) {
				#If it's a value cutoff (eg read_depth, median_quality_score)
				my $local_line = join("\t",@cols);
				$local_line =~ s/OVERALL_PASS/$filter < $filter_set{$filter} FAIL/;			
				$fail_rows{$local_line} = 1;
				next ROW;
			} 
	    }
	
	
		if ($gene) {
			#Additional file cutoffs
			for my $file_filter ( keys %filter_files ) {
				if (exists $filter_files{$file_filter}->{$gene}) {
					$filter_count{$file_filter}{fail}++;
					my $local_line = join("\t",@cols);
					$local_line =~ s/OVERALL_PASS/problem_gene FAIL/;
					$fail_rows{$local_line}  = 1;
					next ROW;
				}
			}  
		}
		
	    my $or_pass = 0;

		#If this field is either NON-SYN or SPLICE then the row is passed
		if ($variant_type eq 'snv' && $cols[$header_col_lookup{snv_exon_type}] ne 'SYN') {
			$or_pass = 1;
		} elsif ($variant_type eq 'indel') {
			$or_pass = 1;
		}
	
	    unless ($or_pass) {
	    	my $local_line = join("\t",@cols);
			$local_line =~ s/OVERALL_PASS/synonomous FAIL/;
	    	$fail_rows{$local_line} = 1;
			next ROW;
	    }
	    
	    #Finally determine whether it's rare allele or full pass
	    if ($cols[$dbsnp_column_number] eq 'RARE_ALLELE') {
			push @rare_allele_rows, \@cols;    	
	    } else {
		    push @pass_rows, \@cols;
	    }
	    
	}
		
	#Now generate the data structure used for custom sorting the passed rows
	my %sorted_pass_rows;
	my %sorted_allele_rows;
	
	if ($variant_type eq 'snv') {
	
		foreach my $pass_row (@pass_rows) {
			my $mutant_depth = $pass_row->[$header_col_lookup{read_depth}];
			
			if ($pass_row->[$header_col_lookup{snv_exon_type}] =~ /SPLICE/) {
		    	$sorted_pass_rows{'SNV'}{1}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{coord}]} = $pass_row;
		    } else {
		    	$sorted_pass_rows{'SNV'}{0}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{coord}]} = $pass_row;
		    }
		}
		
		foreach my $allele_row (@rare_allele_rows) {
			my $mutant_depth = $allele_row->[$header_col_lookup{read_depth}];
			
			if ($allele_row->[$header_col_lookup{snv_exon_type}] =~  /SPLICE/) {
		    	$sorted_allele_rows{'SNV'}{1}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{coord}]} = $allele_row;
		    } else {
				$sorted_allele_rows{'SNV'}{0}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{coord}]} = $allele_row;
		    }
		}
		
		
	} else {
		foreach my $pass_row (@pass_rows) {
			my $indel_type = $pass_row->[$header_col_lookup{var_type}];
			my $mutant_depth = $pass_row->[$header_col_lookup{read_depth}];
			
			if ($pass_row->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
		    	$sorted_pass_rows{$indel_type}{1}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{start_coord}]}{$pass_row->[$header_col_lookup{end_coord}]}{$pass_row->[$header_col_lookup{var_allele}]} = $pass_row;
		    } else {
		    	$sorted_pass_rows{$indel_type}{0}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{start_coord}]}{$pass_row->[$header_col_lookup{end_coord}]}{$pass_row->[$header_col_lookup{var_allele}]} = $pass_row;
		    }
		}
		
		foreach my $allele_row (@rare_allele_rows) {
			my $indel_type = $allele_row->[$header_col_lookup{var_type}];
			my $mutant_depth = $allele_row->[$header_col_lookup{read_depth}];
			
			if ($allele_row->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
		    	$sorted_allele_rows{$indel_type}{1}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{start_coord}]}{$allele_row->[$header_col_lookup{end_coord}]}{$allele_row->[$header_col_lookup{var_allele}]} = $allele_row;
		    } else {
		    	$sorted_allele_rows{$indel_type}{0}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{start_coord}]}{$allele_row->[$header_col_lookup{end_coord}]}{$allele_row->[$header_col_lookup{var_allele}]} = $allele_row;
		    }
		}
		
	}
	
	my @fail_rows = ();
	for my $fail_row (sort {my @afields = split("\t",$a); my @bfields = split ("\t",$b); $afields[0] cmp $bfields[0] || $afields[1] <=> $bfields[1]} keys %fail_rows) {
		my @cols = split ("\t",$fail_row);
		push @fail_rows, \@cols;
	}
	
	my @all_rows = (@pass_rows,@rare_allele_rows,@fail_rows);
	$self->{all_rows} = \@all_rows;
	$self->{pass_rows} = \@pass_rows;
	$self->{allele_rows} = \@rare_allele_rows;
	$self->{fail_rows} = \@fail_rows;
	$self->{filter_count} = \%filter_count;
	$self->{sorted_pass_rows} = \%sorted_pass_rows;
	$self->{sorted_allele_rows} = \%sorted_allele_rows;
	$self->{filtered_genes} = \%filtered_genes;
}

#Get the formatted indel lines
sub _generate_indel_lines {
	my $self = shift;
	my %args = @_;

	my @required_args = (
						'-indel_ids'
			 			);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $indel_ids = $args{-indel_ids};
	my $gene_mapper = $self->{gene_mapper};
	
	#Needed for indel_row entries
	my %indel_lookup = ();
	
	my @filter_names;
    foreach my $filter (@{$self->{filter_info}}){
		push @filter_names, $filter->name;
    }

	my @lines = ();
	foreach my $indel_id (keys %{$indel_ids}) {
		my ($indel) = modules::Adaptors::Variant->search(id => $indel_id);
		
		my $filt_indel
		    = modules::FilteredVariant->new('variant' => $indel,
						    'filter_list_order' => $self->{filter_info});
		

		
	
#		my $tumour_other_alleles = join("\:", @{$pl_tumour->base_frequencies_indel->bases}) . '_' . join("\:", @{$pl_tumour->base_frequencies_indel->counts});
#		
#		my $pl_normal = modules::PileupLine->new();
#		$pl_normal->base_string_indel($indel->normal_base_string);
#	
#		my $normal_alleles = join("\:", @{$pl_normal->base_frequencies_indel->bases}) . '_' . join("\:", @{$pl_normal->base_frequencies_indel->counts});
		
		my $pl_tumour = modules::PileupLine->new();
		$pl_tumour->base_string_indel($indel->tumour_base_string);
		
		my $tumour_other_alleles = "N/A";
		my $ref_allele_count = 0;
		my $ref_allele = 'N/A';
		my $var_allele_count = 0;
		my $var_bases;

		if ($indel->var_type eq 'INS') {
			$var_bases =  '+'. length($indel->inserted_bases) . $indel->inserted_bases;
		} else {
			my $length = $indel->end_coord - $indel->start_coord + 1;
			
			$ref_allele = $self->{adaptor}->fetch_by_region('chromosome',$indel->chr, $indel->start_coord, $indel->end_coord)->seq();;
			$var_bases = '-'. $length .$ref_allele;
		}
		
		my @tumour_bases = @{$pl_tumour->base_frequencies_indel->bases};
		my @tumour_counts = @{$pl_tumour->base_frequencies_indel->counts};
		my @tumour_remaining_alleles = ();
		
		for ( my $count = 0 ; $count < @tumour_bases ; $count++ ) {
			if ($tumour_bases[$count] eq 'REF') {
				$ref_allele_count = $tumour_counts[$count];
			} elsif ($tumour_bases[$count] eq $var_bases) {
				$var_allele_count = $tumour_counts[$count];
			} else {
				push @tumour_remaining_alleles,"$tumour_bases[$count]:$tumour_counts[$count]";
			}
		}
	
		if (@tumour_remaining_alleles) {
			$tumour_other_alleles = join('_', @tumour_remaining_alleles);
		}
	
		my $pl_normal = modules::PileupLine->new();
		$pl_normal->base_string_indel($indel->normal_base_string);
	
		my @normal_bases = @{$pl_normal->base_frequencies_indel->bases};
		my @normal_counts = @{$pl_normal->base_frequencies_indel->counts};
		my @normal_alleles = ();
	
		for ( my $count = 0 ; $count < @normal_bases ; $count++ ) {
			push @normal_alleles, "$normal_bases[$count]:$normal_counts[$count]"
		}
	
		my $normal_alleles = join("_", @normal_alleles);
		
		
		my %filter_values;
		my $exon_overlap;
		my $aa_position = "N/A";
		my $aa_length = "N/A";			
		
		#Get the variant filter info
		foreach my $filter_name (@filter_names) {
			if ($filter_name eq $EXON_FILTER) { 
				if ($filt_indel->variant_filter_hash->{$filter_name} ne 'nodata') {
						if ($filt_indel->variant_filter_hash->{$filter_name}->filtermatch == 1) {
							$exon_overlap = "EXON";
							my $attr = $filt_indel->variant_filter_hash->{$EXON_FILTER}->attribute;
				
							if ($attr =~ /aa_pos\=([^\;]+)/) {
								$aa_position = $1;
							}
							if ($attr =~ /aa_len\=([^\;]+)/) {
								$aa_length = $1;
							}
						}
				}
			} elsif ($filter_name eq $SPLICE_FILTER) { 
				if ($filt_indel->variant_filter_hash->{$filter_name} ne 'nodata') {
						if ($filt_indel->variant_filter_hash->{$filter_name}->filtermatch == 1) {
							my $attr = $filt_indel->variant_filter_hash->{$SPLICE_FILTER}->attribute;
							if ($attr =~ /splice_dist\=([^\;]+)/) {
								$exon_overlap = "SPLICE (".$1.")";
							}							
						}
				}
			} else {
				if ($filt_indel->variant_filter_hash->{$filter_name} ne 'nodata') {
					if ($filt_indel->variant_filter_hash->{$filter_name}->filterpass == 1) {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} else {
					$filter_values{$filter_name} = 'FAIL';
				}
			}
			
			
		}
		
		#Get the allele freq if relevant
		my $read_allele_freq = "reads:R".$indel->ref_base_freq.':V'.$indel->var_base_freq;
		my $dbsnp_allele_freq = "N/A";
		
		if ($filt_indel->variant_filter_hash->{$DBSNP_FILTER_INDEL} ne 'nodata' && $filt_indel->variant_filter_hash->{$DBSNP_FILTER_INDEL}->filtermatch eq '1') {
	    	if ($filt_indel->variant_filter_hash->{$DBSNP_FILTER_INDEL}->attribute) {
		    	my $attr = $filt_indel->variant_filter_hash->{$DBSNP_FILTER_INDEL}->attribute;
				if ($attr =~ /rs=(rs\d+);ref_allele_freq\=([0-9\.]+);var_allele_freq\=([0-9\.]+)/) {
					$dbsnp_allele_freq = 'dbsnp:'. $1 . ':R' . $2 . ':V'. $3;
				} elsif ($attr =~ /rs=(rs\d+)/) {
					$dbsnp_allele_freq = 'dbsnp:'.$1.':No_freq';
				} else {
					$dbsnp_allele_freq = 'DBSNP_Diff_Allele';
				}
	    	}
		}
		
		
		my $cosmic_coord = 'N/A';

		if ($filt_indel->variant_filter_hash->{$COSMIC_FILTER} ne 'nodata' && $filt_indel->variant_filter_hash->{$COSMIC_FILTER}->filtermatch eq '1') {
	    	if ($filt_indel->variant_filter_hash->{$COSMIC_FILTER}->attribute) {
	    		my $attr = $filt_indel->variant_filter_hash->{$COSMIC_FILTER}->attribute;
	    		if ($attr =~ /melanoma/) {
	    			$cosmic_coord = "MELANOMA:".$attr;
	    		} else {
	    			$cosmic_coord = "OTHERCANCER:".$attr;
	    		}	
	    	}
		}
		
		
		my $gene_info = $gene_mapper->get_gene_from_coord(-chr => $indel->chr, -coord => $indel->start_coord);
	
		#If the start coord doesn't overlap use the end coord
		if (!exists $gene_info->{'ensembl'}) {
			$gene_info = $gene_mapper->get_gene_from_coord(-chr => $indel->chr, -coord => $indel->end_coord);
		}			
	
		#Store the annotation_header values for reporting
		my %annotation_headers = ();
	
		$annotation_headers{ensembl} = join(",",keys %{$gene_info->{ensembl}});
	
		for my $annotation_header (@{$self->{annotation_headers}}) {
			$annotation_headers{$annotation_header} = join(",",keys %{$gene_info->{$annotation_header}});
		}
		
	
#		my $ens_name = defined $gene_info->{ensembl}?join(",",keys %{$gene_info->{exon}}):'';
#		
#		my $ccds_name = defined $gene_info->{ccds}?join(",",keys %{$gene_info->{ccds}}):'';
#		my $hgnc = defined $gene_info->{gene_name}?join(",",keys %{$gene_info->{gene_name}}):'';
#		my $gene_desc = defined $gene_info->{gene_desc}?join(",",keys %{$gene_info->{gene_desc}}):'';
#		my $gene_GO = defined $gene_info->{go}?join(",",keys %{$gene_info->{go}}):'';
#		my $uniprot_name = defined $gene_info->{uniprot_name}?join(",",keys %{$gene_info->{uniprot_name}}):'';
#		my $omim = defined $gene_info->{omim}?join(",",keys %{$gene_info->{omim}}):'';
#		my $refseq_name = defined $gene_info->{refseq}?join(",",keys %{$gene_info->{refseq}}):'';
		
		my $var_type;
		my $var_allele;
		my $var_length;
		
		if ($indel->var_type eq 'INS') {
			$var_type = 'INS';
			$var_allele = '+' . $indel->inserted_bases;
			$var_length = length($var_allele) - 1;
		} elsif ($indel->var_type eq 'DEL') {
			$var_type = 'DEL';
			$var_allele = 'N/A'; 
			$var_length = $indel->end_coord - $indel->start_coord + 1;
		} else {
			my $indel_str = $indel->chr.':'. $indel->start_coord.'-' .$indel->end_coord;
			modules::Exception->throw("ERROR: Can't get var_type for indel $indel_str");
		}
		
		$indel_lookup{$indel->chr}{$indel->start_coord}{$indel->end_coord}{$var_type}{$var_allele} = $indel_id;
		
		
		my @entries = ();
		
		#Iterate over the headers; the variables will be the same name except in the special cases; this allows us to have custom reports based on the headers
		for my $header (@{$self->{headers}}) {
			my $variable = '$'.$header;
			my $value;
			
			if (exists $filter_values{$header}) {
				$value = $filter_values{$header};
			} elsif ($header eq 'chr') {
				$value = $indel->chr;
			} elsif ($header eq 'start_coord') {
				$value = $indel->start_coord;
			} elsif ($header eq 'end_coord') {
				$value = $indel->end_coord;
			} elsif ($header eq 'var_score') {
				$value = $indel->var_score;
			} elsif ($header eq 'clr_score') {
				$value = $indel->clr_score;
			} elsif ($header eq 'median_quality_score') {
				$value = $indel->median_quality_score;
			} elsif ($header eq 'read_depth') {
				$value = $indel->read_depth;
			} elsif ($header eq 'final_status') {
		    	#Don't know this yet so use a default of pass
		    	$value = 'OVERALL_PASS'
		    } elsif (exists $annotation_headers{$header}) {
		    	$value = $annotation_headers{$header};
		    } else {
				#Here the variable name matches
				$value = eval($variable);
			}
			
			if ($value !~ /\w/) {
				my $chr = $indel->chr;
				my $coord = $indel->start_coord;
				if ($cosmic_coord eq 'N/A') {
					print "ERROR $header has no value for $chr $coord\n";
				} else {
					print "ERROR $header has no value for $chr $coord\n" unless exists $annotation_headers{$header};
				}
			}
			
			push @entries, $value;
			
			#print "$header $value\n";
		}
		
		my $line_str = join("\t",@entries);
		push @lines,$line_str;
    }
    
    $self->{indel_lookup} = \%indel_lookup;
    
	return @lines;
		
	

}


#Get the formatted summary lines
sub _generate_snv_lines {
	my $self = shift;
	my %args = @_;

	my @required_args = (
						'-snv_ids'
			 			);


    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $snv_ids = $args{-snv_ids};


	my $gene_mapper = $self->{gene_mapper};

	my @filter_names;
    foreach my $filter (@{$self->{filter_info}}){
		push @filter_names, $filter->name;
    }

	#Lookup for creating snv_rows laters (need snv_id from db)
	my %snv_lookup = ();

	my @lines = ();
	foreach my $snv_id (keys %{$snv_ids}) {
		
		my ($snv) = modules::Adaptors::SNV->search(id => $snv_id);
	
		my $filt_snv
		    = modules::FilteredSNV->new('snv' => $snv,
						    'filter_list_order' => $self->{filter_info});
	
		my $pl_tumour = modules::PileupLine->new();
		$pl_tumour->ref_base($snv->ref_base);
		$pl_tumour->base_string($snv->tumour_base_string);
	
		#my $tumour_other_alleles = join("\:", @{$pl_tumour->base_frequencies->bases}) . '_' . join("\:", @{$pl_tumour->base_frequencies->counts});
		my $tumour_other_alleles = "N/A";
		my $ref_allele_count = 0;
		my $var_allele_count = 0;
		
		my @tumour_bases = @{$pl_tumour->base_frequencies->bases};
		my @tumour_counts = @{$pl_tumour->base_frequencies->counts};
		my @tumour_remaining_alleles = ();
		print $snv->chr ." ". $snv->coord."\n";
		print Dumper \@tumour_bases;
		
		for ( my $count = 0 ; $count < @tumour_bases ; $count++ ) {
			if ($tumour_bases[$count] eq $snv->ref_base) {
				$ref_allele_count = $tumour_counts[$count];
			} elsif ($tumour_bases[$count] eq $snv->var_base) {
				$var_allele_count = $tumour_counts[$count];
			} else {
				push @tumour_remaining_alleles,"$tumour_bases[$count]:$tumour_counts[$count]";
			}
		}
	
		if (@tumour_remaining_alleles) {
			$tumour_other_alleles = join('_', @tumour_remaining_alleles);
		}
	
		my $pl_normal = modules::PileupLine->new();
		$pl_normal->ref_base($snv->ref_base);
		$pl_normal->base_string($snv->normal_base_string);
	
		my @normal_bases = @{$pl_normal->base_frequencies->bases};
		my @normal_counts = @{$pl_normal->base_frequencies->counts};
		my @normal_alleles = ();
	
		for ( my $count = 0 ; $count < @normal_bases ; $count++ ) {
			push @normal_alleles, "$normal_bases[$count]:$normal_counts[$count]"
		}
	
		my $normal_alleles = join("_", @normal_alleles);
	
	
		# try and parse the aa change, is possible
		my $aa_change = "N/A";
		my $polyphen_prediction = "N/A";
		my $polyphen_score = "N/A";
		my $sift_prediction = "N/A";
		my $sift_score = "N/A";
					
	    #Check the filter entry exists and that it's passed
	    if ($filt_snv->snv_filter_hash->{$EXON_NS_FILTER} ne 'nodata' && $filt_snv->snv_filter_hash->{$EXON_NS_FILTER}->filtermatch eq '1') {
	    	my $attr = $filt_snv->snv_filter_hash->{$EXON_NS_FILTER}->attribute;
			#aa change is required
			if ($attr =~ /aa_change\=([^\;]+)/) {
				$aa_change = $1;
				if ($attr =~ /combined\=([^\;]+)/) {
					$aa_change .= " (COMBINED:$1)";	
				}
			} else {
				modules::Exception->throw("ERROR: snv filter_exon_ns doesn't have the aa change");
			}
			
			
			#polyphen and sift may be absent
			if ($attr =~ /poly_pred\=([^\;]+)/) {
				$polyphen_prediction = $1;
			}
			if ($attr =~ /poly_score\=([^\;]+)/) {
				$polyphen_score = $1;
			}
			if ($attr =~ /sift_pred\=([^\;]+)/) {
				$sift_prediction = $1;
			}
			if ($attr =~ /sift_score\=([^\;]+)/) {
				$sift_score = $1;
			}
			
		} 

				


		#Get the allele freq if relevant
		my $read_allele_freq = "reads:R".$snv->ref_base_freq.':V'.$snv->var_base_freq;
		my $dbsnp_allele_freq = "N/A";
		
		if ($filt_snv->snv_filter_hash->{$DBSNP_FILTER} ne 'nodata' && $filt_snv->snv_filter_hash->{$DBSNP_FILTER}->filtermatch eq '1') {
	    	if ($filt_snv->snv_filter_hash->{$DBSNP_FILTER}->attribute) {
		    	my $attr = $filt_snv->snv_filter_hash->{$DBSNP_FILTER}->attribute;
				if ($attr =~ /rs=(rs\d+);ref_allele_freq\=([0-9\.]+);var_allele_freq\=([0-9\.]+)/) {
					$dbsnp_allele_freq = 'dbsnp:'. $1 . ':R' . $2 . ':V'. $3;
				} elsif ($attr =~ /rs=(rs\d+)/) {
					$dbsnp_allele_freq = 'dbsnp:'.$1.':No_freq';
				} else {
					$dbsnp_allele_freq = 'DBSNP_Diff_Allele';
				}
	    	}
		}
		
		my $cosmic_coord = 'N/A';

		if ($filt_snv->snv_filter_hash->{$COSMIC_FILTER} ne 'nodata' && $filt_snv->snv_filter_hash->{$COSMIC_FILTER}->filtermatch eq '1') {
	    	if ($filt_snv->snv_filter_hash->{$COSMIC_FILTER}->attribute) {
	    		my $attr = $filt_snv->snv_filter_hash->{$COSMIC_FILTER}->attribute;
	    		if ($attr =~ /melanoma/) {
	    			$cosmic_coord = "MELANOMA:".$attr;
	    		} else {
	    			$cosmic_coord = "OTHERCANCER:".$attr;
	    		}	
	    	}
		}

  		my ($gene_info) = $gene_mapper->get_gene_from_coord(-chr => $snv->chr, -coord => $snv->coord);
	
		#my $ens_name = defined $gene_info->{ensembl}?join(",",keys %{$gene_info->{ensembl}}):'';
#		my $nt_change = $snv->ref_base .'->'. $snv->var_base;
#		
#		if (exists $polyphen{$snv->chr} && exists $polyphen{$snv->chr}{$snv->coord} && exists $polyphen{$snv->chr}{$snv->coord}{$nt_change}) {
#			($polyphen_prediction,$polyphen_score) = split(",",$polyphen{$snv->chr}{$snv->coord}{$nt_change});
#		}
#	
#		if (exists $polyphen_info{$snv->chr} && exists $polyphen_info{$snv->chr}{$snv->coord}) {
#			$polyphen_info = $polyphen_info{$snv->chr}{$snv->coord};
#		}
	
		my $snv_exon_type = 'SYN';
		my %filter_values;
		my $aa_position = "N/A";
		my $aa_length = "N/A";
		
		#Get the snv filter info
		foreach my $filter_name (@filter_names) {
			if ($filter_name eq $EXON_NS_FILTER) { 
				if ($filt_snv->snv_filter_hash->{$EXON_NS_FILTER} ne 'nodata') {
						if ($filt_snv->snv_filter_hash->{$EXON_NS_FILTER}->attribute) {
	    				my $attr = $filt_snv->snv_filter_hash->{$EXON_NS_FILTER}->attribute;
	    				if ($filt_snv->snv_filter_hash->{$EXON_NS_FILTER}->filtermatch == 1 && $attr =~ /Stop/) {
	    					$snv_exon_type = "NONSENSE";
	    				} elsif ($filt_snv->snv_filter_hash->{$EXON_NS_FILTER}->filtermatch == 1) {
							$snv_exon_type = "MISSENSE";
						}
					}
				}
			} elsif ($filter_name eq $EXON_FILTER) {
				if ($filt_snv->snv_filter_hash->{$EXON_FILTER} ne 'nodata') {
					if ($filt_snv->snv_filter_hash->{$EXON_FILTER}->filtermatch eq '1') {
						my $attr = $filt_snv->snv_filter_hash->{$EXON_FILTER}->attribute;
			
						if ($attr =~ /aa_pos\=([^\;]+)/) {
							$aa_position = $1;
						}
						if ($attr =~ /aa_len\=([^\;]+)/) {
							$aa_length = $1;
						}
					}
				}
			
			} elsif ($filter_name eq $SPLICE_FILTER) { 
				if ($filt_snv->snv_filter_hash->{$filter_name} ne 'nodata') {
						if ($filt_snv->snv_filter_hash->{$filter_name}->filtermatch == 1) {
							my $attr = $filt_snv->snv_filter_hash->{$SPLICE_FILTER}->attribute;
							if ($attr =~ /splice_dist\=([^\;]+)/) {
								$snv_exon_type = "SPLICE (".$1.")";
							}							
						}
				}
			} else {
				if ($filt_snv->snv_filter_hash->{$filter_name} ne 'nodata') {
					if ($filt_snv->snv_filter_hash->{$filter_name}->filterpass == 1) {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} else {
					$filter_values{$filter_name} = 'FAIL';
				}
			}
			
			
		}
		
		#Store the annotation_header values for reporting
		my %annotation_headers = ();
	
		$annotation_headers{ensembl} = join(",",keys %{$gene_info->{ensembl}});
	
		for my $annotation_header (@{$self->{annotation_headers}}) {
			$annotation_headers{$annotation_header} = join(",",keys %{$gene_info->{$annotation_header}});
		}
			
#		my $ccds_name = defined $gene_info->{ccds}?join(",",keys %{$gene_info->{ccds}}):'';
#		my $hgnc = defined $gene_info->{gene_name}?join(",",keys %{$gene_info->{gene_name}}):'';
#		my $gene_desc = defined $gene_info->{gene_desc}?join(",",keys %{$gene_info->{gene_desc}}):'';
#		my $gene_GO = defined $gene_info->{go}?join(",",keys %{$gene_info->{go}}):'';
#		my $uniprot_name = defined $gene_info->{uniprot_name}?join(",",keys %{$gene_info->{uniprot_name}}):'';
#		my $omim = defined $gene_info->{omim}?join(",",keys %{$gene_info->{omim}}):'';
#		my $refseq_name = defined $gene_info->{refseq}?join(",",keys %{$gene_info->{refseq}}):'';
#		
		
		#Need this lookup for creating snv_rows later
		$snv_lookup{$snv->chr}{$snv->coord}{$snv->var_base} = $snv_id;
		my @entries = ();
		
		#Iterate over the headers; the variables will be the same name except in the special cases; this allows us to have custom reports based on the headers
		for my $header (@{$self->{headers}}) {
			my $variable = '$'.$header;
			my $value;
			
			if (exists $filter_values{$header}) {
				$value = $filter_values{$header};
			} elsif ($header eq 'chr') {
				$value = $snv->chr;
			} elsif ($header eq 'coord') {
				$value = $snv->coord;
			} elsif ($header eq 'ref_allele') {
				$value = $snv->ref_base;
			} elsif ($header eq 'var_allele') {
				$value = $snv->var_base;
			} elsif ($header eq 'snv_score') {
				$value = $snv->snv_score;
			} elsif ($header eq 'clr_score') {
				$value = $snv->clr_score;
			} elsif ($header eq 'median_quality_score') {
				$value = $snv->median_quality_score;
			} elsif ($header eq 'read_depth') {
				$value = $snv->read_depth;
			} elsif ($header eq 'final_status') {
		    	#Don't know this yet so use a placeholder
		    	$value = 'OVERALL_PASS'
		    } elsif (exists $annotation_headers{$header}) {
		    	$value = $annotation_headers{$header};
		    } else {
				#Here the variable name matches
				$value = eval($variable);
			}
			
			push @entries, $value;
			
			if ($value !~ /\w/) {
				my $chr = $snv->chr;
				my $coord = $snv->coord;
				#Cosmic coords don't have to overlap gene coords
				if ($cosmic_coord eq 'N/A') {
					print "ERROR $header has no value for $chr $coord\n";
				} else {
					print "ERROR $header has no value for $chr $coord\n" unless exists $annotation_headers{$header};
				}
			}
			
			#print "$header $value\n";
		}
	

		my $line_str = join("\t",@entries);
		push @lines,$line_str;
		
		#print "LINE $line_str\n";
		
    }
    
    $self->{snv_lookup} = \%snv_lookup;
    
	return @lines;
}

#Getter for passed rows
sub get_fail_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};		
	
	for my $field_array ( @{$self->{fail_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		}
	}
	
	return \%sorted_rows;
}

#Getter for passed rows
sub get_pass_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{pass_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		}
	}
	
	return \%sorted_rows;
}

#Getter for allele rows
sub get_allele_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{allele_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		}
	}
	
	return \%sorted_rows;
}


#Getter for all the rows
sub get_all_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{all_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		}
	}
	
#	for my $line ( @{$self->{all_rows}}) {
#		my @cols = split("\t",$line);
#		if ($variant_type eq 'snv') {
#			$sorted_rows{$cols[0]}{$cols[1]} = $line;
#		} else {
#			$sorted_rows{$cols[0]}{$cols[1]}{$cols[2]}{$cols[3]}{$cols[4]} = $line;
#		}
#	}
	
	return \%sorted_rows;
}


#Get the snvs to report on from the primary filter field
sub _get_variants {
	my $self = shift;
	my ($chr,$start,$end) = @_;
    my $var_id_iterator;  
    my  %var_ids = ();
    my @primary_filter_ids;
    
    foreach my $primary_filter (@{$self->{primary_filters}}){
		push @primary_filter_ids, $primary_filter->id;
    }

	my $var_count = 0;

	my $variant_type = $self->{variant_type};

    foreach my $primary_filter_id (@primary_filter_ids) {

		if (defined $chr && defined $start && defined $end) {
			if ($variant_type eq 'snv') {
		    	$var_id_iterator = modules::Adaptors::SNV->search_region_by_experiment_and_passed_filter($self->run->id, 
												     $primary_filter_id,
												     $chr, 
												     $start, 
												     $end);
			} else {
				$var_id_iterator = modules::Adaptors::Variant->search_region_by_experiment_and_passed_filter($self->run->id, 
												     $primary_filter_id,
												     $chr, 
												     $start, 
												     $end);
			}
	
		} elsif (defined $chr) {
			if ($variant_type eq 'snv') {
		    	$var_id_iterator = modules::Adaptors::SNV->search_chr_by_experiment_and_passed_filter($self->run->id, 
												     $primary_filter_id,
												     $chr);
			} else {
				$var_id_iterator = modules::Adaptors::Variant->search_chr_by_experiment_and_passed_filter($self->run->id, 
												     $primary_filter_id,
												     $chr);
			}
			
		} else {
			if ($variant_type eq 'snv') {
		    	$var_id_iterator = modules::Adaptors::SNV->search_by_experiment_and_passed_filter($self->run->id, 
												     $primary_filter_id);
			} else {
				$var_id_iterator = modules::Adaptors::Variant->search_by_experiment_and_passed_filter($self->run->id, 
												     $primary_filter_id);
			}
		}
	
		while (my $var_id = $var_id_iterator->next){
			$var_count++;
		    $var_ids{$var_id}++;
		}
    }
    return %var_ids;
}

#Print the summary stats to a file; used in the pipeline for passed snvs
sub summarize {
	my $self = shift;
	my %args = @_;
	
	my @required_args = ('-summary_file');

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $variant_type = $self->{variant_type};
	
	my $summary_file = $args{-summary_file};
	my %filter_set = %{$self->{filter_set}};
	my %report_filters = %{$self->{report_set}};
	my @primary_filters =  @{$self->{primary_filters}};
	
	my %filter_count = %{$self->{filter_count}};
	
	my %filter_files = %{$self->{filter_files}};

	
	
	open(my $SUMM, ">$summary_file")
	    or modules::Exception->throw("Unable to open file [$summary_file] for writing");
	
	  # record details of filters applied
	
	print $SUMM "Filters applied:\n";
	
	foreach my $filter (keys %filter_set){
		if ($filter_set{$filter} =~ /allele_freq=([0-9\.]+)/) {
			#Special dbsnp for human case
			print $SUMM "\t$filter (PASS DEF: no_match OR allele_freq <= $1)\n";
		} elsif ($filter_set{$filter} !~ /^[0-9]+$/) {
			#Non numerical filter
			print $SUMM "\t$filter (PASS DEF: $filter_set{$filter})\n";
		} elsif ($filter_set{$filter} == 1) {
			print $SUMM "\t$filter (PASS DEF: no_match)\n";	
		} else {
	    	print $SUMM "\t$filter (PASS DEF: >= $filter_set{$filter})\n";
		}
	}
	
	for my $file_filter ( keys %filter_files ) {
	    print $SUMM "\t$file_filter (PASS DEF: no_match)\n";
	}
	
	my @required_filters_names;
	for my $primary_filter (@primary_filters) {
		push @required_filters_names, $primary_filter->name;
	}
	my $required_filter = join (" OR ",@required_filters_names);
	print $SUMM "\tMatch EITHER $required_filter\n";
	
	
	my $variant_count;
	if ($variant_type eq 'snv') {
		$variant_count = $filter_count{filter_snv}{pass};
	} else {
		$variant_count = $filter_count{filter_indel}{pass};
		
	}
	print $SUMM "\nTotal $required_filter $variant_type called: $variant_count\n";
	
	
	#Breakdown the snv_type_exon fields
	my @all_rows = @{$self->{all_rows}};
		
	my %header_col_lookup = %{$self->{header_lookup}};
	
	
	my %var_classes = ();
	my %var_type_exons = ();
	
	for my $row ( @all_rows ) {
		if ($variant_type eq 'snv') {
			$var_type_exons{$row->[$header_col_lookup{snv_exon_type}]}++;
		} else {
			$var_classes{$row->[$header_col_lookup{var_type}]}++;
			$var_type_exons{$row->[$header_col_lookup{exon_overlap}]}++;
			
		}
	}
	
	for my $var_class_key (keys %var_classes) {
		my $pass_percent = sprintf("%.2f",$var_classes{$var_class_key} / $variant_count * 100 );
		my $lc_entry = lc($var_class_key);
		print $SUMM "\t$lc_entry MATCH: $var_classes{$var_class_key} ($pass_percent%)\n";
	}
	
	my %splice_entries = ();
	my $splice_total = 0;
	my $exonic_all = 0;	
	for my $var_exon_key (keys %var_type_exons) {
		if ($var_exon_key !~ /SPLICE/) {
			$exonic_all += $var_type_exons{$var_exon_key};
		} else {
			$splice_entries{$var_exon_key} = $var_type_exons{$var_exon_key};
			$splice_total += $var_type_exons{$var_exon_key};
			next;
		}				
		my $pass_percent = sprintf("%.2f",$var_type_exons{$var_exon_key} / $variant_count * 100 );
		my $lc_entry = lc ($var_exon_key);
		print $SUMM "\t$lc_entry MATCH: $var_type_exons{$var_exon_key} ($pass_percent%)\n";	
		
	}
	if ($splice_total) {
		my $pass_percent = sprintf("%.2f",$splice_total / $variant_count * 100 );
	
		print $SUMM "\ttotal splice MATCH: $splice_total ($pass_percent%) ; ";
	
		my $splice_line = '';
		for my $splice_entry (sort keys %splice_entries) {
			$splice_line .= "$splice_entry->$splice_entries{$splice_entry} , ";
		}
		$splice_line =~ s/ , $//;
		print $SUMM "$splice_line\n"; 
	} else {
		print $SUMM "\tfilter_splice MATCH: 0 (0%)\n";
	}
	
	if ($exonic_all > 0) {
		my $pass_percent = sprintf("%.2f",$exonic_all / $variant_count * 100 );
		print $SUMM "\tfilter_exon MATCH: $exonic_all ($pass_percent%)\n";
	} else {
		print $SUMM "\tfilter_exon MATCH: 0 (0%)\n";
	}
		
	#Print the report filters	
	for my $report_filter ( sort keys %report_filters ) {
		if (defined $filter_count{$report_filter}{pass} ) {
	   		my $pass_percent = sprintf("%.2f",$filter_count{$report_filter}{pass} / $variant_count * 100 );
	  		print $SUMM "\t$report_filter PASS: $filter_count{$report_filter}{pass} ($pass_percent%)\n";
		} else {
			print $SUMM "\t$report_filter PASS: 0 (0%)\n";
		}
	}
	 
	my $gene_filtered_count = 0;
	for my $gene ( keys %{$self->{filtered_genes}} ) {
	    $gene_filtered_count += $self->{filtered_genes}{$gene};
	}
	if ($gene_filtered_count > 0) {
		my $pass_percent = sprintf("%.2f",$gene_filtered_count / $variant_count * 100 );
		print $SUMM "\tproblem genes FAIL: $gene_filtered_count ($pass_percent%) (", join(",", keys %{$self->{filtered_genes}}), ")\n";
	} else {
		print $SUMM "\tproblem genes FAIL: 0 (0%)\n";		
	}
	
	my %sorted_pass_rows = %{$self->{sorted_pass_rows}};
	
	# tally a few items of information between rows
	my $ss_passed = 0;
	my $exonic_passed = 0;
	my $heter_count_passed = 0;
	my $homo_count_passed = 0;
	my $other_count_passed = 0;
	my $insertions_passed = 0;
	my $deletions_passed = 0;
	
	my @depth_values;
	my @quality_values;
	my %genes;

	my @pass_rows = @{$self->{pass_rows}};
	
	#Get some meta info from the passed rows
	foreach my $pass_row (@pass_rows) {
		
		if ($variant_type eq 'indel') {
			if ($pass_row->[$header_col_lookup{var_type}] eq 'DEL') {
				$deletions_passed++;
			} elsif ($pass_row->[$header_col_lookup{var_type}] eq 'INS') {
				$insertions_passed++;
			} else {
				modules::Exception->throw("ERROR: INDEL is none of accepted types");
			}
		}
		
		
		if ($variant_type eq 'snv') {
			if ($pass_row->[$header_col_lookup{snv_exon_type}] =~ /SPLICE/) {
				$ss_passed++;
			} else {
				$exonic_passed++;
			}
		} else {
			if ($pass_row->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
				$ss_passed++;
			} else {
				$exonic_passed++;				
			}
		}
		
		
		my $mutant_depth;
		if ($variant_type eq 'snv') {
			$mutant_depth = $pass_row->[$header_col_lookup{read_depth}];
		
		    if ($mutant_depth =~ /[0-9]/) {
				push @depth_values, $mutant_depth;
		    } else {
				modules::Exception->warning("Non-numeric characters in read depth column");
		    }
		    if ($pass_row->[$header_col_lookup{median_quality_score}] =~ /[0-9\.]/) {
				push @quality_values, $pass_row->[$header_col_lookup{median_quality_score}];
		    } else {
				modules::Exception->warning("Non-numeric characters in median base quality column");
		    }
		}
		
		$genes{$pass_row->[$header_col_lookup{$self->{gene_col_name}}]}++;
	}
	
	# print summary output
	
	print $SUMM "\n";
	print $SUMM "Total passed variants: " . scalar @pass_rows . "\n";
	print $SUMM "\tTotal homozygous: $homo_count_passed\n" if $homo_count_passed>0;
	print $SUMM "\tTotal heterozygous: $heter_count_passed\n" if $heter_count_passed>0;
	print $SUMM "\tTotal other: $other_count_passed\n" if $other_count_passed>0;
	print $SUMM "\tTotal deletions: $deletions_passed\n" if $deletions_passed>0;
	print $SUMM "\tTotal insertions: $insertions_passed\n" if $insertions_passed>0;
	
	print $SUMM "\nBreakdown of passed variants:\n\tDistinct genes: " .  scalar (keys %genes) . "\n";
	print $SUMM "\tExonic: $exonic_passed\n";
	print $SUMM "\tSplice-site: $ss_passed\n";
	print $SUMM "\tMedian Depth: " . modules::Utils->median(\@depth_values) . "\n" if @depth_values;
	print $SUMM "\tMedian Quality: " . modules::Utils->median(\@quality_values) . "\n" if @quality_values;
	

	
	print $SUMM "\tGenes with more than one variant: ";
	my $gene_str = 'NONE';
	foreach my $gene (keys %genes){
	    if ($genes{$gene} > 1){
	    	$gene_str .= ",$gene";
	    }
	}
	$gene_str =~ s/NONE,//;
	print $SUMM "($gene_str)\n";
	
	close($SUMM);
		
}

sub insert_SNV_rows {
	my $self = shift;
	my %args = @_;
	my @snv_rows = ();
	my @snv_annotations = ();
	
	my %header_col_lookup = %{$self->{header_lookup}};
	my @headers = @{$self->{headers}};
	my %snv_lookup = %{$self->{snv_lookup}};
	
		
	for my $line ( @{$self->{all_rows}}) {
		#my @cols = split("\t",$line);
		my $chr = $line->[$header_col_lookup{chr}];
		my $coord = $line->[$header_col_lookup{coord}];
		my $var_base = $line->[$header_col_lookup{var_allele}];
		
		if (!exists $snv_lookup{$chr}{$coord}{$var_base}) {
			modules::Exception->throw("ERROR: Can't get snv id for $chr:$coord $var_base");
		}
		
		my $snv_id = $snv_lookup{$chr}{$coord}{$var_base};
		
		my %snv_row = (
						snv_id => $snv_id,
						run_id => $self->run->id
						);
						
#		for my $db_column (@db_columns) {
#			my $column_name = $db_column->{column_name};
#			next if $column_name eq 'snv_id' || $column_name eq 'run_id' || $column_name eq 'id';
#			$snv_row{$column_name} = $line->[$header_col_lookup{$column_name}];
#		}

		push @snv_rows, \%snv_row;
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'snv_rows', -data=>\@snv_rows);
	
	#Get the first column for annotations; snv_exon_type is the last in the snv_row table
	#my $first_col_count = $header_col_lookup{snv_exon_type} + 1;
	my $last_col_count = keys %header_col_lookup;
	$last_col_count--;
	
	my %snv_row_lookup = ();
	
	my $snv_row_iterator = modules::Adaptors::SNV_Row->search(run_id => $self->run->id);	
	
	#Get all the ids for the inserted records
	while (my $snvrow_db = $snv_row_iterator->next){
		my ($snv_obj) = modules::Adaptors::SNV->search(id => $snvrow_db->snv_id);
		$snv_row_lookup{$snv_obj->chr}{$snv_obj->coord}{$snv_obj->var_base} = $snvrow_db->id; 
	}

	
	#Now we need to add the snv_annotations
	for my $line ( @{$self->{all_rows}}) {
		my $chr = $line->[$header_col_lookup{chr}];
		my $coord = $line->[$header_col_lookup{coord}];
		my $var_base = $line->[$header_col_lookup{var_allele}];
		
		if (!exists $snv_row_lookup{$chr}{$coord}{$var_base}) {
			modules::Exception->throw("ERROR: Can't find snv_row entry for $chr:$coord:$var_base");
		}
		
		#Get the snv_row_id using the snv_lookup 
		my $snv_row_id = $snv_row_lookup{$chr}{$coord}{$var_base};

		#This is range of columns not in snv_row; these are added as annotations
		for my $col_count (0..$last_col_count) {
			my $column_value = $line->[$col_count];
			$column_value =~ s/'/"/g if defined $column_value;
			if ($snv_row_id !~ /\d/) {
				print "ERROR: for $line\n";
			}
			my %snv_annotation = (
									column_name => $headers[$col_count],
									column_value => $column_value,
									column_number => $col_count,
									snv_row_id => $snv_row_id
								  );
			push @snv_annotations, \%snv_annotation;
		}
	}
	#print Dumper \@snv_annotations;
	modules::Adaptors::BulkInsert->insert(-table_name=>'snvrow_values', -data=>\@snv_annotations);
	
	
}

sub insert_Indel_rows {
	my $self = shift;
	my %args = @_;
	my @indel_rows = ();
	my @indel_annotations = ();
	
	my %header_col_lookup = %{$self->{header_lookup}};
	my @headers = @{$self->{headers}};
	my %indel_lookup = %{$self->{indel_lookup}};
	
	#print Dumper \%snv_lookup;
	
	for my $line ( @{$self->{all_rows}}) {
		#my @cols = split("\t",$line);
		my @cols = @{$line};
		my $chr = $cols[$header_col_lookup{chr}];
		my $start_coord = $cols[$header_col_lookup{start_coord}];
		my $end_coord = $cols[$header_col_lookup{end_coord}];
		my $var_base = $cols[$header_col_lookup{var_allele}];
		my $var_type = $cols[$header_col_lookup{var_type}];
		
		if (!exists $indel_lookup{$chr}{$start_coord}{$end_coord}{$var_type}{$var_base}) {
			modules::Exception->throw("ERROR: Can't get indel id for $chr:$start_coord $var_base");
		}
		
		my $indel_id = $indel_lookup{$chr}{$start_coord}{$end_coord}{$var_type}{$var_base};
		
		my %indel_row = (
						variant_id => $indel_id,
						run_id => $self->run->id
						);
		
#		for my $db_column (@db_columns) {
#			my $column_name = $db_column->{column_name};
#			next if $column_name eq 'variant_id' || $column_name eq 'run_id' || $column_name eq 'id';
#			$indel_row{$column_name} = $cols[$header_col_lookup{$column_name}];
#		}
		
		push @indel_rows, \%indel_row;
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'variant_rows', -data=>\@indel_rows);
	
	#Get the first column for annotations; snv_exon_type is the last in the snv_row table
	my $first_col_count = 0;
	my $last_col_count = keys %header_col_lookup;
	$last_col_count--;
	
	my %indel_row_lookup = ();
	
	my $indel_row_iterator = modules::Adaptors::Variant_Row->search(run_id => $self->run->id);	
	
	#Get all the ids for the inserted records
	while (my $indelrow_db = $indel_row_iterator->next){
		my ($var_db) = modules::Adaptors::Variant->search(id => $indelrow_db->variant_id);
		my $bases;
		if ($var_db->var_type eq 'DEL') {
			$bases = 'N/A';
		} elsif ($var_db->var_type eq 'INS') {
			$bases = '+'. $var_db->inserted_bases;
		} else {
			modules::Exception->throw("ERROR: Unknown variant type");
		}
		$indel_row_lookup{$var_db->chr}{$var_db->start_coord}{$var_db->end_coord}{$var_db->var_type}{$bases} = $indelrow_db->id; 
	}
	
	#Now we need to add the indel_annotations
	for my $line ( @{$self->{all_rows}}) {
		my $chr = $line->[$header_col_lookup{chr}];
		my $start_coord = $line->[$header_col_lookup{start_coord}];
		my $end_coord = $line->[$header_col_lookup{end_coord}];
		my $var_base = $line->[$header_col_lookup{var_allele}];
		my $indel_type = $line->[$header_col_lookup{var_type}];
		
		if (!exists $indel_row_lookup{$chr}{$start_coord}{$end_coord}{$indel_type}{$var_base}) {
			modules::Exception->throw("ERROR: Can't get indel row id for $chr:$start_coord $var_base");
		}
		
		#Get the snv_row_id using the snv_lookup 
		my $indel_row_id = $indel_row_lookup{$chr}{$start_coord}{$end_coord}{$indel_type}{$var_base};

		#This is range of columns not in snv_row; these are added as annotations
		for my $col_count ($first_col_count..$last_col_count) {
			my $column_value = $line->[$col_count];
			$column_value =~ s/'/"/g if defined $column_value;
			my %indel_annotation = (
									column_name => $headers[$col_count],
									column_value => $column_value,
									column_number => $col_count,
									variant_row_id => $indel_row_id
								  );
			push @indel_annotations, \%indel_annotation;
		}
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'variantrow_values', -data=>\@indel_annotations);
	
	
}

return 1;
