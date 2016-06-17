package modules::Vcf;

use strict;
use Data::Dumper;
use modules::Exception;
use modules::Pipeline;

sub new {
    my ($class, @args) = @_;

	my @required_args = (
			             -sample_type
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }


    my $self = bless {}, $class;
    
	$self->sample_type($args{-sample_type});
    
    return $self;
}

#Get the runid
sub sample_type {
    my ($self, $sample_type) = @_;

    if (defined $sample_type) {
		$self->{'sample_type'} = $sample_type;
    } elsif (! defined $self->{'sample_type'}) {
		modules::Exception->throw("sample_type not set");
    }

    return $self->{'sample_type'};
}

#Parse a vcf file
sub parse_vcf {
    my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    open(VCF,$args{-vcf_file}) || modules::Exception->throw("Can't open file $args{-vcf_file}\n");
    my %variant_data = ();
    
    while (<VCF>) {
    	next if /^#/;
    	
    	if ($_ !~ /^[0-9XYM]/) {
    		modules::Exception->throw("ERROR: Error with VCF format at line $_");
    	}
    	
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,$normal_gt,$tumour_gt) = split("\t");
    	if (!$qual) {
    		modules::Exception->throw("ERROR: Error with VCF format at line $_");
    	}

		my @normal_fields = split(':',$normal_gt);
    	my @tumour_fields = split(':',$tumour_gt) if $tumour_gt;
		
		my $var_type;
		if (@tumour_fields) {
			if ($normal_fields[0] eq '0/0') {
				$var_type =  'SOMATIC';
			} elsif ($normal_fields[0] eq $tumour_fields[0]) {
				$var_type = 'GERMLINE';
			} elsif ($normal_fields[0] eq '0/1' && ($tumour_fields[0] eq '0/0' || $tumour_fields[0] eq '1/1')) {
				$var_type = 'LOH';
			} else {
				$var_type = 'COMPLEX';
			}
		} else {
			#Can't tell at this point when it's not matched
			$var_type = 'UNKNOWN';
		}
		
		
    	#Get the CLR score
    	my $clr_str = 'NO_CLR';
    	if (/CLR=(\d+)/) {
    		
    		if ($normal_fields[0] eq '0/0') {
    			#Check it's the tumour that is variant
    			$clr_str = 'REF_NORM'.$1;
    		} elsif ($tumour_fields[0] eq '0/0') {
	    		#If it's the tumour that matches reference
    			$clr_str = 'REF_TUM'.$1;
    		} else {
				$clr_str = 'REF_NONE'.$1;				
    		}
    	} 
    	
		my @vars = split(",",$var_str);
		
		for my $var ( @vars ) {
			my $variant_type;
			my $start_coord;
			my $end_coord;
			my $bases;
			my $length_ref = length($ref);
			my $length_var = length($var);
			
			if (length($ref) > length($var)) {
				$variant_type = 'DEL';
				$start_coord = $first_coord + 1;
				$end_coord = $start_coord + $length_ref - $length_var - 1;				
				my $del_length = $length_ref - $length_var;
				$bases = substr($ref,1,$del_length);
			} elsif (length($ref) < length($var)) {
				$variant_type = 'INS';
				#Add the ref length and var length difference to the coord
				$start_coord = $end_coord = $first_coord + 1;
				my $ins_length = $length_var - $length_ref;
				$bases = substr($var,1,$ins_length);
			} else {
				$variant_type = 'SNV';
				$start_coord = $end_coord = $first_coord;
				$bases = $ref . '->' .$var;
			}
			my $key = 'Q'.$qual . '^^^'.$clr_str .'^^^'. $var_type;
			$variant_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $key;
		}
    }
    
    $self->{data}{$args{-vcf_file}} = \%variant_data;
}

#Get vcf data for a file
sub get_vcf {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $self->{data}{$args{-vcf_file}} doesn't exist");	
    }
    return $self->{data}{$args{-vcf_file}};
}

#Apply filters to vcf data and return passed lines only
sub filter_vcf {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file,
			             -snv_depth_file,
			             -indel_depth_file,
			             -tumour_flag
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $self->{data}{$args{-vcf_file}} doesn't exist");	
    }
    
    if ( !-e $args{-snv_depth_file} ) {
    	modules::Exception->throw("File $args{-snv_depth_file} doesn't exist");	
    }
    
    if ( !-e $args{-indel_depth_file} ) {
    	modules::Exception->throw("File $args{-indel_depth_file} doesn't exist");	
    }
    
    my $tumour_flag = $args{-tumour_flag};
    
    my %indel_depth = ();
    my %snv_depth = ();
    
    #Get the depth from the merge_vcf files
    open(SNV,"$args{-snv_depth_file}") || modules::Exception->throw("Can't open file $args{-snv_depth_file}\n");
    while (<SNV>) {
    	my @fields = split("\t");
    	my ($depth) = $fields[3] =~ /D(\d+)$/;
    	$snv_depth{$fields[0]}{$fields[1]} = $depth;
    }
    
    open(INDEL,"$args{-indel_depth_file}") || modules::Exception->throw("Can't open file $args{-indel_depth_file}\n");
    while (<INDEL>) {
    	my @fields = split("\t");
    	my ($depth) = $fields[3] =~ /D(\d+)$/;
    	$indel_depth{$fields[0]}{$fields[1]} = $depth;
    }
    
    
    my $pipe_config = modules::Pipeline->get_pipe_conf();
    my $report_config = modules::Pipeline->get_report_conf();
    
    my $min_snv_quality = $pipe_config->read('cutoffs','snv_quality_cutoff');
	my $min_indel_quality = $pipe_config->read('cutoffs','indel_quality_cutoff');
	my $min_clr = $pipe_config->read('cutoffs','clr_cutoff');
	my $min_depth = $pipe_config->read('cutoffs','min_variant_cutoff');
	my $max_depth = $pipe_config->read('cutoffs','max_variant_cutoff');
	my $sample_type = $self->{sample_type};
	my %snv_types = ();
	if ($report_config->exists('snv','sample_types',$sample_type,'db_types')) {
		%snv_types = map {$_ => 1} split(",",$report_config->read('snv','sample_types',$sample_type,'db_types'));
	} else {
		%snv_types = map {$_ => 1} split(",",$report_config->read('snv','common','db_types'));	
	}
	
	
	my %filter_vcf_data = ();
	my $all_vcf_data = $self->get_vcf(-vcf_file => $args{-vcf_file});
	
	if (exists $all_vcf_data->{SNV}) {
		for my $chr (sort keys %{$all_vcf_data->{SNV}}) {
			for my $start_coord (sort {$a<=>$b} keys %{$all_vcf_data->{SNV}{$chr}}) {
				for my $end_coord (keys %{$all_vcf_data->{SNV}{$chr}{$start_coord}}) {
					for my $bases (keys %{$all_vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}}) {
						my $rest = $all_vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}{$bases};
						my $depth = $snv_depth{$chr}{$start_coord};
						my $new_rest = $rest . '^^^D' . $depth;
						#Here we need to consider CLR score as well
						my ($quality_str,$clr_str,$type) = split('\^\^\^',$rest);
						
						if (!exists $snv_types{$type} && !exists $snv_types{'ALL'}) {
							#Only pass SNVs we're interested in
							next;
						}
						my ($quality) = $quality_str =~ /(\d+)/;

						if ($tumour_flag) {
							
							#Now filter CLR score
							next if $clr_str eq 'NO_CLR';
							
							my ($clr_score) = $clr_str =~ /(\d+)/;
													
							if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_snv_quality && $clr_score >= $min_clr ) {
								$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
							}
						} else {
							if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_snv_quality) {
								$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
							}
						}
					}
				}
			}
		}
	}
	
	
	my %indel_types = ();
	if ($report_config->exists('indel','sample_types',$sample_type,'db_types')) {
		%indel_types = map {$_ => 1} split(",",$report_config->read('indel','sample_types',$sample_type,'db_types'));
	} else {
		%indel_types = map {$_ => 1} split(",",$report_config->read('indel','common','db_types'));	
	}
	my @indel_type = qw(DEL INS);
	
	for my $variant_type (@indel_type) {
		if (exists $all_vcf_data->{$variant_type}) {
			for my $chr (sort keys %{$all_vcf_data->{$variant_type}}) {
				for my $start_coord (sort {$a<=>$b} keys %{$all_vcf_data->{$variant_type}{$chr}}) {
					for my $end_coord (keys %{$all_vcf_data->{$variant_type}{$chr}{$start_coord}}) {
						for my $bases (keys %{$all_vcf_data->{$variant_type}{$chr}{$start_coord}{$end_coord}}) {
							my $rest = $all_vcf_data->{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases};
							my $depth = $indel_depth{$chr}{$start_coord};
							my $new_rest = $rest . '^^^D' . $depth;
							#Here we need to consider CLR score as well
							my ($quality_str,$clr_str,$type) = split('\^\^\^',$rest);
					
							if (!exists $indel_types{$type} && !exists $snv_types{'ALL'}) {
								#Only pass SNVs we're interested in
								next;
							}
							my ($quality) = $quality_str =~ /(\d+)/;

							if ($tumour_flag) {
								
								#Now filter CLR score
								next if $clr_str eq 'NO_CLR';
								
								my ($clr_score) = $clr_str =~ /(\d+)/;
										
								if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_indel_quality && $clr_score >= $min_clr ) {
									$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
								}
							} else {
								if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_indel_quality) {
									$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
								}
							}
						}
					}
				}
			}
		}
	}
	
	return \%filter_vcf_data;
}

#Check that a vcf file is 'complete'
sub check_vcf {
	my ($vcf,$chr) = @_;
	my %chr_length = (
						1 => 249250621,
						2 => 243199373, 	
						3 => 198022430, 
						4 => 191154276, 
						5 => 180915260, 
						6 => 171115067, 
						7 => 159138663, 
						8 => 146364022, 
						9 => 141213431, 
						10 => 135534747, 
						11 => 135006516, 
						12 => 133851895, 
						13 => 115169878, 
						14 => 107349540, 
						15 => 102531392, 
						16 => 90354753, 
						17 => 81195210, 
						18 => 78077248, 
						19 => 59128983, 
						20 => 63025520, 
						21 => 48129895, 
						22 => 51304566, 
						'X' => 155270560, 
						'Y' => 59373566,
						);
	
	open(VCF,$vcf) || modules::Exception->throw("ERROR: Can't open vcf file $vcf to check");
	my @vcf_lines = <VCF>;
	if (!@vcf_lines) {
		print STDERR "WARNING: Vcf has no lines and chr $chr length is $chr_length{$chr}\n";
		return 0;
	}

	my $line_count = @vcf_lines;
	my (undef,$last_coord) = split("\t",$vcf_lines[-1]);
	

	my $end_proportion = $last_coord / $chr_length{$chr};
	my $snv_proportion = $line_count / $chr_length{$chr};
	
	#If the ratio of lines to chr_length is < 1:10000 then we likely have a problem
	if ($snv_proportion < 0.0001) {
		print STDERR "WARNING: Vcf has $line_count lines and chr length $chr is $chr_length{$chr} ($snv_proportion < 0.0001)\n";
		return 0;
	}
	
	#If we're not at least 99% of the way to the end of the chromosome then we likely have a problem
	if ($end_proportion < 0.99) {
		print STDERR "WARNING: Vcf has last coord $last_coord and chr $chr length is $chr_length{$chr} ($end_proportion < 0.99)\n";
		return 0;
	}
	
	

	return 1;
}


return 1;