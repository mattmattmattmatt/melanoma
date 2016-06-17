package modules::SetOperations;

use strict;
use modules::Overlap;
use Math::Combinatorics;
use modules::Exception;
use Set::Scalar;
use Data::Dumper;

sub new {
    my ($class,%args) = @_;
    my $self = bless {}, $class;	
    
    #Variable determines which columns will be used for making the overlap key (default = 0 and 1)
    my @key_columns;
    if (defined $args{-key_columns}) {
    	if ($args{-key_columns} =~ /[^0-9,]/) {
    		modules::Exception->throw("Argument -key_columns requires a comma delimited list of integers");
    	} else {
    		$self->{key_columns} = $args{-key_columns};
    	}
    } else {
   		$self->{key_columns} = '0,1';
    }
    
    #Set the delimiter if not the default (\t)
	if (defined $args{-delim}) {
    	$self->{delim} = $args{-delim};
    } else {
   		$self->{delim} = "\t";
    }
    
    return $self;
}

#Load the data into sets
sub load_sets {
	my $self = shift;
	my (%args) = @_;
	
	if (!defined $args{-files}) {
		modules::Exception->throw("Missing required argument -files, a comma separated list of input files");
	} else {
		my @files = @{$args{-files}};
		my @col_counts = split(",",$self->{key_columns});

		my $set_count = 0;
		for my $file ( @files ) {
		    if ( !-e $file ) {
		    	modules::Exception->throw("Input file $file doesn't exist");
		    }
		    #Generate the input set from the file
		    my %set = ();
		    open(FILE,$file) || modules::Exception->throw("Can't open file $file");
		    while (<FILE>) {
		    	chomp;
		    	my @columns = split($self->{delim},$_);
		    	#Check the line has enough columns to make a key
		    	if (!@columns || @columns < $col_counts[-1]) {
		    		next;
		    	}
		    	#Build up the key; can't use array slice because treats as key_columns as string
		    	my $key;
		    	for my $col_count ( sort {$a<=>$b} @col_counts ) {
		    	    $key .= ':'.$columns[$col_count];
		    	}
		    	$key =~ s/^://;
		    	$set{$key}++;
		    }
		    close FILE;
		    my $set_obj = Set::Scalar->new(keys %set);
		    $self->{sets}{$set_count} = $set_obj;
		    $set_count++;
		}
	}
}
  
 
#Master set_operation function takes the input files and the operation to perform
sub set_operation {
	my $self = shift;
	my (%args) = @_;
	
	my @singular_required_args = qw(-intersect -union -min_cutoff -max_cutoff);
	my $required_arg_count = 0;
	
	for my $required_arg ( @singular_required_args ) {
	    if (exists $args{$required_arg}) {
	    	$required_arg_count++;
	    }
	}
	
	my $temp_set;
	if ($required_arg_count > 1) {
		modules::Exception->throw("Problem with arguments: Requires one of ".join(" ",@singular_required_args));
	} elsif ($required_arg_count == 0) {
		#default is intersect all the existing sets
		$temp_set = $self->_intersect_sets(-all=>1);
	} elsif (exists $args{-intersect}) {
		if (exists $args {-setlist}) {
			$temp_set = $self->_intersect_sets(-setlist=>$args{-setlist});		
		} else {
			$temp_set = $self->_intersect_sets(-all=>1);
		}
	} elsif (exists $args{-union}) {
		if (exists $args {-setlist}) {
			$temp_set = $self->_union_sets(-setlist=>$args{-setlist});		
		} else {
			$temp_set = $self->_union_sets(-all=>1);
		}
	} elsif (exists $args{-min_cutoff}) {
		#If we want to run on a subset of the input sets
		if (exists $args{-setlist}) {
			$temp_set = $self->_group_overlap(-cutoff=>$args{-min_cutoff},-setlist=>$args{-setlist},-min=>1);		
		} else {
			#default is to run on all the loaded sets
			$temp_set = $self->_group_overlap(-cutoff=>$args{-min_cutoff},-all=>1,-min=>1)		
		}
	} elsif (exists $args{-max_cutoff}) {
		#To make it <= we need to add one; this is a result of sharing a function 
		my $cutoff = $args{-max_cutoff} + 1;
		#If we want to run on a subset of the input sets
		if (exists $args{-setlist}) {
			$temp_set = $self->_group_overlap(-cutoff=>$cutoff,-setlist=>$args{-setlist},-max=>1);		
		} else {
			#default is to run on all the loaded sets
			$temp_set = $self->_group_overlap(-cutoff=>$cutoff,-all=>1,-max=>1);		
		}
	} else {
		#Shouldn't get here
		modules::Exception->throw("Problem with arguments");
	}
	
	my $count;
	my $return_set;
	
	my $temp_count = $temp_set->size;
	#print "Before $temp_count\n";
	
	#Can be done after any previous set operations
	if (exists $args{-excludelist}) {
		#Need to exclude these sets
		$return_set = $self->_remove_sets(-excludelist=>$args{-excludelist},-current_set=>$temp_set);
	} else {
		$return_set = $temp_set;
	}
	$count = $return_set->size;		
	
	return ($count,$return_set);
	
}

#Intersect sets...
sub _intersect_sets {
	my $self = shift;
	my (%args) = @_;
	
	my $return_set; 

	if (exists $args{-all}) {
		#Here we find the intersection of all the sets
		$return_set = $self->{sets}{0};
		for my $set_count (keys %{$self->{sets}}) {
			next if $set_count == 0;
			$return_set = $return_set * $self->{sets}{$set_count};	
		}		
	} else {
		#Here we find the intersection of some subset of the sets
		my @set_numbers = split(",",$args{-setlist});
		$return_set = $self->{sets}{$set_numbers[0]};
		shift @set_numbers;
		for my $set_num ( @set_numbers ) {
		    $return_set = $return_set * $self->{sets}{$set_num};
		}
	}
	
	return $return_set;
}

#Union sets...
sub _union_sets {
	my $self = shift;
	my (%args) = @_;
	my $return_set; 

	if (exists $args{-all}) {
		#Here we find the union of all the sets
		$return_set = $self->{sets}{0};
		for my $set_count (keys %{$self->{sets}}) {
			next if $set_count == 0;
			$return_set = $return_set + $self->{sets}{$set_count};	
		}		
	} else {
		#Here we find the union of some subset of the sets
		my @set_numbers = split(",",$args{-setlist});
		$return_set = $self->{sets}{$set_numbers[0]};
		shift @set_numbers;
		for my $set_num ( @set_numbers ) {
		    $return_set = $return_set + $self->{sets}{$set_num};
		}
	}
	
	return $return_set;
}

#Remove sets...
sub _remove_sets {
	my $self = shift;
	my (%args) = @_;
	my $return_set = $args{-current_set}; 
	
	my @set_numbers = split(",",$args{-excludelist});
	for my $set_num ( @set_numbers ) {
	    $return_set = $return_set - $self->{sets}{$set_num};
	}
	
	return $return_set;
}

#Find all entries in <= N or >= N number of sets (eg in >= 7 of 10 sets OR in <= 7 of 10 sets)
sub _group_overlap {
	my $self = shift;
	my (%args) = @_;
	my $cutoff = $args{-cutoff};
	my $min = exists $args{-min}?1:0;
		
	my @setlist = ();
	if (exists $args{-setlist}) {
		@setlist = split(",",$args{-setlist});
	} else {
		my $max_count = keys %{$self->{sets}};
		@setlist = (0 .. $max_count-1);
	}
	
	my $return_set = Set::Scalar->new();
	my $combinator = $self->_get_combos(-data=>\@setlist,-count=>$cutoff);
	while (my @combo = sort {$a<=>$b} $combinator->next_combination) {
		my $temp_combo_set = $self->{sets}{$combo[0]};
		shift @combo;
		#Build up an intersection set for the particular combo
		for my $combo_entry ( @combo ) {
		    $temp_combo_set = $temp_combo_set * $self->{sets}{$combo_entry};
		}
		#Add the set of it has enough entries; use union to account for duplicates
		$return_set = $return_set + $temp_combo_set;
	}
	
	if ($min) {
		return $return_set;	
	} else {
		#Here we get the total union set and then subtract the frequently occurring entries
		my $total_set;
		if (exists $args{-setlist}) {
			$total_set =  $self->_union_sets(-setlist=>$args{-setlist});
		} else {
			$total_set = $self->_union_sets(-all=>1);
		}
		my $max_return_set  = $total_set - $return_set;
		return $max_return_set;
	}
	
}

#use combinatorics to get the combos when the min_overlap or max_overlap flag is passed
sub _get_combos {
	my $self = shift;
	my (%args) = @_;
	my $data = $args{-data};
	my $count = $args{-count};
	my $combinator = Math::Combinatorics->new(count=>$count,data=>$data);
	return $combinator;
}

#Write a set to a file; useful for writing intermediate sets from iterative queries
sub set_to_file {
	my $self = shift;
	my (%args) = @_;
	my $set = $args{-set};
	my $file = $args{-file};
	my %existing_entries = ();
	if (defined $args{-append}) {
		#To avoid duplicate entries load existing file to prevent writing duplicates
		if (-e $file) {
			open(TMP,"$file");
			while (<TMP>) {
				chomp;
				$existing_entries{$_} = 1;
			}
			close TMP;
		}
		open(FILE,">>$file") || modules::Exception->throw("Can't open file $file\n");
	} else {
		open(FILE,">$file") || modules::Exception->throw("Can't open file $file\n");
	}
	
	for my $element (sort $set->elements ) {
		(my $print_element = $element) =~ s/:/$self->{delim}/g;
		if (!exists $existing_entries{$print_element}) {
		    print FILE "$print_element\n";
		} 
	}
	
	close FILE;
	
}

1;