package modules::PileupLine;

use strict;
use modules::BaseFrequencies;
use Data::Dumper;

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

    $self->chr($args{'chr'})                       if defined $args{'chr'};
    $self->coord($args{'coord'})                   if defined $args{'coord'};
    $self->ref_base($args{'ref_base'})             if defined $args{'ref_base'};
    $self->var_base($args{'var_base'})             if defined $args{'var_base'};
    $self->base_string($args{'base_string'})       if defined $args{'base_string'};
    $self->quality_string($args{'quality_string'}) if defined $args{'quality_string'};

    return $self;
}

sub chr {
    my ($self, $chr) = @_;

   if (defined $chr) {
       $chr =~ s/chr//i;
       $chr = uc($chr);
       $self->{'chr'} = $chr;
   } else {
       return $self->{'chr'}
   }
}

sub coord {
    my ($self, $coord) = @_;

    if (defined $coord) {
		$self->{'coord'} = $coord;
    } else {
		return $self->{'coord'}
    }
}

sub ref_base {
    my ($self, $ref_base) = @_;

    if (defined $ref_base) {
		$self->{'ref_base'} = uc($ref_base);
    } else {
		return $self->{'ref_base'}
    }
}

sub var_base {
    my ($self, $var_base) = @_;

    if (defined $var_base) {
		$self->{'var_base'} = uc($var_base);
    } else {
		return $self->{'var_base'}
    }
}

sub is_variant {
    my $self = shift;

    if ($self->ref_base eq $self->var_base){
		return 0;
    } else {
		return 1;
    }
}

sub base_string {
    my ($self, $base_string) = @_;

    if (defined $base_string) {
		my $ref_base = $self->ref_base;
		$base_string =~ s/[\.\,]/$ref_base/g;
		$self->{'base_string'} = $base_string;
    } else {
		return $self->{'base_string'}
    }
}

sub base_string_indel {
    my ($self, $base_string) = @_;

    if (defined $base_string) {
		$self->{'base_string'} = $base_string;
    } else {
		return $self->{'base_string'}
    }
}

sub quality_string {
    my ($self, $quality_string) = @_;

    if (defined $quality_string) {
	$self->{'quality_string'} = $quality_string;
    } else {
	return $self->{'quality_string'}
    }
}

sub get_base_array {
    my ($self) = @_;

    my $base_string = $self->base_string;
    #Remove quality strings
    $base_string =~ s/\^.//g;
	
    my @bases = split("",uc($base_string));
    my $count = 0;
    my $base_count = @bases;
    my @final_bases = ();
    
    #Skip over indels entries   
   	while ($count < $base_count) {
	    if (defined $bases[$count+2] && ($bases[$count+1] eq '+' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
	    	#Skip all the indel bases from the array
	    	my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		$count +=  $length_indel+2+length($length_indel);
	    	next;
	    } elsif (defined $bases[$count+2] && ($bases[$count+1] eq '-' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
	    	my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		$count +=  $length_indel+2+length($length_indel);
	    	next;
	    } elsif ($bases[$count] =~ /[ATCGN]/) {
	    	push @final_bases, $bases[$count];
	    } elsif ($bases[$count] eq '*') {
	    	push @final_bases, 'OTHER_DEL';
	    }
	    $count++;
	}
    
    return \@final_bases;
}

#Here we read the indel patterns; for snvs we skip these entries
sub get_base_array_indel {
    my ($self) = @_;

    my $base_string = $self->base_string;
    #Remove quality strings
    $base_string =~ s/\^.//g;
	
    my @bases = split("",uc($base_string));
    my $count = 0;
    my $base_count = @bases;
    my @final_bases = ();
    
    
    #Skip over indels entries   
   	while ($count < $base_count) {
     	if (defined $bases[$count+2] && ($bases[$count+1] eq '+' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
    		my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		my $indel_end = $length_indel+1+length($length_indel)+$count;
    		my $indel_bases = join("",@bases[$count+1..$indel_end]);
    		$count +=  $length_indel+2+length($length_indel);
    		push @final_bases, $indel_bases;
    		next;
    	} elsif (defined $bases[$count+2] && ($bases[$count+1] eq '-' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
    		my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		my $indel_end = $length_indel+1+length($length_indel)+$count;
    		my $indel_bases = join("",@bases[$count+1..$indel_end]);
    		$count +=  $length_indel+2+length($length_indel);
    		push @final_bases, $indel_bases;
    		next;
    		
   		} elsif ($bases[$count] =~ /[ATCGN]/) {
	    	push @final_bases, $bases[$count];
	    } elsif ($bases[$count] eq '.' || $bases[$count] eq ',') {
	    	push @final_bases, 'REF';
	    } elsif ($bases[$count] eq '*') {
	    	push @final_bases, 'OTHER_DEL';
	    }
	    $count++;
	}
    
    return \@final_bases;
}

#Now storing this in the database instead
#sub read_depth {
#    my $self = shift;
#	my @bases =  split("",uc($self->base_string));
#	my $base_count = scalar @bases;
#	my $read_depth = 0;
#	my $count = 0;
#	
#	while ($count < $base_count) {
#		if (defined $bases[$count+2] && ($bases[$count+1] eq '+' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
#	    	#Skip all the indel bases from the array
#	    	$read_depth++;
#	    	$count += 3+$1;
#	    	next;
#	    } elsif (defined $bases[$count+2] && ($bases[$count+1] eq '-' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
#	    	$read_depth++;
#	    	$count += 3+$1;
#	    	next;
#	    } elsif ($bases[$count] =~ /[ATCGN]/) {
#			$read_depth++;
#	    }
#	    $count++;
#	}
#		
#	return $read_depth;
#    #return scalar @{$self->get_base_array};
#}

sub base_frequencies {
    my $self = shift;

    if (! defined $self->{'base_frequencies'}){
		$self->{'base_frequencies'} = modules::BaseFrequencies->new('pileup_line' => $self);
    }

    return $self->{'base_frequencies'};
}

sub base_frequencies_indel {
    my $self = shift;

    if (! defined $self->{'base_frequencies'}){
		$self->{'base_frequencies'} = modules::BaseFrequencies->new('pileup_line' => $self, 'indel' => 1);
    }

    return $self->{'base_frequencies'};
}

sub predominant_base {
    my $self = shift;

    if (defined $self->base_frequencies->bases->[0]){
		return $self->base_frequencies->bases->[0];
    } else {
		return 0;
    }
}

sub predominant_base_freq {
    my $self = shift;

    if (defined $self->base_frequencies->frequencies->[0]){
		return $self->base_frequencies->frequencies->[0];
    } else {
		return 0;
    }
}

sub next_most_common_base {
    my $self = shift;

    if (defined $self->base_frequencies->bases->[1]){
		return $self->base_frequencies->bases->[1];
    } else {
		return 0;
    }
}

sub next_most_common_base_freq {
    my $self = shift;

    if (defined $self->base_frequencies->frequencies->[1]){
		return $self->base_frequencies->frequencies->[1];
    } else {
		return 0;
    }
}

sub print {
    my $self = shift;

    print join("\t",
	       $self->chr,
	       $self->coord,
	       $self->ref_base,
	       $self->var_base,
	       $self->is_variant,
	       $self->read_depth,
	       $self->predominant_base,
	       $self->predominant_base_freq,
	       $self->next_most_common_base,
	       $self->next_most_common_base_freq,
	       $self->base_string,
	       $self->quality_string
	) . "\n";
}

return 1;
