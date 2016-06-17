package modules::BaseFrequencies;

use strict;

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

    $self->pileup_line($args{'pileup_line'}) if defined $args{'pileup_line'};

	if ($args{'indel'}) {
		$self->{indel} = 1;
	}

    return $self;
}


sub pileup_line {
    my ($self, $pl) = @_;

    if (defined $pl) {

		if ($pl->isa("modules::PileupLine")){
		    $self->{'pileup_line'} = $pl;
		} else {
		    warn "Not a modules::PileupLine";
		    return 0;
		}
    } elsif (! defined $self->{'pileup_line'}) {
		warn "modules::PileupLine object is not set";
    }

    return $self->{'pileup_line'};
}


sub _compute {
    my $self = shift;

    die "modules::PileupLine object is not set"
	unless defined $self->pileup_line;

    my @bases = ();
    if ($self->{indel}) {
    	@bases = @{$self->pileup_line->get_base_array_indel};
    } else {
    	@bases = @{$self->pileup_line->get_base_array};
    }
    my %base_counts = ();

    foreach my $base (@bases){
		$base_counts{uc($base)}++;
    }

    my @sorted_bases = sort {$base_counts{$b} <=> $base_counts{$a}} keys %base_counts;

    $self->{'bases'} = \@sorted_bases;

    my %base_frequencies;
    for (my $i = 0; scalar @bases && $i < scalar @sorted_bases; $i++) {
		my $base = $sorted_bases[$i];
		$self->{'frequencies'}->[$i] = $base_counts{$base} / scalar @bases; 
		$self->{'hash_lookup'}->{$base} = $base_counts{$base} / scalar @bases;
		$self->{'counts'}->[$i] = $base_counts{$base};
	 }

    if (! defined $self->{'counts'}) { # if no base counts, set some null values and fail gracefully
		$self->{'hash_lookup'} = \%base_counts;
		$self->{'bases'} = ['N'];
		$self->{'counts'} = [0];
		$self->{'frequencies'} = [0];
    }

    return 1;
}

sub lookup {
    my $self = shift;

    $self->_compute if ! defined $self->{'hash_lookup'};

    return $self->{'hash_lookup'};
}

sub bases {
    my $self = shift;

    $self->_compute if ! defined $self->{'bases'};

    return $self->{'bases'};
}

sub frequencies {
    my $self = shift;

    $self->_compute if ! defined $self->{'frequencies'};

    return $self->{'frequencies'};
}

sub counts {
    my $self = shift;

    $self->_compute if ! defined $self->{'counts'};

    return $self->{'counts'};
}


return 1;
