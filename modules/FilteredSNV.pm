package modules::FilteredSNV;

use strict;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Exception;
use Data::Dumper;


sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

    my @required_args = ('snv', 'filter_list_order');

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

    $self->snv($args{'snv'});
    $self->filter_list_order($args{'filter_list_order'});

    return $self;
}

sub snv {
    my ($self, $snv) = @_;

    if (defined $snv) {
		$self->{'snv'} = $snv;
    }

    return $self->{'snv'};
}

sub filter_list_order {
    my ($self, $filter_list_order) = @_;

    if (defined $filter_list_order){
		if (defined $filter_list_order->[0] 
		    && $filter_list_order->[0]->isa("modules::Adaptors::Filter")) {
		    $self->{'filter_list_order'} = $filter_list_order;
		} else {
		    modules::Exception->throw("Filter list does not contain modules::DBAdaptor::Filter");
		}
    }

    return $self->{'filter_list_order'};
}

sub filter_ids {
    my $self = shift;

    if (! defined $self->{'filter_ids'}){
		my @filter_ids;
		foreach my $filter (@{$self->filter_list_order}){
		    push @filter_ids, $filter->id;
		}
		$self->{'filter_ids'} = \@filter_ids;
    }

    return $self->{'filter_ids'};
}

sub filter_names {
    my $self = shift;

    if (! defined $self->{'filter_names'}){
		my @filter_names;
		foreach my $filter (@{$self->filter_list_order}){
		    push @filter_names, $filter->name;
		}
		$self->{'filter_names'} = \@filter_names;
    }

    return $self->{'filter_names'};
}

sub filter_ids_vs_names {
    my $self = shift;

    if (! defined $self->{'filter_id_vs_names'}){

		foreach my $filter (@{$self->filter_list_order}){
		    $self->{'filter_id_vs_names'}->{$filter->id} = $filter->name;
		}
    }

    return $self->{'filter_id_vs_names'};
}

sub snv_filter_hash {
    my $self = shift;

    if (! defined $self->{'snv_filter_hash'}) {
		my @snv_filters = modules::Adaptors::SNV_Filter->search('snv_id' => $self->snv->id);
		foreach my $snv_filter (@snv_filters){
		    
		    next unless defined $self->filter_ids_vs_names->{$snv_filter->filter_id};
		    my $filter_name = $self->filter_ids_vs_names->{$snv_filter->filter_id};
		    $self->{'snv_filter_hash'}->{$filter_name} = $snv_filter;
		}

	
		foreach my $filter_name (@{$self->filter_names}){
		    if (! defined $self->{'snv_filter_hash'}->{$filter_name}){
				$self->{'snv_filter_hash'}->{$filter_name} = 'nodata';
		    }
		}

    }

    return $self->{'snv_filter_hash'};
}

return 1;
