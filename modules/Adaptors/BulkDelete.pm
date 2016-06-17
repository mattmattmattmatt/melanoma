package modules::Adaptors::BulkDelete;

use strict;
use DBI;
use modules::Adaptors::MelDB;
use modules::Exception;

#Remove db entries in bulk given a table name and an array ref of primary keys
sub delete {
	my $self = shift;
	my %args = @_;
	if (!defined $args{-table_name} || !defined $args{-key_ids}) {
		modules::Exception->throw("delete requires arguments -table_name and -key_ids arrayref");
	}
	my $table_name = $args{-table_name};
	my $keys = $args{-key_ids}; 
	#Assumes the primary key has the format table_name_id
	my $pk_name = 'id';

    my $query = "Delete from $table_name where $pk_name IN (" . join(",", @{$keys}) . ")";
 	#print "$query\n";
 
    my $sth = $self->dbh->prepare($query);
    $sth->execute || modules::Exception->throw("Delete failed for $table_name: $DBI::errstr");

    return 1;
}

#Custom delete designed to remove variant_filter entries
sub delete_indel_filter {
	my $self = shift;
	my %args = @_;
	
	if (!defined $args{-run_id} || (!defined $args{-filter_id} && !defined $args{-filter_name})) {
		modules::Exception->throw("delete requires arguments -run_id and (-filter_id or -filter_name)");
	}
	
	my $run_id = $args{-run_id};
    my $filter_id;
    my $filter_name;
    if (defined $args{-filter_id}) {
		$filter_id= $args{-filter_id};
    } else {
    	$filter_name = $args{-filter_name};
    	my $filterid_sth = $self->dbh->prepare("SELECT id FROM filter WHERE name = ?");
		$filterid_sth->execute($filter_name);

		if ($filterid_sth->rows == 1){
		    ($filter_id) = $filterid_sth->fetchrow_array;
		} else {
		    modules::Exception->throw("Fetching id for filter [$filter_name] returned " 
			. $filterid_sth->rows .  " rows");
		}
    }
		
	# Search for snv ids directly with DBI
	
	my $indelid_sth = $self->dbh->prepare("SELECT id FROM variants WHERE run_id = ?");
	
	$indelid_sth->execute($run_id);
	
	my @indel_ids;
	
	while (my ($indel_id) = $indelid_sth->fetchrow_array){
	    push @indel_ids, $indel_id;
	}
	
	# Delete snv_filter rows
	if (@indel_ids) {
		my $sf_delete_sth = $self->dbh->prepare('DELETE FROM variants_filters WHERE variant_id IN (' . join(',', @indel_ids) . ') AND filter_id = ' . $filter_id . ";\n");
		$sf_delete_sth->execute || modules::Exception->throw("Delete failed for indels: $DBI::errstr");;
	}
}

#Custom delete designed to remove snv_filter entries
sub delete_snv_filter {
	my $self = shift;
	my %args = @_;
	
	if (!defined $args{-run_id} || (!defined $args{-filter_id} && !defined $args{-filter_name})) {
		modules::Exception->throw("delete requires arguments -run_id and (-filter_id or -filter_name)");
	}
	
	my $run_id = $args{-run_id};
    my $filter_id;
    my $filter_name;
    if (defined $args{-filter_id}) {
		$filter_id= $args{-filter_id};
    } else {
    	$filter_name = $args{-filter_name};
    	my $filterid_sth = $self->dbh->prepare("SELECT id FROM filter WHERE name = ?");
		$filterid_sth->execute($filter_name);

		if ($filterid_sth->rows == 1){
		    ($filter_id) = $filterid_sth->fetchrow_array;
		} else {
		    modules::Exception->throw("Fetching id for filter [$filter_name] returned " 
			. $filterid_sth->rows .  " rows");
		}
    }
		
	# Search for snv ids directly with DBI
	
	my $snvid_sth = $self->dbh->prepare("SELECT id FROM snvs WHERE run_id = ?");
	
	$snvid_sth->execute($run_id);
	
	my @snv_ids;
	
	while (my ($snv_id) = $snvid_sth->fetchrow_array){
	    push @snv_ids, $snv_id;
	}
	
	# Delete snv_filter rows
	if (@snv_ids) {
		my $sf_delete_sth = $self->dbh->prepare('DELETE FROM snvs_filters WHERE snv_id IN (' . join(',', @snv_ids) . ') AND filter_id = ' . $filter_id . ";\n");
		$sf_delete_sth->execute || modules::Exception->throw("Delete failed for snvs: $DBI::errstr");;
	}
}

sub dbh {
    return DBI->connect(modules::Adaptors::MelDB->getConfig)
	or modules::Exception->throw("Unable to connect to database: $DBI::errstr");
}

return 1;
