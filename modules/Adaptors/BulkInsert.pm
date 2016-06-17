package modules::Adaptors::BulkInsert;

use strict;
use DBI;
use modules::Adaptors::MelDB;
use modules::Exception;
use Data::Dumper;

sub insert 
{
    my $self = shift;
    my %args = @_;
    if (!defined $args{-table_name} || !defined $args{-data}) {
         modules::Exception->throw("insert requires arguments -table_name and -data arrayref");
    }
    my $ra_row_data = $args{-data};
    my $table_name = $args{-table_name};
    my @fields = keys %{$ra_row_data->[0]};

    my $query = "INSERT INTO $table_name (" . join(",", @fields) . ")";
 
    my @rows;
    foreach my $rh_row (@{$ra_row_data}){
                my @row_data;
                foreach my $field (@fields){
                    push @row_data, defined $rh_row->{$field} ? $rh_row->{$field} : 'NULL';
		    		$row_data[-1] =~ s/\"/\\\"/g; # because pileup format allows ^' and ^" in the base string (argh!)
                }
                push @rows, "\'" . join("\',\'", @row_data) . "\'";
    }

    $query .= ' VALUES (' . join('),(', @rows) . ')';
    my $db = defined $args{-db}?$args{-db}:'MelDB';
    my $sth = $self->dbh(-db=>$db)->prepare($query);
    #print "Q $query\n";
    $sth->execute
        or modules::Exception->throw("Insert failed: $DBI::errstr");

    return 1;
}

sub dbh {
	my $self = shift;
	my %args = @_;
    if (!defined $args{-db}) {
         modules::Exception->throw("dbh requires arguments -db");
    }
    
    my $db = $args{-db};
    
    if ($args{-db} eq 'MelDB') {
    	return DBI->connect(modules::Adaptors::MelDB->getConfig) || modules::Exception->throw("Unable to connect to database: $DBI::errstr");
    } else {
    	modules::Exception->throw("ERROR: db $db isn't a known database");
    }
    
    
}

return 1;
