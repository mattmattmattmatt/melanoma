package modules::Log;

use strict;
use modules::Exception;
use Data::Dumper;

sub new {
    my ($class, $logfile,$append) = @_;

    my $self = bless {}, $class;
    
    $self->{'logfile'} = $logfile;
    if ($append && -e $logfile) {
    	open(FILE,">>","$logfile") || modules::Exception->throw("Can't open log file for writing; Check permission for logfile directory");
    } else {
		open(FILE,">","$logfile") || modules::Exception->throw("Can't open log file for writing; Check permission for logfile directory");
    }
	$self->{'fh'} = \*FILE;

    return $self;
}

sub append {
	my ($self, $message) = @_;
	my $fh = $self->{'fh'};
	chomp $message;
	print $fh "$message\n";
}

sub close {
	my ($self) = @_;
	close $self->{'fh'};
}

1;