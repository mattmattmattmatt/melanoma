package modules::QualityEncoding;
use modules::Exception;
use vars qw(%SIG);
use strict;

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    return $self;
}

sub encoding {
	my %args = @_;
	if (!-e $args{-readfile}) {
		modules::Exception->throw("Problem with argument -readfile $args{-readfile}");
	}
	my $reads_to_check = defined $args{-reads_to_check}?$args{-reads_to_check}:10000;

	my $encoding = 'phred64';
	$SIG{PIPE} = '_IGNORE_'; #Doesn't give the annoying zcat closed pipe error
	#Parse the first reads
	my $error_line_count = 4;
	my $line_num = 1;
	my $reads_checked = 0;
	my $fh;
	
	if ($args{-readfile} =~ /.gz$/) {
		open($fh,"gunzip -c $args{-readfile} |") || modules::Exception->throw("Cannot open read file $args{-readfile}");
	} else {
		open($fh,$args{-readfile}) || modules::Exception->throw("Cannot open read file $args{-readfile}");
	}
	
	if (eof $fh) {
		modules::Exception->throw("ERROR: Empty filehandle for $args{-readfile}");
	}
	
	while(<$fh>) {
		if ($reads_checked == $reads_to_check) {
			last;
		}
		if ($line_num % $error_line_count == 0) {
			chomp;
			my @quals = split("",$_);
			for my $qual ( @quals ) {
			    my $phred = ord($qual);
			    #Check the range that only exists for phred33
				if ($phred > 33 && $phred < 64) {
					#close $fh || die "Can't close $fh";
					return 'phred33';
				}
			}
			$reads_checked++;
		}
		$line_num++;
	}
	
	
	close $fh;
	return $encoding;
}

1;