package modules::CLI_Options;
use Data::Dumper;
use strict;

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    return $self;
}

sub get_option {
	my $self = shift;
	my %args = @_;
    my @required_args = (
			             -name,
			             -default,
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
  
    my $error;
    my $variable = $args{-default};
    my $name = $args{-name};
    my $question = defined $args{-yes}?"Is it":"What is the";
    my @possibilities = ();
    if (defined $args{-possibilities}) {
		@possibilities = @{$args{-possibilities}};
    } elsif (defined $args{-yes}) {
		@possibilities = qw(yes no);
    }
    #If it's a retry then we die instead of trying again
    my $retry = 0;
    if (defined $args{-retry}) {
    	$retry = 1;
    }
    
    my $format_str = defined $args{-format}?$args{-format}:'';
    my $regex = defined $args{-regex}?$args{-regex}:'';
    my $not_empty = defined $args{-not_empty}?1:0;
    
    #If passed on command line
    if (defined $args{-opt_arg}) {
		$variable = $args{-opt_args};
	} else {
		#Otherwise enter the value for STDIN
		if ($variable) {
			print "$question $name $format_str [default=$variable]?";
		} else {
			print "$question $name $format_str [default=NONE]?";			
		}
		my $user_input = <STDIN>;
		chomp $user_input;
		if ($user_input) {
			$variable = $user_input;
		} 
		if (!$variable && $not_empty) {
			$error = "ERROR: $name cannot be empty";
		}
	}
	
	if ($args{-directory}) {
		if ( !-d $variable ) {
			$error = "ERROR: Directory $variable is not a directory";	
		}
	}
	
	
	if (@possibilities) {
		my $acceptable_value = 0;
		for my $acceptable ( @possibilities ) {
		  	if ($acceptable eq $variable) {
		  		$acceptable_value = 1;
		  	}
		}
		if (!$acceptable_value) {
			my $acceptable_str = join(" OR ",@possibilities);
			$error = "ERROR for input: $name must be $acceptable_str";
		}
	}
	
	if ($regex) {
		my @entries = split(",",$variable);
		for my $entry ( @entries ) {
		    if ($entry !~ /^$regex$/) {
				$error = "$entry does not match acceptable format ($regex)";
		    }
		}
	}
	
	#If there is any error
	if ($retry && $error =~ /\w/) {
		modules::Exception->throw("ERROR: $error\n");
	} elsif ($error =~ /\w/) {
		return ($variable,1,$error);
	} else {
		return ($variable,0,$error);
	}
	
	
	
}

1;