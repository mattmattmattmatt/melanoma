package modules::AnnoVar;

use strict;
use modules::Exception;
use modules::SystemCall;

sub new {
    my ($class, @args) = @_;

    my %args = @args;

    my $self = bless \%args, $class;
    
    return $self;
}

sub executable_path {
    my ($self, $executable_path) = @_;
    if (defined $executable_path 
	&& -e $executable_path){
		$self->{'executable_path'} = $executable_path;
    } elsif (! defined $self->{'executable_path'}) {
		$self->{'executable_path'} = 'annotate_variation.pl';
    }

    return $self->{'executable_path'};
}

sub build_version {
    my ($self, $build_version) = @_;

    if (defined $build_version){
		$self->{'build_version'} = $build_version;
    } elsif (! defined $self->{'build_version'}) {
		modules::Exception->warn("Build version is not set");
    }

    return $self->{'build_version'};
}

sub dbtype {
    my ($self, $dbtype) = @_;

    if (defined $dbtype){
		$self->{'dbtype'} = $dbtype;
    } elsif (! defined $self->{'dbtype'}) {
		$self->{'dbtype'} = 'ccdsGene';
    }

    return $self->{'dbtype'};
}

sub db_dir {
    my ($self, $db_dir) = @_;

    if (defined $db_dir){
		$self->{'db_dir'} = $db_dir;
    } elsif (! defined $self->{'db_dir'}) {
		modules::Exception->warn("Annovar database directory not set");
    }

    return $self->{'db_dir'};
}

sub working_dir {
    my ($self, $working_dir) = @_;

    if (defined $working_dir){
		$self->{'working_dir'} = $working_dir;
    } elsif (! defined $self->{'working_dir'}) {
		modules::Exception->warn("Annovar working directory not set");
    }

    return $self->{'working_dir'};
}

sub outfile_stub {
    my ($self) = @_;

    my @path = split /\//,  $self->input_file;

    return $self->{'working_dir'} . '/' . $path[-1] . '_annovar_output';
}

sub input_file {
    my ($self, $input_file) = @_;

    if (defined $input_file){
		$self->{'input_file'} = $input_file;
    } elsif (! defined $self->{'input_file'}) {
		modules::Exception->warn("No input file set");
    }

    return $self->{'input_file'};
}

sub run {
    my ($self) = @_;

    my $command 
	= $self->executable_path 
	. ' -dbtype ' . $self->dbtype 
	. ' -buildver ' . $self->build_version 
	. ' -geneanno ' . $self->input_file
	. ' -outfile ' . $self->outfile_stub 
	. ' ' . $self->db_dir;

    print STDERR "Running command: $command\n";

    if (system($command)){
		modules::Exception->throw("Command: $command\nError: Non-zero exit status.")
    }

    return 1;
}

sub parse_result {
    my ($self) = @_;

    open(my $OUTPUT, $self->outfile_stub . '.exonic_variant_function');

    my @parsed_result;

    while (<$OUTPUT>){
		chomp;
	
		my (undef, $result, $attributes, $chr, $start, $end) = split /\t/;
	
		my (undef, $gene, $exon, undef, $aa_change) = split /\:/, $attributes;
	
		unless ($aa_change =~ /p\.(\D)\d+([^\,])\,/){
		    modules::Exception->throw("Unable to parse amino acid identities from file [$aa_change]");
		}
	
		my $aa_ref = $1;
		my $aa_var = $2;
	
		$aa_ref = 'Stop' if ($aa_ref eq 'X');
		$aa_var = 'Stop' if ($aa_var eq 'X');
	
		push @parsed_result, [$result, $chr, $start, $gene, $exon, $aa_ref, $aa_var];

    }

    close($OUTPUT);

    return \@parsed_result;
}

return 1;
