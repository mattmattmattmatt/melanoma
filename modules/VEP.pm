package modules::VEP;

use strict;
use modules::Exception;
use modules::SystemCall;

#./variant_effect_predictor/variant_effect_predictor.pl -i vep.in --poly b --sift b --canonical --cache --dir ~/Desktop/ --offline -o STDOUT --coding_only

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
		$self->{'executable_path'} = 'variant_effect_predictor.pl';
    }

    return $self->{'executable_path'};
}

sub organism {
    my ($self, $organism) = @_;

    if (defined $organism){
		$self->{'organism'} = $organism;
    } elsif (! defined $self->{'organism'}) {
		$self->{'organism'} = 'human';
    }

    return $self->{'organism'};
}

sub db_dir {
    my ($self, $db_dir) = @_;

    if (defined $db_dir){
		$self->{'db_dir'} = $db_dir;
    } elsif (! defined $self->{'db_dir'}) {
		modules::Exception->warn("variant_predictor database directory not set");
    }

    return $self->{'db_dir'};
}

sub working_dir {
    my ($self, $working_dir) = @_;

    if (defined $working_dir){
		$self->{'working_dir'} = $working_dir;
    } elsif (! defined $self->{'working_dir'}) {
		modules::Exception->warn("variant_predictor working directory not set");
    }

    return $self->{'working_dir'};
}

sub output_file {
    my ($self) = @_;

    my @path = split /\//,  $self->input_file;

    return $self->{'working_dir'} . '/vep.out';
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

	#Only run polyphen/sift for human
	my $args = $self->organism eq 'human'?' --force_overwrite --poly b --sift b --canonical --cache --offline --coding_only --no_intergenic':' --force_overwrite --canonical --cache --offline --coding_only --no_intergenic';


    my $command 
	= $self->executable_path 
	#Get sift and polyphen scores, canonical state, use the cache file(no db), and only find coding exon entries
	. $args
	. ' --species ' . $self->organism
	. ' --dir ' . $self->db_dir
	. ' --o ' . $self->output_file
	. ' --i ' . $self->input_file;

    print STDERR "Running command: $command\n";

    if (system($command)){
		modules::Exception->throw("Command: $command\nError: Non-zero exit status.")
    }

    return 1;
}


sub parse_result {
    my ($self) = @_;

    open(my $OUTPUT, $self->output_file) || modules::Exception->throw("ERROR: Can't open output file $self->output_file");

	my %grouping_data = ();

    while (<$OUTPUT>){
		chomp;
		#Skip headers
		next if $_ =~ /^#/;
	
		my ($identifier, $coord_str, $ref_base, $ens_gene, $ens_transcript, undef, $aa_type, undef, undef, undef, $aa_change, $codon_change, $strand, $attribute_str ) = split /\t/;
		
		my $aa_ref = my $aa_var;
		
		if ($aa_change =~ /([A-Z\*])\/([A-Z\*])/){
			$aa_ref = $1;
			$aa_var = $2;
		} else {
		    next;
		}
		
		#Only record results from canonical transcripts and skip synonomous mutants
		next unless $attribute_str =~ /CANONICAL/; 
		my $chr = my $start = my $end;
		if ($coord_str =~ /([0-9XYMT]+):(\d+)\-(\d+)/) {
			$chr = $1;
			$start = $2;
			$end = $3;
		} else {
			($chr,$start) = $coord_str =~ /([0-9XYMT]+):(\d+)/;
			$end = $start;
		}
		
		my $polyphen_pred = "N/A";
		my $polyphen_score = "N/A";
		my $sift_pred = "N/A";	
		my $sift_score = "N/A";
	
		my @attribute_pairs = split /;/, $attribute_str;
		#Sample attribute line: 
		#PolyPhen=possibly_damaging(0.593);CANONICAL=YES;SIFT=deleterious(0)
		for my $attribute_pair ( @attribute_pairs ) {
		    next if $attribute_pair =~ /CANONICAL/;
		    
		    if ($attribute_pair =~ /PolyPhen/) {
		    	($polyphen_pred, $polyphen_score) = $attribute_pair =~ /PolyPhen=(.+)\(([0-9\.]+)\)/;
		    } elsif ($attribute_pair =~ /SIFT/) {
		    	($sift_pred, $sift_score) = $attribute_pair =~ /SIFT=(.+)\(([0-9\.]+)\)/;
		    }
		}
	
		
	
	
		$aa_ref = 'Stop' if ($aa_ref eq '*');
		$aa_var = 'Stop' if ($aa_var eq '*');
	
		#Just keep one copy of redundant entries
		$grouping_data{$chr}{$start}{$end}{"$aa_ref->$aa_var"} = [$aa_type, $chr, $start, $end, $ens_gene, $ens_transcript, $aa_ref, $aa_var, $polyphen_pred, $polyphen_score, $sift_pred, $sift_score];
	
		#push @parsed_result, [$aa_type, $chr, $start, $ens_gene, $ens_transcript, $aa_ref, $aa_var, $polyphen_pred, $polyphen_score, $sift_pred, $sift_score];
    }

    close($OUTPUT);

	my @parsed_result;

	for my $chr ( keys %grouping_data ) {
	    for my $start_coord ( keys %{$grouping_data{$chr}} ) {
	    	for my $end_coord (keys %{$grouping_data{$chr}{$start_coord}}) {
	    		for my $aa_change (keys %{$grouping_data{$chr}{$start_coord}{$end_coord}}) {
	    			push @parsed_result, $grouping_data{$chr}{$start_coord}{$end_coord}{$aa_change};
	    		}
	    	}
		}
	}
    return \@parsed_result;
}

return 1;
