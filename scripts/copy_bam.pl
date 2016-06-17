#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"bam=s",
	   	"patient_id=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{bam} || !$OPT{patient_id});

	   
=pod

=head1 SYNOPSIS

copy_bam.pl -bam bam_file -patient_id patient_id(needed_for_destdir) [options]

Required flags: -bam -patient_id

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

copy_bam.pl -> Add lanes to existing sample

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

copy_bam.pl 

=cut

my $bam_file = $OPT{bam};
my $bam_index = $bam_file . '.bai';

if (!-e $bam_file) {
	modules::Exception->throw("ERROR: File $bam_file doesn't exist");
}

if (!-e $bam_index) {
	modules::Exception->throw("ERROR: File $bam_index doesn't exist");
}

my $cluster_config = modules::Pipeline->get_cluster_conf();
my $destdir = $cluster_config->read('base_directories','base_results_directory') . '/' . $OPT{patient_id};

#Add the patient_id directory if required
mkdir($destdir) if !-d $destdir;

my $command = "cp $bam_file $bam_index $destdir";
print STDERR "$command\n";
my $sys_call = modules::SystemCall->new();
$sys_call->run($command);










