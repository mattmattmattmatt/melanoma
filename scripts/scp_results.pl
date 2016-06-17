#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use modules::Adaptors::Release_File;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"files=s",
	   	"runid=i",
	   	"release_files",
	   	"workdir=s",
	   	"patient_id=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{files});

	   
=pod

=head1 SYNOPSIS

scp_results.pl -files files_to_copy -runid run_id [options]

Required flags: -files

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

scp_results.pl -> Wrapper to copy tarballs to destination directory either on local or remote host

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

scp_results.pl 

=cut

my $cluster_config = modules::Pipeline->get_cluster_conf();
my $archive_command = $cluster_config->read('archive','archive_command');


my $sys_call = modules::SystemCall->new();


my $copy_flag = my $scp_flag = 0;
if ($archive_command eq 'cp') {
	$copy_flag = 1;
} else {
	$scp_flag = 1;
}

my @files = ();
push @files, $OPT{files};


if ($OPT{release_files}) {
	if (!$OPT{workdir}) {
		modules::Exception->throw("ERROR: Need workdir arg with release_files arg");
	}
	my $workdir = $OPT{workdir};
	my @release_files_obj = modules::Adaptors::Release_File->search_release_files_run($OPT{runid});
	for my $release_file_obj (@release_files_obj) {
		my $release_file_name = $release_file_obj->file_name;
		next if $release_file_name =~ /bam/; #These are already handled in copy_bam step
		if ($release_file_name =~ /vcf/) {
        	push @files, "$workdir/vcf/$release_file_name";
      	} elsif ($release_file_name =~ /summary/) {
        	push @files, "$workdir/summary/$release_file_name";
      	} else {
        	push @files, "$workdir/$release_file_name";
    	}
	}
}

my $files = join(" ",@files);

if ($OPT{runid}) {
	$files =~ s/RUN/$OPT{runid}/g;
}

#Error check if we're scp'ing or cp'ing
if ($scp_flag || $copy_flag) {
	if ( $files =~ /.gz$/ && !-e $files ) {
		modules::Exception->throw("File $files doesn't exist");	
	}	
	my $destdir = $cluster_config->read('base_directories','base_results_directory') . '/'. $OPT{patient_id};

	#This directory should have been created during copy_bam step
	if ($copy_flag && !-d $destdir) {
		modules::Exception->throw("Need to define the source directory for output");
	} 
	if ($scp_flag && $cluster_config->empty('archive','file_server')) {
		modules::Exception->throw("Need to define the fileServer directory for scp'ing");
	}
	my $copy_command;
	
	
	if ($copy_flag) {
		$copy_command = "cp $files $destdir";			
	} elsif ($scp_flag) {
		my $user = `whoami`;
		chomp $user;		
		my $scp_server = $cluster_config->empty('archive','file_server');	
		$copy_command = "scp $files $user\@$scp_server:$destdir";
	}
	#First we run the copy command
	$sys_call->run("$copy_command");
	#Then we remove the local tarball
	if ( $files =~ /.gz$/) {
		$sys_call->run("rm $files");
	}
}


