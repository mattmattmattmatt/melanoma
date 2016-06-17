#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::Adaptors::Run;
use modules::Adaptors::Sample;
use modules::Adaptors::Pipeline_Step;
use modules::Adaptors::Pipeline_Step_Run;
use modules::Adaptors::Syscall;
use modules::Adaptors::BulkDelete;
use modules::Vcf;
use modules::Pipeline;
use File::Basename;
use Pod::Usage;
use Cwd;
use Cwd 'abs_path';
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"chr=s",
		"vcf=s",
		"runid=i",
		"step_name=s",
		"qsub_file=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf} || !$OPT{chr} || !$OPT{runid} || !$OPT{qsub_file} || !$OPT{step_name});

	   
=pod

=head1 SYNOPSIS

record_snv_call.pl -step_name step_name_to_record -runid runid -chr chr -vcf vcf_file -qsub_file qsub_file 

Required flags: -chr -runid -qsub_file -vcf -step_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

call_snvs.pl -> Submits all the snv calling jobs for each chromosome

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

call_snvs.pl 

=cut

my $vcf_file = $OPT{vcf};
my $chr = $OPT{chr};
my $vcf_step_name = $OPT{step_name};
my $runid = $OPT{runid};
my $qsub_file = $OPT{qsub_file};

if (!-e $qsub_file) {
	modules::Exception->throw("ERROR: Qsub file $qsub_file doesn't exist");
}

my ($step_obj) = modules::Adaptors::Pipeline_Step->search('name' => $vcf_step_name);

if (!defined $step_obj) {
	modules::Exception->throw("ERROR: No step db object for $vcf_step_name");
}
my $step_id = $step_obj->id;


my ($run_obj) = modules::Adaptors::Run->search('id' => $runid);

if (!defined $run_obj) {
	modules::Exception->throw("ERROR: No run db object with runid $runid");
}

if (! -e $vcf_file) {
	modules::Exception->throw("ERROR: vcf file $vcf_file doesn't exist");
} elsif (-s $vcf_file == 0) {
	modules::Exception->throw("ERROR: vcf file $vcf_file is empty");
}

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

#Get the xml files to get variables we need first
my $cluster_config = modules::Pipeline->get_cluster_conf();
my $pipe_config = modules::Pipeline->get_pipe_conf();


#Check if the vcf ran properly
if (!modules::Vcf::check_vcf($vcf_file,$chr)) {
	modules::Exception->throw("ERROR: $vcf_file fails for chr $chr\n") unless $chr eq 'Y';
}

#Check if the pipeline step exists first; get all the steps for all chromosomes first 
my @singlevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$runid,pipeline_step_id=>$step_id);
								
#Get the syscalls for the single_vcfs
my %chr_done = ();
						
for my $singlevcf_step ( @singlevcf_steps ) {
	my ($syscall_obj) = modules::Adaptors::Syscall->search('pipeline_steps_run_id'=>$singlevcf_step->id);
   	my $pipe_step_id = $singlevcf_step->id;
    if (!defined $syscall_obj) {
    	#Delete orphaned pipeline_steps_runs entries (no syscall)
    	modules::Adaptors::Pipeline_Step_Run->search(id=>$pipe_step_id)->delete_all;
    } else {
		push @{$chr_done{$syscall_obj->chr}},$pipe_step_id;	
    }
}


#Delete previous steps if they exist
if (exists $chr_done{$chr}) {
	modules::Adaptors::BulkDelete->delete(-table_name=>'pipeline_steps_runs',-key_ids=>\@{$chr_done{$chr}});
}

#Get the commands run from the vcf
my $commands = modules::Pipeline::get_commands($qsub_file);
my $command = join(" ; ",@{$commands});

#Finally insert the pipeline_steps_runs into the db
my $pipeline_step_id = modules::Pipeline::Insert_Step_Command($vcf_step_name,$command,$runid,$chr);
print STDERR "Insert pipeline step $pipeline_step_id for chr $chr\n";







