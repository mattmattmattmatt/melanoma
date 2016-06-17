#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::Adaptors::Align_Lane;
use modules::Adaptors::Lane;
use modules::Adaptors::Syscall;
use modules::Adaptors::Lane_Run;
use modules::Adaptors::Sample;
use modules::ConfigXML;
use modules::Pipeline;
use File::Basename;
use Pod::Usage;
use Cwd;
use Cwd 'abs_path';
use vars qw(%OPT);

GetOptions(\%OPT, 
		   	"help|h",
		   	"man|m",
		   	"submit"
		   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

commence_run.pl -submit submit_jobs

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

commence_run.pl -> Create runs when individual lanes are aligned

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

commence_run.pl 

=cut  

#Keep the svndir 
my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

#Get the xml files and create the pipeline object
my $cluster_xml = "$svndir/conf/cluster.xml";

if ( !-e $cluster_xml ) {
	modules::Exception->throw("File $cluster_xml doesn't exist");	
}

my $cluster_config = modules::ConfigXML->new($cluster_xml);

my @commands = ();

#First retrieve all the samples
my @samples_obj = modules::Adaptors::Sample->search_all();
my $unaligned_lanes = 0;


for my $sample_obj ( @samples_obj ) {
    my $sample_name = $sample_obj->sample_name;
    
    #Get the total lanes variable
    my $total_lanes = $sample_obj->total_lanes;

	#Get the number of lanes db entries and ensure they match    
    my @lanes_objs = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);
    
    if (@lanes_objs != $total_lanes) {
    	my $lane_obj_count = @lanes_objs;
    	modules::Exception->throw("ERROR: Sample $sample_name has 'total_lanes' of $total_lanes and $lane_obj_count db lanes entries; These should be the same");
    }
    
    #Check if each lane has a lane_align entry
    my $align_lane_total = 0;
    my @unaligned_lanes = ();
    for my $lane_obj ( @lanes_objs ) {
        my $lane_name = $lane_obj->lane_name;
        my @align_lane = modules::Adaptors::Align_Lane->search(lane_id=>$lane_obj->id);
        if (@align_lane > 1) {
        	modules::Exception->throw("ERROR: Lane number $lane_name has more than one align_lane entry");
        } elsif (@align_lane == 1) {
        	$align_lane_total++;
        } else {
        	push @unaligned_lanes, $lane_name;
        }
    }
    
    if ($align_lane_total == $total_lanes) {
	    #Check if all these lanes are part of existing run
	    #We need to trigger new run if there exists at least one lane not used in previous run (either completely new or new lane added)  
	    my $run_lane_total = 0;
	    for my $lane_obj ( @lanes_objs ) {
			my @lane_run_objs = modules::Adaptors::Lane_Run->search(lane_id=>$lane_obj->id);
			#TODO: Fix for multiple runs
			if (@lane_run_objs) {
				$run_lane_total++;
			}
		}
		
		#If we need to create a new run
		if ($run_lane_total != $total_lanes) {
			print STDERR "Begin new run for sample $sample_name\n";
			my $qsub_base = $cluster_config->read('base_directories','base_qsub_directory');
			my ($patient_name) = modules::Pipeline::get_patient_id($sample_name);
			my $qsub_dir = $qsub_base . '/' . $patient_name . '/runs/';
			my $qsub_wrapper = $qsub_dir. $sample_name . '.pipe1.wrapper.qsub';
			if (!-e $qsub_wrapper) {
				modules::Exception->throw("ERROR: No qsub file for $sample_name (file = $qsub_wrapper)");
			}
			push @commands, "cd $qsub_dir; bash $qsub_wrapper";
			
		} else {
			#print STDERR "Skip sample $sample_name: Run with these lanes already exists\n";
		}
	
	
    	
    } else {
    	my $lane_count = @unaligned_lanes;
    	my $lane_str = join(", ",@unaligned_lanes);
    	print STDERR "Skip sample $sample_name: There are $lane_count lanes unaligned: ($lane_str)\n";
    	$unaligned_lanes += $lane_count;
    }
    
	
	    
    
}

print STDERR "\nRun the following commands:\n\n" if @commands;
for my $command ( @commands ) {
   print STDERR "$command\n";
   system("$command") if $OPT{submit};
}

print "\nTotal unaligned lanes: $unaligned_lanes\n";










