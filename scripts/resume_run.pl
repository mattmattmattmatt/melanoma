#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::Adaptors::Lane;
use modules::Adaptors::Syscall;
use modules::Adaptors::Sample;
use modules::Adaptors::Sample_Group;
use modules::Adaptors::Release_File;
use modules::Adaptors::Pipeline_Step;
use modules::ConfigXML;
use modules::Pipeline;
use modules::Cluster;
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

resume_run.pl -submit submit_jobs -user cluster_user(default=mxf221)

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

resume_run.pl -> Resume existing runs when we reach various states

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

resume_run.pl 

=cut

#Keep the svndir 
my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}


#If this step is run the sample run is complete and ready for snv_call
my $end_pipe1_step = 'bam_stats';
my $submit_snvs_step = 'submit_snvs';
my $merge_vcf_step  = 'merge_vcf';
my $single_vcf_step = 'single_vcf';
my $end_pipe3_step_tumour = 'scp_results';
my $merge_bam_step = 'merge_bam';

my %step_ids = ();

my @step_objs = modules::Adaptors::Pipeline_Step->search_all();

for my $step_obj (@step_objs) {
	my $step = $step_obj->name;
	$step_ids{$step} = $step_obj->id;
}

#Get the xml files and create the pipeline object
my $pipe_config = modules::Pipeline->get_pipe_conf();
my $cluster_config = modules::Pipeline->get_cluster_conf();
my %chromosomes = map {$_ => 1} split(" ",$pipe_config->read('annotation_version','chr'));

my @commands = ();

#First retrieve all the sample groups
my @sample_group_objs = modules::Adaptors::Sample_Group->search_all();

#Keep track on normal samples that don't have bams yet; these are required for tumour snv calling
my %missing_normals = ();

for my $sample_group_obj ( @sample_group_objs) {
	my $sample_total_field = $sample_group_obj->total_samples;
	my $sample_group_name = $sample_group_obj->sample_group_name;
	
	
	#Check if the normal sample is ready for snv calling
	my @normal_samples = modules::Adaptors::Sample->search(sample_group_id=>$sample_group_obj->id,tumour=>0);	
	
	if (!@normal_samples) {
		modules::Exception->throw("ERROR: Sample group $sample_group_name has no normal sample");
	} elsif (@normal_samples > 1) {
		modules::Exception->throw("ERROR: Sample group $sample_group_name has >1 normal samples");
	} 
		    
		    
	#my $normal_sample_id = $normal_samples[0]->id;
	
	my $normal_sample_name = $normal_samples[0]->sample_name;
	my $normal_total_lanes = $normal_samples[0]->total_lanes;
	my ($patient_name) = modules::Pipeline::get_patient_id($normal_sample_name);

	#Now check if the tumour sample is ready for snv calling 
	my @tumour_samples = modules::Adaptors::Sample->search(sample_group_id=>$sample_group_obj->id,tumour=>1);
	#The +1 is for the normal sample
	my $sample_total_db = 1 + @tumour_samples; 


	#Check the sample numbers agree
	if ($sample_total_field != $sample_total_db) {
		modules::Exception->throw("ERROR: Sample group $sample_group_name expects $sample_total_field samples and there are $sample_total_db samples in the database\n");
	}
	
	my $sample_ready_count = 0;
	my $complete = 0;
	
	
	#now check the runs are complete
	for my $sample ( @normal_samples,@tumour_samples ) {
		
		my $db_total_lanes;
		my $tumour_flag;
		#Tumour uses normal lanes as well
		if ($sample->tumour) {
			$db_total_lanes = $normal_total_lanes  + $sample->total_lanes;
			$tumour_flag = 1;
		} else {
			$db_total_lanes = $normal_total_lanes;
			$tumour_flag = 0;
		}
		
	    my $sample_id = $sample->id;
	    my $sample_name = $sample->sample_name;
	    #Get the latest run if any
	    my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_id);
	    
	    #Check if bam_stats has been run
	    if (defined $run_obj) {
	    	my $run_id = $run_obj->id;
	    	my @bamstats_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$end_pipe1_step});	
	    	if (@bamstats_steps == 1) {
	    		#Don't proceed past here with normal samples...
	    		if (!$tumour_flag) {
		    		$complete++;
		    		next;				
	    		}
	    		
				#Now check if submit_snvs step has been run
				my @callsnv_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$submit_snvs_step});
				
				#If snvs not submitted yet
				if (!@callsnv_steps) {
					#Tumour snv calling requires the normal bam file so skip if it doesn't exist
					if (exists $missing_normals{$sample_group_name}) {
						print "\n";
						next;
					}
					print STDERR "Sample $sample_name can submit snvs with $db_total_lanes lanes..\n";
					my $qsub_dir = $cluster_config->read('base_directories','base_qsub_directory') . '/'. $patient_name .'/runs';
					my $qsub_file = $sample_name . '.pipe2.'.$submit_snvs_step.'.qsub'; 
					my $qsub_snvcalls = $cluster_config->read('base_directories','base_qsub_directory') . '/'. $patient_name .'/snvcalls';
					push @commands, "grep submit $qsub_dir/$qsub_file | bash; cd $qsub_snvcalls; for f in $sample_name*qsub; do qsub \$f; done";
						
				} else {
					#Here jobs have been submitted so check if there is a complete run
					my @mergevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$merge_vcf_step});
					
					
					if (@mergevcf_steps) {
						#Get the release file and check the total lane numbers
						my @total_lanes = modules::Adaptors::Release_File->search_total_lanes($run_id,$step_ids{$merge_vcf_step});
						if (! @total_lanes) {
							modules::Exception->throw("ERROR: Sample group $sample_group_name sample $sample_name has run $merge_vcf_step but no release file");
						} elsif (@total_lanes > 1) {
							modules::Exception->throw("ERROR: Sample group $sample_group_name sample $sample_name has multiple release files for $merge_vcf_step");
						}
						
						my $release_file_lane_total = $total_lanes[0]->total_lanes;
											
						#If the lane number is different then submit jobs 
						if ($release_file_lane_total != $db_total_lanes) {
							print STDERR "Sample $sample_name can submit snvs as it has different lane number ($db_total_lanes) than existing run ($release_file_lane_total) for run $run_id\n";
							my $qsub_dir = $cluster_config->read('base_directories','base_qsub_directory') . '/'. $patient_name .'/runs';
							my $qsub_file = $sample_name . '.pipe2.'.$submit_snvs_step.'.qsub'; 
							my $qsub_snvcalls = $cluster_config->read('base_directories','base_qsub_directory') . '/'. $patient_name .'/snvcalls';
							push @commands, "grep submit $qsub_dir/$qsub_file | bash; cd $qsub_snvcalls; for f in $sample_name*qsub; do qsub \$f; done";
						} else {
							#Check the last step is run
							my ($last_step_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$end_pipe3_step_tumour});
							if (defined $last_step_obj) {
								$complete++;
								#print STDERR "\tSample $sample_name with $db_total_lanes lanes has finished step $end_pipe3_step and is complete..\n";					
							} else {
								#Last step not complete; find the latest step run
								my $last_step_run = modules::Pipeline::get_latest_step($run_id);
								print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
							}
						}
					} else {
						#Check how many (if any) single_vcf jobs are done
						my @singlevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$single_vcf_step});
						
						#Get the syscalls for the single_vcfs
						my %chr_done = ();
						
						for my $singlevcf_step_obj ( @singlevcf_steps ) {
						    my ($syscall_obj) = modules::Adaptors::Syscall->search('pipeline_steps_run_id'=>$singlevcf_step_obj->id);
						    my $pipe_step_id = $singlevcf_step_obj->id;
						    
						    if (!defined $syscall_obj) {
    							#Delete orphaned pipeline_steps_runs entries (no syscall)
    							modules::Adaptors::Pipeline_Step_Run->search(id=>$pipe_step_id)->delete_all;
    						} else {
						    	$chr_done{$syscall_obj->chr} = 1 if defined $syscall_obj;
    						}
						}
						
						if (keys %chr_done == 24) {
							#If run for each chr we start running pipe3 from merge_vcf to the end unless job already started
							my $qsub_dir = $cluster_config->read('base_directories','base_qsub_directory') . '/'. $patient_name .'/runs';
							my $qsub_file = $sample_name . '.pipe3.wrapper.qsub';
							 
							#Check the job hasn't been started yet
							#TODO: Test for SGE
							my $scheduler = $cluster_config->read('scheduler');
							my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);		
							my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
							
							if ($job_id) {
								print STDERR "Sample $sample_name has completed all 24 $single_vcf_step but $job_id has started running..\n";
							} else {
								print STDERR "Sample $sample_name has completed all 24 $single_vcf_step so begin pipe3 ..\n";											
								push @commands, "cd $qsub_dir; bash $qsub_dir/$qsub_file";
							}
							
						} else {
							my @missing_chr = ();
							for my $chr (keys %chromosomes) {
								if (!exists $chr_done{$chr}) {
									push @missing_chr, $chr;
								}
							}
							my $chr_missing = join(",",sort @missing_chr);
							my $missing_count = @missing_chr;
							print STDERR "Sample $sample_name has submitted jobs but $missing_count chromosome are not complete yet ($chr_missing)..\n";	
															
						}	
					} 
				}
	    	} elsif (@bamstats_steps == 0) {
	    		$missing_normals{$sample_group_name}++ if $tumour_flag == 0;
	    		my $last_step_run = modules::Pipeline::get_latest_step($run_id);
				print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
	    	} else {
	    		modules::Exception->throw("ERROR: Sample $sample_name has run $end_pipe1_step multiple times");
	    	}
	    } else {
	    	$missing_normals{$sample_group_name}++ if $tumour_flag == 0;
	    	print STDERR "Sample group $sample_group_name sample $sample_name has no runs created..\n";
	    }
	}
	my $total_samples = @normal_samples + @tumour_samples;
    print "\n" unless $total_samples == $complete;
}

print STDERR "\nRun the following commands:\n\n" if @commands;
for my $command ( @commands ) {
   print STDERR "$command\n";
   system("$command") if $OPT{submit};
}











