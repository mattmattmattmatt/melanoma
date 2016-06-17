#! /usr/bin/perl -w

use strict;
use modules::ConfigXML;
use modules::SystemCall;
use modules::Exception;
use modules::Log;
use modules::Pipeline;
use modules::Adaptors::Release_File;
use modules::Adaptors::Run;
use modules::Adaptors::Sample;
use modules::Adaptors::Pipeline_Step;
use modules::Adaptors::Pipeline_Step_Run;
use modules::Adaptors::SNV;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Filter;
use modules::Adaptors::BulkDelete;
use modules::Adaptors::Variant;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::SNV_Row;
use modules::Adaptors::Variant_Row;
use modules::Adaptors::Syscall;
use modules::Adaptors::Lane_Run;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Sys::Hostname;
use Benchmark;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "xml=s",
		   "resume=s",
		   "latest_sample=s",
		   "clean_only",
		   "last_step=s",
		   "no_cleanup_files",
		   "no_cleanup_db",
		   "create_run",
		   "no_run",
		   "email"
		   )  || modules::Exception->throw("Invalid command-line option for script\n");
	   	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{xml});

	   
=pod

=head1 SYNOPSIS

run_pipeline_steps.pl -xml xml_config_file -no_run do_everything_except_the_step -no_cleanup_db no_cleanup_db(except pipeline_steps_run AND syscalls) -no_cleanup_files no_cleanup_files -email email_when_done -latest_sample get_the_latest_run_for_this_sample -create_run just_create_run_and_quit -resume resume_at_specific_step(stepname_from_xml) -last_step last_step_to_run -clean_only only_clean_steps [options]

Required flags: -xml

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

run_pipeline_steps.pl -> Wrapper script for running the pipeline

=head1 DESCRIPTION

Feb 16, 2011

a script that runs the individual steps in the pipeline using an xml config file

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./run_pipeline_steps.pl -xml file.xml -last_step sort_bam
./run_pipeline_steps.pl -xml file.xml -resume filter_basic
./run_pipeline_steps.pl -xml file.xml -latest_sample A15_sg1_tumour2

=cut
my $xml_config_file = $OPT{xml};
my ($xml_full_path) = `readlink -f $xml_config_file`;
chomp $xml_full_path;
my $cleanup_files = defined $OPT{no_cleanup_files}?0:1;
my $cleanup_db = defined $OPT{no_cleanup_db}?0:1;
my $step_resume = defined $OPT{resume}?$OPT{resume}:0;
my $last_step = defined $OPT{last_step}?$OPT{last_step}:0;
my $email = defined $OPT{email}?1:0;
my $latestlib_flag = defined $OPT{latest_sample}?1:0;
my $finish_flag = 0; #Flag to tell whether the last step has been run

if (!-e $xml_config_file) {
	modules::Exception->throw("ERROR:  XML file $xml_config_file doesn't exist $!\n");
}

if (!exists $ENV{'SVNDIR'}) {
	print "ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n";
	exit;
}
my $svndir = $ENV{'SVNDIR'};

my $start = new Benchmark;
my $run_config = modules::ConfigXML->new($xml_config_file);
my $sys_call = modules::SystemCall->new();

my $pipe_config = modules::Pipeline->get_pipe_conf();
my $cluster_config = modules::Pipeline->get_cluster_conf();

#The base to the run directory (/home/matt/runs/)
my $rundir_base = $cluster_config->read('base_directories','base_run_directory');

#The sample name (eg A15_sg1_tumour2); get from config file or the latest_sample flag
my $sample_name;
if ($latestlib_flag) {
	$sample_name = $OPT{latest_sample};
} else {
	$sample_name = $run_config->read('sample_name');
}

my ($patient_name) = modules::Pipeline::get_patient_id($sample_name);

#Get the sample and confirm it exists
my ($sample) = modules::Adaptors::Sample->search('sample_name' => $sample_name);
if (!defined $sample) {
	modules::Exception->throw("Sample $sample_name not in database");
}

#Need the sample id a few places
my $sample_id = $sample->id;

#The run id (eg 12)
my $runid;

#The name of the run directory (eg A15_sg1_tumour2_12)
my $rundir_name;

#The full directory (eg /home/matt/runs/A15_sg1_tumour2_12)
my $workdir;

#The overlap directory (eg /home/matt/runs/A15_sg1_tumour2_12/overlap)
my $overlapdir;

#The bam dir
my $bamdir;

#The log dir
my $logdir;

#The summary dir
my $summarydir;

#The conf dir
my $confdir;

#The polyphen dir
my $polydir;

#The vcf dir
my $vcfdir;

my %steps_run = ();
my %steps_to_clean = ();

my $steps = $run_config->read('steps_order', 'step');
my @steps;
my %step_mapping = ();

#Get the steps to run	
if (ref($steps) eq 'ARRAY'){ # Cope with steps being a single step or an array of steps
    @steps = @$steps;
} else {
    @steps = ($steps);
}

#Get the step order to check last_step is after resume_step
my $step_count = 0;
for my $step ( @steps ) {
    $step_mapping{$step} = $step_count;
    $step_count++;
}

#Hold the db lanes involved in the run; create lanes_runs if it's a new run
my @lanes_objs;

#Check the step exists if -resume is used
if ($step_resume) {
	if (!&Step_Exists($step_resume)) {
		modules::Exception->throw("ERROR: Resume step $step_resume isn't contained in the xml step list");
	}
} else {
	#If it's the first step check the align_lanes entries exist and match the total_lanes number; if not we should be starting a run
	if ($run_config->exists('steps', 'step', $steps[0],'first_step')) {
		my $total_lanes = $sample->total_lanes;

		#Get the number of lanes db entries and ensure they match    
	    @lanes_objs = modules::Adaptors::Lane->search(sample_id=>$sample->id);
	    
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
	        } 
	    }
	    
	    if ($align_lane_total != $total_lanes) {
	    	modules::Exception->throw("ERROR: Sample $sample_name expect $total_lanes total_lanes but there are $align_lane_total align_lanes db entries");
	    }   
	}
}

if ($last_step) {
	if (!&Step_Exists($last_step)) {
		modules::Exception->throw("ERROR: Last step $last_step isn't contained in the xml step list");
	}
}

if ($last_step && $step_resume) {
	if (!&Check_Step_Order($step_resume,$last_step)) {
		modules::Exception->throw("ERROR: resume step $step_resume is after last step $last_step");
	}
}

my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime time;
my $timestr = sprintf( "%02d-%02d-%d_%02d:%02d",
						   $mday, 
						   $mon + 1, 
						   $year + 1900,
						   $hour, 
						   $min );

my $sqltime = sprintf( "%d-%02d-%02d",$year+1900,$mon+1,$mday);
my $logfile;

#log object
my $log; 
my $host = hostname;

# Check all error combos of -resume, -latest_sample and directory and db status; set the directory paths
if ( $step_resume || $latestlib_flag) {
		
	if ($latestlib_flag) {
		#Get the runid from the database with this flag
		my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_id);
		if (!defined $run_obj) {
			modules::Exception->throw("Cannot retrieve any runs for $sample_name");
		}
		$runid = $run_obj->id;
		
		#Get the step to resume at from the database if needed
		if (!$step_resume) {
			$step_resume = &Get_Step_Resume();
			if (!$step_resume) {
				modules::Exception->throw("No steps available to run for $sample_name (run: $runid)");
			}
		}
	} else {
		#Here we need to have passed in the runid if we didn't use -latest_sample
		if (!defined $OPT{run_id}) {
			modules::Exception->throw("ERROR: Cannot Run $runid doesn't exist in the database and -resume flag used\n");
		} else {
			$runid = $OPT{run_id};
		}
	}
	&Set_Paths ($runid);
	if (! -d $workdir) {
		modules::Exception->throw("ERROR:  Directory $workdir doesn't exist and -resume or -latest_sample flag used\n");
	}
	
	#Can't resume from first step 
	if ( $step_mapping{$step_resume} == 0) {
		#Bit of hack to present resuming at first step
		modules::Exception->throw("Can't resume from the step $step_resume");
	}
	
	#Set which steps to run, clean up, etc...
	&Resume();
	
	
} elsif ($run_config->exists('steps', 'step', $steps[0],'first_step')) {
	&Initialize(\%ENV);
} else {
	modules::Exception->throw("ERROR: Expecting first_step flag in xml but it is missing");
}

if ($OPT{create_run}) {
	exit;
}

#Now run the common functions
&Run_Pipeline_Steps();
&Finish_Run() unless $OPT{clean_only};

#Checks if a step exists in the config file when -resume flag used
sub Step_Exists {
	my ($check_step) = @_;
	for my $step ( @steps ) {
	    if ($check_step eq $step) {
	    	return 1;
	    }
	}
	return 0;
}

#Check the order of -resume and -last_step makes sense
sub Check_Step_Order {
	my ($first_step,$second_step) = @_;
	if ($step_mapping{$first_step} <= $step_mapping{$second_step}) {
		return 1;
	} else {
		return 0;
	}
}

#Get the step to resume at when using -latest_sample flag and no -resume flag
sub Get_Step_Resume {
	#Iterate over the steps and find the first one that doesn't exist in the database
	for (my $i = 0; $i < scalar @steps; $i++){
		my $step = $steps[$i];
		#Here we are using the latest_sample flag without the -resume flag; use the db to see if the step has been run
		my ($stepobj) = modules::Adaptors::Pipeline_Step->search('name'=>$step);
		my $stepid = $stepobj->id;
	
		#Check if the step has been already run
		my ($runstep) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$runid,pipeline_step_id=>$stepid);
		if (! defined $runstep) {
			#If we haven't run this step yet we'll start here
			return $step;
		}
	}
	
	return 0;
}

#Resuming a run that was stopped for some reason
#CAVEAT: This flag assumes steps are being rerun in the same order as the previous run; won't run properly if this is not the case
sub Resume
{	
	my $step_found = 0;
	#Use the step name to tell what to run
	for (my $i = 0; $i < scalar @steps; $i++){
		my $step = $steps[$i];
		
		#Find the step to start at
		if (!$step_found) {
			#Here we're using the -resume flag or -latest_flag
			if ($step eq $step_resume) {
				#Once we find the step flag everything after for cleanup
				$step_found = 1;			
			} else {
				#Here we confirm there is a db entry if we plan to skip a step
				my ($stepobj) = modules::Adaptors::Pipeline_Step->search('name'=>$step);
				if (!defined $stepobj) {
					modules::Exception->throw("ERROR: Can't find step for $step");
				}
				my $stepid = $stepobj->id;
				#Check if the step has been already run
				my ($runstep) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$runid,pipeline_step_id=>$stepid);
				if (!defined $runstep) {
					modules::Exception->throw("ERROR: Can't resume at step $step_resume because haven't run step $step in the database");
				}
			}

		}
		
		if ($step_found) {
			#Steps after the -resume step potentially needing cleaning
			$steps_to_clean{$step} = 1;
		} else {
			#flag the previous steps as run to avoid re-running
			$steps_run{$step} = 1;
		}
	}
	
	if (!$step_found) {
		modules::Exception->throw("Cannot find step $step_resume in the xml so cannot resume from here");		
	}
	
	#Move the files back to the work directory so we can rerun steps
	modules::Pipeline::move_run_files_back($workdir);
	
	#Resume pipeline run; open the log file for appending
	$log = modules::Log->new($logfile,1);	
	$log->append("Resuming pipeline run $rundir_name at $timestr on $host");
}


#Start the new run; create db entries, record variables used, etc
sub Initialize
{
	my ($LOCAL_ENV) = @_;
	
	
	#Create the run in the db
	my %run = (
				'status' => 'inprocess',
				'production' => 0,
				'commenced' => $sqltime,
				'sample_id' => $sample_id
			   );
	$runid = modules::Adaptors::Run->insert(\%run);
	 
	print STDERR "Create run $runid\n";
	 
	#Now link the lanes to the run by inserting into db; this is used for checking when to start new runs in case new lanes are added to samples
	for my $lane_obj ( @lanes_objs ) {
        my $lane_id = $lane_obj->id;
        my %lanes_runs = (
        				lane_id=>$lane_id,
        				run_id=>$runid
        				);
        modules::Adaptors::Lane_Run->insert(\%lanes_runs);
    } 
	 
	&Set_Paths($runid);

	if (-d $workdir) {
		modules::Exception->throw("ERROR:  directory $workdir already exists: Run with -resume flag and set the runid\n");
	}
	#Get the run object so we can update the working directory
	my ($runobj) = modules::Adaptors::Run->search(id=>$runid);
	$runobj->run_directory($workdir);
	$runobj->update();
	 
	#Create the directories
	&Create_Dirs($workdir);
		
	my $xml_file_dest = $confdir.'/'.$rundir_name.'.xml';
	$sys_call->run("cp $xml_config_file $xml_file_dest");
	
	#Start a new pipeline run; open a new log file
	$log = modules::Log->new($logfile,0);
	$log->append("Starting pipeline run $rundir_name at $timestr on $host");

	#Capture the svn revision number
	my $svn_output = `svn info $ENV{'SVNDIR'} 2>&1`;
	my ($svn_revision_num) = $svn_output =~ /Revision:\s+(\d+)/;
	
	if ( $svn_revision_num !~ /\d/ ) {		
		modules::Exception->throw("Can't obtain svn reversion number");
	}
	
	$log->append("Svn revision: $svn_revision_num");
	$log->append('');

	#Capture the environment variables
	$log->append("ENV VARIABLES:");
	my @env_var = qw(SVNDIR BIOPERL ENSEMBLLIB PERL5LIB ENSEMBL_REGISTRY);
	for my $env ( @env_var ) {
	    $log->append("$env: $LOCAL_ENV->{$env}") if exists $LOCAL_ENV->{$env};
	}
	$log->append('');
	

	
	#Capture the versions of the binaries from the config file
	$log->append("BINARIES:");
	my $bin_struct = $pipe_config->read('binaries');
	
	for my $bin ( keys %{$bin_struct} ) {
		my $version = $bin_struct->{$bin}{version};	    
	    $log->append("$bin: $version");		    
	}
	$log->append('');
}

#Run the actual pipeline steps
sub Run_Pipeline_Steps {
	my %work_files;
		
	#Flag telling us whether to keep running on not (used with -last_step option)
	my $keep_running = 1;	
		
	for (my $i = 0; $i < scalar @steps; $i++) {
		
		last if !$keep_running;
		
		my $step = $steps[$i];
		
		#Set the last step flag
		if ($last_step eq $step) {
			$keep_running = 0;
		}
		
		
		
	    # get command structure
	    my $command = $run_config->read('steps', 'step', $step, 'command');
	    
	    #Special handling for SNV and indel calling steps...
	    my $deletesnv = my $deleteindel = my $deletesnvrow = my $deleteindelrow = 0;
	    
	    #Special handling for filter_basic; need to delete snv entries and not just snv_filter entries
	    if ($run_config->exists('steps', 'step', $step,'deleteSnp')) {
			$deletesnv = 1;	
	    }
	    
	    if ($run_config->exists('steps', 'step', $step,'deleteVar')) {
			$deleteindel = 1;	
	    }
	    
	    if ($run_config->exists('steps', 'step', $step,'deleteSnpRow')) {
			$deletesnvrow = 1;	
	    }
	    
	    if ($run_config->exists('steps', 'step', $step,'deleteVarRow')) {
			$deleteindelrow = 1;	
	    }
	    
	    
	    
	    # figure out input and output files
		my $input_step = $run_config->read('steps', 'step', $step, 'input_from');
		my $overlap_input_step = $run_config->read('steps', 'step', $step, 'input_overlap');
		
	    my $workdir_input =  $workdir .'/' . $rundir_name . '.' .$input_step . '.out';
	    (my $workdir_output = $workdir_input) =~ s/$input_step/$step/;
	    
		my $overlap_input = $workdir .'/'.$rundir_name . '.' . $overlap_input_step . '.match';
		(my $overlap_output = $overlap_input) =~ s/$overlap_input_step/$step/;
		
		my $summary = $workdir . '/summary/' . $rundir_name;

		my $outbam = $workdir . '/bam/' . $rundir_name;

		my $outvcf = $workdir . '/vcf/' . $rundir_name;
		

		#any argument requiring runid as an argument
	    if ($run_config->exists('steps', 'step', $step, 'runid_args')) {
	    	$command =~ s/(\.pl\s+)/$1-runid $runid /g;
	    }
	   	
	   	#Add files to be deleted whether the step is run or not
	    if ( $run_config->exists('steps', 'step', $step, 'delete_files') && $run_config->read('steps', 'step', $step, 'delete_files')) {
			$work_files{$workdir_output.'*'} = 1;
		}	
	   	   	
	   	if ($step_resume && exists $steps_run{$step}) {
	   		$log->append("SKIPPING PREVIOUSLY COMPLETED STEP: $step\n");
	   		print STDERR "SKIPPING STEP: $step\n";
	   		next;
	   	} elsif (exists $steps_to_clean{$step}) {
	   		$log->append("CLEANING UP FILES AND DB FROM PREVIOUS RUN: $step\n");
	   		print STDERR "CLEANING UP STEP: $step\n";
	   		&CleanUp_Step($step,$deletesnv,$deleteindel,$deletesnvrow,$deleteindelrow);
	   		
	   		#Special handling to clean up any existing single_vcf steps when resubmitting snv jobs
	   		if ($step eq 'submit_snvs') {
				&CleanUp_Step('single_vcf',$deletesnv,$deleteindel,$deletesnvrow,$deleteindelrow);   			
	   		}
	   	}
	
	    if ($OPT{clean_only}) {
	    	next;
	    }
	    
	    
	    
		$log->append("RUNNING STEP: $step\n");
		
	    print STDERR "\nRUNNING STEP : $step\n\n";
	
	    #print $workdir_input . "\n";
	    #print $workdir_output . "\n";
	
		
	
		#Replace all the variables with the specifics
	    $command =~ s/WORKDIR/$workdir/g;
		$command =~ s/SUMMARY/$summary/g;
	    $command =~ s/INPUT/$workdir_input/g;
	    $command =~ s/OUTPUT/$workdir_output/g;
	    $command =~ s/OVERLAPIN/$overlap_input/g;
	    $command =~ s/XML/$xml_full_path/g;
	    $command =~ s/OUTBAM/$outbam/g;
	   	$command =~ s/OUTVCF/$outvcf/g;
	    
	    
	    if ($command =~ /OVERLAPOUT/) {
			#Check if the file already exists; error if it does as it should have been previously cleaned up
			if (-e $overlap_output) {
				modules::Exception->throw("ERROR: Step $step has file $overlap_output already existing\n");
			}
		    $command =~ s/OVERLAPOUT/$overlap_output/g;
	    }
	    
	    
	    # run step
	    if ($run_config->read('steps', 'step', $step, 'by_chr')) {
	    	my @chromosomes = split(" ",$pipe_config->read('annotation_version','chr'));
			foreach my $chr (@chromosomes) {
		    	my $per_chr_command = $command;
		    	$per_chr_command =~ s/CHR/$chr/g;
		    	my $t1 = new Benchmark;				
		    	$sys_call->run($per_chr_command) unless $OPT{no_run};
				my $t2 = new Benchmark;				
				my $td = Benchmark::timediff($t2,$t1);
  				my $perchr_command_str = "$per_chr_command took ".Benchmark::timestr($td)."\n";
  				$log->append($perchr_command_str);
		    	#Insert into the database
		    	modules::Pipeline::Insert_Step_Command($step,$per_chr_command,$runid,$chr);
			}
	    } else {
			my $t1 = new Benchmark;				
			$sys_call->run($command) unless $OPT{no_run};
			my $t2 = new Benchmark;				
			my $td = Benchmark::timediff($t2,$t1);
  			my $command_str = "$command took ".Benchmark::timestr($td)."\n";
			$log->append($command_str);
			#Insert into the database
			my $pipeline_step_run_id = modules::Pipeline::Insert_Step_Command($step,$command,$runid);
			if ($run_config->exists('steps', 'step',$step,'release_file')) {
	    		my $release_suffix = $run_config->read('steps', 'step',$step,'release_file');
				my $total_lanes;
				my $file_type;
				
				if ($release_suffix =~ /bam/) {
					$file_type = 'bam';
					$total_lanes = $sample->total_lanes;
				} elsif ($release_suffix =~ /summary/) {
					$file_type = 'summary';
					$total_lanes = $sample->total_lanes;
				} elsif ($release_suffix =~ /vcf/) {
					$file_type = 'vcf';
					
					if ($sample->tumour) {
						#Here we need to count the normal lanes as well
						my ($sample_group_obj)  = modules::Adaptors::Sample_Group->search('id'=>$sample->sample_group_id);
						my ($normal_sample) = modules::Adaptors::Sample->search('sample_group_id'=>$sample->sample_group_id,'tumour'=>0);
						if (!defined $normal_sample) {
							modules::Exception->throw("ERROR: Cannot retrieve normal sample for sample $sample_name");
						}
						$total_lanes = $sample->total_lanes + $normal_sample->total_lanes;
					} else {
						$total_lanes =  $sample->total_lanes;
					}
					
				}  else {
					modules::Exception->throw("ERROR: File type for release file must be bam, vcf, or summary");
				}
				
	    		my $release_file_name = $rundir_name.$release_suffix;
	    		
	    		
	    		#TODO: Copy to central repository; perhaps cronjob after run is complete
	    		my %release_file = (
	    							file_name=>$release_file_name,
	    							file_type=>$file_type,
	    							pipeline_steps_run_id=>$pipeline_step_run_id,
	    							total_lanes=>$total_lanes
	    							);
	    							
	    		modules::Adaptors::Release_File->insert(\%release_file);
	    	}
	    }
	    
		$log->append("FINISHED STEP: $step\n\n");
	    print STDERR "\nFINISHED STEP : $step\n\n";			
	    
	    #If we've run the last pipeline step update the complete flag
	    if ($run_config->exists('steps', 'step', $step,'last_step')) {
	    	$finish_flag = 1;
	    }
	    
	}
	
	
	
	# Delete temporary files
	foreach my $work_file ( keys %work_files) {
		#Check the file was created
		if (`ls $work_file 2>/dev/null`) {		
			$sys_call->run("rm $work_file");
		} 
	}
}

#Clean up the files and db for a specific command before rerunning
sub CleanUp_Step {
	
	my ($stepname,$deletesnv,$deletevar,$deletesnvrow,$deletevarrow) = @_;
	#First clean up the files
	my $workdir_files_to_remove = $workdir .'/*'. $stepname .'*';
	my $overlap_files_to_remove = $overlapdir .'/*'. $stepname .'*';
	my $bam_files_to_remove = $bamdir .'/*' . $stepname . '*';
	
	#delete files (default)
	if ($cleanup_files) { 
		#Don't use sys_call as we don't want to die if there are no files to remove
		system("rm $workdir_files_to_remove 2>/dev/null");
		system("rm $overlap_files_to_remove 2>/dev/null");
		system("rm $bam_files_to_remove 2>/dev/null");
	}
	
	#First get the step_id
	my ($stepobj) = modules::Adaptors::Pipeline_Step->search('name'=>$stepname);
	my $stepid = $stepobj->id;
	
	#In all cases we remove the record of the unsuccessful step
	modules::Adaptors::Pipeline_Step_Run->search(run_id=>$runid,pipeline_step_id=>$stepid)->delete_all;

	#If it's a command that doesn't add db entries
	if (!$run_config->exists('steps', 'step', $stepname, 'runid_args')) {
		return;
	}

	#Now clean up the database (default)
	
	if ($cleanup_db) {
	
		#Now delete the snp entries if it's filter_basic; otherwise just remove the snp filters
		if ($deletesnv) {
			#Cascading delete will take care of the snv_filters
			my @snv_ids = modules::Adaptors::SNV->search_snv_ids($runid);
			modules::Adaptors::BulkDelete->delete(-table_name=>'snvs',-key_ids=>\@snv_ids) if @snv_ids;
		} 
		if ($deletevar) {
			#Cascading delete will take care of the variant_filters; need to remove small insertions and small deletions here
			my @del_ids = modules::Adaptors::Variant->search_variant_ids_type($runid,'DEL');
			modules::Adaptors::BulkDelete->delete(-table_name=>'variants',-key_ids=>\@del_ids) if @del_ids;
			my @ins_ids = modules::Adaptors::Variant->search_variant_ids_type($runid,'INS');
			modules::Adaptors::BulkDelete->delete(-table_name=>'variants',-key_ids=>\@ins_ids) if @ins_ids;
		} 
		if ($deletesnvrow) {
			#here we need to delete the snv_rows; cascading delete takes care of snp_annotations
			my @snv_row_ids = modules::Adaptors::SNV_Row->search(run_id=>$runid);
			modules::Adaptors::BulkDelete->delete(-table_name=>'snv_rows',-key_ids=>\@snv_row_ids) if @snv_row_ids;
		} 
		if ($deletevarrow) {
			#here we need to delete the indel_rows; cascading delete takes care of indel_annotations
			my @indel_row_ids = modules::Adaptors::Variant_Row->search(run_id=>$runid);
			modules::Adaptors::BulkDelete->delete(-table_name=>'variant_rows',-key_ids=>\@indel_row_ids) if @indel_row_ids;	
		} 
		
		#Here we need to remove the snp_filters or variant filters; get all the filters linked to that step
		
		my @filters = modules::Adaptors::Filter->search(pipeline_step_id=>$stepid);
		
		#If it's a snv only step; like filter_dbsnp_snv
		if ($stepname =~ /snv$/) {
			for my $filterid ( @filters ) {
				modules::Adaptors::BulkDelete->delete_snv_filter(-run_id=>$runid,-filter_id=>$filterid);
			}
		} elsif ($stepname =~ /indel$/) {
			#If it's an indel only step; like report_indel
			for my $filterid ( @filters ) {
				modules::Adaptors::BulkDelete->delete_indel_filter(-run_id=>$runid,-filter_id=>$filterid);
			}
		} else {
			#Otherwise remove both filter_snp and filter_variant; like filter_exon
			for my $filterid ( @filters ) {
				modules::Adaptors::BulkDelete->delete_snv_filter(-run_id=>$runid,-filter_id=>$filterid);
				modules::Adaptors::BulkDelete->delete_indel_filter(-run_id=>$runid,-filter_id=>$filterid);
			}
		}
	}
}


sub Finish_Run {
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
	  localtime time;
	my $timestr = sprintf( "%02d-%02d-%d_%02d:%02d",
						   $mday, $mon + 1, $year + 1900,
						   $hour, $min );

	my $host = hostname;
	$log->append("Finished pipeline block $rundir_name at $timestr on $host");
	#Update the run entry status when we reach the last step
	&Set_Production_Flag() if $finish_flag;
	#Move the files to their final subdirectories
	modules::Pipeline::move_run_files($workdir);
	my $end = new Benchmark;				
	my $td = Benchmark::timediff($end,$start);
  	my $final_time_str = "Total run time: ".Benchmark::timestr($td)."\n";
  	$log->append($final_time_str);
	$log->close();
}


#Set the production flag of all the other runs for the sample to false
sub Set_Production_Flag {
	my @runs = modules::Adaptors::Run->search(sample_id=>$sample_id);
	for my $run_obj (@runs) {
		#Update the current run
		if ($run_obj->id == $runid) {
			$run_obj->production(1);		
			$run_obj->status('complete');					
		} else {
			#Change all other runs to non-production
			$run_obj->production(0);
		}
		$run_obj->update();
	}
	
}

#Create the subdirectories
sub Create_Dirs {
    my ( $base_dir ) = @_;
    $sys_call->run("mkdir $base_dir");
    my @subdirs = split(",",$pipe_config->read('run_subdirs'));
    for my $subdir ( @subdirs ) {
        $sys_call->run("mkdir $base_dir/$subdir");
    }
    
}

#Set the local dir variables used throughout from the runid
sub Set_Paths
{
	my ($local_runid) = @_;
	$rundir_name = $sample_name . '_' . $local_runid;
	$workdir = join('/',$rundir_base,$patient_name,$sample_name.'_runs',$rundir_name);
	$overlapdir = $workdir."/overlap";
	$bamdir = $workdir."/bam";
	$summarydir = $workdir."/summary";
	$logdir = $workdir."/log";
	$confdir = $workdir."/conf";
	$vcfdir = $workdir."/vcf";
	$logfile = $logdir.'/'.$rundir_name.'.log';
}
