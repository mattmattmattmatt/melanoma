package modules::Cluster;

use strict;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;

sub new {
	my ($class, @args) = @_;
    my %args = @args;
    my $self = bless {}, $class;
	my @required_args = (
			             -scheduler,
			             -svndir
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
    }

	my $scheduler = $args{-scheduler};
	my $svndir = $self->svndir($args{-svndir});

    if ($scheduler eq 'PBS') {
		$self->{pbs} = 1;	
    } elsif ($scheduler eq 'SGE') {
    	$self->{sge} = 1;
    } else {
		modules::Exception->throw("Only PBS and SGE support at this time");
    }

	return $self;
}

#Subroutine to create the qsub file
sub create_single_qsub {
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -qsub_dir,
			             -qsub_file,
			             -commands
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $qsub_dir = $args{-qsub_dir};
	if ( !-d $qsub_dir ) {
		system("mkdir -p $qsub_dir");
	}
	my $qsub_file = $args{-qsub_file};
	my $qsub_out = $qsub_dir . '/' . $qsub_file;
	
	
	my $bash = defined $args{-bash}?$args{-bash}:"#!/bin/bash";
	#Add shebang if not included
	if ($bash !~ /^#/) {
		$bash = "#!".$bash;
	}
	my $queue = defined $args{-queue}?$args{-queue}:"normal";
	my $project = defined $args{-project}?$args{-project}:"u86";
	
	my $walltime = defined $args{-walltime}?$args{-walltime}:"24:00:00";
	if ($walltime !~ /\d+:\d\d:\d\d/) {
		modules::Exception->throw("ERROR: Walltime arguments must be in the format [H]H:MM:SS\n");
	}
	my $jobfs = defined $args{-jobfs}?$args{-jobfs}:"10GB";
	my $mem = defined $args{-mem}?$args{-mem}:"4GB";
	$mem = uc($mem);
	
	if ($mem !~ /GB$/ && $mem !~ /MB$/) {
		modules::Exception->throw("ERROR: Memory must end in GB or MB\n");
	}
		
	my $cpus = defined $args{-cpus}?$args{-cpus}:1;
	if ($cpus !~ /^\d+$/) {
		modules::Exception->throw("ERROR: cpus must be an integer\n");
	}
	
	my $modules = defined $args{-modules}?$args{-modules}:"";
	my $commands = $args{-commands};
	if (ref($commands) ne 'ARRAY') {
		modules::Exception->throw("ERROR: commands args must be array ref");
	} 
	my $command_str = join("\n",@{$commands});

	my $source_line = 'source '. $self->svndir(). '/conf/export_env.txt '.$self->svndir();

	my $resource_line;

	open(QSUB,">$qsub_out") || modules::Exception->throw("ERROR: Cannot open qsub directory $!");
	print QSUB "$bash\n";

	#Now create the qsub file
	if ($self->{pbs}) {
		$resource_line = "walltime=$walltime,vmem=$mem,jobfs=$jobfs,ncpus=$cpus";
		print QSUB "#PBS -P $project\n";
		print QSUB "#PBS -q $queue\n";
		print QSUB "#PBS -l $resource_line\n";
		print QSUB "$source_line\n";
		if ($modules =~ /module/) {
			print QSUB "$modules\n";
		}
		print QSUB "$command_str\n";			
	} elsif ($self->{sge}) {
		$resource_line = "h_rt=$walltime,virtual_free=$mem\n#\$ -pe threads $cpus";
		print QSUB "#\$ -cwd\n";
		print QSUB "#\$ -S $bash\n";
		print QSUB "#\$ -q $queue\n";
		print QSUB "#\$ -l $resource_line\n";
		print QSUB "$source_line\n";
		print QSUB "$command_str\n";			
	}
	close QSUB;
}

#Subroutine to create the bash wrapper file with dependent job submission
sub create_qsub_wrapper {
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -qsub_dir,
			             -qsub_file,
			             -jobs
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $qsub_dir = $args{-qsub_dir};
	if ( !-d $qsub_dir ) {
		system("mkdir -p $qsub_dir");
	}
	my $qsub_file = $args{-qsub_file};
	my $qsub_out = $qsub_dir . '/' . $qsub_file;
	my $jobs = $args{-jobs};
	if (ref($jobs) ne 'ARRAY') {
		modules::Exception->throw("ERROR: jobs args must be array ref");
	} 
	my @full_jobs = map {$_ = $qsub_dir .'/'.$_} @{$jobs};
	my $first_job = shift @full_jobs;
	my $dependent_jobs = join(" ",@full_jobs);
		
	
	my $bash = defined $args{-bash}?$args{-bash}:"#!/bin/bash";
	#Add shebang if not included
	if ($bash !~ /^#/) {
		$bash = "#!".$bash;
	}
	
	#Now create the qsub wrapper file
	open(WRAPPER,">$qsub_out") || modules::Exception->throw("ERROR: Cannot open qsub directory $!");
	print WRAPPER "$bash\n";
	print WRAPPER "jobid=`qsub $first_job`\n";
	print WRAPPER "regex=\"([0-9]+)\"\n";
	print WRAPPER "if [[ \$jobid =~ \$regex ]]\nthen\n";
	print WRAPPER "\tfor job in $dependent_jobs\n";	
	if ($self->{pbs}) {
		print WRAPPER "\t\n\tdo\n\t\tif [[ \$jobid =~ \$regex ]]\n\t\tthen\n\t\t\tjobid=`qsub -W depend=afterok:\${BASH_REMATCH[1]} \$job`\n\t\telse\n\t\t\techo 'ERROR: Later submission did not return jobid' \n\t\tfi\n\tdone\nelse\n\techo 'ERROR: First submission did not return jobid'\nfi\n";
	} elsif ($self->{sge}) {
		print WRAPPER "\t\n\tdo\n\t\tif [[ \$jobid =~ \$regex ]]\n\t\tthen\n\t\t\tjobid=`qsub -hold_jid \${BASH_REMATCH[1]} \$job`\n\t\telse\n\t\t\techo 'ERROR: Later submission did not return jobid' \n\t\tfi\n\tdone\nelse\n\techo 'ERROR: First submission did not return jobid'\nfi\n";
	}
	close WRAPPER;
	
}

#Set svndir
sub svndir
{
    my ($self, $svndir) = @_;

    if (defined $svndir) {
		$self->{'svndir'} = $svndir;
    } elsif (! defined $self->{'svndir'}) {
		modules::Exception->throw("svndir not set");
    }

    return $self->{'svndir'};
}

#Checks whether merge_vcf job is running already or not; required for resume_run to avoid submitting same job twice
sub check_job_running
{
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -sample_name
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $sample_name = $args{-sample_name};
	if ($self->{pbs}) {
		my @qstat_lines = split("\n",`qstat -f 2>&1`);
		my $job_id;
		my $job_state;
		my $job_name;
	
	    for my $qstat_line (@qstat_lines) {
	        if ($qstat_line =~ /Job Id: (\S+)/) {
	        	$job_id = $1;
	        } elsif ($qstat_line =~ /Job_Name = (\S+)/) {
	        	$job_name = $1;
	        } elsif ($qstat_line =~ /job_state = (\S+)/) {
	        	$job_state = $1;
	        	if ($job_state ne 'H' && $job_name =~ /$sample_name/) {
	        		return "$job_id:$job_name";
	        	}
	        	#Reset the values
	        	$job_name = '';
	        	$job_id = '';
	        	$job_state = '';
	        }
		}
	} elsif ($self->{sge}) {
		#TODO: Write for sge
	}
	return 0;
	
	
	
}


# THIS IS INITIAL MODULE
#sub new {
#	my ($class, @args) = @_;
#
#    my %args = @args;
#
#    my $self = bless {}, $class;
#
#	my @required_args = (
#			             -jobs
#						 );
#
#    foreach my $required_arg (@required_args){
#
#		if (! defined $args{$required_arg}){
#		    modules::Exception->throw("Required argument [$required_arg] not set");
#		}	
#		
#		
#			
#		$self->{jobs} = $args{$required_arg};	
#    }
#    
#	#By default don't resubmit failed jobs
#	my $max_resubmit = defined $args{-max_resubmit}?$args{-max_resubmit}:0;
#	$self->{max_resubmit} = $max_resubmit;
#	
#	if (defined $args{-resubmit_changes}) {
#		$self->{resubmit_changes} = $args{-resubmit_changes};
#	}
#	
#	my %status_map = (
#						Q => 'QUEUE',
#						H => 'HELD',
#						E => 'COMPLETE',
#						R => 'RUNNING',
#						T => 'MOVING',
#						W => 'WAITING',
#						S => 'SUSPENDED'
#					);
#	
#	$self->{status_map} = \%status_map;
#	
#    return $self;
#}
#
##Return the max resubmit count
#sub get_max_resubmit {
#	my $self = shift;
#	return $self->{max_resubmit};
#}
#
##Subroutine to run jobs in the data structure
#sub submit_jobs {
#	my $self = shift;
#	
#	
#	for my $job_count (keys %{$self->{jobs}}) {
#		if ($self->{jobs}{$job_count}{job_status} eq 'READY') {
#			my $qsub_file = $self->{jobs}{$job_count}{qsub_file};
#			my $submit_command =  "qsub $qsub_file";
#			
#			my $job_id = `$submit_command`;
#		
#			
#			#Check we got a job id returned
#			if (!$job_id) {
#				modules::Exception->throw("ERROR: job submission failed for $submit_command; nothing returned");				
#			} elsif ($job_id !~ /^\d+/) {
#				modules::Exception->throw("ERROR: job submission failed for $submit_command; non numeric job id $job_id");	
#			} elsif ($job_id =~ /^(\d+)/) {
#				$self->{jobs}{$job_count}{job_id} = $1;
#				$self->{jobs}{$job_count}{job_status} = 'SUBMITTED';
#			}
#		}
#	}
#}
#

#
##Subroutine to create new qsub file with different resources; used with -change_resources flag
#sub _create_new_qsub {
#	my ($self, @args) = @_;
#	
#	my %args = @args;
#
#	my @required_args = (
#			             -input_dir,
#			             -job_count
#						 );
#
#    foreach my $required_arg (@required_args){
#		if (! defined $args{$required_arg}){
#		    modules::Exception->throw("Required argument [$required_arg] not set");
#		}	
#	}
#	
#	my $job_count = $args{-job_count};
#	my $resubmit_count = $self->{jobs}{$job_count}{resubmits};
#	
#	my $file = $args{-input_dir} . '/polyphen.qsub.'.$job_count.'.'.$resubmit_count;
#	$self->{jobs}{$job_count}{qsub_file} = $file;
#	my @new_resources = ();
#	
#	my $resource_changes = $self->{resubmit_changes};
#	
#	for my $resource_request ( @{$self->{jobs}{$job_count}{resources}} ) {
#	    my ($field,$value) = split("=",$resource_request);
#	    if (exists $resource_changes->{$field}) {
#	    	my $new_value_str;
#	    	#Here we increase the value
#	    	if ($value =~ /GB/i) {
#	    		#Special case for vmem or jobfs
#	    		my ($number) = $value =~ /(\d+)/;
#	    		my $new_value = $number * $resource_changes->{$field};
#	    		$new_value_str = $new_value. 'GB';
#	    	} else {
#	    		
#	    		#This will work for ncpus field
#	    		$new_value_str = $value * $resource_changes->{$field};
#	    	}
#
#	    	my $new_resource = $field.'='.$new_value_str;
#	    	push @new_resources, $new_resource;
#	    } else {
#	    	#Here the value stays the same
#	    	push @new_resources, $resource_request;
#	    }
#	}
#	
#	$self->{jobs}{$job_count}{resources} = \@new_resources;
#	my $resources = join(",",@new_resources);
#	open(NEWQSUB,">$file") || modules::Exception->throw("Can't open qsub file $file");
#	print NEWQSUB "#!/bin/bash\n";	
#	print NEWQSUB "#PBS -P u86\n";
#	print NEWQSUB "#PBS -q $self->{jobs}{$job_count}{queue}\n";
#	print NEWQSUB "#PBS -l $resources\n";
#	#Need to load the java module here for polyphen
#	if (defined $self->{jobs}{$job_count}{modules} && @{$self->{jobs}{$job_count}{modules}}) {
#		for my $module ( @{$self->{jobs}{$job_count}{modules}} ) {
#			print NEWQSUB "module load $module\n";			    
#		}
#	}
#	
#	print NEWQSUB join("\n", @{$self->{jobs}{$job_count}{commands}}), "\n";	
#	
#}
#
##Subroutine to update the status of jobs
#sub update_jobs {
#	my $self = shift;
#	
#	for my $job_count (keys %{$self->{jobs}}) {
#		
#		#If the job has been submitted and isn't complete or error update the status
#		if ($self->{jobs}{$job_count}{job_id} != 0 && $self->{jobs}{$job_count}{job_status} ne 'COMPLETE' && $self->{jobs}{$job_count}{job_status} ne 'ERROR') {
#			my $job_status = $self->check_job(-job_id=>$self->{jobs}{$job_count}{job_id},-job_count=>$job_count);
#			
#			#update status if it's changed
#			if ($job_status ne $self->{jobs}{$job_count}{job_status}) {
#				$self->{jobs}{$job_count}{job_status} = $job_status;	
#			}
#		}
#	}
#	
#}
#
##Subroutine to check job status of a particular job
#sub check_job {
#	my ($self, @args) = @_;
#	
#	my %args = @args;
#
#	my @required_args = (
#						 -job_count,
#			             -job_id
#						 );
#
#    foreach my $required_arg (@required_args){
#		if (! defined $args{$required_arg}){
#		    modules::Exception->throw("Required argument [$required_arg] not set");
#		}	
#	}
#	
#	my $job_id = $args{-job_id};
#	my $job_count = $args{-job_count};
#	
#	my @qstat_lines = split("\n",`qstat $job_id 2>&1`);
#	
#	my @qstat_args = split(" ",$qstat_lines[-1]);
#
#	#This means job is done when we receive
#	if ($qstat_lines[0] =~ /qstat: Unknown Job Id/) {
#		#Check the job generated output if we're expecting it
#		my $output = $self->{jobs}{$job_count}{output};
#		if (exists $self->{jobs}{$job_count}{output} && -s $self->{jobs}{$job_count}{output}) {
#			return 'COMPLETE';
#		} else {
#			return 'ERROR';
#		}
#	} else {
#		#Return the mapped value of the letter pbs returns
#		return $self->{status_map}{$qstat_args[4]};	
#	}
#
#}
#
##Checks if there are queued jobs
#sub has_jobs_queued {
#	my $self = shift;
#	
#	for my $job_count (keys %{$self->{jobs}}) {
#		
#		if ($self->{jobs}{$job_count}{job_status} eq 'QUEUED') {
#			return 1;
#		}
#	}
#	return 0;
#}
#
##Checks if there are running jobs
#sub has_jobs_running {
#	my $self = shift;
#	
#	for my $job_count (keys %{$self->{jobs}}) {
#		
#		if ($self->{jobs}{$job_count}{job_status} ne 'COMPLETE' || $self->{jobs}{$job_count}{job_status} ne 'COMPLETE') {
#			return 1;
#		}
#	}
#	return 0;
#	
#}
#
##Checks all jobs are done
#sub has_working_jobs {
#	my $self = shift;
#	
#	for my $job_count (keys %{$self->{jobs}}) {
#		
#		#If the job isn't complete or in error
#		if ($self->{jobs}{$job_count}{job_status} ne 'COMPLETE' && $self->{jobs}{$job_count}{job_status} ne 'FAIL') {
#			return 1;
#		}
#		
#	}
#	return 0;
#}
#
##Give a summary of all jobs
#sub summary {
#	my $self = shift;
#	my %job_status_summary = ();
#	my $total_jobs = keys %{$self->{jobs}};
#	for my $job_count (keys %{$self->{jobs}}) {
#		$job_status_summary{$self->{jobs}{$job_count}{job_status}}++;
#	}
#	
#	print "Total jobs: $total_jobs\n";
#	
#	for my $job_type ( sort {$job_status_summary{$a} <=> $job_status_summary{$b}} keys %job_status_summary ) {
#	    print "\t$job_type: $job_status_summary{$job_type}\n";
#	}
#	
#	my @headers = qw(job_count job_id job_status resubmit);
#	
#	print "\n". join("\t",@headers) . "\n";
#	
#	for my $job_count (sort {$a<=>$b} keys %{$self->{jobs}}) {
#		print join("\t", 
#					$job_count,
#					$self->{jobs}{$job_count}{job_id},
#					$self->{jobs}{$job_count}{job_status},
#					$self->{jobs}{$job_count}{resubmits}
#					) ."\n";
#	}
#	
#	
#}
#
##Get the jobs in state of error
#sub get_error_jobs {
#	my $self = shift;
#	my %error_jobs = ();
#	for my $job_count (keys %{$self->{jobs}}) {
#		if ($self->{jobs}{$job_count}{job_status} eq 'ERROR') {
#			$error_jobs{$job_count} = $self->{jobs}{$job_count};
#		}
#	}
#	return \%error_jobs;
#}
#
##Set the jobs in state of error once it been resubmitted enough times
#sub set_fail_job {
#	my ($self, @args) = @_;
#	
#	my %args = @args;
#
#	my @required_args = (
#						 -job_count,
#						 );
#
#    foreach my $required_arg (@required_args){
#		if (! defined $args{$required_arg}){
#		    modules::Exception->throw("Required argument [$required_arg] not set");
#		}	
#	}
#	
#	my $job_count = $args{-job_count};
#	$self->{jobs}{$job_count}{job_status} = 'FAIL';
#	
#}
#
##Resubmit a job
#sub set_job_for_resubmit {
#	my ($self, @args) = @_;
#	
#	my %args = @args;
#
#	my @required_args = (
#						 -job_count,
#						 -input_dir
#						 );
#						 
#	foreach my $required_arg (@required_args){
#		if (! defined $args{$required_arg}){
#		    modules::Exception->throw("Required argument [$required_arg] not set");
#		}	
#	}
#	
#	my $job_count = $args{-job_count};
#	my $input_dir = $args{-input_dir};
#	
#	$self->{jobs}{$job_count}{resubmits}++;
#	$self->{jobs}{$job_count}{job_status} = 'READY';
#	if (defined $self->{resubmit_changes}) {
#		$self->_create_new_qsub(-job_count=>$job_count,-input_dir=>$input_dir);
#	}
#}


return 1;
