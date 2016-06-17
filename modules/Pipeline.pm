package modules::Pipeline;

use strict;
use modules::QualityEncoding;
use modules::Exception;
use Data::Dumper;
use modules::Cluster;
use modules::ConfigXML;
use modules::Adaptors::Syscall;
use modules::Adaptors::Pipeline_Step;

sub new  
{
	my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

	#xml objects that contain the config data we need to access
    my @required_args = (
			         	'-patient_id',
			         	'-cluster_obj'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    #Set the global variables
    $self->pipe_xml(get_pipe_conf());
	$self->cluster_xml(get_cluster_conf());    
    $self->{patient_id} = $args{'-patient_id'};
    $self->{cluster_obj} = $args{-cluster_obj};
    return $self;
}

#Set steps_xml
sub pipe_xml
{
    my ($self, $pipe_xml) = @_;

    if (defined $pipe_xml) {
		$self->{'pipe_xml'} = $pipe_xml;
    } elsif (! defined $self->{'pipe_xml'}) {
		modules::Exception->throw("pipe_xml not set");
    }

    return $self->{'pipe_xml'};
}

#Set steps_xml
sub cluster_xml
{
    my ($self, $cluster_xml) = @_;

    if (defined $cluster_xml) {
		$self->{'cluster_xml'} = $cluster_xml;
    } elsif (! defined $self->{'cluster_xml'}) {
		modules::Exception->throw("cluster_xml not set");
    }

    return $self->{'cluster_xml'};
}

sub svndir
{
	my ($self) = @_;
	return $self->cluster_xml->read('svn','svndir');
}


sub threadnum
{
	my ($self) = @_;
	return $self->cluster_xml->read('qsub_vars','thread_num');
}

sub samtools
{
	my ($self) = @_;
	return $self->pipe_xml->read('binaries','samtools','binary');
}

sub bwa
{
	my ($self) = @_;
	return $self->pipe_xml->read('binaries','bwa','binary');
}

sub qsub_dir
{
	my ($self) = @_;
	return $self->cluster_xml->read('base_directories','base_qsub_directory');
}

sub create_sample_xml 
{
	my $self = shift;
    my %args = @_;
	
    my @required_args = (
			             '-sample_name',
			             '-steps_xml_template',
			             '-xml_out',
			             '-encoding',
			             '-bams'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $sample_name = $args{-sample_name};
    my $steps_xml_template = $args{-steps_xml_template};
    my $xml_out = $args{-xml_out};
    my $encoding = $args{-encoding};
    my $patient_name = $self->{patient_id};
    my $sample_group_name = modules::Pipeline::get_sample_group_name($sample_name);
    my $bam_files = join(" ",@{$args{-bams}});
	open(TEMPLATE, $steps_xml_template) || modules::Exception->throw("ERROR: Can't open file $steps_xml_template");
	open(XMLOUT, ">$xml_out")  || modules::Exception->throw("ERROR: Can't open file for writing $xml_out");   
	
	while (<TEMPLATE>) {
		if (/COMMANDS/) {
			#Here we write out the sample specific variables
			print XMLOUT "<!ENTITY nameSample \"$sample_name\">\n";
			print XMLOUT "<!ENTITY nameSampleGroup \"$sample_group_name\">\n";
			print XMLOUT "<!ENTITY patientID \"$patient_name\">\n";
			print XMLOUT "<!ENTITY encoding \"$encoding\">\n\n";
			print XMLOUT $_;
		} elsif (/(.*ENTITY.*)\"(.*)->\((.*)\)/) {
			my $start = $1;
			my $file = $2;
			my @xml_args = split(" ",$3);
			my $value = $self->parse_xml($file,\@xml_args);
			print XMLOUT "$start \"$value\">\n";
		} elsif (/BAMFILES/) {
			#Need to get the input bam files from the lanes into the xml
			$_ =~ s/BAMFILES/$bam_files/g;
			print XMLOUT $_;
		} elsif (/SYM_BAM/) {
			my $bamlinks_dir = "$patient_name/bam_links/";
			#Need to have standard names for bam files for snv_calling later
			$_ =~ s/SYM_BAM/$bamlinks_dir$sample_name.bam/g;
			print XMLOUT $_;
		} else {
			print XMLOUT $_;
		}
	}
	close XMLOUT;
}

#Get patient id from sample_name
sub get_patient_id
{
	my ($sample_name) = @_;
	if ($sample_name =~ /^(.+)_sg[0-9]_\w+$/) {
		return $1;
	} elsif ($sample_name =~ /^(.+)_sg[0-9]$/) {
		return $1;
	} else {
		modules::Exception->throw("ERROR: Can't get patient name from $sample_name");
	}
}

#Get patient id from sample_name
sub get_sample_group_name
{
	my ($sample_name) = @_;
	if ($sample_name =~ /^(.+_sg[0-9])_\w+$/) {
		return $1;
	} elsif ($sample_name =~ /^(.+_sg[0-9])$/) {
		return $1;
	} else {
		modules::Exception->throw("ERROR: Can't get patient name from $sample_name");
	}	
}

#Get the latest stepname for a run
sub get_latest_step {
	my ($run_id) = @_;
	my $last_step_run = 'no_steps_recorded';
	my @pipeline_steps_obj = modules::Adaptors::Pipeline_Step->search_all();
								
	for my $pipeline_step_obj (@pipeline_steps_obj) {
		my ($pipeline_step_run_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$pipeline_step_obj->id);
		if (defined $pipeline_step_run_obj) {
			$last_step_run = $pipeline_step_obj->name;
		}
	}
	return $last_step_run;
}


#Generic parser for converting steps.xml template file into xml used by pipeline
sub parse_xml 
{
	my ($self,$file,$args) = @_;
	my $xml_value;
	if ($file eq 'cluster') {
		$xml_value = $self->cluster_xml->read(@{$args});
	} elsif ($file eq 'pipe') {
		$xml_value = $self->pipe_xml->read(@{$args});
	} else {
		modules::Exception->throw("ERROR: File argument must be pipe or cluster\n");
	}
	return $xml_value;
}

#Create the snv_calling qsubs
sub create_snvcall_qsubs
{
	my $self = shift;
    my %args = @_;
	
    my @required_args = (
			             '-normal_sample',
			             '-tumour_samples',
			             '-sample_group_name',
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }	
    my $normal_sample = $args{-normal_sample};
	my @tumour_samples = @{$args{-tumour_samples}};
	
	my $sg_name = $args{-sample_group_name};
	my $patient_name = $self->{patient_id};

	#These are from the xml files
    my $svndir = $self->svndir();
 	 	
    my $qsub_dir = $self->qsub_dir().'/'.$patient_name.'/snvcalls';
    
   	if ( !-d $qsub_dir ) {
   		system("mkdir -p $qsub_dir");	
   	}
   	
   	my $samtools_bin = $self->samtools();
   	my $bcftools_bin = $self->pipe_xml->read('binaries','bcftools','binary');
   	
   	my $rundir = $self->cluster_xml->read('base_directories','base_run_directory');
   	my $fasta_ref = $self->cluster_xml->read('svn','fasta_file');
   	my @chrs = split(" ",$self->pipe_xml->read('annotation_version','chr'));
   	my $mpileup_args = $self->pipe_xml->read('binaries','samtools','mpileup','args');
   	my $bcftools_tumour_args = $self->pipe_xml->read('binaries','bcftools','view','tumour_args');
   	my $bcftools_normal_args = $self->pipe_xml->read('binaries','bcftools','view','normal_args');
 	my $record_snv_bin = "$svndir/scripts/record_snv_call.pl";  	
 
 	#If we're adding a new sample group only don't bother with this step -> DON'T DO THIS ANYMORE AS HAPPENS WHEN SNVS CALLED SIMULTANEOUSLY
 	unless ($args{-no_normal}) {
	   	#Create the normal sample files for normal snv calls
	   	for my $chr ( @chrs ) {
	   			my $vcf_file = "$rundir/$patient_name/$sg_name"."_snvcalls/$normal_sample.$chr.vcf";
	   			my $snv_call_qsub_tmp = $normal_sample.".snvcall.$chr.qsub.tmp"; #tmp b/c we need to add the run id later
	   			my $snv_call_qsub = $qsub_dir."/".$normal_sample.".snvcall.$chr.qsub"; #pass this to file
	   			#Mpileup command to create vcf
	   			my $pileup_command = "$samtools_bin $mpileup_args $fasta_ref -r $chr $rundir/$patient_name/bam_links/$normal_sample.bam | $bcftools_bin $bcftools_normal_args - | grep -v NNNNN | grep -v ^# > $vcf_file";
	   			#Create a record of vcf being created; need to allow mpileup jobs to finish before proceeding to next pipeline steps
	   			my $record_command = "$record_snv_bin -vcf $vcf_file -chr $chr -step_name single_vcf -qsub_file $snv_call_qsub -runid RUNID"; #RUNID gets replaced by submit_snvs.pl
	   			my @commands = ($pileup_command,$record_command);
	    		$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$snv_call_qsub_tmp,-commands=>\@commands,-cpus=>1, -mem=>'4Gb',-walltime=>'48:00:00');
	   	}
 	}
 
 	#Now create all the tumour sample files (paired with normal) for each tumour snv call
   	for my $tumour_sample ( @tumour_samples ) {
   		for my $chr ( @chrs ) {
   			my $vcf_file = "$rundir/$patient_name/".$sg_name."_snvcalls/$tumour_sample.$chr.vcf";
   			my $snv_call_qsub_tmp = $tumour_sample.".snvcall.$chr.qsub.tmp"; #tmp b/c we need to add the run id later
   			my $snv_call_qsub = $qsub_dir."/".$tumour_sample.".snvcall.$chr.qsub"; #pass this to file
   			#Mpileup command to create vcf
   			my $pileup_command = "$samtools_bin $mpileup_args $fasta_ref -r $chr $rundir/$patient_name/bam_links/$normal_sample.bam $rundir/$patient_name/bam_links/$tumour_sample.bam | $bcftools_bin $bcftools_tumour_args - | grep -v NNNNN | grep -v ^# > $vcf_file";
   			#Create a record of vcf being created; need to allow mpileup jobs to finish before proceeding to next pipeline steps
   			my $record_command = "$record_snv_bin -vcf $vcf_file -chr $chr -step_name single_vcf -qsub_file $snv_call_qsub -runid RUNID"; #RUNID gets replaced by submit_snvs.pl
   			my @commands = ($pileup_command,$record_command);
    		$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$snv_call_qsub_tmp,-commands=>\@commands,-cpus=>1, -mem=>'4Gb',-walltime=>'48:00:00');
   		}
   	}
}

#Generate the bwa align lane to sort bam commands
sub align_lane_qsub 
{
	my $self = shift;
    my %args = @_;
	
    my @required_args = (
			             '-read1',
			             '-read2',
			             '-lane_name',
						 '-lane_id',
						 '-outdir'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $read1 = $args{-read1};
    my $read2 = $args{-read2};
    my $lane_name = $args{-lane_name};
    my $lane_id = $args{-lane_id};
    my $outdir = $args{-outdir};
    
    #These are from the xml files
    my $svndir = $self->svndir();
 	my $threadnum = $self->threadnum();
 	my $max_memory = $self->cluster_xml->read('qsub_vars','max_memory');
 	
 	#Account for samtools poor memory estimation
 	my ($mem_number) = $max_memory =~ /(\d+)/;
 	my $max_sam = $mem_number * 1000000000 * 0.45; 
 	
    my $qsub_dir = defined $args{-qsub_dir}?$args{-qsub_dir}:$self->qsub_dir().'/'.$self->{patient_id}.'/lanes';
   	if ( !-d $qsub_dir ) {
   		system("mkdir -p $qsub_dir");	
   	}
   	my $index = $self->cluster_xml->read('svn','bwa_index');
    

    my %commands = ();
    #If the files need decompression first (only bz2 need this so shouldn't happen anymore)
    if (defined $args{-commands}) {
    	%commands = %{$args{-commands}};
    	$read1 =~ s/.bz2//; #Change the read files to be the uncompressed versions
    	$read2 =~ s/.bz2//;
    }
    
    #Create the file name; all match lane name
    my $read1_sai = $lane_name.'_r1.sai';
    my $read2_sai = $lane_name.'_r2.sai';
    my $sam = $lane_name.'.sam';
	(my $bam = $sam) =~ s/.sam//; 
	my $aln_args;
	my $encoding;
	
	if ($args{-skip_quality}) {
		$encoding = 'phred33';	
	} else {
    	$encoding = modules::QualityEncoding::encoding(-readfile=>"$read1", -reads_to_check=>100000);
	}
	
	if ($encoding eq 'phred33') {
		$aln_args = "aln";
	} else {
		$aln_args = "aln -I";
	}

    my $bwa_qsub_file = $lane_name.'.bwa.qsub';
    my $sam_qsub_file =  $lane_name.'.sam.qsub';
	my $bam_qsub_file =  $lane_name.'.bam.qsub';
	my $bam_stats_qsub_file = $lane_name.'.bam_stats.qsub';
	my $align_lane_qsub_file = $lane_name.'.align_lane.qsub';
	my $qsub_wrapper_file =  $lane_name.'.wrapper.qsub';
	my @jobs = ($bwa_qsub_file,$sam_qsub_file,$bam_qsub_file,$bam_stats_qsub_file,$align_lane_qsub_file);

	my $bwa_full_qsub = $qsub_dir . '/' . $bwa_qsub_file;
	my $sam_full_qsub = $qsub_dir . '/' . $sam_qsub_file;
	my $bam_full_qsub = $qsub_dir . '/' . $bam_qsub_file;
	my $bam_stats_full_qsub = $qsub_dir .'/'.$bam_stats_qsub_file;
	my $samtools_bin = $self->samtools();
	my $bwa_bin = $self->bwa();

	push @{$commands{bwa}}, "$bwa_bin $aln_args -t $threadnum $index $read1 > $outdir/$read1_sai";
	push @{$commands{bwa}}, "$bwa_bin $aln_args -t $threadnum $index $read2 > $outdir/$read2_sai";
	push @{$commands{sam}}, "$bwa_bin sampe $index $outdir/$read1_sai $outdir/$read2_sai $read1 $read2 > $outdir/$sam";
	push @{$commands{bam}}, "grep -v XT:A:R $outdir/$sam | $samtools_bin view -S -b - | $samtools_bin sort -m$max_sam - $outdir/$bam";
	push @{$commands{bam}}, "$samtools_bin index $outdir/$bam.bam";
	push @{$commands{bam_stats}}, "$svndir/scripts/get_bam_stats.pl -bam $outdir/$bam.bam -output $outdir/$bam.stats ";
	push @{$commands{align_lane}}, "$svndir/scripts/align_lanes_db.pl -lane_id $lane_id -phred_quality $encoding -qsub_bam_stats $bam_stats_full_qsub -qsub_bwa $bwa_full_qsub -qsub_sam $sam_full_qsub -qsub_bam $bam_full_qsub";
	
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$bwa_qsub_file,-commands=>$commands{bwa},-cpus=>$threadnum, -mem=>'8Gb', -walltime=>'24:00:00');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$sam_qsub_file,-commands=>$commands{sam},-cpus=>1, -mem=>'8Gb', -walltime=>'24:00:00');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$bam_qsub_file,-commands=>$commands{bam},-cpus=>1, -mem=>$max_memory, -walltime=>'48:00:00');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$bam_stats_qsub_file,-commands=>$commands{bam_stats},-cpus=>1, -mem=>'1Gb');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$align_lane_qsub_file,-commands=>$commands{align_lane},-cpus=>1, -mem=>'1Gb');
	$self->{cluster_obj}->create_qsub_wrapper(-qsub_dir=>$qsub_dir,-qsub_file=>$qsub_wrapper_file,-jobs=>\@jobs);
    
    return ($encoding,"$outdir/$bam.bam");
}

#Create all the individual qsub files and a wrapper for a range of steps
sub create_run_qsubs
{
	my $self = shift;
    my %args = @_;
    my @required_args = (
			            -steps_xml,
			            -start_step,
			            -sample_name,
			            -pipe_block
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $svndir = $self->svndir();
	my $pipeline_script = $svndir . '/scripts/run_pipeline_steps.pl';
    my $qsub_dir = defined $args{-qsub_dir}?$args{-qsub_dir}:$self->qsub_dir().'/'.$self->{patient_id}.'/runs';
    my $new = '';
	my $start_step = $args{-start_step};
	my $sample_name = $args{-sample_name};
	my $pipe_block = $args{-pipe_block};
	
	
	#Get the steps xml
	my $steps_xml = $args{-steps_xml};
	if ( !-e $steps_xml ) {
		modules::Exception->throw("File $steps_xml doesn't exist");	
	}
	my $steps_config = modules::ConfigXML->new($steps_xml);
	my $steps = $steps_config->read('steps_order', 'step');
	my @steps = ();
	

	#Get the steps list	
	if (ref($steps) eq 'ARRAY'){ # Cope with steps being a single step or an array of steps
    	@steps = @$steps;
	} else {
    	@steps = ($steps);
	}
	
	my %step_numbers = my %step_lookup = ();
	#Get the step order to check last_step is after resume_step
	my $step_count = 0;
	for my $step ( @steps ) {
    	$step_numbers{$step} = $step_count;
    	$step_lookup{$step_count} = $step;
    	$step_count++;
	}
	
	#Check the step arguments exist in the xml
	if (!$self->_step_exists($start_step,\@steps)) {
		modules::Exception->throw("ERROR: Step $start_step doesn't exist in xml file $steps_xml");
	}
	
	my $end_step;
	if (defined $args{-end_step}) {
		$end_step = $args{-end_step};
		if (!$self->_step_exists($end_step,\@steps)) {
			modules::Exception->throw("ERROR: Step $end_step doesn't exist in xml file $steps_xml");
		}
	} else {
		$end_step = $start_step;
	}
	
	my $first_step_count = $step_numbers{$start_step};
	my $last_step_count = $step_numbers{$end_step};
	
	my @dependent_jobs = ();
	
	
	for ( my $count = $first_step_count ; $count <= $last_step_count ; $count++ ) {
		my $local_step = $step_lookup{$count};
		my $qsub_file = $sample_name.'.pipe'.$pipe_block.'.'.$local_step.'.qsub';
		
		
		#Check if we're dealing with the first step
		my $first_step = $steps_config->exists('steps', 'step', $local_step, 'first_step')?1:0;
		my $pipe_args = "-xml $steps_xml";	
	    if ($first_step) {
			$pipe_args .= " -last_step $local_step";
	    } else {
	    	$pipe_args .= " -resume $local_step -latest_sample $sample_name -last_step $local_step";
	    }
	    
	   	if ( !-d $qsub_dir ) {
	   		system("mkdir -p $qsub_dir");	
	   	}
	   	
	   	my @commands = ("$pipeline_script $pipe_args");
		
		#Get the resources required for the step in question from the xml
		my $mem = $steps_config->exists('steps', 'step', $local_step, 'mem')? $steps_config->read('steps', 'step', $local_step, 'mem'):"4GB";
		my $cpus = $steps_config->exists('steps', 'step', $local_step, 'cpus')? $steps_config->read('steps', 'step', $local_step, 'cpus'):"1";
		my $walltime = $steps_config->exists('steps', 'step', $local_step, 'walltime')? $steps_config->read('steps', 'step', $local_step, 'walltime'):"24:00:00";	
	    if ($steps_config->exists('steps','step',$local_step,'module')) {
			my $module = $steps_config->read('steps','step',$local_step,'module');
			$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$qsub_file,-commands=>\@commands,-cpus=>$cpus,-mem=>$mem,-walltime=>$walltime,-modules=>$module);
	    } else {
		    $self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$qsub_file,-commands=>\@commands,-cpus=>$cpus,-mem=>$mem,-walltime=>$walltime);
	    }
		push @dependent_jobs, $qsub_file;
		last if $local_step eq $end_step;
	}
	
	#Now create the wrapper if it's more than one command
	my $qsub_wrapper_file = $sample_name.'.pipe'.$pipe_block.'.wrapper.qsub';
	$self->{cluster_obj}->create_qsub_wrapper(-qsub_dir=>$qsub_dir,-qsub_file=>$qsub_wrapper_file,-jobs=>\@dependent_jobs) if $end_step ne $start_step;
}

#Checks if a step exists in the config file when -resume flag used
sub _step_exists {
	my ($self,$check_step,$steps) = @_;
	
	for my $step ( @{$steps} ) {
	    if ($check_step eq $step) {
	    	return 1;
	    }
	}
	return 0;
}

#Get the commands from a qsub file
sub get_commands {
	my ($qsub_file) = @_;
	
	if ( !-e $qsub_file ) {
		modules::Exception->throw("File $qsub_file doesn't exist");	
	}
	
	my @commands;
	
	open(FILE,"$qsub_file") || modules::Exception->throw("Can't open file $qsub_file\n");
	my $command_flag = 0;
	while (<FILE>) {
		chomp;
		push @commands, $_ if $command_flag;
		if (/source/) {
			$command_flag = 1;
		}
	}
	return \@commands;
}

#Move the run files to subdirectories from workdir
sub move_run_files {
	my ($workdir) = @_;
   	my $overlap_files = $workdir .'/*sg*';
   	my $overlapdir = $workdir.'/overlap';
   	my $overlap_flag = `ls $overlap_files 2>/dev/null`;
   	system("mv $overlap_files $overlapdir") if $overlap_flag;
}


#Move the run files back to workdir from subdirectories
sub move_run_files_back {
	my ($workdir) = @_;
   	my $overlap_files = $workdir .'/overlap/*sg*'; #Remainder goes to overlap
   	my $overlap_flag = `ls $overlap_files 2>/dev/null`;
   	system("mv $overlap_files $workdir") if $overlap_flag;
}


#Insert the step_command and bin_vers entry if required
sub Insert_Step_Command {
	my ($stepname,$command,$runid,$chr) = @_;
	
	
	
	#Now create the pipeline_step_run entry
	my ($stepobj) = modules::Adaptors::Pipeline_Step->search('name'=>$stepname);
	if (! defined $stepobj) {
		modules::Exception->throw("Step $stepname doesn't exist in the database");
	}
	my $stepid = $stepobj->id;
	#Insert the step_command entry
    my %step_command = (
    					'pipeline_step_id' => $stepid,
    					'run_id' => $runid
    					);
    my $pipeline_step_run_id = modules::Adaptors::Pipeline_Step_Run->insert(\%step_command);		
    
    #First create the syscall entry
	my %syscall_entry = (
						'command' => $command,
						'level' => 'sample',
						'pipeline_steps_run_id'=>$pipeline_step_run_id,
						'analysis_step_name'=>$stepname
						);
	
	#If it's a chromosome specific command					
	if (defined $chr) {
		$syscall_entry{'chr'} = $chr;
	}	
					
	my $syscall_id = modules::Adaptors::Syscall->insert(\%syscall_entry);	
	
    return $pipeline_step_run_id;
}

sub get_cluster_conf {
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
	
	return $cluster_config;
}

sub get_pipe_conf {
	my $svndir;
	if (!exists $ENV{'SVNDIR'}) {
		modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
	} else {
		$svndir = $ENV{'SVNDIR'};
	}
	
	#Get the xml files and create the pipeline object
	my $pipe_xml = "$svndir/conf/pipe.xml";

	if ( !-e $pipe_xml ) {
		modules::Exception->throw("File $pipe_xml doesn't exist");	
	}
	
	my $pipe_config = modules::ConfigXML->new($pipe_xml);
	
	return $pipe_config;
}

sub get_report_conf {
	my $svndir;
	if (!exists $ENV{'SVNDIR'}) {
		modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
	} else {
		$svndir = $ENV{'SVNDIR'};
	}
	
	#Get the xml files and create the pipeline object
	my $report_xml = "$svndir/conf/report.xml";

	if ( !-e $report_xml ) {
		modules::Exception->throw("File $report_xml doesn't exist");	
	}
	
	my $report_config = modules::ConfigXML->new($report_xml);
	
	return $report_config;
}

sub get_tumour_flag {
	my $self = shift;
    my %args = @_;
	
    my @required_args = (
			             '-sample_name'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $sample_name = $args{-sample_name};
	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);
	if (!$sample_obj) {
		modules::Exception->throw("ERROR: Cannot retrieve sample object for $sample_name");
	}
	
	if ($sample_obj->tumour) {
		return 1;	
	} else {
		return 0;
	}
	
	
}


return 1;
