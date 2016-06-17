#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Cluster;
use modules::Pipeline;
use modules::Adaptors::Patient;
use modules::Adaptors::Sample_Group;
use modules::Adaptors::Sample;
use modules::Adaptors::Lane;
use modules::Adaptors::Read_File;
use modules::Adaptors::Run;
use modules::Adaptors::Sequencing_Centre;
use Pod::Usage;
use Cwd;

use vars qw(%OPT);


GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"researcher_id=s",
	   		"normal_bpa_ids=s",
	   		"tumour_bpa_ids=s",
	   		"normal_readdirs=s",
	   		"tumour_readdirs=s",
	   		"normal_sample_types=s",
	   		"tumour_sample_types=s",
	   		"sequencing_centre=s",
	   		"tumour_types=s",
	   		"tumour_tissues=s",
	   		"read1_pattern=s",
	   		"cluster_xml=s",
	   		"steps_xml=s",
	   		"skip_quality",
	   		"new_sample_group",
	   		"test",
	   		"only_qsubs",
	   		"sample_group_number=i"
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{researcher_id} || !$OPT{normal_bpa_ids} || !$OPT{tumour_bpa_ids} || !$OPT{normal_readdirs} || !$OPT{tumour_readdirs} || !$OPT{normal_sample_types} || !$OPT{tumour_sample_types} || !$OPT{tumour_types});

	   
=pod

=head1 SYNOPSIS

add_cancer_sample.pl -normal_bpa_id bpa_barcode_normal -new_sample_group new_sample_group_to_existing_patient(default=new_patient) -tumour_bpa_id bpa_barcode_tumour -researcher_id researcher_id(e.g.tb412) -read1_pattern read1_regex_pattern(default=R1) -steps_xml xml_steps_file(default=../conf/steps/xml) -cluster_xml xml_cluster_file(default=../conf/cluster.xml) -sequencing_center sequencing_centre(option=BRF,AGRF,RAM;default=AGRF) -normal_sample_types comma_delim_sample_type -tumour_sample_types comma_delim_sample_type -tumour_types comma_delim_tumour_types(options=primary,metatstatic)  -normal_readdirs comma_delim_fq_directories -tumour_readdirs comma_delim_fq_directories -skip_quality Skip_read_check_and_assume_phred33_qualities -only_qsubs only_create_qsubs -sample_group_number sample_group_number_for_only_qsubs_flag [options]

Required flags: -normal_bpa_id -tumour_bpa_id -researcher_id -normal_readdirs -tumour_readdirs -normal_sample_types -tumour_sample_types -tumour_types

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

add_cancer_sample.pl -> Add patients, sample_groups, and sample to db

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./add_cancer_sample.pl -tumour_readdirs MASCRI-434_27604_Primary,MASCRI-434_27604_Metastasis -tumour_types primary,metastatic -sequencing_centre AGRF -read1_pattern R1 -tumour_sample_types tumour_primary_patient,tumour_metastatic_patient -researcher_id MASCRI-434_27604 -normal_readdirs MASCRI-434_27604_Blood -normal_sample_types normal_patient -normal_bpa_ids 102.100.100/7701 -tumour_bpa_ids 102.100.100/7699,102.100.100/7700
./add_cancer_sample.pl -tumour_readdirs A15_metastatic_2,A15_primary -tumour_types metastatic,metastatic -sequencing_centre AGRF -read1_pattern R1 -tumour_sample_types tumour_metastatic_cell_line,tumour_metastatic_cell_line -researcher_id A15 -normal_readdirs A15_LCL_2 -normal_sample_types normal_cell_line -normal_bpa_ids 102.100.100/7710 -tumour_bpa_ids 102.100.100/7709,102.100.100/7709 -skip_quality -only_qsubs -sample_group_number 2 -new_sample_group -test
./add_cancer_sample.pl -tumour_tissues brain -tumour_readdirs tumour -tumour_types metastatic -sequencing_centre BRF -read1_pattern R1 -tumour_sample_types tumour_metastatic_patient -researcher_id SCC11 -normal_readdirs normal -normal_sample_types normal_patient -normal_bpa_ids 102.100.100/7731 -tumour_bpa_ids 102.102.100/7730 -skip_quality -only_qsubs -test -read1_pattern read1

=cut


my @sample_types = modules::Adaptors::Sample->search_sample_type();

#ugly way to get enum values 
my $text = $sample_types[0]->{column_type};
$text =~ s/enum\(//;
$text =~ s/\)//;

my %possible_sample_types;
for my $field (split(",",$text)) {
	$field =~ s/\'//g;
	$possible_sample_types{$field} = 1;
}

my $writedb = defined $OPT{only_qsubs}?0:1;

#Get the arguments
my $researcher_id = $OPT{researcher_id};

#First create the patient entry
my %patient_info = (
					patient_external_name=>$researcher_id
					);

my $patient_db_id;
my ($patient_obj) = modules::Adaptors::Patient->search('patient_external_name'=>$researcher_id);
if ($OPT{new_sample_group}) {
	$patient_db_id = $patient_obj->id;
	print STDERR "Add sample group to patient $researcher_id with id $patient_db_id\n";
} elsif ($writedb) {
	#normal case; writing out new objects
	if (defined $patient_obj) {
		modules::Exception->throw("ERROR: A patient already exists with this researcher_id ($researcher_id). Must delete patient first");
	}
	#normal case here; insert objects into database
	$patient_db_id = modules::Adaptors::Patient->insert(\%patient_info);
	($patient_obj) = modules::Adaptors::Patient->search('id'=>$patient_db_id);
	print STDERR "Create patient $researcher_id with id $patient_db_id\n";
} else {
	$patient_db_id = $patient_obj->id;
}



if (!defined $patient_obj) {
	modules::Exception->throw("ERROR: Can't retrieve patient object for $researcher_id");
}

my $sample_group_count;

#If we're adding a sample_group to existing patient
if ($OPT{new_sample_group}) {
	if ($writedb) {
		my @sample_group_objs = modules::Adaptors::Sample_Group->search(patient_id=>$patient_obj->id);
		if (!@sample_group_objs) {
			modules::Exception->throw("ERROR: No existing sample_groups for $researcher_id");
		}
		#Get the highest sample group
		my $max_sg = 1;
		for my $sample_group_obj ( @sample_group_objs ) {
		    if ($sample_group_obj->sample_group_number > $max_sg) {
		    	$max_sg = $sample_group_obj->sample_group_number;
		    }
		}
		$sample_group_count = $max_sg  + 1;
	} else {
		if (! defined $OPT{sample_group_number}) {
			modules::Exception->throw("ERROR: Need sample_group_number input for this sample");
		} else {
			$sample_group_count = $OPT{sample_group_number};
		}
	}
} elsif ($writedb) {
	
	
	$sample_group_count = 1;
} else {
	#Writing out qsubs in patient; confirm there aren't multiple sample_groups
	my @sample_group_objs = modules::Adaptors::Sample_Group->search(patient_id=>$patient_obj->id);
	if (@sample_group_objs > 1) {
		modules::Exception->throw("ERROR: Multiple sample groups for patient $researcher_id");
	}
	$sample_group_count = 1;
}

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

if ($svndir =~ /trunk/) {
	modules::Exception->throw("ERROR: Do not use trunk for svn directory") unless $OPT{test};
}

my $sequencing_centre = defined $OPT{sequencing_centre}?$OPT{sequencing_centre}:'AGRF';
if ($sequencing_centre ne 'AGRF' && $sequencing_centre ne 'BRF' && $sequencing_centre ne 'RAM') {
	modules::Exception->throw("ERROR: Sequencing centre must be BRF, AGRF, or RAM");
}

my ($sequencing_centre_obj) = modules::Adaptors::Sequencing_Centre->search(sequencing_centre_name=>$sequencing_centre);
if ( !$sequencing_centre_obj ) {
	modules::Exception->throw("ERROR: Can't retrieve db entry for sequencing center $sequencing_centre");
}
my $sequencing_centre_db_id = $sequencing_centre_obj->id;

my @normal_sample_types = split(",",$OPT{normal_sample_types});

if (@normal_sample_types > 1){
	modules::Exception->throw("ERROR: Can't have more than 1 normal sample");
}

my @tumour_sample_types = split(",",$OPT{tumour_sample_types});

my @normal_bpa_ids = split(",",$OPT{normal_bpa_ids});
my @tumour_bpa_ids = split(",",$OPT{tumour_bpa_ids});

my @normal_readdirs = split(",",$OPT{normal_readdirs});
my @tumour_readdirs = split(",",$OPT{tumour_readdirs});

my $normal_num = @normal_readdirs;
my $tumour_num = @tumour_readdirs;

my @tumour_types = split(",",$OPT{tumour_types});
my @tumour_tissues = ();

if (defined $OPT{tumour_tissues}) {
	@tumour_tissues = split(",",$OPT{tumour_tissues});
} else {
	#Default is skin
	@tumour_tissues = (('skin') x $tumour_num);
}

#Check all the argument numbers match
if (@tumour_bpa_ids != @tumour_readdirs || @tumour_sample_types != @tumour_readdirs || @tumour_sample_types != @tumour_types || @tumour_sample_types != @tumour_tissues) {
	modules::Exception->throw("ERROR: Number of tumour arguments must match");
}

if (@normal_sample_types != @normal_readdirs) {
	modules::Exception->throw("ERROR: Number of normal arguments must match");
}



#Check the sample_types are ok
for my $sample_type ( @normal_sample_types ) {
	if (!exists $possible_sample_types{$sample_type}) {
		modules::Exception->throw("ERROR: $sample_type is not a correct sample_type option");
	}    
}
for my $sample_type ( @tumour_sample_types ) {
	if (!exists $possible_sample_types{$sample_type}) {
		modules::Exception->throw("ERROR: $sample_type is not a correct sample_type option");
	}    
}


my $pipe1_start = 'merge_bam';
my $pipe1_end = 'bam_stats';
my $pipe2_start = 'submit_snvs';
my $pipe2_end = 'submit_snvs';
my $pipe3_start = 'merge_vcf';
my $pipe3_end = 'scp_results';
my $pipe3_end_normal = 'filter_vcf'; #Don't record these snvs in the database as they will bog down performance



#Next create the sample group
my $sample_total = @normal_sample_types + @tumour_sample_types;

#This is the first sample_group
my $sample_group_name = $researcher_id.'_sg'.$sample_group_count;


my $samplegroup_db_id;
if ($writedb) {
	my %sample_group_info = (
							total_samples=>$sample_total,
							sample_group_number=>$sample_group_count,
							sample_group_name=>$sample_group_name,
							patient_id=>$patient_db_id
							);
	$samplegroup_db_id = modules::Adaptors::Sample_Group->insert(\%sample_group_info);
	print STDERR "\tCreate sample_group $sample_group_name with id $samplegroup_db_id\n";
} else {
	my ($sample_group_obj) = modules::Adaptors::Sample_Group->search(sample_group_name=>$sample_group_name,patient_id=>$patient_db_id);
	if (!defined $sample_group_obj) {
		modules::Exception->throw("No sample_group object for $sample_group_name and patient id $patient_db_id");
	}
	$samplegroup_db_id = $sample_group_obj->id;
}

#Create the pipeline object
#Get the xml files to get variables we need first
my $cluster_config = modules::Pipeline->get_cluster_conf();
my $pipe_config = modules::Pipeline->get_pipe_conf();

my $steps_xml = defined $OPT{steps_xml}?$OPT{steps_xml}:"$svndir/conf/steps.xml";

if ( !-e $steps_xml ) {
	modules::Exception->throw("File $steps_xml doesn't exist");	
}

#Read the cluster specific variables 
my $qsub_dir = $cluster_config->read('base_directories','base_qsub_directory').'/'.$researcher_id;
system("mkdir $qsub_dir") if (!-d $qsub_dir);
my $read_directory = $cluster_config->read('base_directories','base_read_directory').'/'.$researcher_id;
if ( !-d $read_directory ) {
	modules::Exception->throw("File $read_directory doesn't exist");	
}
my $outdir = $cluster_config->read('base_directories','base_run_directory').'/'.$researcher_id;
system("mkdir -p $outdir") if (!-d $outdir);
my $snv_call_dir = $outdir .'/' .$sample_group_name.'_snvcalls';
system("mkdir -p $snv_call_dir") if (!-d $snv_call_dir);

my $bam_link_dir = $outdir . '/bam_links';
system("mkdir -p $bam_link_dir") if (!-d $bam_link_dir);


my $threadnum = $cluster_config->read('qsub_vars','thread_num');
my $scheduler = $cluster_config->read('scheduler');
my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
my $pipeline = modules::Pipeline->new(-patient_id=>$researcher_id,-cluster_obj=>$cluster_obj);

#Next create the samples

#Acceptable suffices
my %suffices = (
				gz=>1,
				bz2=>1
				);

my %qsub_commands = ();

#Now create the samples..
#First normal
my $normal_sample;
&Create_Samples('normal',$normal_num);
#Next tumour
my @tumour_samples = ();
&Create_Samples('tumour',$tumour_num);

#Create the snv_calling qsubs # germline calls are now obtained as subset of pair calls
$pipeline->create_snvcall_qsubs(-normal_sample=>$normal_sample,-tumour_samples=>\@tumour_samples,-sample_group_name=>$sample_group_name,-no_normal=>1);

if (keys %qsub_commands) {
	print STDERR "Run the following commands...\n";
	for my $command ( keys %qsub_commands ) {
	    print "\t$command\n";
	}
	print "\n";
}

#Does all the work creating samples, lanes, and read_files
sub Create_Samples () {
	my ($tumour_normal,$sample_total) = @_;	
	my $tumour_flag = $tumour_normal eq 'tumour'?1:0;


	for ( my $sample_count = 0 ; $sample_count < $sample_total ; $sample_count++ ) {
	    my $cell_line = 0;
	    if ($tumour_flag) {
	    	if ($tumour_sample_types[$sample_count] =~ /cell_line/) {
		    	$cell_line = 1;
		    }
	    } else {
		    if ($normal_sample_types[$sample_count] =~ /cell_line/) {
		    	$cell_line = 1;
		    }
	    }
		my $sample_number = $sample_count+1;
	
		my $sample_name = $sample_group_name.'_'.$tumour_normal.$sample_number;
	
		#Figure out how many lanes there are
		my $lane_data;
		
		#Get the local sample_directory for getting reads later
		my $sample_readdir;
			
		
		
		if ($tumour_flag) {
			push @tumour_samples, $sample_name;
			$sample_readdir = $tumour_readdirs[$sample_count];
			($lane_data) = &Get_Lane_Info($tumour_readdirs[$sample_count]);
		} else {
			$normal_sample = $sample_name;
			$sample_readdir = $normal_readdirs[$sample_count];
			($lane_data) = &Get_Lane_Info($normal_readdirs[$sample_count]);
		}
	    my $lane_total = keys %{$lane_data};
	    	    
	    	    
	    
	    my $sample_db_id;
	    if ($writedb) {
		    my %sample_info = (
		    					sample_number=>$sample_number,
		    					tumour => $tumour_flag,
		    					cell_line => $cell_line,
		    					total_lanes => $lane_total,
		    					sample_name => $sample_name,
		    					sample_group_id=>$samplegroup_db_id 					
		    					);
		    
		    
		    if ($tumour_flag) {
				$sample_info{tissue}= $tumour_tissues[$sample_count];
				$sample_info{tumour_type} = $tumour_types[$sample_count];
		    	$sample_info{sample_type} = $tumour_sample_types[$sample_count];
		    	$sample_info{bpa_id} = $tumour_bpa_ids[$sample_count];	    						    	
		    } else {
		    	$sample_info{tissue}='blood';
		    	$sample_info{sample_type} = $normal_sample_types[$sample_count];
		    	$sample_info{bpa_id} = $normal_bpa_ids[$sample_count];	    					
		    }
	    	$sample_db_id =  modules::Adaptors::Sample->insert(\%sample_info);
		    print STDERR "\t\tCreate sample $sample_name with id $sample_db_id\n";   
		    
		    
		     	
	    } else {
	    	#Get the existing sample
	    	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name,sample_group_id=>$samplegroup_db_id);
	    	print "Sample $sample_name sample_groupdb $samplegroup_db_id\n";
	    	$sample_db_id = $sample_obj->id;
	    }				
	    
	    
  
	    
	    my @bams = ();			
	    my $lanes_outdir = $outdir.'/'.$sample_name.'_lanes';
	    my $runs_outdir = $outdir.'/'.$sample_name.'_runs';
	    
	    my $encoding;
		#Create the lanes for each sample
	    for ( my $lane_count = 0 ; $lane_count < $lane_total ; $lane_count++ ) {
	    	my $lane_number = $lane_count;
	    	$lane_number++;
	    	my $lane_name = $sample_name.'_l'.$lane_number;
	        
	        #Get the read_file info
	        my $read1_compressed = 0;
	        my $read2_compressed = 0;
	        my $read1_suffix;
	        my $read2_suffix;
			my $read1_name = $lane_name.'_r1';
			my $read2_name = $lane_name.'_r2';
			
			
			my $read1_actual = $lane_data->{$lane_count}{read1};
			my $read2_actual = $lane_data->{$lane_count}{read2};
			
			my $read1_full_file = $read_directory .'/'. $sample_readdir. '/'. $lane_data->{$lane_count}{read1};
			my $read2_full_file = $read_directory .'/' . $sample_readdir. '/'. $lane_data->{$lane_count}{read2};
			my $read1_symlink = $read_directory .'/' . $sample_readdir. '/'. $read1_name;
			my $read2_symlink = $read_directory .'/' . $sample_readdir. '/'. $read2_name;
			
			
			my $read_directory_full = $read_directory .'/'. $sample_readdir;
			
			
			#First check the read files exist
			if ( !-e $read1_full_file ) {
				modules::Exception->throw("File $read1_full_file doesn't exist");	
			}	
			if ( !-e $read2_full_file ) {
				modules::Exception->throw("File $read2_full_file doesn't exist");	
			}	
			
			#Create symlinks to the files if they don't exist
			if (!-e $read1_symlink) {
				system("cd $read_directory/$sample_readdir; ln -s $read1_actual $read1_name");
			}
			
			if (!-e $read2_symlink) {
				system("cd $read_directory/$sample_readdir; ln -s $read2_actual $read2_name");
			}
	
			my %commands = ();
	
			#Check if data is compressed        
	        for my $suffix ( keys %suffices ) {
	            if ($lane_data->{$lane_count}{read1} =~ /$suffix$/) {
	            	$read1_compressed = 1;
	            	$read1_suffix = $suffix;
	            	#Add decompression if required
	            	if ($suffix eq 'bz2') {
						push @{$commands{bwa}},"$svndir/scripts/compress.pl -keep -suffix $suffix -threadNum $threadnum -files $read1_full_file";
					}

	            }
	            
	            if ($lane_data->{$lane_count}{read2} =~ /$suffix$/) {
	            	$read2_compressed = 1;
	            	$read2_suffix = $suffix;
	            	if ($suffix eq 'bz2') {
						push @{$commands{bwa}}, "$svndir/scripts/compress.pl -keep -suffix $suffix -threadNum $threadnum -files $read2_full_file";
					}
	            }
	        }
	        
	        
			
			
			#Create the aligning lane qsub files
			my $lane_bam;
			system("mkdir -p $lanes_outdir") if (!-d $lanes_outdir);
			
			#By default don't check the quality (assume phred33)
			my $quality = defined $OPT{skip_quality}?1:0;
			
	        
	        
	        my $lane_db_id;
	        if ($writedb) {
		        my %lane_info = (
	        				lane_number=>$lane_number,
	        				lane_name=>$lane_name,
	        				sequencing_centre_id=>$sequencing_centre_db_id,
	        				sample_id=>$sample_db_id
	        				);
	        	$lane_db_id = modules::Adaptors::Lane->insert(\%lane_info);
	        	print STDERR "\t\t\tCreate lane $lane_name with id $lane_db_id\n";
	        	
	        	#Create the read_files
				my %read1_info = (
								file_name => $lane_data->{$lane_count}{read1},
								is_compressed=>$read1_compressed,
								compression_suffix=>$read1_suffix,
								read_file_number=>1,
								read_directory=>$read_directory_full,
								read_file_name=>$read1_name,
								lane_id=>$lane_db_id
								);
								
				my %read2_info = (
								file_name => $lane_data->{$lane_count}{read2},
								is_compressed=>$read2_compressed,
								compression_suffix=>$read2_suffix,
								read_file_number=>2,
								read_directory=>$read_directory_full,
								read_file_name=>$read2_name,
								lane_id=>$lane_db_id
								);        
			
				modules::Adaptors::Read_File->insert(\%read1_info);
				modules::Adaptors::Read_File->insert(\%read2_info);
			
				print STDERR "\t\t\t\tCreate two read_files $read1_name and $read2_name\n";
	        	
	        	
	        } else {
	        	my ($lane_obj) = modules::Adaptors::Lane->search(sample_id=>$sample_db_id,lane_name=>$lane_name);
	        	$lane_db_id = $lane_obj->id;
	        }
	        
			if (keys %commands) {
				($encoding,$lane_bam) = $pipeline->align_lane_qsub(-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id,-commands=>\%commands);
			} else {
				($encoding,$lane_bam) = $pipeline->align_lane_qsub(-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id);
			}
			push @bams, $lane_bam;
			print "\n";
	    }
	    
	    $qsub_commands{"cd $qsub_dir/lanes/; for f in $sample_name*wrapper.qsub; do bash \$f; done"}++;
		system("mkdir -p $runs_outdir") if (!-d $runs_outdir);
	    
	   	#Now create the config xml for the runs
	    my $steps_xml_dir = $qsub_dir . '/runs/';
	    if (!-d $steps_xml_dir) {
	    	system("mkdir -p $steps_xml_dir");
	    }
	    
	    my $steps_xml_out = $steps_xml_dir.$sample_name.'.xml';
	    #Assumes same encoding for a sample but likely ok
		$pipeline->create_sample_xml(-sample_name=>$sample_name,-steps_xml_template=>$steps_xml,-xml_out=>$steps_xml_out,-encoding=>$encoding,-bams=>\@bams);

	    #Now create the run qsubs
		$pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe1_start,-end_step=>$pipe1_end,-pipe_block=>1);	  
	    
	    
	    #Generate the snv processing part of the pipeline run
	    if ($tumour_flag == 1) {
		    $pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe2_start,-end_step=>$pipe2_end,-pipe_block=>2);
		    $pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe3_start,-end_step=>$pipe3_end,-pipe_block=>3);
	    } else {
	    	#germline calls done simultaneously
	    	#$pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe2_start,-end_step=>$pipe2_end,-pipe_block=>2);
	    	#$pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe3_start,-end_step=>$pipe3_end_normal,-pipe_block=>3);
	    }
	}
}

#Get lane info from a directory
sub Get_Lane_Info {
	my ($readdir) = @_;
	my %lane_info;
	my $read1_regex = my $read2_regex;

	if (defined $OPT{read1_pattern}) {
		$read1_regex = $OPT{read1_pattern};
		($read2_regex = $read1_regex) =~ s/1/2/;
	} else {
		$read1_regex = 'R1';
		$read2_regex = 'R2';
	}
	
	opendir(DIR,"$read_directory/$readdir") || modules::Exception->throw("ERROR: Can't open readdir $read_directory/$readdir $!\n");
	my @files = readdir DIR;
	my @reads1_all = sort grep {/$read1_regex/} @files;
	my @reads2_all = sort grep {/$read2_regex/} @files;
	closedir(DIR);
	
	if (!@reads1_all) {
		modules::Exception->throw("ERROR: No reads match regex $read1_regex");
	} elsif (!@reads2_all) {
		modules::Exception->throw("ERROR: No reads match regex $read2_regex");
	} elsif (@reads1_all != @reads2_all) {
		modules::Exception->throw("ERROR: Read numbers don't match");
	}
	
	for ( my $count = 0 ; $count < @reads1_all; $count++ ) {
	    $lane_info{$count}{read1} = $reads1_all[$count];
	    $lane_info{$count}{read2} = $reads2_all[$count];
	}
	return (\%lane_info);
}









