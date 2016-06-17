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
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"sample_group_name=s",
	   	"tumour_bpa_id=s",
   		"tumour_readdir=s",
	   	"tumour_sample_type=s",
	   	"sequencing_centre=s",
	   	"tumour_type=s",
	   	"tumour_tissue=s",
	   	"read1_pattern=s",
	   	"cluster_xml=s",
	   	"steps_xml=s",
	   	"skip_quality",
	   	"test"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{tumour_bpa_id} || !$OPT{tumour_readdir} || !$OPT{tumour_sample_type} || !$OPT{tumour_type});

	   
=pod

=head1 SYNOPSIS

add_sample_to_samplegroup.pl -tumour_bpa_id bpa_barcode_tumour -sample_group_name sample_group_name(e.g. A15_sg1) -read1_pattern read1_regex_pattern(default=R1) -steps_xml xml_steps_file(default=../conf/steps/xml) -cluster_xml xml_cluster_file(default=../conf/cluster.xml) -sequencing_centre sequencing_centre(option=BRF,AGRF,RAM;default=AGRF) -tumour_sample_type comma_delim_sample_type -tumour_type comma_delim_tumour_types(options=primary,metatstatic) -tumour_readdir comma_delim_fq_directories -skip_quality Skip_read_check_and_assume_phred33_qualities [options]

Required flags: -tumour_bpa_id -tumour_readdir -tumour_sample_type -tumour_type

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

add_sample_to_samplegroup.pl -> Add sample to existing sample_group

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

add_sample_to_samplegroup.pl 

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


my $tumour_sample_type = $OPT{tumour_sample_type};

my $tumour_bpa_id = $OPT{tumour_bpa_id};

my $tumour_readdir = $OPT{tumour_readdir};

my $tumour_type = $OPT{tumour_type};
my $tumour_tissue;

if (defined $OPT{tumour_tissue}) {
	$tumour_tissue = $OPT{tumour_tissue};
} else {
	#Default is skin
	$tumour_tissue = 'skin';
}

if (!exists $possible_sample_types{$tumour_sample_type}) {
	modules::Exception->throw("ERROR: $tumour_sample_type is not a correct sample_type option");
}    


#Get the arguments
my $sample_group_name = $OPT{sample_group_name};


my ($samplegroup_obj) = modules::Adaptors::Sample_Group->search(sample_group_name=>$sample_group_name);

if (!defined $samplegroup_obj) {
	modules::Exception->throw("ERROR: There is no sample group for sample_group $sample_group_name");
}

my $samplegroup_db_id = $samplegroup_obj->id;

my $pipe1_start = 'merge_bam';
my $pipe1_end = 'bam_stats';
my $pipe2_start = 'submit_snvs';
my $pipe2_end = 'submit_snvs';
my $pipe3_start = 'merge_vcf';
my $pipe3_end = 'scp_results';


#Create the pipeline object
#Get the xml files to get variables we need first
my $cluster_config = modules::Pipeline->get_cluster_conf();
my $pipe_config = modules::Pipeline->get_pipe_conf();

my $steps_xml = defined $OPT{steps_xml}?$OPT{steps_xml}:"$svndir/conf/steps.xml";

if ( !-e $steps_xml ) {
	modules::Exception->throw("File $steps_xml doesn't exist");	
}

#Read the cluster specific variables 
(my $patient_name = $sample_group_name) =~ s/_sg[0-9]//;
my $qsub_dir = $cluster_config->read('base_directories','base_qsub_directory').'/'.$patient_name;
system("mkdir $qsub_dir") if (!-d $qsub_dir);
my $read_directory = $cluster_config->read('base_directories','base_read_directory').'/'.$patient_name;
if ( !-d $read_directory ) {
	modules::Exception->throw("File $read_directory doesn't exist");	
}
my $outdir = $cluster_config->read('base_directories','base_run_directory').'/'.$patient_name;
system("mkdir -p $outdir") if (!-d $outdir);
my $snv_call_dir = $outdir .'/' .$sample_group_name.'_snvcalls';
system("mkdir -p $snv_call_dir") if (!-d $snv_call_dir);

my $bam_link_dir = $outdir . '/bam_links';
system("mkdir -p $bam_link_dir") if (!-d $bam_link_dir);


my $threadnum = $cluster_config->read('qsub_vars','thread_num');
my $scheduler = $cluster_config->read('scheduler');
my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
my $pipeline = modules::Pipeline->new(-patient_id=>$patient_name,-cluster_obj=>$cluster_obj);

#Next create the samples

#Acceptable suffices
my %suffices = (
				gz=>1,
				bz2=>1
				);

my %qsub_commands = ();

my @sample_objs = modules::Adaptors::Sample->search(sample_group_id=>$samplegroup_obj->id);
if (!@sample_objs) {
	modules::Exception->throw("ERROR: No existing samples for $sample_group_name");
}
#Get the highest sample group
my $max_sample = 1;
for my $sample_obj ( @sample_objs ) {
    if ($sample_obj->sample_number > $max_sample) {
    	$max_sample = $sample_obj->sample_number;
    }
}
my $sample_number = $max_sample  + 1;
my $new_sample_name = $sample_group_name . '_tumour'.$sample_number;

#Check this sample doesn't already belong to any other sample_group
my $sample_out = `ls $read_directory/$tumour_readdir/*r1 2>/dev/null`;
if ($sample_out) {
	modules::Exception->throw("ERROR: Symlinks in $read_directory/$tumour_readdir already exist indicating sample already part of sample group");
} else {
	#Now create the sample..
	&Create_Sample($new_sample_name);
}



my $normal_sample = $sample_group_name . '_normal1';

#Create the snv_calling qsubs
my @tumour_samples = ($new_sample_name);
$pipeline->create_snvcall_qsubs(-normal_sample=>$normal_sample,-tumour_samples=>\@tumour_samples,-sample_group_name=>$sample_group_name,-no_normal=>1);

#Increase the sample group count
#Finally update the sample to add a new lane
my $current_samples = $samplegroup_obj->total_samples;
my $new_total_samples = $current_samples + 1;
$samplegroup_obj->total_samples($new_total_samples);
$samplegroup_obj->update();

if (keys %qsub_commands) {
	print STDERR "Run the following commands...\n";
	for my $command ( keys %qsub_commands ) {
	    print "\t$command\n";
	}
	print "\n";
}

#Does all the work creating samples, lanes, and read_files
sub Create_Sample () {
		my ($sample_name) = @_;
	    my $cell_line = $tumour_sample_type =~ /cell_line/?1:0;
	
		#Figure out how many lanes there are
		my $lane_data;
		
		#Get the local sample_directory for getting reads later
		my $sample_readdir = $tumour_readdir;
			
		($lane_data) = &Get_Lane_Info($tumour_readdir);
	
		
	    my $lane_total = keys %{$lane_data};
	    	    
	    	    
	    my %sample_info = (
	    					sample_number=>$sample_number,
	    					tumour => 1,
	    					cell_line => $cell_line,
	    					total_lanes => $lane_total,
	    					sample_name => $sample_name,
	    					sample_group_id=>$samplegroup_db_id,
	    					tissue=>$tumour_tissue,
	    					tumour_type=>$tumour_type,
	    					sample_type=>$tumour_sample_type,
	    					bpa_id=>$tumour_bpa_id				
	    					);
	       				
	    my $sample_db_id =  modules::Adaptors::Sample->insert(\%sample_info);
	    
	    print STDERR "\t\tCreate sample $sample_name with id $sample_db_id\n";
	    
  
	    
	    my @bams = ();			
	    my $lanes_outdir = $outdir.'/'.$sample_name.'_lanes';
	    my $runs_outdir = $outdir.'/'.$sample_name.'_runs';
	    
	    my $encoding;
		#Create the lanes for each sample
	    for ( my $lane_count = 0 ; $lane_count < $lane_total ; $lane_count++ ) {
	    	my $lane_number = $lane_count;
	    	$lane_number++;
	    	my $lane_name = $sample_name.'_l'.$lane_number;
	        
	        my %lane_info = (
	        				lane_number=>$lane_number,
	        				lane_name=>$lane_name,
	        				sequencing_centre_id=>$sequencing_centre_db_id,
	        				sample_id=>$sample_db_id
	        				);
	        
	        
	        my $lane_db_id = modules::Adaptors::Lane->insert(\%lane_info);
	        
	        print STDERR "\t\t\tCreate lane $lane_name with id $lane_db_id\n";
	        
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
			
			
			#Create the aligning lane qsub files
			my $lane_bam;
			system("mkdir -p $lanes_outdir") if (!-d $lanes_outdir);
			
			#By default don't check the quality (assume phred33)
			my $quality = defined $OPT{skip_quality}?1:0;
			
			if (keys %commands) {
				($encoding,$lane_bam) = $pipeline->align_lane_qsub(-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id,-commands=>\%commands);
			} else {
				($encoding,$lane_bam) = $pipeline->align_lane_qsub(-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id);
			}
			push @bams, $lane_bam;
			print "\n";
	    
	    
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
	    
	    #Generate the snv calling part of the pipeline run
	    $pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe2_start,-end_step=>$pipe2_end,-pipe_block=>2);
	    
	    #Generate the snv processing part of the pipeline run
		$pipeline->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$pipe3_start,-end_step=>$pipe3_end,-pipe_block=>3);
	   
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


