#! /usr/bin/perl -w

use strict;
use modules::Overlap;
use modules::Exception;
use modules::ConfigXML;
use modules::SystemCall;
use modules::Vcf;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "sample_name=s",
		   "output=s",
		   "outvcf=s",
		   "ref=s",
		   "bam=s",
		   "rmbam=s",
		   "no_cleanup_files"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sample_name} || !$OPT{output}  || !$OPT{ref} || !$OPT{bam} || !$OPT{outvcf} || !$OPT{rmbam});

	   
=pod

=head1 SYNOPSIS

merge_vcf.pl -sample_name sample_name -output output_file -ref fasta_ref -bam input_bam -rmbam merge_bam_to_replace_with_symlink -outvcf outvcf_file -no_cleanup_file leave_sam_sai_files(default=delete) [options]

Required flags: -sample_name -output -ref -bam -outvcf

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

merge_vcf.pl -> Script to create single vcf file 

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./merge_vcf.pl  

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $outfile = $OPT{output};
my $sys_call = modules::SystemCall->new();

my $bam = $OPT{bam};
if ( !-e $bam ) {
	modules::Exception->throw("File $bam doesn't exist");	
}

my $ref = $OPT{ref};
if ( !-e $ref ) {
	modules::Exception->throw("File $ref doesn't exist");	
}

my $sample_name = $OPT{sample_name};


#Get the xml files and create the pipeline object
my $pipe_config = modules::Pipeline->get_pipe_conf();
my $cluster_config = modules::Pipeline->get_cluster_conf();
my $samtools_bin = $pipe_config->read('binaries','samtools','binary');

#Get the normal bam file from the 
my $rundir = $cluster_config->read('base_directories','base_run_directory');
#We want the pileup from the normal bam as well

my $tumour_flag = modules::Pipeline->get_tumour_flag(-sample_name=>$sample_name);

my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);
my $sample_type;
if (!$sample_obj) {
	modules::Exception->throw("ERROR: Cannot retrieve sample object for $sample_name");
}
$sample_type = $sample_obj->sample_type;

my ($patient_name) = modules::Pipeline::get_patient_id($sample_name);
my ($sample_group_name) = modules::Pipeline::get_sample_group_name($sample_name);

my $normal_bam;

#Get the bam file for mpileup and remove the merge_bam files and replace with symlinks
if ($tumour_flag) {
	#Need the normal pileup as well for tumour samples
	my $bam_norm = $sample_group_name . '_normal1.bam';
	$normal_bam = $rundir . '/' . $patient_name . '/bam_links/'. $bam_norm;
	if (!-e $normal_bam ) {
		modules::Exception->throw("File $normal_bam doesn't exist");	
	}
}

#Remove the merge_bam files from the run directories and replace them with symlinks to the results files
my $destdir = $cluster_config->read('base_directories','base_results_directory') . '/'. $patient_name;

#Cleanup up the tumour merge_bam files and replace with symlinks; clean up lanes .sam and .sai files 
unless ($OPT{no_cleanup_files}) {
	my $rmbam = $OPT{rmbam};
	my $rmbam_index = $rmbam .'.bai';
	if ( !-e $rmbam ) {
		modules::Exception->throw("File $rmbam doesn't exist");	
	}
	
	my $rmbam_just_file = basename($rmbam);
	
	my $bam_mass = $destdir . '/' . $rmbam_just_file;
	my $bam_mass_index = $bam_mass . '.bai';
	
	if ( !-e $bam_mass ) {
		modules::Exception->throw("File $bam_mass doesn't exist");	
	}
	
	$sys_call->run("rm $rmbam; rm $rmbam_index");
	$sys_call->run("ln -s $bam_mass $rmbam; ln -s $bam_mass_index $rmbam_index");

	my $lane_dir = $rundir . '/' . $patient_name . '/' . $sample_name . '_lanes';
	opendir(LANES,$lane_dir) || modules::Exception->throw("ERROR: Can't open lanes dir $lane_dir");
	my @delete_files = grep {/sa[im]$/} readdir LANES;
	closedir LANES;
	
	for my $file ( @delete_files ) {
	    my $full_delete_file = $lane_dir . '/' . $file;
	    if ( !-e $full_delete_file ) {
	    	modules::Exception->throw("File $full_delete_file doesn't exist");	
	    }
	    $sys_call->run("rm $full_delete_file") unless $OPT{no_cleanup_files};
	}


	#Remove the merge_bam files and replace with symlinks
	if ($tumour_flag) {
		#Clean up the normal bam files as well
		#Get the normal bam merge_bam step file by reading the symlink and changing the rmdup file to merge_bam file
		my $rmbam_normal = readlink($normal_bam);
		$rmbam_normal =~ s/remove_duplicates/merge_bam/;
		my $rmbam_normal_index = $rmbam_normal .'.bai';
		if ( !-e $rmbam_normal ) {
			modules::Exception->throw("File $rmbam_normal doesn't exist");	
		}
		
	
		my $rmbam_normal_just_file = basename($rmbam_normal);
		
		my $bam_normal_mass = $destdir . '/' . $rmbam_normal_just_file;
		my $bam_normal_mass_index = $bam_normal_mass . '.bai';
		if ( !-e $bam_normal_mass ) {
			modules::Exception->throw("File $bam_normal_mass doesn't exist");	
		}
		
		$sys_call->run("rm $rmbam_normal; rm $rmbam_normal_index");
		$sys_call->run("ln -s $bam_normal_mass $rmbam_normal; ln -s $bam_normal_mass_index $rmbam_normal_index");
		
		#Clean up the normal lanes as well
		my $normal_lane_dir = $rundir . '/' . $patient_name . '/' . $sample_group_name . '_normal1_lanes';
		opendir(NLANES,$normal_lane_dir) || modules::Exception->throw("ERROR: Can't open lanes dir $normal_lane_dir");
		my @normal_delete_files = grep {/sa[im]$/} readdir NLANES;
		closedir NLANES;
	
		for my $file ( @normal_delete_files ) {
	    	my $full_delete_file = $normal_lane_dir . '/' . $file;
	    	if ( !-e $full_delete_file ) {
	    		modules::Exception->throw("File $full_delete_file doesn't exist");	
	    	}
	    	$sys_call->run("rm $full_delete_file") unless $OPT{no_cleanup_files};
		}	
		
	} 
}

my %chromosomes = map {$_ => 1} split(" ",$pipe_config->read('annotation_version','chr'));

#run the commands to generate the vcf
my $vcf = modules::Vcf->new(-sample_type=>$sample_type);
my $all_vcf_file = $OPT{outvcf};

#Remove file if exists as we're concatonating
if (-e $all_vcf_file) {
	$sys_call->run("rm $all_vcf_file");
}

(my $snv_vcf_file = $all_vcf_file) =~ s/all/snv/;
(my $indel_vcf_file = $all_vcf_file) =~ s/all/indel/;

my $snv_out = $outfile.'.snv';
my $indel_out = $outfile .'.indel';


my $outdir = $cluster_config->read('base_directories','base_run_directory');

my $vcf_dir = $outdir.'/'.$patient_name.'/'.$sample_group_name.'_snvcalls';
if (!-d $vcf_dir) {
	modules::Exception->throw("ERROR: Vcf dir $vcf_dir doesn't exist");
}

my @single_vcfs = ();
for my $chr (sort keys %chromosomes) {
	my $vcf_file = $sample_name . '.'. $chr . '.vcf';
	my $full_vcf = $vcf_dir . '/' . $vcf_file;
	if ( !-e $full_vcf ) {
		modules::Exception->throw("File $full_vcf doesn't exist");	
	}
	push @single_vcfs, $full_vcf;
}

my $file_str = join(" ",@single_vcfs);
my $merge_command = "cat $file_str >> $all_vcf_file";
print STDERR "$merge_command\n";
$sys_call->run($merge_command);

#Create vcfs for the different variant types
my $snv_command = "grep -v INDEL $all_vcf_file > $snv_vcf_file";
print STDERR "$snv_command\n";
$sys_call->run($snv_command);

my $indel_command = "grep INDEL $all_vcf_file > $indel_vcf_file";
print STDERR "$indel_command\n";
$sys_call->run($indel_command);

#Parse the unfiltered vcf
$vcf->parse_vcf(-vcf_file => $all_vcf_file);
my $vcf_data = $vcf->get_vcf(-vcf_file => $all_vcf_file);



my $pileup_indel_file = $indel_out . '.pileupcoord';
my $pileup_snv_file = $snv_out . '.pileupcoord';

#Open the coord file for writing; needed for hom/het calls later
open(SNVCOORD,">$pileup_snv_file") || modules::Exception->throw("Can't open file to write $pileup_snv_file\n");
if (exists $vcf_data->{SNV}) {

	for my $chr (sort keys %{$vcf_data->{SNV}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{SNV}{$chr}}) {
			print SNVCOORD "$chr\t$start_coord\n";
		}
	}
}

close SNVCOORD;


#Open the coord file for writing; needed for hom/het calls later
open(INDELCOORD,">$pileup_indel_file") || modules::Exception->throw("Can't open file to write $pileup_indel_file\n");
if (exists $vcf_data->{DEL}) {

	for my $chr (sort keys %{$vcf_data->{DEL}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{DEL}{$chr}}) {
			$start_coord--;
			print INDELCOORD "$chr\t$start_coord\n";
		}
	}
}

if (exists $vcf_data->{INS}) {

	for my $chr (sort keys %{$vcf_data->{INS}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{INS}{$chr}}) {
			$start_coord--;
			print INDELCOORD "$chr\t$start_coord\n";
		}
	}
}

close INDELCOORD;

#Now generate the mpileup files
my $mpileup_snv_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_snv_file $bam  > $snv_out.pileup";
print STDERR "$mpileup_snv_command\n";
$sys_call->run($mpileup_snv_command);
	
my $mpileup_indel_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_indel_file $bam  > $indel_out.pileup";
print STDERR "$mpileup_indel_command\n";
$sys_call->run($mpileup_indel_command);

#Need normal mpileups for each tumour
if ($tumour_flag) {
	my $mpileup_normal_snv_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_snv_file $normal_bam  > $snv_out.normal.pileup";
	print STDERR "$mpileup_normal_snv_command\n";
	$sys_call->run($mpileup_normal_snv_command);
	
	my $mpileup_normal_indel_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_indel_file $normal_bam  > $indel_out.normal.pileup";
	print STDERR "$mpileup_normal_indel_command\n";
	$sys_call->run($mpileup_normal_indel_command);
} 


#Get the depths from the pileup; this is needed for filtering in the next step
my %snv_depths = ();
my %indel_depths = ();

open(SNVPILEUP,"$snv_out.pileup") || modules::Exception->throw("ERROR: Can't open file $snv_out.pileup");
while (<SNVPILEUP>) {
	my @fields = split("\t");
	$snv_depths{$fields[0]}{$fields[1]} = $fields[3];
}

open(INDELPILEUP,"$indel_out.pileup") || modules::Exception->throw("ERROR: Can't open file $indel_out.pileup");
while (<INDELPILEUP>) {
	my @fields = split("\t");
	my $indel_start = $fields[1]  + 1;
	$indel_depths{$fields[0]}{$indel_start} = $fields[3];
}

open(INDEL,">$indel_out") || modules::Exception->throw("Can't open file to write $indel_out\n");
open(SNV,">$snv_out") || modules::Exception->throw("Can't open file to write $snv_out\n");


#Write out unfiltered variants to text files
if (exists $vcf_data->{DEL}) {
	for my $chr (sort keys %{$vcf_data->{DEL}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{DEL}{$chr}}) {
			for my $end_coord (keys %{$vcf_data->{DEL}{$chr}{$start_coord}}) {
				for my $bases (keys %{$vcf_data->{DEL}{$chr}{$start_coord}{$end_coord}}) {
					my $depth = 0;
					#If not defined then tumour has no reads
					if (defined $indel_depths{$chr}{$start_coord}) {
						$depth = $indel_depths{$chr}{$start_coord};
					} 
					
					my $rest = 'DEL^^^' . $bases . '^^^' .$vcf_data->{DEL}{$chr}{$start_coord}{$end_coord}{$bases}.'^^^D'.$depth;
					print INDEL join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							$rest
		  							) ."\n";
				}
			}
		}
	}
}

if (exists $vcf_data->{INS}) {
	for my $chr (sort keys %{$vcf_data->{INS}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{INS}{$chr}}) {
			for my $end_coord (keys %{$vcf_data->{INS}{$chr}{$start_coord}}) {
				for my $bases (keys %{$vcf_data->{INS}{$chr}{$start_coord}{$end_coord}}) {
					my $depth = 0;
					if (defined $indel_depths{$chr}{$start_coord}) {
						$depth = $indel_depths{$chr}{$start_coord};
					} 
					my $rest = 'INS^^^' . $bases . '^^^' .$vcf_data->{INS}{$chr}{$start_coord}{$end_coord}{$bases}.'^^^D'.$depth;
					print INDEL join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							$rest
		  							) ."\n";
				}
			}
		}
	}
}

close INDEL;

if (exists $vcf_data->{SNV}) {
	for my $chr (sort keys %{$vcf_data->{SNV}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{SNV}{$chr}}) {
			for my $end_coord (keys %{$vcf_data->{SNV}{$chr}{$start_coord}}) {
				for my $base_change (keys %{$vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}}) {
					my $depth = 0;
					if (defined $snv_depths{$chr}{$start_coord}) {
						$depth = $snv_depths{$chr}{$start_coord};
					} 
					my $rest = 'SNV^^^'.$base_change .'^^^' . $vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}{$base_change} . '^^^D'.$depth;
					print SNV join("\t",
									$chr,
									$start_coord,
									$end_coord,
									$rest
									)."\n";
				}
			}
		}
	}
}

close SNV;


#Don't need this file
$sys_call->run("rm $pileup_indel_file");
$sys_call->run("rm $pileup_snv_file");


