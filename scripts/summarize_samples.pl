#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;
use modules::Pipeline;
use modules::Adaptors::Sample;
use modules::Adaptors::Patient;
use modules::Adaptors::Sample_Group;
use modules::Adaptors::Release_File;
use modules::Adaptors::Lane;
use modules::Adaptors::Sequencing_Centre;
use modules::Utils;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "exclude_sample=s",
	   "include_sample=s",
	   "sample_type=s",
	   "single_sample_group",
	   "fail",
	   "debug"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

summarize_samples.pl -include_sample comma_delim_list_of_samples_to_include(default=all) -exclude_sample comma_delim_list_of_samples_to_ignore(default=none) -fail include_fail (default=report_passed_only) -sample_type only_report_this_sample_type(default=all) -single_sample_group only_report_first_sample_group [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

summarize_samples.pl -> generate summary for multiple samples using files

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./summarize_samples.pl 

=cut



if ($OPT{exclude_sample} && $OPT{include_sample}) {
	modules::Exception->throw("Can't call with -include_sample and -exclude_sample; use one or the other");
}
my %excludes = map {$_=>1} split(",",$OPT{exclude_sample}) if $OPT{exclude_sample};
my %includes = map {$_=>1} split(",",$OPT{include_sample}) if $OPT{include_sample};
my $include_fails = defined $OPT{fail} ? 1 : 0;

my $clus_conf = modules::Pipeline->get_cluster_conf();
my $report_conf = modules::Pipeline->get_report_conf();


my $resultdir = $clus_conf->read('base_directories','base_results_directory');
my $rundir = $clus_conf->read('base_directories','base_run_directory');


if (!-d $resultdir) {
	IGLPipeline::Exception->throw("Directory $resultdir problem");
}

opendir(DIR,$resultdir) || modules::Exception->throw("Can't open directory $resultdir");
my @patients = grep {/\w/} readdir DIR;

my %sample_info = ();

my %required_coverages = (
							'normal_cell_line'=>40,
							'normal_patient'=>40,
							'tumour_primary_cell_line'=>40,
							'tumour_primary_patient'=>60,
							'tumour_metastatic_cell_line'=>40,
							'tumour_metastatic_patient'=>60
							);

my ($date_stamp) = split('_',modules::Utils->GetTime());
my $coverage_file = "coverage_".$date_stamp.'.tsv';
my $snv_file = "snv_summary_".$date_stamp.'.tsv';
my $indel_file = "indel_summary_".$date_stamp.'.tsv';
open(COVER,">$coverage_file") || modules::Exception->throw("Can't open file to write $coverage_file\n");
my %reported_samples = ();


#Get the sample_types from the database; this gives us the columns
for my $patient (@patients) {
	if (keys %excludes) {
		next if exists $excludes{$patient};
	} elsif (keys %includes) {
		next unless $includes{$patient};
	}
	my ($patient_obj) = modules::Adaptors::Patient->search(patient_external_name=>$patient);
	if (! defined $patient_obj) {
		next;
		#modules::Exception->throw("ERROR: Can't retrieve patient with name $patient");
	}
	my @sample_group_objs = modules::Adaptors::Sample_Group->search(patient_id=>$patient_obj->id);
	
	if (! @sample_group_objs ) {
		modules::Exception->throw("ERROR: Can't retrieve sample_group with patient id $patient_obj->id");
	}
	
	for my $sample_group_obj (@sample_group_objs) {
		if ($OPT{single_sample_group}) {
			next unless $sample_group_obj->sample_group_name =~ /_sg1/;
		}
		
		my @sample_objs = modules::Adaptors::Sample->search(sample_group_id=>$sample_group_obj->id);
		if (! @sample_objs ) {
			modules::Exception->throw("ERROR: Can't retrieve sample with sample_group id $sample_group_obj->id");
		}
		
		for my $sample_obj (@sample_objs) {
			#Get the sequencing centre from a lane obj
			my ($lane_obj) = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);
			my ($seq_centre_obj) = modules::Adaptors::Sequencing_Centre->search_sequencing_centre($lane_obj->id);
			my $seq_centre_name = $seq_centre_obj->sequencing_centre_name;
			
			
			my $sample_name = $sample_obj->sample_name;
			if ($OPT{sample_type}) {
				next unless $OPT{sample_type} eq $sample_obj->sample_type;
			}
			my @release_file_objs = modules::Adaptors::Release_File->search_release_files_sample($sample_obj->id);
			if (! @release_file_objs ) {
				next;
			}
			for my $release_file_obj (@release_file_objs) {
				my ($sample_and_run) = split('\.',$release_file_obj->file_name);
				my $local_rundir = $rundir  . '/' . $patient . '/' . $sample_name . '_runs/' . $sample_and_run;

				if ($release_file_obj->file_name =~ /summary.tsv$/) {
					$sample_info{$patient}{$sample_name}{type} = $sample_obj->sample_type;
					my $full_file = $resultdir . '/' . $patient . '/' . $release_file_obj->file_name;
					if ( !-e $full_file ) {
						modules::Exception->throw("File $full_file doesn't exist for $sample_name");	
					}
					$reported_samples{$sample_name}++;
					if ($full_file =~ /snv/) {
						$sample_info{$patient}{$sample_name}{snv_file} =  $full_file;	
					} else {
						$sample_info{$patient}{$sample_name}{indel_file} =  $full_file;	
					}
				} elsif ($release_file_obj->file_name =~  /bam$/) {
					#Get the bam_stats file for coverage stats
					my $read_report_file = $local_rundir . '/summary/' . $sample_and_run . '.readReport.summary.txt';
					if ( !-e $read_report_file ) {
						modules::Exception->throw("File $read_report_file doesn't exist");	
					}
					open(FILE,$read_report_file) || modules::Exception->throw("ERROR: Can't open $read_report_file");
					while (<FILE>) {
						if (/^aligned: (\d+)/) {
							my $coverage = sprintf ("%.2f",100 * $1 / 2897310462);
							my $status = "COVER_OK";
							
							if ($coverage < $required_coverages{$sample_obj->sample_type}) {
								$status = "COVER_LOW(".$required_coverages{$sample_obj->sample_type}.")";
							}
							
							
							print COVER join("\t",
												$sample_name,
												$seq_centre_name,
												$coverage . 'X',
												$status
												) . "\n";
						}
					}
					close FILE;		 			
					
				}
			}
		}
		
		
	}
	
	
	
}

close COVER;


my %snv_data = my %indel_data = ();
my %snv_counts = my %indel_counts = ();
my %snv_headers = my %indel_headers = ();
my $first = 1;
my $merge_distance = 10; #Used to merge indels of the same type within 10bp of each other


for my $patient ( keys %sample_info ) {
	for my $sample ( keys %{$sample_info{$patient}} ) {
		my $sample_type = $sample_info{$patient}{$sample}{type};
		
		open(SNV_REPORT,$sample_info{$patient}{$sample}{snv_file}) || modules::Exception->throw("ERROR: Can't open snv file $patient $sample $sample_info{$patient}{$sample}{snv_file}");
		if (!$report_conf->exists('snv','sample_types',$sample_type)) {
    		modules::Exception->throw("ERROR: Cannot get annotation columns for snv $sample_type");
    	}
		
		my $fail_flag_snv = 0;
		my %snv_header_lookup = ();
		my @snv_report_headers = ();
		
		while (<SNV_REPORT>) {
			chomp;
			if (/^chr/) {
				@snv_report_headers = split("\t");
				for (my $i = 0; $i < scalar @snv_report_headers; $i++){
    				$snv_header_lookup{$i} = $snv_report_headers[$i];
				}
				#Set the master ordering from the first file we see
				if ($first) {
					for (my $i = 0; $i < scalar @snv_report_headers; $i++){
    					$snv_headers{$snv_report_headers[$i]} = $i;
					}
				}
				
				next;
			}
			
			if (/---/) {
				if (/FAIL/) {
					$fail_flag_snv = 1;
					next;
				}
				if (/PASS/) {
					next;
				}
				next;
			}
			

			if ($fail_flag_snv) {
				next unless $include_fails;
			}

			my @values = split("\t");
			my ($chr,$coord) = @values;
			$snv_counts{"$chr:$coord"}++;
			
			$snv_data{$chr}{$coord}{"0:samples"}{$sample}++;
			
			for (my $i = 1; $i < scalar @values; $i++){
				next if $snv_header_lookup{$i} eq 'chr' || $snv_header_lookup{$i} eq 'coord';
				my $header_name = $snv_header_lookup{$i};
				my $index = $snv_headers{$header_name}; #Get the index from the first file seen as that represents master order
				$snv_data{$chr}{$coord}{"$index:$snv_header_lookup{$i}"}{$values[$i]}++;
			}
			
			
		}
		close SNV_REPORT;
		
		
		
		my %indel_header_lookup = ();		
		my @indel_report_headers = ();
		my $fail_flag_indel;
		my %local_indels = ();
		
		open(INDEL_REPORT,$sample_info{$patient}{$sample}{indel_file}) || modules::Exception->throw("ERROR: Can't open indel file $sample_info{$patient}{$sample}{indel_file}");
		while (<INDEL_REPORT>) {
			chomp;
			if (/^chr/) {
				@indel_report_headers = split("\t");
				for (my $i = 0; $i < scalar @indel_report_headers; $i++){
    				$indel_header_lookup{$i} = $indel_report_headers[$i];
				}
				if ($first) {
					for (my $i = 0; $i < scalar @indel_report_headers; $i++){
    					$indel_headers{$indel_report_headers[$i]} = $i;
					}
					
					$first = 0;
				}
				next;
			}
			
			if (/---/) {
				if (/FAIL/) {
					$fail_flag_indel = 1;
					next;
				}
				if (/PASS/) {
					next;
				}
				next;
			}


			if ($fail_flag_indel) {
				next unless $include_fails;
			}

			my @values = split("\t");
			my ($chr,$start_coord,$end_coord,undef,undef,$type,undef,$var_base) = @values;
			
			my $min_merge = $start_coord - $merge_distance;
			my $max_merge = $start_coord + $merge_distance;
			my $increase_count = 1;
			
			#Check if there is an existing entry to merge with
			for my $count ( $min_merge .. $max_merge ) {
				if (exists $local_indels{"$chr:$count:$type"}) {
					#Here we should merge the events
					$start_coord = $count; #Set the lookup coordinate to be the same; real start coord will be recorded in indel_data
					$increase_count = 0; #don't increase the count in these cases
				}			    
			}
			
			$local_indels{"$chr:$start_coord:$type"}++;
			$indel_counts{"$chr:$start_coord:$type"}++ if $increase_count;
			
			$indel_data{$chr}{$start_coord}{$type}{"0:samples"}{$sample}++;
			for (my $i = 0; $i < scalar @values; $i++){
				#next if $indel_header_lookup{$i} eq 'chr' || $indel_header_lookup{$i} eq 'start_coord';
				my $header_name = $indel_header_lookup{$i};
				my $index = $indel_headers{$header_name}; #Get the index from the first file seen as that represents master order
				$indel_data{$chr}{$start_coord}{$type}{"$index:$indel_header_lookup{$i}"}{$values[$i]}++;
			}
		}
		close INDEL_REPORT;
		
		
	}
}




#Just get the uniq counts for iterating over hash
my %snv_uniq_counts = map {$_ => 1} values %snv_counts;
my %indel_uniq_counts = map {$_ => 1} values %indel_counts;

open(SNV,">$snv_file") || modules::Exception->throw("Can't open file to write $snv_file\n");
my $sample_count = keys %reported_samples;
print SNV "ANALYSED $sample_count SAMPLES on $date_stamp:\n\n".join("\n",sort keys %reported_samples)."\n\n";

print SNV join("\t",'Occurences','Samples','chr','coord');
my $snv_headers_printed = 0;

my %snv_lines = ();

for my $count (sort {$b<=>$a} keys %snv_uniq_counts) {
	for my $coord_str (sort {my ($a_chr,$a_coord) = split(':',$a); my ($b_chr,$b_coord) = split(':',$b); $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %snv_counts) {
		next unless $snv_counts{$coord_str} == $count;  #Only report snvs that have that count
		my ($chr,$coord) = split(':',$coord_str);
		my $sample_str;
		my @ordered_values = ();
		my $collapse_count = 0;
		for my $col_str (sort {my ($a_index) = split(':',$a); my ($b_index) = split(':',$b); $a_index <=> $b_index} keys %{$snv_data{$chr}{$coord}}) {
			
			my ($index,$col_name) = split(':',$col_str);
			if ($col_name eq 'samples') {
				#Combine the primary/meta tumours for count stats
				my %short_sample_names = ();
				for my $full_sample (keys %{$snv_data{$chr}{$coord}{$col_str}}) {
					(my $short_name = $full_sample) =~ s/_sg\d+_tumour\d+//;
					$short_sample_names{$short_name}++;
				}
				$collapse_count = keys %short_sample_names;
				$sample_str = join("; ",sort keys %{$snv_data{$chr}{$coord}{$col_str}});
				next;
			}
			
			if (!$snv_headers_printed) {
				print SNV "\t$col_name";
			}
			
			push @ordered_values, join("; ",keys %{$snv_data{$chr}{$coord}{$col_str}});
		
		}
		
		if (!$snv_headers_printed) {
			print SNV "\n\n";
			$snv_headers_printed = 1;
		}
		
		push @{$snv_lines{$collapse_count}}, join("\t",
													$collapse_count,
													$sample_str,
													$chr,
													$coord,
													@ordered_values
													);
		
	}
}

for my $snv_count (sort {$b<=>$a} keys %snv_lines) {
	print SNV join("\n",@{$snv_lines{$snv_count}});
	print SNV "\n\n";
}

close SNV;

open(INDEL,">$indel_file") || modules::Exception->throw("Can't open file to write $indel_file\n");
print INDEL "ANALYSED $sample_count SAMPLES on $date_stamp:\n\n".join("\n",sort keys %reported_samples)."\n\n";

print INDEL join("\t",'Occurences','Samples');
my $indel_headers_printed = 0;
my %indel_lines = ();

for my $count (sort {$b<=>$a} keys %indel_uniq_counts) {
	for my $coord_str (sort {my ($a_chr,$a_coord) = split(':',$a); my ($b_chr,$b_coord) = split(':',$b); $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %indel_counts) {
		next unless $indel_counts{$coord_str} == $count;
		my ($chr,$start_coord,$type) = split(':',$coord_str);
		my $sample_str;
		my @ordered_values = ();
		my $collapse_count;
		
		
		for my $col_str (sort {my ($a_index) = split(':',$a); my ($b_index) = split(':',$b); $a_index <=> $b_index} keys %{$indel_data{$chr}{$start_coord}{$type}}) {
			my ($index,$col_name) = split(':',$col_str);
			if ($col_name eq 'samples') {
				#Combine the primary/meta tumours for count stats
				my %short_sample_names = ();
				for my $full_sample (keys %{$indel_data{$chr}{$start_coord}{$type}{$col_str}}) {
					(my $short_name = $full_sample) =~ s/_sg\d+_tumour\d+//;
					$short_sample_names{$short_name}++;
				}
				$collapse_count = keys %short_sample_names;
				$sample_str = join("; ",sort keys %{$indel_data{$chr}{$start_coord}{$type}{$col_str}});
				next;
			}
			
			if (!$indel_headers_printed) {
				print INDEL "\t$col_name";
			}
			
			push @ordered_values, join("; ",keys %{$indel_data{$chr}{$start_coord}{$type}{$col_str}});
		
		}
		if (!$indel_headers_printed) {
			print INDEL "\n\n";
			$indel_headers_printed = 1;
		}
		
		push @{$indel_lines{$collapse_count}}, join("\t",
													$collapse_count,
													$sample_str,
													@ordered_values
													);
		
			
	}
}


for my $indel_count (sort {$b<=>$a} keys %indel_lines) {
	print INDEL join("\n",@{$indel_lines{$indel_count}});
	print INDEL "\n\n";
}

close INDEL;
#print Dumper \%indel_data;
#print Dumper \%snv_data;
