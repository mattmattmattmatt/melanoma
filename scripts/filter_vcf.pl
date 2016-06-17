#! /usr/bin/perl -w

use strict;
use modules::Adaptors::Variant;
use modules::Adaptors::Filter;
use modules::Adaptors::Variant_Filter;
use modules::Overlap;
use modules::Exception;
use modules::Adaptors::BulkInsert;
use modules::ConfigXML;
use modules::SystemCall;
use modules::Vcf;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "overlap_outfile_snv=s",
		   "overlap_outfile_indel=s",
		   "vcf_file=s",
		   "snv_infile=s",
		   "indel_infile=s",
		   "sample_name=s"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf_file} || !$OPT{overlap_outfile_snv} || !$OPT{overlap_outfile_indel} || !$OPT{snv_infile} || !$OPT{indel_infile} || !$OPT{sample_name});

	   
=pod

=head1 SYNOPSIS

filter_vcf.pl -overlap_outfile_snv output_snv_overlap_file -overlap_outfile_indel output_indel_overlap_file -vcf_file vcf_file -snv_infile  [options]

Required flags: -overlap_outfile_snv -overlap_outfile_indel -vcf_file -snv_infile merge_vcf_snv_file -indel_infile merge_vcf_indel_file -sample_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_vcf.pl -> Script to parse and filter a vcf file

=head1 DESCRIPTION

Feb 7, 2012

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_vcf.pl

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $sys_call = modules::SystemCall->new();

#Get the flags for the particular library
my $outfile_indel = $OPT{overlap_outfile_indel};
my $outfile_snv = $OPT{overlap_outfile_snv};

my $vcf_file = $OPT{vcf_file};

my $sample_name = $OPT{sample_name};
my $tumour_flag = modules::Pipeline->get_tumour_flag(-sample_name=>$sample_name);

my $sample_type;
my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);

if (!$sample_obj) {
	modules::Exception->throw("ERROR: Cannot retrieve sample object for $sample_name");
}
$sample_type = $sample_obj->sample_type;


my $vcf_obj = modules::Vcf->new(-sample_type=>$sample_type);

#Parse the filtered vcf
$vcf_obj->parse_vcf(-vcf_file => $vcf_file);
my $vcf_data = $vcf_obj->get_vcf(-vcf_file => $vcf_file);
my $filter_vcf_data = $vcf_obj->filter_vcf(-vcf_file=>$vcf_file, -snv_depth_file=>$OPT{snv_infile},-indel_depth_file=>$OPT{indel_infile},-tumour_flag=>$tumour_flag);

open(VAR,">$outfile_indel") || modules::Exception->throw("Can't open file $outfile_indel\n");

#Parse the filtered deletions; these will be added to the database
if (exists $filter_vcf_data->{DEL}) {
	for my $chr (sort keys %{$filter_vcf_data->{DEL}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$filter_vcf_data->{DEL}{$chr}}) {
			for my $end_coord (keys %{$filter_vcf_data->{DEL}{$chr}{$start_coord}}) {
				for my $bases (keys %{$filter_vcf_data->{DEL}{$chr}{$start_coord}{$end_coord}}) {
					my $rest = $filter_vcf_data->{DEL}{$chr}{$start_coord}{$end_coord}{$bases};
									
		  			print VAR join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							"DEL^^^$bases^^^$rest"
		  							) . "\n";
				}
			}
		}
	}
}

#Parse the insertions
if (exists $filter_vcf_data->{INS}) {
	for my $chr (sort sort keys %{$filter_vcf_data->{INS}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$filter_vcf_data->{INS}{$chr}}) {
			for my $end_coord (keys %{$filter_vcf_data->{INS}{$chr}{$start_coord}}) {
				for my $bases (keys %{$filter_vcf_data->{INS}{$chr}{$start_coord}{$end_coord}}) {
					my $rest = $filter_vcf_data->{INS}{$chr}{$start_coord}{$end_coord}{$bases};
					
		  			 print VAR join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							"INS^^^$bases^^^$rest"
		  							) . "\n";
		  			 
				}
			}
		}
	}
}




#Parse the vcf and generate the pileup coords file
open(SNV,">$outfile_snv") || modules::Exception->throw("Can't open file $outfile_snv\n");

#Parse the filtered deletions; these will be added to the database
if (exists $filter_vcf_data->{SNV}) {
	for my $chr (sort keys %{$filter_vcf_data->{SNV}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$filter_vcf_data->{SNV}{$chr}}) {
			for my $end_coord (keys %{$filter_vcf_data->{SNV}{$chr}{$start_coord}}) {
				for my $bases (keys %{$filter_vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}}) {
					my $rest = $filter_vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}{$bases};
					
		  			print SNV join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							"SNV^^^$bases^^^$rest"
		  							) . "\n";
				}
			}
		}
	}
}





