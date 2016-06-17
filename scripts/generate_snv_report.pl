#! /usr/bin/perl -w

use strict;
use modules::Report;
use modules::Adaptors::Filter;
use modules::Adaptors::BulkInsert;
use modules::Annotation;
use modules::Pipeline;
use modules::ConfigXML;
use Getopt::Long;
use modules::Exception;
use modules::Utils;
use File::Basename;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

GetOptions(\%OPT, 
		    "help|h",
		    "man|m",
		 	"runid=i",
		    "chr=s",
		 	"start=i",
		 	"end=i",
		 	"ref=s",
		   	"report_xml=s",
		 	"filter_summary_file=s",
		 	"annotation_file=s",
		 	"exon_coord_file=s",
		    "tsv_file=s",
		    "pass_summary_file=s",
		    "pass_file=s",
		    "genefilter_files=s",
		    "ignore_filter=s",
		    "polyphen_file=s",
		    "pass_filter=s",
		    "rareallele_filter=s",
		    "writeDB=i",
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{report_xml} || !$OPT{tsv_file} || !$OPT{runid} || !$OPT{filter_summary_file} || !$OPT{annotation_file} || !$OPT{exon_coord_file} || !$OPT{pass_summary_file} || !$OPT{pass_file});

	   
=pod

=head1 SYNOPSIS

generate_snv_eport.pl -runid runid -rareallele_filter filter_rare_allele(default=filter_rareallele_snv) -report_xml report_xml_file -polyphen_file polyphen_score_file -pass_filter filter_report_name(default=filter_pass_snv) -ignore_filter don't_require_these_filter_to_pass -ref ref_name(default=NCBIM37) -annotation_file annotation_file -exon_coord_file exon_coord_file -debug print_stderr -filter_summary_file filter_summary_output_file -tsv_file output_tsv_file -chr specific_chr(default=all) -start chr_start_coord(default=all) -end chr_end_coord(default=all) -pass_summary_file pass_rows_summary [options]

Required flags: -runid -report_xml -tsv_file -filter_summary_file -annotation_file -exon_coord_file -pass_summary_file -pass_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

generate_snv_report.pl -> Script to generate a snv report

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Dan Andrews

=head1 EXAMPLE

generate_snv_report.pl

=cut

my $pipe_conf = modules::Pipeline->get_pipe_conf();

#Lots of variables to load...
my $report_file = $OPT{report_file};
my $output_file = $OPT{output_file};
my $summary_file = $OPT{pass_summary_file};
my $pass_file = $OPT{pass_file};
my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:1;
my $pass_filter = defined $OPT{pass_filter}?$OPT{pass_filter}:'filter_pass_snv';
my $rareallele_filter = defined $OPT{rareallele_filter}?$OPT{rareallele_filter}:'filter_rareallele_snv';


my ($filter) = modules::Adaptors::Filter->search('name' => $pass_filter);
my ($allele_filter) = modules::Adaptors::Filter->search('name' => $rareallele_filter);

unless ( defined $filter && defined $allele_filter ) {
    modules::Exception->throw("Unable to retrieve filter object from database for $pass_filter or $rareallele_filter ");
}

my @filter_files = ();
if (defined $OPT{genefilter_files}) {
	@filter_files = split(",",$OPT{genefilter_files});
}

#Flag tells us whether to report polyphen values
my $polyphen_score = my $polyphen_info = 0;
my $polyphen_score_file;
my $polyphen_info_file;
if (defined $OPT{polyphen_file}) {
	$polyphen_score_file = $OPT{polyphen_file};
	($polyphen_info_file = $polyphen_score_file) =~ s/match/nomatch/;
	#One of the two must exist...
	if (!-e $polyphen_score_file && !-e $polyphen_info_file) {
		modules::Exception->throw("ERROR: Polyphen files $polyphen_score_file and $polyphen_info_file don't exist");
	} 
	
	if (-e $polyphen_score_file) {
		$polyphen_score = 1;
	}
	
	if (-e $polyphen_info_file) {
		$polyphen_info = 1;
	}
	
} 

my $runid = $OPT{runid};
my $tsv_file = $OPT{tsv_file};
my $chr = defined $OPT{chr}?$OPT{chr}:'';
my $start = defined $OPT{start}?$OPT{start}:'';
my $end = defined $OPT{end}?$OPT{end}:'';
my $annotation_file = $OPT{annotation_file};
my $ref = defined $OPT{ref}?$OPT{ref}:'GRCh37';
my $splice_size = $pipe_conf->read('cutoffs','splice_length_cutoff');
my $exon_coord_file = $OPT{exon_coord_file};
my $report_xml = $OPT{report_xml};

if ( !-e $report_xml ) {
	modules::Exception->throw("File $report_xml doesn't exist");	
}

if ( !-e $annotation_file ) {
	modules::Exception->throw("File $annotation_file doesn't exist");	
}


if ( !-e $exon_coord_file ) {
	modules::Exception->throw("File $exon_coord_file doesn't exist");	
}

my ($run) = modules::Adaptors::Run->search('id' => $runid);
unless (defined $run) {
    modules::Exception->throw("Did not get all necessary things from database [run_id:$runid]");
}

my $sample_id = $run->sample_id;
my ($sample_obj) = modules::Adaptors::Sample->search('id'=>$sample_id);
unless (defined $sample_obj) {
	modules::Exception->throw("Could not get sample from db using runid $runid");
}
unless (defined $sample_obj->sample_type) {
	modules::Exception->throw("Could not get sample_type from sample");
}
my $sample_type = $sample_obj->sample_type;

#Get the gene column name
my $config = modules::ConfigXML->new($report_xml);
if (!$config->exists('snv','sample_types',$sample_type)) {
   	modules::Exception->throw("ERROR: Cannot get annotation columns for snv $sample_type");
}
 
my $gene_col_name;
if ($config->exists('snv','sample_types',$sample_type,'annotations')) {
	($gene_col_name) = split(",",$config->read('snv','sample_types',$sample_type,'annotations'));
} else {
	($gene_col_name) = split(",",$config->read('snv','common','annotations'));
} 
    

my %ignore_filters = ();
if (defined $OPT{ignore_filter}) {
	%ignore_filters = map {$_ => 1} split(",",$OPT{ignore_filter});
}

#generate the gene mapper 
my $annotation = modules::Annotation->new(-annotation_file=>$annotation_file,-exon_coord_file=>$exon_coord_file,-splice_size=>$splice_size,-report_xml=>$report_xml,-sample_type=>$sample_type,-mutant_type=>'snv');

#Create the snv report params

my %snv_report_param = ( -run => $run,
	                     -gene_mapper => $annotation,
	                     -sample_type => $sample_type,
	                     -gene_col_name => $gene_col_name
	                     );


if (@filter_files) {
	$snv_report_param{-filter_files} = \@filter_files
}

if ($polyphen_score) {
	$snv_report_param{-polyphen_file} = $polyphen_score_file;
}

if ($polyphen_info) {
	$snv_report_param{-polyphen_info_file} = $polyphen_info_file;
}


#Check the filters exist where ignoring
foreach my $ignore_filter (keys %ignore_filters) {
	my ($filter) = modules::Adaptors::Filter->search('name' => $ignore_filter);

	unless ( defined $filter ) {
    	modules::Exception->throw("Unable to retrieve filter object from database for ignore_filter $ignore_filter");
	}
}

my $snv_report = modules::Report->new(%snv_report_param);


#Load the filters, this is specific to each sample group
if (keys %ignore_filters) {
	$snv_report->load(-report_xml=>$report_xml,-ignore_filters=>\%ignore_filters);
} else {
	$snv_report->load(-report_xml=>$report_xml);
}

if ($chr && $start) {
	$snv_report->generate_pass_fail(-chr=>$chr, -start=>$start, -end=>$end);
} elsif ($chr) {
	$snv_report->generate_pass_fail(-chr=>$chr);
} else {
	$snv_report->generate_pass_fail();
}

$snv_report->print_to_files(-tsv_file=>$tsv_file,-pass_file=>$pass_file);
$snv_report->summarize(-summary_file=>$summary_file);


my $overlap = $OPT{filter_summary_file};
my $overlapdir = dirname($overlap);
$overlapdir =~ s/summary//;
opendir(DIR,$overlapdir) || die "Can't open directory $overlapdir\n";
my @local_filter_files = grep {/filter_/} readdir DIR;
closedir (DIR);

my %all_snvs = ();
my %bp_lookup = ();

for my $filter_file ( @local_filter_files ) {
	next if $filter_file =~ /indel/;
	#This filter is covered by filter_basic
    open(FILTER,"$overlapdir/$filter_file") || modules::Exception->throw("Can't open file $overlapdir/$filter_file\n");
    while (<FILTER>) {
        my ($chr,$start,$end,$rest) = split;
        if ($rest =~ /([AGTCN]->[ATGC])/) {
        	$bp_lookup{$chr}{$start} = $1;
        }
        my ($filter) = $filter_file =~ /filter_(.*)$/;
        $filter =~ s/\.snv//;
        push @{$all_snvs{$chr}{$start}},$filter;
    }
    close FILTER;
}

my @overlap_lines = ();
open(OVERLAP,">$overlap") || die "Can't open file to write $overlap\n";
for my $chr (sort keys %all_snvs) {
	for my $coord ( sort {$a<=>$b} keys %{$all_snvs{$chr}} ) {
		my $filter_str = join(",",sort @{$all_snvs{$chr}{$coord}});
		my $bp_change = $bp_lookup{$chr}{$coord};		
	    print OVERLAP "$chr $coord $coord $bp_change $filter_str\n";
	}
}

close OVERLAP;

#Get the snv filter
my @snv_filter_inserts;

#First the passed rows
my $pass_row_data = $snv_report->get_pass_rows();
my $rare_allele_rows = $snv_report->get_allele_rows();
my $fail_row_data = $snv_report->get_fail_rows();


my %rows = (
			'pass'=>$pass_row_data,
			'allele'=>$rare_allele_rows,
			'fail'=>$fail_row_data
			);


for my $row_type (keys %rows) {
	my $row_count = 0;
	for my $chr ( sort keys %{$rows{$row_type}} ) {
		    for my $coord (sort {$a<=>$b} keys %{$rows{$row_type}->{$chr}}) {
		    	
		    	my %search_params = ('run_id' => $runid);
				$search_params{'chr'} = $chr;
				$search_params{'coord'} = $coord;		
				my @snvs = modules::Adaptors::SNV->search(\%search_params);
				if (@snvs > 1) {
					modules::Exception->throw("ERROR: Multiple snvs for chr $chr coord $coord");
				}
				my %snv_filter = (
								 'snv_id'    => $snvs[0]->id,
				  				 'filter_id' => $filter->id,
				  				 'filtermatch' => 1,
				  				 'filterpass' => 1
				  			 );
	
				push @snv_filter_inserts, \%snv_filter;
				$row_count++;
		    }
	}
	print STDERR "Reporting $row_count $row_type rows\n";
}

modules::Adaptors::BulkInsert->insert(-table_name=>'snvs_filters', -data=>\@snv_filter_inserts) if $write_to_db && @snv_filter_inserts;

#Finally insert the snv_
$snv_report->insert_SNV_rows();


