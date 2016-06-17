#! /usr/bin/perl -w

use strict;
use modules::Polyphen;
use modules::Adaptors::BulkInsert;
use modules::Adaptors::Run;
use modules::Adaptors::SNV;
use modules::Adaptors::Sample;
use modules::Adaptors::Filter;
use modules::Report;
use modules::Overlap;
use modules::Utils;
use modules::Pipeline;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "overlap_outfile=s",
	   "tmpdir=s",
	   "runid=i",
	   "writeDB=i",
	   "annotation_file=s",
	   "ccds_coord_file=s",
	   "ignore_filter=s",
	   "polyphen_filter=s",
	   "report_xml=s"
    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{report_xml} || !$OPT{overlap_outfile} || !$OPT{tmpdir} || !$OPT{annotation_file} || !$OPT{ccds_coord_file});


=pod

=head1 SYNOPSIS

filter_polyphen.pl -report_xml report_xml_file_for_passed definitions -polyphen_filter polyphen_filter_name(default=failter_polyphen) -overlap_outfile <overlap_outfile> -writeDB 1|0 -tmpdir scratch_dir -ignore_filter don't_report_these_filters -annotation_file annotation_file -ccds_coord_file ccds_coord_file [options]

Required flags: -overlap_outfile -tmpdir -annotation_file -ccds_coord_file -runid

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_polyphen.pl -> Script to drive polyphen and obtain the polyphen score

=head1 DESCRIPTION

Mar 30, 2011

Single processor polyphen wrapper

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_polyphen.pl -overlap_outfile test_polyphen/filter_polyphen.match -tmpdir test_polyphen/ -ccds_coord_file /home/matt/work/pipeline/conf/mouse/NCBIM37/ccds/040411/NCBIM37.ccdsGene.txt -annotation_file /home/matt/work/pipeline/conf/mouse/NCBIM37/gene/160611/NCBIM37.BiomartMapping.txt -runid 15

=cut

# Put command line options into the right places

my $outfile = $OPT{overlap_outfile};
(my $polyphen_info_file = $outfile) =~ s/match/nomatch/;

# Put command line options into the right places
my $pipe_conf = modules::Pipeline->get_pipe_conf();
my $clus_conf = modules::Pipeline->get_cluster_conf();
my $polyphen_executable = $pipe_conf->read('binaries','polyphen','binary');
my $polyphen_conf = $pipe_conf->read('binaries','polyphen','polyphen_conf');
my $conf_base = $clus_conf->read('svn','conf_dir') . '/polyphen/' . $pipe_conf->read('binaries','polyphen','version');

#Scratch directory for the program; create if needed, otherwise reuse alignments already there
if (!-d $OPT{tmpdir}) {
	mkdir($OPT{tmpdir});
}
my $working_dir = `readlink -f $OPT{tmpdir}`;
chomp $working_dir;

if (!-d $working_dir) {
	modules::Exception->throw("Directory $working_dir doesn't exist");
} else {
    #Remove locks from previous run if they exist
    my $lockdir = $working_dir . '/lock';
    if (-d $lockdir) {
        my $lock_files = $lockdir .'/*lock';
        my $lock_flag = `ls $lock_files 2>/dev/null`;
        if ($lock_flag) {
            print STDERR "Removing locks...\n";
            system("rm $lockdir/*");
        }
    }
}

if (!-d $polyphen_conf) {
	modules::Exception->throw("Polyphen conf $polyphen_conf doesn't exist");
}

if (!-x $polyphen_executable) {
	modules::Exception->throw("Polyphen executable $polyphen_executable problems");
}

my $report_xml = $OPT{report_xml};

if ( !-e $report_xml ) {
	modules::Exception->throw("File $report_xml doesn't exist");	
}

my $runid = $OPT{runid};
my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;
my $annotation_file = $OPT{annotation_file};
my $ref = defined $OPT{ref}?$OPT{ref}:'GRCh37';

my $ccds_coord_file = $OPT{ccds_coord_file};

my $splice_length = $pipe_conf->read('cutoffs','splice_length_cutoff');

my $polyphen_filter_name = defined $OPT{polyphen_filter_name}?$OPT{polyphen_filter_name}:'filter_polyphen';

if ( !-e $annotation_file ) {
	modules::Exception->throw("File $annotation_file doesn't exist");	
}


if ( !-e $ccds_coord_file ) {
	modules::Exception->throw("File $ccds_coord_file doesn't exist");	
}

my ($run) = modules::Adaptors::Run->search('id' => $runid);
my ($polyphen_filter) = modules::Adaptors::Filter->search('name' => $polyphen_filter_name);

unless (defined $run && defined $polyphen_filter) {
    modules::Exception->throw("Did not get all necessary things from database [run_id:$runid] [polyphen_filter:$polyphen_filter_name]");
}

#Get the sample group
my $sample_id = $run->sample_id;
my ($sample_obj) = modules::Adaptors::Sample->search('id'=>$sample_id);
unless (defined $sample_obj) {
	modules::Exception->throw("Could not get sample from db using runid $runid");
}
unless (defined $sample_obj->sample_type) {
	modules::Exception->throw("Could not get sample_type from sample");
}
my $sample_type = $sample_obj->sample_type;

#Get a list of passed snvs to run polyphen on

my %ignore_filters = ();
if (defined $OPT{ignore_filter}) {
	%ignore_filters = map {$_ => 1} split(",",$OPT{ignore_filter});
}


#generate the gene mapper 
my $annotation = modules::Annotation->new(-annotation_file=>$annotation_file,-ccds_coord_file=>$ccds_coord_file,-splice_size=>$splice_length);
#Always report for exons and splice
my $snv_report = modules::Report->new(
											 -run        => $run,
                                             -gene_mapper => $annotation,
                     						 -sample_type => $sample_type
                                             );
                                             
#Check the filters exist where ignoring
foreach my $ignore_filter (keys %ignore_filters) {
	my ($filter) = modules::Adaptors::Filter->search('name' => $ignore_filter);

	unless ( defined $filter ) {
    	modules::Exception->throw("Unable to retrieve filter object from database for ignore_filter $ignore_filter");
	}
}


#Load the filters, this is specific to each sample group
if (keys %ignore_filters) {
	$snv_report->load(-report_xml=>$report_xml,-ignore_filters=>\%ignore_filters);
} else {
	$snv_report->load(-report_xml=>$report_xml);
}

$snv_report->generate_pass_fail();

#Now generate the polyphen input files
my $polyphen_input_mapped = $working_dir . '/polyphen.mapped.input.txt';
my $polyphen_input_unmapped = $working_dir. '/polyphen.unmapped.input.txt';
my $fasta = $working_dir . '/polyphen.fa';

open(my $ph_mapped, ">$polyphen_input_mapped")
    or modules::Exception->throw("Unable to open file for writing polyphen input format [$polyphen_input_mapped]");

open(my $ph_unmapped, ">$polyphen_input_unmapped")
    or modules::Exception->throw("Unable to open file for writing polyphen input format [$polyphen_input_unmapped]");

open(my $fh_fasta, ">$fasta")
    or modules::Exception->throw("Unable to open file for writing polyphen input format [$fasta]");


# Create the polyphen object for both mapped and unmapped cases
my $polyphen
    = modules::Polyphen->new(
    							-input_mapped_file    => $polyphen_input_mapped,
    							-input_unmapped_file    => $polyphen_input_unmapped,
								-executable_path => $polyphen_executable,
								-config_path => $polyphen_conf,
								-working_dir     => $working_dir,
								-fasta => $fasta					
								);
								
#Get the passed cases from the snvreport object
my ($polyphen_data,$polyphen_count) = $snv_report->generate_polyphen_input();
print STDERR "Analysing $polyphen_count entries for polyphen...\n";

#Create the overlap structure for generating the mappings
my $overlap = modules::Overlap->new();
my $mapped_found = my $unmapped_found = 0;

for my $chr ( sort keys %{$polyphen_data} ) {
	#Open an tmp file for overlapping
	my $chr_ref_file = $polyphen_input_mapped.".$chr";
	open(CHR,">$chr_ref_file") || modules::Exception->throw("Can't open file $chr_ref_file\n");
	
	for my $coord ( sort keys %{$polyphen_data->{$chr}} ) {
		#1	3062302	3062302	C G QU43432 V->L
		#1	3083416	3083416	A T QZ85904 D->S
		
		my $refbase = $polyphen_data->{$chr}{$coord}{ref}; 
		my $varbase = $polyphen_data->{$chr}{$coord}{var};
		my $uniprot = $polyphen_data->{$chr}{$coord}{uniprot};
		my $aa_change = $polyphen_data->{$chr}{$coord}{aachange};
			
		my $snv_name = "$chr:$coord:$refbase->$varbase:$aa_change:$uniprot";
		print CHR "$chr\t$coord\t$coord\t$snv_name\n";
   	}
   	close CHR;

	my $chr_mapped_file = $conf_base. '/' . $ref . '.uniprot_mapping.' . $chr;
	my $chr_unmapped_file = $conf_base. '/' . $ref . '.uniprot_nomapping.' . $chr;
	
	#First check if it's mapped
	my ($overlaps_mapped,$mapped_flag) =  $overlap->overlap(-ref=>$chr_ref_file,-coord=>$chr_mapped_file,-silent=>1);
   	if ($mapped_flag) {
   		$mapped_found  = 1;
		#Get it from the uniprot mapping file
		$polyphen->generate_input(-overlap=>$overlaps_mapped,-mapped=>1,-fh=>$ph_mapped);
   	} 

   	my ($overlaps_unmapped,$unmapped_flag) =  $overlap->overlap(-ref=>$chr_ref_file,-coord=>$chr_unmapped_file, -silent=>1);
   	if ($unmapped_flag) {
   		#Get it from the ccds file
   		$unmapped_found = 1;
   		$polyphen->generate_input(-overlap=>$overlaps_unmapped,-mapped=>0,-fh=>$ph_unmapped,-fh_fasta=>$fh_fasta);
   	}
	system("rm $chr_ref_file");
}

close($fasta);
close($ph_mapped);
close($ph_unmapped);

#Now generate the polyphen_info file; this occurs for all summary entries allowing users to run polyphen if required
open(INFO,">$polyphen_info_file") || modules::Exception->throw("Can't open file $polyphen_info_file\n");
my $row_data = $snv_report->get_all_rows();

for my $chr ( sort keys %{$row_data} ) {
	#print STDERR "CHR $chr\n";
	my $chr_info_file = $polyphen_input_mapped.".info.".$chr;
	open(CHR,">$chr_info_file") || modules::Exception->throw("Can't open file $chr_info_file\n");
    for my $coord (sort {$a<=>$b} keys %{$row_data->{$chr}}) {
    	print CHR "$chr $coord $coord\n";	
    }
     
    close CHR;
    
    my $chr_mapped_file = $conf_base. '/' . $ref . '.uniprot_mapping.' . $chr;
	my $chr_unmapped_file = $conf_base. '/' . $ref . '.uniprot_nomapping.' . $chr;
	
	my ($overlaps_mapped,$mapped_flag) =  $overlap->overlap(-coord=>$chr_info_file,-ref=>$chr_mapped_file,-silent=>1);
   	if ($mapped_flag) {
   		for my $chr (keys %{$overlaps_mapped->{PASS}}) {
			for my $overlap_coord (keys %{$overlaps_mapped->{PASS}{$chr}}) {
				for my $snv_str (keys %{$overlaps_mapped->{PASS}{$chr}{$overlap_coord}{$overlap_coord}}) {
					print INFO join("\t",
									$chr,
									$overlap_coord,
									$overlap_coord,
									$snv_str
									) . "\n";
				}
			}
   		}
   	} 

   	my ($overlaps_unmapped,$unmapped_flag) =  $overlap->overlap(-coord=>$chr_info_file,-ref=>$chr_unmapped_file, -silent=>1);
   	if ($unmapped_flag) {
   		for my $chr (keys %{$overlaps_unmapped->{PASS}}) {
			for my $overlap_coord (keys %{$overlaps_unmapped->{PASS}{$chr}}) {
				for my $snv_str (keys %{$overlaps_unmapped->{PASS}{$chr}{$overlap_coord}{$overlap_coord}}) {
					print INFO join("\t",
									$chr,
									$overlap_coord,
									$overlap_coord,
									$snv_str
									) . "\n";
				}
			}
   		}
   	}
    
    system("rm $chr_info_file");
    
}

#The two output files
my $polyphen_output_mapped = $working_dir . '/polyphen.mapped.output.txt';
my $polyphen_output_unmapped = $working_dir . '/polyphen.unmapped.output.txt';

if ($mapped_found) {
	$polyphen->run(-output=>$polyphen_output_mapped,-mapped=>1);
	$polyphen->parse_result(-output=>$polyphen_output_mapped);
}

if ($unmapped_found) {
	$polyphen->run(-output=>$polyphen_output_unmapped,-mapped=>0);
	$polyphen->parse_result(-output=>$polyphen_output_unmapped);
}

open(FILE,">$outfile") || modules::Exception->throw("Can't open file to write $outfile\n");

my $combined_results = $polyphen->get_results();

my @polyphen_inserts = ();

for my $chr (sort keys %{$combined_results}) {
	
	my %chr_snvs;

    my $snv_iterator
	= modules::Adaptors::SNV->search(
										  'run_id' => $runid,
					    				  'chr' => $chr
					    				  );

	#Get a list of snvs for the chromosome
    while (my $snv = $snv_iterator->next){
		$chr_snvs{$snv->coord} = $snv->id;
    }
    
    
	
	for my $coord (sort {$a<=>$b} keys %{$combined_results->{$chr}}) {
		for my $nt_change (keys %{$combined_results->{$chr}{$coord}}) {
			my @fields = split(':',$combined_results->{$chr}{$coord}{$nt_change}{name});
			
			#Create the snv filter
			my %snv_filter_exon = ('filter_id' => $polyphen_filter->id,
							       'snv_id'    => $chr_snvs{$coord},
							       'filtermatch' => 1,
							       'filterpass' => 1,
							       'attribute' => 'poly_pred=' . $combined_results->{$chr}{$coord}{$nt_change}{prediction} . ';poly_score=' . $combined_results->{$chr}{$coord}{$nt_change}{score});
			
			push @polyphen_inserts, \%snv_filter_exon;     
			
			print FILE join("\t",
								$chr,
								$coord,
								$coord,
								$combined_results->{$chr}{$coord}{$nt_change}{score},
								$combined_results->{$chr}{$coord}{$nt_change}{prediction},
								$combined_results->{$chr}{$coord}{$nt_change}{name}
								). "\n";
			
		}
	}
}

close FILE;
#Do the insert if flag set
modules::Adaptors::BulkInsert->insert(-table_name=>'snvs_filters', -data=>\@polyphen_inserts) if $write_to_db && @polyphen_inserts;








