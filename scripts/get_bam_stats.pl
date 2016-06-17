#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Set::IntSpan;
use modules::Utils;
use modules::Exception;
use modules::SystemCall;
use modules::Pipeline;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "bam=s",
	   "exon=s",
	   "output=s",
	   "debug"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{bam});



=pod

=head1 SYNOPSIS

get_bam_stats.pl -bam bamfile -exon exonfile_if_want_exon_stats -output output_file -samtoolsbin samtools_binary [options]

Required flags: -bam -exon

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

get_bam_stats.pl -> Generate alignment statistics for a run

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./get_bam_stats.pl -bam /home/matt/work/analysis/pipeline/runs/IGL00054_11/overlap/IGL00054_11_sort_bam_match.bam -exon /home/matt/work/pipeline/conf/mouse/NCBIM37/exon/040411/NCBIM37.exon -samtoolsbin ../bin/samtools

=cut

my $output;
my $tmp_out;
if (defined $OPT{output}) {
	$output = $OPT{output};
	my $basedir = dirname($output);
	my $randnumber = int(rand(10000000));
	$tmp_out = $basedir.'/tmp'.$randnumber;
} else {
	$output = "summary.txt";
	$tmp_out = "out.tmp";
}

my $sys_call = modules::SystemCall->new();

open(FILE,">$output") || modules::Exception->throw("Can't open file to write $output\n");

my $bam = $OPT{bam};
if (!-e $bam) {
	modules::Exception->throw("ERROR: Cannot open bam file $OPT{bam}");
}

#Get the xml files and create the pipeline object
my $pipe_config = modules::Pipeline->get_pipe_conf();
my $samtools_bin = $pipe_config->read('binaries','samtools','binary');

if (!-x $samtools_bin) {
	modules::Exception->throw("ERROR: Cannot execute samtools");
}

#First get the overall read stats
my %run_stats = &Read_Stats();

for my $field ( keys %run_stats ) {
    print FILE $field,': ',$run_stats{$field},"\n";
    print STDERR $field,': ',$run_stats{$field},"\n" if $OPT{debug};
}

#Generate the overlap alignment stats
sub Read_Stats {
    my %stats = ();
    my $command = "$samtools_bin flagstat $bam";
    my @output = split("\n",`$command`);
	for my $line ( @output ) {
	    if ($line =~ /^(\d+).*in\s+total/) {
	    	$stats{total} = $1;
	    } elsif ($line =~ /^(\d+).*mapped\s+\(\d+\.\d+%/) {
	    	$stats{aligned} = $1;
	    } elsif ($line =~ /^(\d+).*singletons\s+\(\d+\.\d+%/) {
	    	$stats{singles} = $1;
	    } elsif ($line =~ /^(\d+).*with itself and mate mapped/) {
	    	$stats{paired} = $1;
	    } elsif ($line =~ /^(\d+).*with mate mapped to a different chr$/) {
	    	$stats{mispaired} = $1;
	    } 
	}
	$stats{unaligned} = $stats{total} - $stats{aligned};
	
	if ($OPT{exon}) {
		my $exon = $OPT{exon};
		if (!-e $exon) {
			modules::Exception->throw("ERROR: Cannot open exon file $OPT{exon}");
		}
		open(EXON,"$exon") || die "Can't open file $exon\n";
	
		
		my $count = 0;
		my $samtools_index = 0;
		my %exon_group = ();
		
		#Bin size for exon lookup index -> don't make much higher as there is a length limit on unix commands 
		my $samtools_binsize = 500;
		
		while (<EXON>) {
			if ($count % $samtools_binsize == 0) {
				$samtools_index++;
			}
					
			my @fields = split;
					
			#Group together these coordinates for samtools arguments
			my $chrcoord = "$fields[0]:$fields[1]-$fields[2]";
			push @{$exon_group{$samtools_index}}, $chrcoord;
			$count++;
		}
			
		#file shouldn't be here but will be removed in case to prevent appending to existing file
		if ( -e $tmp_out ) {
			$sys_call->run("rm $tmp_out");
		}
	
		
		#Run the command with a sub set of coordinates at a time passed in; this makes it run in minutes vs hours with one coordinate at a time
		for my $key ( sort {$a<=>$b} keys %exon_group ) {
		    my $chr_str = join(" ",@{$exon_group{$key}});
			my $command = "$samtools_bin view $bam $chr_str >> $tmp_out";
			system($command); 
		}
		
		#Get the read name and the coord to make a unique key;  Accounts for unaligned pairs where unaligned reads lists aligned reads coords
		my $sort_command = "awk '{print ".'$1,$4}\''." $tmp_out | sort -S 1500M | uniq |wc -l";
		my $read_to_exon_count = `$sort_command`;
		chomp $read_to_exon_count;
		system("rm $tmp_out");
		my $percent;
		if ($stats{paired}) {
			$percent = sprintf("%.2f",$read_to_exon_count/$stats{paired} *100);
		} else {
			$percent = sprintf("%.2f",$read_to_exon_count/$stats{aligned} *100);
		}
		$stats{exon_aligned} = "$read_to_exon_count ($percent% of total aligned paired reads)"; 
	}	
	return %stats;
}
