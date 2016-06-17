#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Set::IntSpan;
use modules::Exception;
use modules::SystemCall;
use modules::Pipeline;
use File::Basename;
use modules::ConfigXML;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "files=s",
	   "dir=s",
	   "threadNum=i",
	   "suffix=s",
	   "compress",
	   "tarball=s",
	   "keep",
	   "workdir",
	   "run=i",
	   "runid=i",
	   "bam",
	   "mv",
	   "debug"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || (!$OPT{files} && !$OPT{dir}));



=pod

=head1 SYNOPSIS

compress.pl -mv move_file_when_done -runid runid_to_add_to_tarball_name -files comma_delim_list_of_files_without_compression_suffix -dir directory_to_compress -keep keep_input_files -suffix compression_suffix(bz2,zip,gzip; default=pbzip2) -threadNum run_with_multiple_threads_if_possible(default=1) -compress compress_files(default=uncompress) -tarball name_for_tarball(default=tarball.tar.gz) -workdir special_handling_for_pipeline_archiving -bam create_bam_archive_for_pipeline(default=all_but_bam) [options]

Required flags: -files || -dir

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

compress.pl -> compress or decompress files

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./compress.pl -files /home/matt/read1_mouse54.fq,/home/matt/read2_mouse54.fq  -compress -suffix tar.gz -tarball test.tar.gz -keep
./compress.pl -files /home/matt/read1_mouse54.fq,/home/matt/read2_mouse54.fq -threadNum 5 -compress
./compress.pl -files /home/matt/read1_mouse54.fq,/home/matt/read2_mouse54.fq -threadNum 5  
./compress.pl -files /home/matt/test -suffix tar.gz 
./compress.pl -dir /home/matt/testdir -compress -suffix tar.gz -tarball test.tar.gz


=cut

my $threadNum = defined $OPT{threadNum}?$OPT{threadNum}:1;
my $tarball_name = defined $OPT{tarball}?$OPT{tarball}:"tarball.tar.gz";
#Add the runid to the final tarballs
if ($OPT{runid}) {
	$tarball_name =~ s/RUN/$OPT{runid}/;
}
#Indicates whether we're compressing or uncompressing
my $compress = defined $OPT{compress}?'compress':'uncompress';
my $suffix = defined $OPT{suffix}?$OPT{suffix}:'bz2';
my $run = defined $OPT{run}?$OPT{run}:1;

#Dummy flag to have the step 'run' in the pipeline
if (!$run) {
	exit;
}

if ($OPT{dir} && $OPT{files}) {
	modules::Exception->throw("ERROR: Cannot run script with both -files and -dir");
}

my $sys_call = modules::SystemCall->new();


my %suffix_mapper = (
						bzip2 => 'bzip2',
						bz2 => 'bzip2',
						gzip => 'gzip',
						gz => 'gzip',
						'tar.gz' => 'tar',
						tgz => 'tar',
						tar => 'tar',
					);
my @files;
my $file_str;
#my $outdir = `pwd`;
#chomp $outdir;
my $rundir;
my $basedir;

if ($OPT{dir}) {
	if ( $suffix_mapper{$suffix} ne 'tar' || $compress eq 'uncompress') {
		modules::Exception->throw("ERROR: -dir argument only works with tar and uncompress");
	}
	if (!-d $OPT{dir}) {
		modules::Exception->throw("ERROR: Directory $OPT{dir} doesn't exist");
	}
	
	if ($OPT{workdir}) {
		#Move the files to final subdirectories if we're archiving a run directory
		(my $fulldir = $OPT{dir}) =~ s/\/$//;
		($rundir,$basedir) = fileparse($fulldir); 
		if (-e "$fulldir/$tarball_name") {
			$sys_call->run("rm $fulldir/$tarball_name");
		}
		modules::Pipeline::move_run_files($fulldir);
		
		#Here we only tar up the specific output directories
		if ($OPT{bam}) {
			$file_str = "$rundir/bam";
		} else {
			#Here we get the sub directories to compress from the config file
			my $pipe_config = modules::Pipeline->get_pipe_conf();
			my @run_subdirs = split(",",$pipe_config->read('run_subdirs'));
			for my $subdir (@run_subdirs) {
				next if $subdir =~ /bam/; #Bams don't need to be compressed
				$file_str .= "$rundir/$subdir ";
			}
		}
	} else {
		#Here we just tar up the directory as per normal
		$file_str = $OPT{dir};
	}
	
} else {
	($file_str = $OPT{files}) =~ s/,/ /g;
	@files = split(",",$OPT{files});
	#Add the appropriate suffix if were uncompressing
	if ($compress eq 'uncompress') {
		$file_str =~ s/ /\.$suffix /g;
		$file_str = $file_str . '.' . $suffix;
	}
}


#if ( !-d $outdir ) {
#	modules::Exception->throw("Tarball directory $outdir doesn't exist");	
#}


my %compression_binaries = (
								bzip2 => {
										compress => 'bzip2',
										compress_args => "-k -z $file_str",
										compress_suffixes => ['bz2','bzip2'],
										uncompress => 'bzip2',
										uncompress_args => "-k -d $file_str"
										},
								pbzip2 => {
										compress => 'pbzip2',
										compress_args => "-k -z -p$threadNum $file_str",
										compress_suffixes => ['bz2','bzip2'],
										uncompress => 'pbzip2',
										uncompress_args => "-k -d -p$threadNum $file_str"
										},
								gzip => {
										compress => 'gzip',
										compress_args => "$file_str",
										compress_suffixes => ['gz','gzip'],
										uncompress => 'gzip',
										uncompress_args => "-d $file_str"
										},
								pigz => {
										compress => 'pigz',
										compress_args => "-k -p $threadNum $file_str",
										compress_suffixes => ['gz','gzip'],
										uncompress => 'gzip',
										uncompress_args => "-k -d -p $threadNum $file_str"	
										},
								tar => {
										compress => 'tar',
										compress_args => "czvf $rundir/$tarball_name $file_str",
										compress_suffixes => ['tgz','tar.gz'],
										uncompress => 'tar',
										uncompress_args => "xzvPf $file_str"
										}
							);

my $type;
my $binary;

#Check the type is handled
if (exists $suffix_mapper{$suffix}) {
	$type = $suffix_mapper{$suffix};
	$binary = $compression_binaries{$type}{$compress};
	
	if ($type eq 'bzip2' && $threadNum > 1) {
		$type = 'pbzip2';
		#Hack to get pbzip2 to work on nodes where it's not installed
		if (&pbzipInstalled()) {
			$binary = 'pbzip2';
		} else {
			$binary = $ENV{SVNDIR}.'/ext/bin/pbzip2';
		}
	}
	
	if ($type eq 'gzip' && $threadNum > 1) {
		$type = 'pigz';
		#Hack to get pigz to work on nodes where it's not installed
		if (&pigzInstalled()) {
			$binary = 'pigz';
		} else {
			$binary = $ENV{SVNDIR}.'/ext/bin/pigz';
		}
	}
	
} else {
	my $suffix_str = join(",",keys %suffix_mapper);
	modules::Exception->throw("Suffix must be one of $suffix_str");
}

#Get the arguments
my $compress_args_key = $compress."_args";
my $compress_command = $binary . ' ' .$compression_binaries{$type}{$compress_args_key};

#Add the proper suffixes if we're uncompressing and check file doesn't exist if we're compressing
for my $file ( @files ) {
	
	if ($compress eq 'uncompress') {
		$file .= '.' . $suffix;
	}
	
	if ( !-e $file ) {
		modules::Exception->throw("File $file doesn't exist");	
	}
	
	if ($compress eq 'compress') {
		for my $suffix ( @{$compression_binaries{$type}{compress_suffixes}} ) {
		    if ($file =~ /$suffix$/) {
				modules::Exception->throw("File $file already has compression suffix $suffix");	
		    }
		}
	}	    
}

#Need to change to the relative directory to avoid having big ugly paths in our tarballs
if ($OPT{workdir}) {
	$compress_command = "cd $basedir;" . $compress_command;
}

$sys_call->run($compress_command);

#After the compression/uncompression remove the input files if they still exist
unless ($OPT{dir} || $OPT{keep}) {
	for my $file ( @files ) {
		if ( -e $file ) {
			$sys_call->run("rm $file");
		}
	}
}

if ($OPT{mv} && -e $tarball_name) {
	`mv $tarball_name $OPT{dir} 2>/dev/null`;
}

sub pbzipInstalled {
	my $pbout = `which pbzip2 2>/dev/null`;
	if ($pbout) {
		return 1;
	} else {
		return 0;
	}
}

sub pigzInstalled {
	my $pgout = `which pigz 2>/dev/null`;
	if ($pgout) {
		return 1;
	} else {
		return 0;
	}
}


							