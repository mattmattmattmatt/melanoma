#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "input=s",
	   "debug=s",
	   "outdir=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{input});



=pod

=head1 SYNOPSIS

parse_ucsc_dbsnp.pl -input dbsnp_file -outdir output_directory(default=pwd) [options]

Required flags: -input

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_ucsc_dbsnp.pl -> Script to parse ucsc dbsnp data

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./parse_ucsc_dbsnp.pl -input /home/matt/work/data/seqdata/mouse/mm9/snps/dbsnp/snp128.txt -outdir /home/matt/work/pipeline/conf/mouse/NCBIM37/dbsnp/128
cd /home/matt/work/pipeline/conf/mouse/NCBIM37/dbsnp/128
for f in *; do for chr in {1..19} X Y; do grep -w ^$chr $f > NCBIM37.$f.$chr; done; done

=cut

my $dbsnp = $OPT{input};
my $outdir = defined $OPT{outdir}?`readlink -f $OPT{outdir}`:`pwd`;
chomp $outdir;
my $count = 0;
my %final_snps = ();

open(SNV,">$outdir/dbsnp.snv") || die "Can't open file to write $outdir/snv\n";
open(INDEL,">$outdir/dbsnp.indel") || die "Can't open file to write $outdir/indel\n";

my @fh = (\*SNV,\*INDEL);

#Global variables
my %snv_entries = ();
my %allele_count = ();
my $end_coord;
my $local_chr;

#now parse the entries
open(DBSNP,"$dbsnp") || die "Can't open file $dbsnp\n";
while (<DBSNP>) {
	chomp;
	if ($count %1000000 == 0) {
		print "$count $_\n";
	}
	$count++;
    if ($OPT{debug}) {
		next unless $_ =~ /$OPT{debug}/;
		print "$_\n";
	}
	#print "$_\n";
	
	my %line_allele_count = ();
	
	
    #next unless $_ =~ /chrY/;
    my @fields = split("\t");
    my $rs = $fields[4];
    (my $chr = $fields[1]) =~ s/chr//;
    next if $chr =~ /random/;
    my $start = $fields[2];
    my $end = $fields[3];
    
#	if ($chr != 1) {
#		last;
#	}
	
#	if ($end_coord > 31733) {
#		last;
#	}

    my $strand = $fields[6];
    my $refbase = $fields[8];
    #Skip the really large cases
    next if $fields[9] =~ /\(/;
    next if $fields[9] =~ /lengthTooLong/;
    next if $fields[8] =~ /\(/;
    
    #reset the coord count when hit a new chromosome
    if ($chr ne $local_chr) {
		$end_coord = 0;
    }
    
    if ($start >= $end_coord) {
		#Don't put snvs in data structure as becomes too large; need to collapse duplicate entries 
       	#Write new entries when we move to a new coordinate; this makes sure we have all the allele info
		my $fh_snv = $fh[0];
		for my $localchr ( keys %snv_entries ) {
		    for my $local_coord ( sort {$a<=>$b} keys %{$snv_entries{$localchr}} ) {
				for my $local_ref (keys %{$snv_entries{$localchr}{$local_coord}}) {
					my $ref_count = defined $allele_count{$local_ref}?$allele_count{$local_ref}:0;
					for my $local_var (keys %{$snv_entries{$localchr}{$local_coord}{$local_ref}}) {
						next if $local_var =~ /^total/;
						my $allele_change = "$local_ref(n/a)->$local_var(n/a)";
						#print "$localchr $local_coord $local_ref $local_var ", Dumper \%snv_entries;
						my $rs_string = join(",",@{$snv_entries{$localchr}{$local_coord}{$local_ref}{$local_var}{rs}});
						
						if (exists $snv_entries{$localchr}{$local_coord}{$local_ref}{$local_var}{varcount}) {
							my $total_varcount = &Get_Total_Vars($snv_entries{$localchr}{$local_coord}{$local_ref});
							my $total = $ref_count + $total_varcount;
							my $varcount = $snv_entries{$localchr}{$local_coord}{$local_ref}{$local_var}{varcount};
							
						
							
							my $ref_freq =  $ref_count==0?0:sprintf("%.3f",$ref_count/$total);
							my $var_freq =  $varcount==0?0:sprintf("%.3f",$varcount/$total);
							$allele_change = "$local_ref($ref_freq)->$local_var($var_freq)";
#							if ($local_coord == 12049 || $local_coord == 39879 ) {
#								print Dumper \%allele_count;
#								print Dumper \%snv_entries;
#								print "VAR $varcount REF $ref_count TOTAL_VAR $total_varcount TOTAL $total\n";				
#								print "$localchr $local_coord $local_coord SNV^^^$rs_string^^^$allele_change\n";
#							}
						}
						#print "$localchr $local_coord $local_coord SNV^^^$rs_string^^^$allele_change\n";
						print $fh_snv "$localchr $local_coord $local_coord SNV^^^$rs_string^^^$allele_change\n";
					}
				}
		    }
		}
    	%allele_count = ();
    	%snv_entries = ();
	}
    
    
    
    my $freq_flag_local = 0;
    
    #If there's allele frequency info
    if ($fields[22] =~ /,/) {
    	$freq_flag_local = 1;
    	my @freq_alleles = ();
    	if ( $strand eq '-' ) {
	    	@freq_alleles = &Complement($fields[22],',');
	    } else {
	    	@freq_alleles = split(",",$fields[22]);
	    }

		#print Dumper \@freq_alleles;

		my @counts = split(",",$fields[23]);
    	for ( my $tmp = 0 ; $tmp < @freq_alleles ; $tmp++ ) {
    	    $allele_count{$freq_alleles[$tmp]} += $counts[$tmp];
    	    $line_allele_count{$freq_alleles[$tmp]} += $counts[$tmp];
    	}
    }
    
    #Handle the separate indel events
    my @alleles = ();

    #Put all entries on positive strand
    if ( $strand eq '-' ) {
    	@alleles = &Complement($fields[9],'/');
    } else {
    	@alleles = split("/",$fields[9]);
    }
    
    
    #Deal with the remaining cases one by one
    for my $variant (@alleles) {
       next if $refbase eq $variant;
                        
             
       if ($variant =~ /^[ACGT]$/ && $refbase =~ /^[ACGT]$/) {
       		#Regular SNV like A->T
       		#print "SIMPLE SNV $refbase->$variant $start $end LINE $fields[9] $chr $start $rs\n";       		

			push @{$snv_entries{$chr}{$end}{$refbase}{$variant}{rs}},$rs;
				
			if ($freq_flag_local) {
				$snv_entries{$chr}{$end}{$refbase}{$variant}{varcount} = defined $allele_count{$variant}?$allele_count{$variant}:0;
			}
       }  elsif ($refbase eq '-') {
       		#Simple insertion like '-'->A
       		#print "SIMPLE INS $refbase->$variant $start $end LINE $fields[9] $chr $start $rs\n";
       		
       		next if $variant eq '-';
	    	push @{$final_snps{$chr}{$start}{$end}{INS}{$variant}{rs}},$rs;

			if ($freq_flag_local) {
				$final_snps{$chr}{$start}{$end}{INS}{$variant}{varcount} = defined $line_allele_count{$variant}?$line_allele_count{$variant}:0;
				$final_snps{$chr}{$start}{$end}{INS}{$variant}{refcount} = defined $line_allele_count{$refbase}?$line_allele_count{$refbase}:0;
			} 
			
       } elsif ($variant eq '-') {
       		#Simple deletion like A->'-'
       		#print "SIMPLE DEL $refbase->$variant $start $end LINE $fields[9] $chr $start $rs\n";
       		
       		#Few of these cases which don't make sense
       		if ($start == $end) {
       			next;
       		}
       		
			push @{$final_snps{$chr}{$start+1}{$end}{DEL}{$refbase}{rs}},$rs;
       		
			if ($freq_flag_local) {
		    	$final_snps{$chr}{$start+1}{$end}{DEL}{$refbase}{varcount} = defined $line_allele_count{$variant}?$line_allele_count{$variant}:0;
		    	$final_snps{$chr}{$start+1}{$end}{DEL}{$refbase}{refcount} = defined $line_allele_count{$refbase}?$line_allele_count{$refbase}:0;
			}
	   } elsif (length($variant) < length($refbase) && $refbase =~ /^$variant/) {
	   		#Complex deletion like TA -> T; bases must agree
	   		#print "COMPLEX DEL $refbase->$variant $start $end LINE $fields[9] $chr $start $rs\n";
	   		my $new_ref = substr($refbase,length($variant),length($refbase));
	   		my $new_start = $start + length($variant);
			push @{$final_snps{$chr}{$new_start+1}{$end}{DEL}{$new_ref}{rs}},$rs;

			if ($freq_flag_local) {
				$final_snps{$chr}{$new_start+1}{$end}{DEL}{$new_ref}{varcount} = defined $line_allele_count{$variant}?$line_allele_count{$variant}:0;
				$final_snps{$chr}{$new_start+1}{$end}{DEL}{$new_ref}{refcount} = defined $line_allele_count{$refbase}?$line_allele_count{$refbase}:0;
			} 
	   		
	   		
	   } elsif (length($variant) > length($refbase) && $variant =~ /^$refbase/) {
       		#Complex insertion like A -> AT; bases must agree
       		#print "COMPLEX INS $refbase->$variant LINE $fields[9] $chr $start $end $rs\n";
       		
			#Here we adjust the inserted bases to account for common bases between insertion
			my $new_variant = substr($variant,length($refbase),length($variant));
			my $new_start = $start + length($refbase);
			push @{$final_snps{$chr}{$new_start}{$end}{INS}{$new_variant}{rs}},$rs;

			if ($freq_flag_local) {
				$final_snps{$chr}{$new_start}{$end}{INS}{$new_variant}{varcount} = defined $line_allele_count{$variant}?$line_allele_count{$variant}:0;
				$final_snps{$chr}{$new_start}{$end}{INS}{$new_variant}{refcount} = defined $line_allele_count{$refbase}?$line_allele_count{$refbase}:0;
			} 
			  		
       } elsif (length($variant) == length($refbase) && $refbase ne $variant) {
       		#Complex SNV AA->TG
       		#print "SNV COMPLEX $refbase->$variant LINE $fields[9] $chr $start $end $rs\n";
       		
       		#Few of these cases which don't make sense
       		if ($start == $end) {
       			next;
       		}
       		
       		my $fh = $fh[0];
       		my @ref_alleles = split("",$refbase);
			my @new_alleles = split("",$variant);
			
			my $index = 0;
			for ( my $newstart = $start+1 ; $newstart <= $end ; $newstart++ ) {
				if ($ref_alleles[$index] ne $new_alleles[$index]) {
					#print "DIFF NS $newstart E $end I $index REF $ref_alleles[$index] VAR $new_alleles[$index]\n";
					push @{$snv_entries{$chr}{$newstart}{$ref_alleles[$index]}{$new_alleles[$index]}{rs}},$rs;
					if ($freq_flag_local) {
						$snv_entries{$chr}{$newstart}{$ref_alleles[$index]}{$new_alleles[$index]}{varcount} = defined $allele_count{$new_alleles[$index]}?$allele_count{$new_alleles[$index]}:0;
					}					
				}
				$index++;				
			}
			
			
       } else {
       		#print "NOT PROCESS $refbase->$variant $fields[9] $chr $start $rs\n";
       }
    }
    
 
    $end_coord = $end if $end >  $end_coord;
    $local_chr = $chr;
    print Dumper \%final_snps if $OPT{debug};
}

close $fh[0];
#print Dumper $final_snps{1}{42777};
#print Dumper \%final_snps;

#Report the events to files
my $fh_indel = $fh[1];
for my $chr ( sort keys %final_snps ) {
    for my $start ( sort {$a<=>$b} keys %{$final_snps{$chr}} ) {
    	for my $end ( sort {$a<=>$b} keys %{$final_snps{$chr}{$start}} ) {
    		if (exists $final_snps{$chr}{$start}{$end}{INS}) {
	    		
				for my $inserted_bases (keys %{$final_snps{$chr}{$start}{$end}{INS}}) {
		    		my $ref_count = defined $final_snps{$chr}{$start}{$end}{INS}{$inserted_bases}{refcount}?$final_snps{$chr}{$start}{$end}{INS}{$inserted_bases}{refcount}:0;
					my $rs_string = join(",",@{$final_snps{$chr}{$start}{$end}{INS}{$inserted_bases}{rs}});
				    my $variant_str = $inserted_bases.'(n/a)';
				    	
				    if (exists $final_snps{$chr}{$start}{$end}{INS}{$inserted_bases}{varcount}) {
				    	my $total_varcount = &Get_Total_Vars($final_snps{$chr}{$start}{$end}{INS});
				    	my $total = $total_varcount + $ref_count;
				    	my $varcount = $final_snps{$chr}{$start}{$end}{INS}{$inserted_bases}{varcount};
				    	my $ref_freq =  $ref_count==0?0:sprintf("%.3f",$ref_count/$total);
						my $var_freq =  $varcount==0?0:sprintf("%.3f",$varcount/$total);
						$variant_str = $inserted_bases . '('. $ref_freq . ':'. $var_freq .')';
				    }
				    
				    #print "$chr $start $end INS^^^$rs_string^^^$variant_str\n";
		       		print $fh_indel "$chr $start $end INS^^^$rs_string^^^$variant_str\n";
		        	
	    		}
    		
    		}
    		
    		if (exists $final_snps{$chr}{$start}{$end}{DEL}) {
    			for my $deleted_bases (keys %{$final_snps{$chr}{$start}{$end}{DEL}}) {
		    		my $ref_count = defined $final_snps{$chr}{$start}{$end}{DEL}{$deleted_bases}{refcount}?$final_snps{$chr}{$start}{$end}{DEL}{$deleted_bases}{refcount}:0;
					#print "$chr $start $end \n", Dumper $final_snps{$chr}{$start}{$end}{DEL};
					my $rs_string = join(",",@{$final_snps{$chr}{$start}{$end}{DEL}{$deleted_bases}{rs}});
				    my $variant_str = $deleted_bases.'(n/a)';
				    	
				    if (exists $final_snps{$chr}{$start}{$end}{DEL}{$deleted_bases}{varcount}) {
				    	my $total_varcount = &Get_Total_Vars($final_snps{$chr}{$start}{$end}{DEL});
				    	my $total = $total_varcount + $ref_count;
				    	my $varcount = $final_snps{$chr}{$start}{$end}{DEL}{$deleted_bases}{varcount};
				    	my $ref_freq =  $ref_count==0?0:sprintf("%.3f",$ref_count/$total);
						my $var_freq =  $varcount==0?0:sprintf("%.3f",$varcount/$total);
						$variant_str = $deleted_bases . '('. $ref_freq . ':'. $var_freq .')';
				    }
				    
				    #print "$chr $start $end DEL^^^$rs_string^^^$variant_str\n";
		       		print $fh_indel "$chr $start $end DEL^^^$rs_string^^^$variant_str\n";
		        	
	    		}
			    
    		}
    		
    	}
	}
}



close $fh[1];

sub Get_Total_Vars {
	my ($data) = @_;
	my $total = 0;
	for my $var_allele (keys %{$data}) {
		$total += $data->{$var_allele}{varcount};
	}
	return $total;
}


#Reverse Complement the bases in an array
sub Complement {
    my ($bases,$delim) = @_;
    my @old_alleles = split($delim,$bases);
    my @new_alleles = ();
    for my $old_allele ( @old_alleles ) {
    	my $rev_allele = reverse $old_allele;
    	(my $new_allele = $rev_allele) =~ tr/ACTG/TGAC/;
    	push @new_alleles, $new_allele;        
    }
    return @new_alleles;
}




