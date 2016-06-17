#Create the pipeline steps
Insert into pipeline_steps VALUES ('10','merge_bam','Merge the bam lane files');
Insert into pipeline_steps VALUES ('20','remove_duplicates','Remove duplicate reads from aligned BAM file');
Insert into pipeline_steps VALUES ('25','copy_bam','Remove duplicate reads from aligned BAM file');
Insert into pipeline_steps VALUES ('30','bam_stats','Generate read stats regarding the alignments');
Insert into pipeline_steps VALUES ('40','exon_depth','Report coverage stats for all the exon regions');
Insert into pipeline_steps VALUES ('50','submit_snvs','Submits all the individual chromosome snv calling jobs');
Insert into pipeline_steps VALUES ('60','single_vcf','Job that checks the vcf and records as done in db');
Insert into pipeline_steps VALUES ('70','merge_vcf','Merge the vcf files from the chr by chr results');
Insert into pipeline_steps VALUES ('80','filter_vcf','Filter the snv and indels based on cutoffs');
Insert into pipeline_steps VALUES ('90','db_variants','Classify the variant calls into heterozygous, homozygous, or other and insert into the database');
Insert into pipeline_steps VALUES ('100','filter_exon','Flag snps matching exons');
Insert into pipeline_steps VALUES ('110','filter_exon_ns','Flag snps matching non-synonomous snps changes in the protein');
Insert into pipeline_steps VALUES ('120','filter_dbsnp_snv','Flag snps matching dbsnp entries');
Insert into pipeline_steps VALUES ('130','filter_dbsnp_indel','Flag indels matching dbsnp entries');
Insert into pipeline_steps VALUES ('140','filter_splicesite','Flag snps matching splice sites of exons');
Insert into pipeline_steps VALUES ('150','filter_polyphen','Run polyphen on passed snps');
Insert into pipeline_steps VALUES ('160','filter_cosmic','Overlap COSMIC cancer mutant');
Insert into pipeline_steps VALUES ('170','report_snv','Generate a full snp summary for each snp reported');
Insert into pipeline_steps VALUES ('180','report_indel','Generate a full indel summary for each indel reported');
Insert into pipeline_steps VALUES ('190','archive','Archive the output files');
Insert into pipeline_steps VALUES ('200','scp_results','Scp resulting tarball to a remote server');

#Create the filters
Insert into filters VALUES ('10','filter_deletion','Identifies a deletion for the sample','90');
Insert into filters VALUES ('20','filter_insertion','Identifies an insertion for the sample','90');
Insert into filters VALUES ('30','filter_snv','Identifies a snv for a sample','90');
Insert into filters VALUES ('40','filter_exon','Snps that overlap exons','100');
Insert into filters VALUES ('50','filter_exon_ns','Snps that overlap exons and cause non-synonomous changes','110');
Insert into filters VALUES ('60','filter_dbsnp_snv','Snps that match dbsnp snv entries','120');
Insert into filters VALUES ('70','filter_dbsnp_indel','Snps that match indel dbsnp entries','130');
Insert into filters VALUES ('80','filter_splicesite','Snps that match splice sites','140');
Insert into filters VALUES ('90','filter_polyphen','Snps that are run through polyphen','150');
Insert into filters VALUES ('100','filter_cosmic','SNPs that overlap COSMIC database','160');
Insert into filters VALUES ('110','filter_pass_snv','Snps that pass according to run parameters','170');
Insert into filters VALUES ('120','filter_causal_snv','Identifies the causal snp for the sample','170');
Insert into filters VALUES ('130','filter_rareallele_snv','Identifies a human snv that is a rare allele in dbsnp','170');
Insert into filters VALUES ('140','filter_pass_indel','Indels that pass according to run parameters','180');
Insert into filters VALUES ('150','filter_causal_indel','Identifies the causal indel for the sample','180');
Insert into filters VALUES ('160','filter_rareallele_indel','Identifies a human indel that is a rare allele in dbsnp','180');


#Create the sequencing centre ids
Insert into sequencing_centres VALUES ('1','BRF');
Insert into sequencing_centres VALUES ('2','AGRF');
Insert into sequencing_centres VALUES ('3','RAM');

