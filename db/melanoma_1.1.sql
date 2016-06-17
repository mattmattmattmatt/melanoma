CREATE TABLE patients (
id int unsigned AUTO_INCREMENT,
patient_external_name varchar(45), -- e.g. tb412
KEY id (id), 
PRIMARY KEY(id)
) ENGINE=InnoDB;

CREATE TABLE sample_groups (
id int unsigned AUTO_INCREMENT,
total_samples int unsigned, -- total number of samples in the sample group
sample_group_number int unsigned, -- number for this sample group; usually only one per patient but pilot has three
sample_group_name varchar(45), -- e.g tb412_sg1, A15_sg2
patient_id int unsigned,
KEY id (id), 
FOREIGN KEY (patient_id) REFERENCES patients(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY(id)
) ENGINE=InnoDB;

CREATE TABLE samples (
id int unsigned AUTO_INCREMENT,
sample_number int unsigned, -- Internal counter within sample group
tumour bool, -- Flag indicating if it's a tumour sample
cell_line bool, -- Flag indicating if it's a cell_line
tumour_type ENUM('primary','metastatic'), -- Tumour type
total_lanes int unsigned, -- Total lanes for sample; this can change
tissue varchar(45), -- Input tissue for sample
sample_type ENUM('normal_cell_line','normal_patient','tumour_primary_cell_line','tumour_primary_patient','tumour_metastatic_cell_line','tumour_metastatic_patient'),
sample_name varchar(45), -- e.g tb412_sg1_tumour3_s1, A15_sg2_normal1_s1
bpa_id varchar(100),
sample_group_id int unsigned,
KEY id (id), 
FOREIGN KEY (sample_group_id) REFERENCES sample_groups(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY(id)
) ENGINE=InnoDB;


CREATE TABLE sequencing_centres (
id int unsigned AUTO_INCREMENT,
sequencing_centre_name ENUM('BRF','AGRF','RAM'),
KEY id (id), 
PRIMARY KEY(id)
) ENGINE=InnoDB;

CREATE TABLE lanes (
id int unsigned AUTO_INCREMENT,
lane_number int unsigned,
lane_name varchar(100), -- e.g. tb412_sg1_tumour3_s1_l2
sequencing_centre_id int unsigned,
sample_id int unsigned,
KEY id (id),
FOREIGN KEY (sequencing_centre_id) REFERENCES sequencing_centres(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (sample_id) REFERENCES samples(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY(id)
) ENGINE=InnoDB;

CREATE TABLE read_files (
id int unsigned AUTO_INCREMENT,
file_name  varchar(100),
is_compressed bool,
compression_suffix varchar(40), 
read_file_number int unsigned,
read_directory varchar(500),
read_file_name varchar(45), -- e.g. tb412_sg1_tumour3_s1_l2_r2
lane_id int unsigned,
FOREIGN KEY (lane_id) REFERENCES lanes(id) ON DELETE CASCADE ON UPDATE CASCADE,
KEY id (id),
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE align_lanes (
id int unsigned AUTO_INCREMENT,
quality_encoding ENUM('phred33','phred64'),
lane_id int unsigned,
KEY id (id),
FOREIGN KEY (lane_id) REFERENCES lanes(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE runs (
id int unsigned AUTO_INCREMENT,
run_directory varchar(300),
status ENUM('complete','inprocess','stopped'),
commenced DATETIME, -- populate this with the mysql NOW() command
production bool, -- flag indicates whether this is the latest production run to use for stats, etc
sample_id int unsigned,
KEY id (id),
FOREIGN KEY (sample_id) REFERENCES samples(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE lanes_runs (
id int unsigned AUTO_INCREMENT,
run_id int unsigned,
lane_id int unsigned,
KEY id (id),
FOREIGN KEY (lane_id) REFERENCES lanes(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE pipeline_steps (
id int unsigned AUTO_INCREMENT,
name varchar(40),
description varchar(300),
KEY id (id),
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE pipeline_steps_runs (
id int unsigned AUTO_INCREMENT,
run_id int unsigned,
pipeline_step_id int unsigned,
KEY id (id),
FOREIGN KEY (pipeline_step_id) REFERENCES pipeline_steps(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE syscalls (
id int unsigned AUTO_INCREMENT,
command varchar(1000),
level ENUM('lane','sample','sample_group'),
analysis_step_name varchar(45),
align_lane_id int unsigned,
pipeline_steps_run_id int unsigned,
chr varchar(45),
KEY id (id),
FOREIGN KEY (align_lane_id) REFERENCES align_lanes(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (pipeline_steps_run_id) REFERENCES pipeline_steps_runs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY(id)
) ENGINE=InnoDB;

CREATE TABLE release_files (
id int unsigned AUTO_INCREMENT,
file_name varchar(100),
file_type ENUM('vcf','bam','summary'),
total_lanes int unsigned, -- Way to allow for running separate analyses when additional lanes are added
pipeline_steps_run_id int unsigned,
KEY id (id),
FOREIGN KEY (pipeline_steps_run_id) REFERENCES pipeline_steps_runs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE filters (
id int unsigned AUTO_INCREMENT,
name varchar(45),
description varchar(100),
pipeline_step_id int unsigned, 
KEY id (id),
FOREIGN KEY (pipeline_step_id) REFERENCES pipeline_steps(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE snvs (
id int unsigned AUTO_INCREMENT,
chr varchar(40),
coord int unsigned,
ref_base varchar(4),
var_base varchar(4),
ref_base_freq float(2,2),
var_base_freq float(2,2),
median_quality_score float(4,2),
tumour_base_string varchar(10000),
normal_base_string varchar(10000),
snv_score float(5,2),
clr_score float(5,2),
read_depth int unsigned,
snv_class ENUM('SOMATIC','GERMLINE','LOH','COMPLEX','UNKNOWN') default 'SOMATIC',
run_id int unsigned,
KEY id (id),
FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE snvs_filters (
id int unsigned AUTO_INCREMENT,
filtermatch bool, -- flag indicates whether filter overlaps with filter
filterpass bool, -- flag indicates whether this is considered a fail or a pass; this is filter specific
attribute varchar(1000),
filter_id int unsigned,
snv_id int unsigned,
KEY id (id),
FOREIGN KEY (filter_id) REFERENCES filters(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (snv_id) REFERENCES snvs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE variants (
id int unsigned AUTO_INCREMENT,
chr varchar(40),
start_coord int unsigned,
end_coord int unsigned,
var_type ENUM('DEL','INS'),
inserted_bases varchar(10000),
var_score float(5,2),
clr_score float(5,2),
read_depth int unsigned,
ref_base_freq float(2,2),
var_base_freq float(2,2),
tumour_base_string varchar(10000),
normal_base_string varchar(10000),
median_quality_score float(4,2),
var_caller ENUM('dindel','samtools_indel'),
variant_class ENUM('SOMATIC','GERMLINE','LOH','COMPLEX','UNKNOWN') default 'SOMATIC',
run_id int unsigned,
KEY id (id),
FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE variants_filters (
id int unsigned AUTO_INCREMENT,
filtermatch bool, -- flag indicates whether filter overlaps with filter
filterpass bool, -- flag indicates whether this is considered a fail or a pass; this is filter specific
attribute varchar(1000),
filter_id int unsigned,
variant_id int unsigned,
KEY id (id),
FOREIGN KEY (filter_id) REFERENCES filters(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (variant_id) REFERENCES variants(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE snv_rows (
id int unsigned AUTO_INCREMENT,
snv_id int unsigned,
run_id int unsigned,
KEY id(id),
FOREIGN KEY (snv_id) REFERENCES snvs(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE, 
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE snvrow_values (
id int unsigned AUTO_INCREMENT,
column_number int unsigned,
column_name varchar(40),
column_value varchar(1000),
snv_row_id int unsigned not null,
KEY id(id),
FOREIGN KEY (snv_row_id) REFERENCES snv_rows(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE variant_rows (
id int unsigned AUTO_INCREMENT,
variant_id int unsigned,
run_id int unsigned,
KEY id(id),
FOREIGN KEY (variant_id) REFERENCES variants(id) ON DELETE CASCADE ON UPDATE CASCADE,
FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE ON UPDATE CASCADE, 
PRIMARY KEY (id)
) ENGINE=InnoDB;

CREATE TABLE variantrow_values (
id int unsigned AUTO_INCREMENT,
column_number int unsigned,
column_name varchar(40),
column_value varchar(1000),
variant_row_id int unsigned,
KEY id(id),
FOREIGN KEY (variant_row_id) REFERENCES variant_rows(id) ON DELETE CASCADE ON UPDATE CASCADE,
PRIMARY KEY (id)
) ENGINE=InnoDB;


CREATE INDEX snv_run_chr_coord_idx ON snvs (run_id, chr, coord);
CREATE INDEX snvfilter_snpid_idx ON snvs_filters (snv_id);
CREATE INDEX snvfilter_filterid_idx ON snvs_filters (filter_id);
CREATE INDEX variantfilter_variantid_idx ON variants_filters (variant_id);
CREATE INDEX variantfilter_filterid_idx ON variants_filters (filter_id);

