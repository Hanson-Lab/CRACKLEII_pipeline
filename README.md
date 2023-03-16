# CRACKLE_II data processing and resistance analysis pipeline

## Version - 1.3

## Program requirements
The programs required for this pipeline, and their respective installation instructions can be found here:
* Java - https://java.com/en/ 
	* version - 1.8.0 or newer 
		* Pipeline tested with version 1.8.0_252
* Python - https://www.python.org
	* version - 3.8.5
* Perl - https://www.perl.org
	* version - 5.26 or newer
		* Pipeline tested with version 5.30.2
	* packages - Moo, List::MoreUtils, JSON, LWP::Simple, Bio::Perl, Path::Tiny
* R - https://www.r-project.org
	* version - 3.6.3
	* packages - XML, dplyr, biostrings, ape
	* note: this is not included in the gitlab as there are binaries for multiple linux distributions. You can download R from the R-project website (https://www.r-project.org), make sure it is version 3.6.3. 
* Samtools - http://www.htslib.org 
	* version - 1.10
* bcftools - http://www.htslib.org
	* version - 1.10.2
* EMBOSS - http://emboss.sourceforge.net/download/
	* version - 6.6.0
* Seqtk - https://github.com/lh3/seqtk
	* version - 1.3 (r106)
* NCBI-blast+ - https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
	* version - 2.10.1+
* Raspberry - https://github.com/CEG-ICRISAT/Raspberry
	* version - 0.3
* Quast - http://quast.sourceforge.net
	* version - 5.0.2
	* database(s) - ss_db_w32_4324 (update 10-feb-2016) & ss_db_w16 (update 21-feb-2017)
* Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic
	* version - 0.39
* StrainSeeker - http://bioinfo.ut.ee/strainseeker/index.php?r=site/page&view=downloadable
	* version - 1.5
* Spades - http://cab.spbu.ru/software/spades/
	* version - 3.14.0
* Bowtie 2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* version - 2.4.1
* MLST - https://github.com/tseemann/mlst
	* version - 2.19.0
	* database(s) - internal blast and pubmlst databases updated 2020-02-23
* Kleborate - https://github.com/katholt/Kleborate
	* version - 1.0.0
* Abricate - https://github.com/tseemann/abricate
	* version - 1.0.1


## Pipeline schematic

![Pipeline schematic](CRACKLE_II_pipeline.png)


## Pipeline installation notes
**A few programs used within this pipeline require configuration and installation. Follow the instructions below before filling out the configuration file**

SAMTOOLS

```

cd /path/to/samtools-1.8
./configure
make
```

BCFTOOLS
```

cd /path/to/bcftools-1.8
./configure
make
```

SEQTK
Note: this may give a warning that variable lc is not set. This is safe to ignore. 
```

cd /path/to/seqtk
make
```

## Pipeline update notes (only needed when creating a new version release)
**The in_house database needs to be formatted for abricate to access it. This must be done after abricate is installed. There are detailed instructions on the abricate github page, but the following instructions were used during pipeline development:**
```
cd /path/to/abricate/db
mkdir in_house_db
cd in_house_db
cp /path/to/in_house_db.fasta sequences
abricate --setupdb
#Test if install worked - check if in_house_db is present in list of abricate databases
abricate --list
```

The complete lists of genes for card, plasmidfinder, and in_house_db for the XML_generator were created using the following code:

NOTE: These are pre-populated for each pipeline iteration and do NOT need to be created.

NOTE: These piped commands will not work on MacOS, only on Linux OS

CARD
```
cd /path/to/abricate/db/card/
grep ">" sequences | cut -d'~' -f4,7 | tr '~' '\t' > ../card_db.tsv
```

PlasmidFinder
```
cd /path/to/abricate/db/plasmidfinder/
grep ">" sequences | cut -d'~' -f4,7 | tr '~' '\t' > ../plasmidfinder_db.tsv
```

In_house_db 
```
cd /path/to/abricate/db/in_house
grep ">" sequences > ../in_house_db.tsv
```

## Input

**Usage: bash /path/to/CRACKLE2_WGS_SingleIsolate.sh ConfigFile Sample_Name Number_Threads Fastq_Read1 Fastq_Read2 OutputLocation**

**Config file - CRACKLEII_configfile.txt**

This file contains site specific details like the site and library prep method used, as well as details on the pipeline, paths of programs used, and database and list paths.

The full script contains the following information: 
```
#CONFIG FILE - CRACKLEII Pipeline v1.01

#Site specific details
SITE="Site Name"
LIBRARY_PREP="PREP KIT"

#Pipeline details
VERSION=1.1
DATE=`date +%Y-%m-%d`

#Program and script paths
RASPBERRY_PATH=/home/bhanson/CRACKLE_II_Updated/programs/Raspberry-0.3/bin/raspberry
TRIMMOMATIC_PATH=/home/bhanson/CRACKLE_II_Updated/programs/Trimmomatic-0.36/trimmomatic-0.36.jar
STRAINSEEKER_PATH=/home/bhanson/CRACKLE_II_Updated/programs/strainseeker/seeker.pl
SPADES_PATH=/home/bhanson/CRACKLE_II_Updated/programs/SPAdes-3.11.1-Linux/bin/spades.py
REMOVESMALLS_PATH=/home/bhanson/CRACKLE_II_Updated/removesmalls.pl
QUAST_PATH=/home/bhanson/CRACKLE_II_Updated/programs/quast-4.6.2/quast.py
MLST_PATH=/home/bhanson/CRACKLE_II_Updated/programs/mlst-2.10/bin/mlst
BLAST_PATH=/home/bhanson/databases/ncbi-blast-2.7.1+/bin/
KLEBORATE_PATH=/home/bhanson/CRACKLE_II_Updated/programs/Kleborate/kleborate-runner.py
ABRICATE_PATH=/home/bhanson/CRACKLE_II_Updated/programs/abricate-0.7/bin/abricate
EMBOSS_PATH=/home/bhanson/CRACKLE_II_Updated/programs/EMBOSS-6.6.0/emboss/
BOWTIE2_PATH=/home/bhanson/ProgramTest/bowtie2-2.3.4.1-linux-x86_64/
SAMTOOLS_PATH=/home/bhanson/ProgramTest/samtools-1.8/samtools
BCFTOOLS_PATH=/home/bhanson/ProgramTest/bcftools-1.8
SEQTK_PATH=/home/bhanson/ProgramTest/seqtk-1.2/seqtk

XML_GENERATOR_PATH=/home/bhanson/CRACKLE_II_Updated/XML_generator.R
XMLTOGENEHITLIST=/home/bhanson/CRACKLE_II_Updated/XML_to_gene_hit_list.R
XML_SHORTREAD_PATH=/home/bhanson/CRACKLE_II_Updated/XML_ShortRead.R

#Database paths
ILLUMINAADAPTERS=/home/bhanson/CRACKLE_II_Updated/programs/Trimmomatic-0.36/adapters
STRAINSEEKERDB=/home/bhanson/CRACKLE_II_Updated/programs/strainseeker/ss_db_w32_4324

#Resistance database list document paths
CARDLIST=/home/bhanson/CRACKLE_II_Updated/resistanceLists/card_db.tsv
PLASMIDLIST=/home/bhanson/CRACKLE_II_Updated/resistanceLists/plasmidfinder_db.tsv
INHOUSELIST=/home/bhanson/CRACKLE_II_Updated/resistanceLists/in_house_db.tsv

###Program specific parameters
###The following section of the config file set specific parameters used for 
###different programs within this pipeline. Each program will be identified
###and specific details 

#Trimmomatic - please specify the library prep kit used (assumes paired end)
#Nextera Paired End:
ADAPTER=NexteraPE-PE.fa
#TruSeq Paired End:
#ADAPTER=TruSeq3-PE-2.fa

#Trimmomatic - other parameters. Each set to default from trimmomatic
#Parameter for illumina read clipping
CLIPPING=2:30:10
#Remove leading low quality bases below defined quality score
LEADING=3
#Remove trailing low quality bases below defined quality score
TRAILING=3
#Define the sliding window size and average quality score for trimming
WINDOW=`4:15`
#Set minimum length for reads following trimming
MIN=36

#Removal of small contigs following assembly (default is 500 bases)
MINCONTIG=500

#Abricate - setting minimum ID % and minimum coverage % for a match
#Minimum ID (pipeline default is 95% ID, use 95)
MINID=95
#Minimum coverage (pipeline default is to include everything, use 0)
MINCOV=0
```

**Number_Threads - This is the number of threads you want the pipeline to use on your computer**

**Fastq_Read1 - The R1 fastq file**

**Fastq_Read2 - The R2 fastq file**

**OutputLocation - The folder you want the output put into**

Example: 

```
for F in *R1*; do bash /path/to/CRACKLE_II/CRACKLE_II_WGS_SingleIsolate.sh /path/to/CRACKLE_II/CRACKLEII_configfile.txt 12 $F ${F%R1_001.fastq.gz}R2_001.fastq.gz ${F%_S*_R1*.fastq.gz}; done
```

## Output

Within the OutputLocation folder, you will find the following files, each with a prefix of the specified sample name:

* abricate_card.tab
	* Output from abricate for CARD database
* abricate_inhouse.tab
	* Output from abricate for in_house database
* abricate_plasmid.tab
	* Output from abricate for PlasmidFinder database
* SAMPLENAME_CRACKLE_V1.0.xml
	* Final XML formatted output for pipeline
* SAMPLENAME_R1.paired.fastq.gz
	* Trimmed fastq file for R1 used in SPAdes assembly
* SAMPLENAME_R2.paired.fastq.gz
	* Trimmed fastq file for R2 used in SPAdes assembly
* contigs.fasta.gz
	* Assembled contigs from SPAdes that are greater than 500 bases
* crackle.err
	* Error log for pipeline
* crackle.log
	* Standard output log for pipeline
* kleborate.tab
	* Output from Kleborate
* mlst.tab
	* Output from MLST program
* raspberry_raw_R1.txt
	* Raspberry report for raw R1 
* raspberry_raw_R2.txt
	* Raspberry report for raw R2
* raspberry_trim_R1.txt
	* Raspberry report for trimmed R1
* raspberry_trim_R2.txt
	* Raspberry report for trimmed R2
* report.tsv
	* Quast report for SPAdes assembly
* strainseeker_R1.csv
	* StrainSeeker output for trimmed R1
* version.csv
	* Versions of all software used in pipeline


## Additional Scripts

**RegenerateXML.sh script**

This script works on the output folder from the CRACKLE_II_WGS_SingleIsolate.sh script, reading in all of the step specific files and recreating the XML output. This must be run from the folder containing your pipeline output folder. 

Note: Make sure you move or remove the XML file in the output folder you are replacing.

Syntax:

```
bash /path/to/CRACKLE_II/CRACKLE_II_RegenerateXML.sh /path/to/CRACKLE_II/CRACKLEII_configfile.txt PIPELINE_OUTPUT_FOLDER/
```

Example:
```
bash /path/to/CRACKLE_II/CRACKLE_II_RegenerateXML.sh /path/to/CRACKLE_II/CRACKLEII_configfile.txt C3828/
```

**CRACKLE_II_Runstrainseeker_RegenerateXML.sh script**

This script works on the output folder from the CRACKLE_II_WGS_SingleIsolate.sh script, reading in all of the step specific files, rerunning the strainseeker file on the assembly file, and recreating the XML output. This must be run from the folder containing your pipeline output folder. 

Note: Make sure you move or remove the XML file in the output folder you are replacing.

Syntax:

```
bash /path/to/CRACKLE_II/CRACKLE_II_Runstrainseeker_RegenerateXML.sh /path/to/CRACKLE_II/CRACKLEII_configfile.txt PIPELINE_OUTPUT_FOLDER/
```

Example:
```
bash /path/to/CRACKLE_II/CRACKLE_II_Runstrainseeker_RegenerateXML.sh /path/to/CRACKLE_II/CRACKLEII_configfile.txt C3828/
```

**XML_to_multiple_csv.R**

This script takes in a list of XML files and writes out a set of CSV files (one for each step in the pipeline) where all of the samples are represented as an individual row. 

Examples:

If all XML files are in same directory:

```
cd /directory/of/XML/files
ls -1 *xml > ListFiles.txt
Rscript /path/to/CRACKLE_II/XML_to_multiple_csv.R ListFiles.txt /directory/of/XML/files 95 90 /path/to/desired/final/output/location
```

If all XML files in nested output folders:
```
cd /directory/of/pipeline/output
ls -1 */*xml > ListFiles.txt
Rscript /path/to/CRACKLE_II/XML_to_multiple_csv.R ListFiles.txt /directory/of/pipeline/output 95 90 /path/to/desired/final/output/location
```

**ShortReadGeneAlignment.sh**

This script takes the genes identified by abricate in the main pipeline and re-aligns the trimmed short reads to the gene sequences in the abricate references to confirm if a gene is present or interrupted, as well as where SNPs and Indels are in the sequence. 

```
bash /path/to/CRACKLE_II/ShortReadGeneAlignment.sh /path/to/CRACKLE_II/CRACKLEII_configfile.txt /path/to/PIPELINE_OUTPUT_FOLDER /path/to/PIPELINE_OUTPUT_FOLDER/XML_FILE /path/to/desired/final/output/location

```

