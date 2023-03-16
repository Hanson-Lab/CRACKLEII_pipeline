#!/bin/bash

usage() {
  echo "Usage: $0 -c <config_file> -s <sample_name> -t <threads> [-a <assembly_file> | -1 <read1_file> -2 <read2_file>] -o <output_dir>"
  exit 1
}

CONFIGFILE=""
SAMPLE_NAME=""
THREADS=""
ASSEMBLY_FILE=""
FQ1=""
FQ2=""
OUTPUTLOCATION=""

options=$(getopt -o c:s:t:a:1:2:o: --long config_file:,sample_name:,threads:,assembly_file:,read1_file:,read2_file:,output_dir: -- "$@")

if [ $? -ne 0 ]; then
  usage
fi

eval set -- "$options"

while true; do
  case "$1" in
    -c|--config_file)
      CONFIGFILE="$2"
      shift 2 ;;
    -s|--sample_name)
      SAMPLE_NAME="$2"
      shift 2 ;;
    -t|--threads)
      THREADS="$2"
      shift 2 ;;
    -a|--assembly_file)
      ASSEMBLY_FILE="$2"
      shift 2 ;;
    -1|--read1_file)
      FQ1="$2"
      shift 2 ;;
    -2|--read2_file)
      FQ2="$2"
      shift 2 ;;
    -o|--output_dir)
      OUTPUTLOCATION="$2"
      shift 2 ;;
    --)
      shift
      break ;;
    *)
      usage ;;
  esac
done

# Check if all the required options are set
if [[ -z "$CONFIGFILE" || -z "$SAMPLE_NAME" || -z "$THREADS" || -z "$OUTPUTLOCATION" ]]; then
  usage
fi

#Check if both -1 and -2 are used if one of the options is used
if [ -n "$FQ1" ] && [ -z "$FQ2" ]; then
    echo "Error: Both -1 and -2 must be provided."
    exit 1
elif [ -n "$FQ2" ] && [ -z "$FQ1" ]; then
    echo "Error: Both -1 and -2 must be provided."
    exit 1
fi

#Check if the user submitted both the -a and -1 or -2 variables and if so generate error since both assembly and short reads provided
if [[ -n "$ASSEMBLY_FILE" && ( -n "$FQ1" || -n "$FQ2" ) ]]; then
  echo "Error: Please provide both -1 and -2 if you wish to generate an assembly with -a"
  usage
fi

#Check if the user did not submit either -1-2 or -a variables and if so provide error message
if [[ -z "$ASSEMBLY_FILE" && ( -z "$FQ1" || -z "$FQ2" ) ]]; then
  echo "Error: Please provide both -1 and -2 if you wish to generate an assembly with -a or provide only -a if you wish to use a pre-generated assembly."
  usage
fi

#Echo user submitted variables
echo "The selected configuration file is: " $CONFIGFILE
echo "The user-defined sample name is: " $SAMPLE_NAME
echo "The number of threads requested is: " $THREADS
echo "The selected fastq read1 is: " $FQ1
echo "The selected fastq read2 is: " $FQ2
echo "The user-defined output location is: " $OUTPUTLOCATION
echo "The selected assembly fasta is: " $ASSEMBLY_FILE

###Visual spacing for output
echo ""
echo ""

###Source config file
echo "Reading config file"
source $CONFIGFILE

###Set variables
RAW_FQ1=$FQ1
RAW_FQ2=$FQ2
USER_ASSEM=$FROMASSEMBLY
LOG=${OUTPUTLOCATION}/${SAMPLE_NAME}_crackle.log
ERR=${OUTPUTLOCATION}/${SAMPLE_NAME}_crackle.err
VERSIONDOC=${OUTPUTLOCATION}/${SAMPLE_NAME}_pipeline_version.csv
RASPBERRY_RAW_R1=${OUTPUTLOCATION}/${SAMPLE_NAME}_raspberry_raw_R1.txt
RASPBERRY_RAW_R2=${OUTPUTLOCATION}/${SAMPLE_NAME}_raspberry_raw_R2.txt
TRIM_PAIRED_R1=${OUTPUTLOCATION}/${SAMPLE_NAME}_R1.paired.fastq
TRIM_PAIRED_R2=${OUTPUTLOCATION}/${SAMPLE_NAME}_R2.paired.fastq
TRIM_UNPAIRED_R1=${OUTPUTLOCATION}/${SAMPLE_NAME}_R1.unpaired.fastq
TRIM_UNPAIRED_R2=${OUTPUTLOCATION}/${SAMPLE_NAME}_R2.unpaired.fastq
RASPBERRY_TRIM_R1=${OUTPUTLOCATION}/${SAMPLE_NAME}_raspberry_trim_R1.txt
RASPBERRY_TRIM_R2=${OUTPUTLOCATION}/${SAMPLE_NAME}_raspberry_trim_R2.txt
STRAINSEEKER_OUT=${OUTPUTLOCATION}/${SAMPLE_NAME}_strainseeker.csv
SPADESDIR=${OUTPUTLOCATION}/spades
SPADES_ASSEM=${SPADESDIR}/contigs.fasta
SPADES_ASSEM_FILTER=${SPADESDIR}/contigs_filter.fasta
QUAST_OUT=${OUTPUTLOCATION}/report.tsv
QUAST_OUT_FINAL=${OUTPUTLOCATION}/${SAMPLE_NAME}_report.tsv
MLST_OUT=${OUTPUTLOCATION}/${SAMPLE_NAME}_mlst.tab
KLEBORATE_OUT=${OUTPUTLOCATION}/${SAMPLE_NAME}_kleborate.tab
ABRICATE_CARD=${OUTPUTLOCATION}/${SAMPLE_NAME}_abricate_card.tab
ABRICATE_PLASMID=${OUTPUTLOCATION}/${SAMPLE_NAME}_abricate_plasmid.tab
ABRICATE_INHOUSE=${OUTPUTLOCATION}/${SAMPLE_NAME}_abricate_inhouse.tab
OUTPUT_XML_TEMP=${OUTPUTLOCATION}/${SAMPLE_NAME}_${DATE}_CRACKLE_v${VERSION}_temp.xml
OUTPUT_XML=${OUTPUTLOCATION}/${SAMPLE_NAME}_${DATE}_CRACKLE_v${VERSION}.xml


###Create directories
echo "Creating directories"
mkdir -p $OUTPUTLOCATION
mkdir -p $SPADESDIR

#Error logging
exec > >(tee -a "${LOG}")
exec 2> >(tee -a "${ERR}")

###Echo config file into standard out
echo "Listing subset of config variables"
echo $SITE
echo $LIBRARY_PREP
echo $VERSION
echo $DATE
echo "Printing config file"
cat $CONFIGFILE

###Add blast and emboss executables and to path
PATH=$MASH_PATH:$ANY2FASTA_PATH:$BLAST_PATH:$EMBOSS_PATH:$PATH

###Get program version numbers
#Note: StrainSeeker, Spades, and Kleborate output version to stderr, which is rerouted to stdout
echo "Getting program versions"
touch $VERSIONDOC
echo "Raspberry,"$($RASPBERRY_PATH -v | head -1) >> $VERSIONDOC
echo "Trimmomatic,"$(java -Xmx72g -jar $TRIMMOMATIC_PATH -version) >> $VERSIONDOC
echo "StrainSeeker,"$(perl $STRAINSEEKER_PATH -v 2>&1 | cut -d' ' -f2) >> $VERSIONDOC
echo "Spades,"$($SPADES_PATH -v 2>&1 | cut -d' ' -f4) >> $VERSIONDOC
echo "Quast,"$($QUAST_PATH -v | cut -d' ' -f2 | cut -d',' -f1) >> $VERSIONDOC
echo "MLST,"$($MLST_PATH -v | cut -d' ' -f2) >> $VERSIONDOC
echo "Kleborate,"$($KLEBORATE_PATH --version 2>&1 | cut -d' ' -f2) >> $VERSIONDOC
echo "Abricate,"$($ABRICATE_PATH -v | cut -d' ' -f2) >> $VERSIONDOC

###Check if the from_assembly flag was used
#Marks that no flag was used
if [ -z "$ASSEMBLY_FILE" ]
then
        echo "Read trimming and assembly will be conducted"

###Raspberry for QC of raw untrimmed reads
echo "Running Raspberry on raw untrimmed reads"
$RASPBERRY_PATH $RAW_FQ1 > $RASPBERRY_RAW_R1
$RASPBERRY_PATH $RAW_FQ2 > $RASPBERRY_RAW_R2

###Trimmomatic filtering
echo "Trimming raw reads"
java -Xmx72g -jar $TRIMMOMATIC_PATH PE -phred33 -threads $THREADS $RAW_FQ1 $RAW_FQ2 $TRIM_PAIRED_R1 $TRIM_UNPAIRED_R1 $TRIM_PAIRED_R2 $TRIM_UNPAIRED_R2 ILLUMINACLIP:${ILLUMINAADAPTERS}/${ADAPTER}:${CLIPPING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${WINDOW} MINLEN:${MIN}

###Raspberry for QC of trimmed reads
echo "Running Raspberry on trimmed reads"
$RASPBERRY_PATH $TRIM_PAIRED_R1 > $RASPBERRY_TRIM_R1
$RASPBERRY_PATH $TRIM_PAIRED_R2 > $RASPBERRY_TRIM_R2

###Spades assembly
echo "Assembling with SPAdes"
$SPADES_PATH -t $THREADS --isolate --only-assembler --cov-cutoff auto --pe1-1 $TRIM_PAIRED_R1 --pe1-2 $TRIM_PAIRED_R2 -o $SPADESDIR

#set else statement for situation when from_assembly flag was used
else
        echo "User indicated assembly is provided, skipping trimming and assembly"
        RASPBERRY_RAW_R1=$RASPBERRY_DUMMY
        RASPBERRY_RAW_R2=$RASPBERRY_DUMMY
        RASPBERRY_TRIM_R1=$RASPBERRY_DUMMY
        RASPBERRY_TRIM_R2=$RASPBERRY_DUMMY
        SPADES_ASSEM=$ASSEMBLY_FILE
        SPADES_ASSEM_FILTER=${SAMPLE_NAME}_filtered.contigs.fasta

fi
###Removecontigs smaller than 500 bases
echo "Removing small contigs"
python3 $REMOVESMALLS_PATH -i $SPADES_ASSEM -o $SPADES_ASSEM_FILTER -n $MINCONTIG

###Quast assembly check
echo "Checking assembly quality"
$QUAST_PATH -t $THREADS -o $OUTPUTLOCATION $SPADES_ASSEM_FILTER

###Check contamination and species for the assembly
echo "Checking for contamination"
perl $STRAINSEEKER_PATH -d $STRAINSEEKERDB -i $SPADES_ASSEM_FILTER -o $STRAINSEEKER_OUT

###MLST
echo "Performing in silico MLST"
$MLST_PATH $SPADES_ASSEM_FILTER > $MLST_OUT

###Kleborate
echo "Running Kleborate"
$KLEBORATE_PATH -o $KLEBORATE_OUT -r -a $SPADES_ASSEM_FILTER

###Abricate resistance identification
echo "Running Abricate - CARD database"
$ABRICATE_PATH --minid $MINID --mincov $MINCOV --nopath --db card $SPADES_ASSEM_FILTER > $ABRICATE_CARD

###Plasmid finder
echo "Running Abricate - PlasmidFinder database"
$ABRICATE_PATH --minid $MINID --mincov $MINCOV --nopath --db plasmidfinder $SPADES_ASSEM_FILTER > $ABRICATE_PLASMID

###Abricate with in-house database
echo "Running Abricate - In House database"
$ABRICATE_PATH --minid $MINID --mincov $MINCOV --nopath --db in_house_db $SPADES_ASSEM_FILTER > $ABRICATE_INHOUSE

###R program to generate xml files
echo "Generating XML file"
echo "Command used: Rscript $XML_GENERATOR_PATH $SAMPLE_NAME $SITE $LIBRARY_PREP $DATE $VERSION $VERSIONDOC $RASPBERRY_RAW_R1 $RASPBERRY_RAW_R2 $RASPBERRY_TRIM_R1 $RASPBERRY_TRIM_R2 $QUAST_OUT $STRAINSEEKER_OUT $MLST_OUT $KLEBORATE_OUT $ABRICATE_CARD $ABRICATE_PLASMID $ABRICATE_INHOUSE $CARDLIST $PLASMIDLIST $INHOUSELIST $OUTPUT_XML_TEMP"
Rscript $XML_GENERATOR_PATH $SAMPLE_NAME $SITE $LIBRARY_PREP $DATE $VERSION $VERSIONDOC $RASPBERRY_RAW_R1 $RASPBERRY_RAW_R2 $RASPBERRY_TRIM_R1 $RASPBERRY_TRIM_R2 $QUAST_OUT $STRAINSEEKER_OUT $MLST_OUT $KLEBORATE_OUT $ABRICATE_CARD $ABRICATE_PLASMID $ABRICATE_INHOUSE $CARDLIST $PLASMIDLIST $INHOUSELIST $OUTPUT_XML_TEMP
#Takes flat xml file and restores tree structure
xmllint --format $OUTPUT_XML_TEMP > $OUTPUT_XML

###List all files in Output folder
ls -lR $OUTPUTLOCATION

###Clean up
echo "Cleaning up"
#Clean up raspberry intermediate files
rm ${RAW_FQ1%_S*_R1*.fastq.gz}*rlen
#Clean up output folder
rm -r ${OUTPUTLOCATION}/basic_stats
rm -r ${OUTPUTLOCATION}/icarus_viewers
rm $OUTPUT_XML_TEMP ${OUTPUTLOCATION}/*rlen ${OUTPUTLOCATION}/*unpaired* ${OUTPUTLOCATION}/icarus.html ${OUTPUTLOCATION}/quast.log ${OUTPUTLOCATION}/report.html ${OUTPUTLOCATION}/report.pdf ${OUTPUTLOCATION}/report.tex ${OUTPUTLOCATION}/report.txt ${OUTPUTLOCATION}/transposed_report.tex ${OUTPUTLOCATION}/transposed_report.txt ${OUTPUTLOCATION}/transposed_report.tsv
#Move spades assembly contigs.fasta to output location
mv $SPADES_ASSEM_FILTER ${OUTPUTLOCATION}/${SAMPLE_NAME}_contigs.fasta
#Move quast report to final output location name
mv $QUAST_OUT $QUAST_OUT_FINAL
#Remove spades folder
rm -r ${OUTPUTLOCATION}/spades
#Gzip fastq and fasta files in output
gzip ${OUTPUTLOCATION}/*fastq
gzip ${OUTPUTLOCATION}/*fasta
echo "All Finished - Goodbye"
