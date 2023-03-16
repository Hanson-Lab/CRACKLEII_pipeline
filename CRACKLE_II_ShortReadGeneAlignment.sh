#!/bin/bash

###############
USAGE="
Usage: bash $0 ConfigFile Sample_Name Folder XMLFile OutputLocation"

if [[ $# -ne 5 ]]; then
    echo "$USAGE"
else
    CONFIGFILE=$1
    SAMPLE_NAME=$2
    FOLDER=$3
    XMLFILE=$4
    OUTPUTLOCATION=$5
fi
###############

###Source config file
echo "Reading config file"
source $CONFIGFILE

###Set variables
LOG=${OUTPUTLOCATION}/crackle_shortread.log
ERR=${OUTPUTLOCATION}/crackle_shortread.err
VERSIONDOC=${OUTPUTLOCATION}/short_read_version.csv
TRIMMEDFASTQ_R1=$(basename "$FOLDER")_R1.paired.fastq.gz
TRIMMEDFASTQ_R2=$(basename "$FOLDER")_R2.paired.fastq.gz
CARDDB=${ABRICATE_PATH%/bin/abricate}/db/card/sequences
PLASMIDFINDERDB=${ABRICATE_PATH%/bin/abricate}/db/plasmidfinder/sequences
INHOUSEDB=${ABRICATE_PATH%/bin/abricate}/db/in_house_db/sequences
CARDTEMP=temp/card_sequences.fasta
PLASMIDFINDERTEMP=temp/plasmidfinder_sequences.fasta
INHOUSETEMP=temp/in_house_sequences.fasta
CARDHIT=temp/Card_ResistanceGeneList.csv
PLASMIDHIT=temp/Plasmid_ResistanceGeneList.csv
INHOUSEHIT=temp/InHouse_ResistanceGeneList.csv
CARDHITTEMP=temp/Card_genes.tmp
PLASMIDHITTEMP=temp/Plasmid_genes.tmp
INHOUSEHITTEMP=temp/InHouse_genes.tmp
CARDHITFASTA=${OUTPUTLOCATION}/Card_genes.fasta
PLASMIDHITFASTA=${OUTPUTLOCATION}/Plasmid_genes.fasta
INHOUSEHITFASTA=${OUTPUTLOCATION}/InHouse_genes.fasta
CARDBOWTIE=temp/CARD_bt2 
PLASMIDBOWTIE=temp/Plasmid_bt2 
INHOUSEBOWTIE=temp/InHouse_bt2 
CARDSAM=temp/CARD_alignments.sam 
PLASMIDSAM=temp/Plasmid_alignments.sam 
INHOUSESAM=tem/InHouse_alignments.sam 
CARDBAM=temp/CARD_alignments.bam
PLASMIDBAM=temp/Plasmid_alignments.bam
INHOUSEBAM=temp/InHouse_alignments.bam
CARDSORT=temp/CARD_alignments.bai
PLASMIDSORT=temp/Plasmid_alignments.bai
INHOUSESORT=temp/InHouse_alignments.bai
CARDDEPTH=temp/CARD_depth.tsv
PLASMIDDEPTH=temp/Plasmid_depth.tsv
INHOUSEDEPTH=temp/InHouse_depth.tsv
CARDCONSENSUSFQ=temp/CARD_consensus.fastq
PLASMIDCONSENSUSFQ=temp/Plasmid_consensus.fastq
INHOUSECONSENSUSFQ=temp/InHouse_consensus.fastq
CARDCONSENSUSFA=${OUTPUTLOCATION}/CARD_consensus.fasta
PLASMIDCONSENSUSFA=${OUTPUTLOCATION}/Plasmid_consensus.fasta
INHOUSECONSENSUSFA=${OUTPUTLOCATION}/InHouse_consensus.fasta
CARDSANITYBOWTIE=temp/CARD_Sanity_bt2
PLASMIDSANITYBOWTIE=temp/Plasmid_Sanity_bt2
INHOUSESANITYBOWTIE=temp/InHouse_Sanity_bt2
CARDSANITYSAM=temp/CARD_alignments_sanity.sam 
PLASMIDSANITYSAM=temp/Plasmid_alignments_sanity.sam 
INHOUSESANITYSAM=temp/InHouse_alignments_sanity.sam 
CARDSANITYBAM=temp/CARD_alignments_sanity.bam
PLASMIDSANITYBAM=temp/Plasmid_alignments_sanity.bam
INHOUSESANITYBAM=temp/InHouse_alignments_sanity.bam
CARDSANITYSORT=temp/CARD_alignments_sanity.bai
PLASMIDSANITYSORT=temp/Plasmid_alignments_sanity.bai
INHOUSESANITYSORT=temp/InHouse_alignments_sanity.bai
CARDSANITYCONSENSUSFQ=temp/CARD_consensus_sanity.fastq
PLASMIDSANITYCONSENSUSFQ=temp/Plasmid_consensus_sanity.fastq
INHOUSESANITYCONSENSUSFQ=temp/InHouse_consensus_sanity.fastq
CARDSANITYCONSENSUSFA=${OUTPUTLOCATION}/CARD_consensus_sanity.fasta
PLASMIDSANITYCONSENSUSFA=${OUTPUTLOCATION}/Plasmid_consensus_sanity.fasta
INHOUSESANITYCONSENSUSFA=${OUTPUTLOCATION}/InHouse_consensus_sanity.fasta
OUTPUTXMLTEMP=${OUTPUTLOCATION}/$(SAMPLE_NAME)_${DATE}_ShortReadAligned_ResistanceGenes_v1.0_temp.xml
OUTPUTXML=${OUTPUTLOCATION}/$(SAMPLE_NAME)_${DATE}_ShortReadAligned_ResistanceGenes_v1.0.xml

###Program paths
BOWTIE2_BUILD_PATH=${BOWTIE2_PATH}/bowtie2-build
BOWTIE2_COMMAND_PATH=${BOWTIE2_PATH}/bowtie2
BCFTOOLS=${BCFTOOLS_PATH}/bcftools
VCFTOOLS=${BCFTOOLS_PATH}/misc/vcfutils.pl

###cd to pipeline output directory
cd $FOLDER

###Create directories
echo "Creating directories"
mkdir -p $OUTPUTLOCATION
mkdir -p temp

###Error logging
exec > >(tee -a "${LOG}")
exec 2> >(tee -a "${ERR}")

###Echo config file into standard out
echo "Listing subset of config variables"
echo $DATE
echo "Printing config file"
cat $CONFIGFILE

###Get program version numbers
#Note: StrainSeeker, Spades, and Kleborate output version to stderr, which is rerouted to stdout
echo "Getting program versions"
touch $VERSIONDOC
echo "Bowtie2,"$($BOWTIE2_COMMAND_PATH --version | head -1 | cut -d' ' -f3) >> $VERSIONDOC
echo "Samtools,"$($SAMTOOLS_PATH --version | head -1 | cut -d' ' -f2) >> $VERSIONDOC
echo "Bcftools,"$($BCFTOOLS --version | head -1 | cut -d' ' -f2) >> $VERSIONDOC
echo "Seqtk,"$($SEQTK_PATH 2>&1 >/dev/null | grep "Version" | cut -d' ' -f2) >> $VERSIONDOC

###Create singleline files for databases
echo "Converting database multi-line fasta files to single-line fasta files"
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $CARDDB > $CARDTEMP
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $PLASMIDFINDERDB > $PLASMIDFINDERTEMP
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $INHOUSEDB > $INHOUSETEMP

###For each resistome section of XML, pull out list of hits 
echo "Getting resistance gene hits from XML files"
Rscript $XMLTOGENEHITLIST $XMLFILE CARD $CARDHIT
Rscript $XMLTOGENEHITLIST $XMLFILE PLASMIDFINDER $PLASMIDHIT
Rscript $XMLTOGENEHITLIST $XMLFILE INHOUSE $INHOUSEHIT

###If loop to check if any genes were identified in XML file
echo "Working on CARD reference genes"
#CARD
if [[ $(wc -l <$CARDHIT) -ge 2 ]]
then
	#Create files with one gene name on each line
	echo "Creating fasta file of identified CARD reference genes"
	cut -d "," -f1 $CARDHIT | tail -n +2 | tr -d '"' > $CARDHITTEMP

	#For each hit, pull gene sequence from DB
	while read LINE; do grep -A1 ~~~${LINE}~~~ $CARDTEMP; done < $CARDHITTEMP > $CARDHITFASTA

	#Create index of genes for all hits for bowtie2 
	echo "Creating bowtie2 reference for CARD identified genes"
	$BOWTIE2_BUILD_PATH -f $CARDHITFASTA $CARDBOWTIE

	#Use bowtie2 to align paired end reads to indexed reference
	echo "Aligning trimmed, paired-end reads to CARD reference"
	$BOWTIE2_COMMAND_PATH -x $CARDBOWTIE --very-sensitive-local -1 $TRIMMEDFASTQ_R1 -2 $TRIMMEDFASTQ_R2 -S $CARDSAM

	#Convert sam files generated by bowtie2 to bam files and sort
	echo "Converting sam file to bam file"
	$SAMTOOLS_PATH view -bS $CARDSAM | $SAMTOOLS_PATH sort -o $CARDBAM -

	#Index bam files
	echo "Indexing bam file"
	$SAMTOOLS_PATH index $CARDBAM $CARDSORT

	#Get coverage at each position of genes 
	echo "Getting coverage at each base position in identified CARD genes"
	$SAMTOOLS_PATH depth $CARDBAM > $CARDDEPTH

	#Use samtools mpileup to generate consensus and sequencing details for each gene 
	echo "Getting consensus fastq of alignment"
	$SAMTOOLS_PATH mpileup -uf $CARDHITFASTA $CARDBAM | $BCFTOOLS call -c --ploidy 1 - | $VCFTOOLS vcf2fq > $CARDCONSENSUSFQ

	#Convert fastq to fasta
	echo "Converting consensus fastq to fasta"
	$SEQTK_PATH seq -a $CARDCONSENSUSFQ > $CARDCONSENSUSFA

	#Sanity check section where short reads are realigned to consensus fasta
	echo "Realigning short reads to consensus sequences as sanity check - this may take a minute, it repeats all previous steps for CARD alignment"
	$BOWTIE2_BUILD_PATH -f $CARDCONSENSUSFA $CARDSANITYBOWTIE
	$BOWTIE2_COMMAND_PATH -x $CARDSANITYBOWTIE --very-sensitive-local -1 $TRIMMEDFASTQ_R1 -2 $TRIMMEDFASTQ_R2 -S $CARDSANITYSAM
	$SAMTOOLS_PATH view -bS $CARDSANITYSAM | $SAMTOOLS_PATH sort -o $CARDSANITYBAM -
	$SAMTOOLS_PATH index $CARDSANITYBAM $CARDSANITYSORT
	$SAMTOOLS_PATH mpileup -uf $CARDCONSENSUSFA $CARDSANITYBAM | $BCFTOOLS call -c --ploidy 1 - | $VCFTOOLS vcf2fq > $CARDSANITYCONSENSUSFQ
	$SEQTK_PATH seq -a $CARDSANITYCONSENSUSFQ > $CARDSANITYCONSENSUSFA
else
	echo "No CARD hits identified"
fi

#PlasmidFinder
if [[ $(wc -l <$PLASMIDHIT) -ge 2 ]]
then
	#Create files with one gene name on each line
	echo "Creating fasta file of identified PlasmidFinder reference genes"
	cut -d "," -f1 $PLASMIDHIT | tail -n +2 | tr -d '"' > $PLASMIDHITTEMP

	#For each hit, pull gene sequence from DB
	while read LINE; do grep -A1 ~~~${LINE}~~~ $PLASMIDFINDERTEMP; done < $PLASMIDHITTEMP > $PLASMIDHITFASTA

	#Create index of genes for all hits for bowtie2 
	echo "Creating bowtie2 reference for PlasmidFinder identified genes"
	$BOWTIE2_BUILD_PATH -f $PLASMIDHITFASTA $PLASMIDBOWTIE

	#Use bowtie2 to align paired end reads to indexed reference
	echo "Aligning trimmed, paired-end reads to PlasmidFinder reference"
	$BOWTIE2_COMMAND_PATH -x $PLASMIDBOWTIE --very-sensitive-local -1 $TRIMMEDFASTQ_R1 -2 $TRIMMEDFASTQ_R2 -S $PLASMIDSAM

	#Convert sam files generated by bowtie2 to bam files and sort
	echo "Converting sam file to bam file"
	$SAMTOOLS_PATH view -bS $PLASMIDSAM | $SAMTOOLS_PATH sort -o $PLASMIDBAM -

	#Index bam files
	echo "Indexing bam file"
	$SAMTOOLS_PATH index $PLASMIDBAM $PLASMIDSORT

	#Get coverage at each position of genes 
	echo "Getting coverage at each base position in identified PlasmidFinder genes"
	$SAMTOOLS_PATH depth $PLASMIDBAM > $PLASMIDDEPTH

	#Use samtools mpileup to generate consensus and sequencing details for each gene
        echo "Getting consensus fastq of alignment"	
	$SAMTOOLS_PATH mpileup -uf $PLASMIDHITFASTA $PLASMIDBAM | $BCFTOOLS call -c --ploidy 1 - | $VCFTOOLS vcf2fq > $PLASMIDCONSENSUSFQ

	#Convert fastq to fasta
	echo "Converting consensus fastq to fasta"
	$SEQTK_PATH seq -a $PLASMIDCONSENSUSFQ > $PLASMIDCONSENSUSFA

	#Sanity check section where short reads are realigned to consensus fasta
	echo "Realigning sholrt reads to consensus sequences as sanity check - this may take a minute, it repeats all previous steps for PlasmidFinder alignment"
	$BOWTIE2_BUILD_PATH -f $PLASMIDCONSENSUSFA $PLASMIDSANITYBOWTIE
	$BOWTIE2_COMMAND_PATH -x $PLASMIDSANITYBOWTIE --very-sensitive-local -1 $TRIMMEDFASTQ_R1 -2 $TRIMMEDFASTQ_R2 -S $PLASMIDSANITYSAM
	$SAMTOOLS_PATH view -bS $PLASMIDSANITYSAM | $SAMTOOLS_PATH sort -o $PLASMIDSANITYBAM -
	$SAMTOOLS_PATH index $PLASMIDSANITYBAM $PLASMIDSANITYSORT
	$SAMTOOLS_PATH mpileup -uf $PLASMIDCONSENSUSFA $PLASMIDSANITYBAM | $BCFTOOLS call -c --ploidy 1 - | $VCFTOOLS vcf2fq > $PLASMIDSANITYCONSENSUSFQ
	$SEQTK_PATH seq -a $PLASMIDSANITYCONSENSUSFQ > $PLASMIDSANITYCONSENSUSFA
else
	echo "No PlasmidFinder hits identified"
fi

#InHouseDB
if [[ $(wc -l <$INHOUSEHIT) -ge 2 ]]
then
	#Create files with one gene name on each line
	echo "Creating fasta file of identified InHouse reference genes"
	cut -d "," -f1 $INHOUSEHIT | tail -n +2 | tr -d '"' | sed 's/^/>/' > $INHOUSEHITTEMP

	#For each hit, pull gene sequence from DB
	grep -w -A1 --no-group-separator -f $INHOUSEHITTEMP $INHOUSETEMP > $INHOUSEHITFASTA

	#Create index of genes for all hits for bowtie2 
	echo "Creating bowtie2 reference for InHouse identified genes"
	$BOWTIE2_BUILD_PATH -f $INHOUSEHITFASTA $INHOUSEBOWTIE

	#Use bowtie2 to align paired end reads to indexed reference
	echo "Aligning trimmed, paired-end reads to InHouse reference"
	$BOWTIE2_COMMAND_PATH -x $INHOUSEBOWTIE --very-sensitive-local -1 $TRIMMEDFASTQ_R1 -2 $TRIMMEDFASTQ_R2 -S $INHOUSESAM

	#Convert sam files generated by bowtie2 to bam files and sort
	echo "Converting sam file to bam file"
	$SAMTOOLS_PATH view -bS $INHOUSESAM | $SAMTOOLS_PATH sort -o $INHOUSEBAM -

	#Index bam files
	echo "Indexing bam file"
	$SAMTOOLS_PATH index $INHOUSEBAM $INHOUSESORT

	#Get coverage at each position of genes
        echo "Getting coverage at each base position in identified InHouse genes"	
	$SAMTOOLS_PATH depth $INHOUSEBAM > $INHOUSEDEPTH

	#Use samtools mpileup to generate consensus and sequencing details for each gene 
	echo "Getting consensus fastq of alignment"
	$SAMTOOLS_PATH mpileup -uf $INHOUSEHITFASTA $INHOUSEBAM | $BCFTOOLS call -c --ploidy 1 - | $VCFTOOLS vcf2fq > $INHOUSECONSENSUSFQ

	#Convert fastq to fasta
	echo "Converting consensus fastq to fasta"
	$SEQTK_PATH seq -a $INHOUSECONSENSUSFQ > $INHOUSECONSENSUSFA

	#Sanity check section where short reads are realigned to consensus fasta
	echo "Realigning short reads to consensus sequences as sanity check - this may take a minute, it repeats all previous steps for InHouse alignment"
	$BOWTIE2_BUILD_PATH -f $INHOUSECONSENSUSFA $INHOUSESANITYBOWTIE
	$BOWTIE2_COMMAND_PATH -x $INHOUSESANITYBOWTIE --very-sensitive-local -1 $TRIMMEDFASTQ_R1 -2 $TRIMMEDFASTQ_R2 -S $INHOUSESANITYSAM
	$SAMTOOLS_PATH view -bS $INHOUSESANITYSAM | $SAMTOOLS_PATH sort -o $INHOUSESANITYBAM -
	$SAMTOOLS_PATH index $INHOUSESANITYBAM $INHOUSESANITYSORT
	$SAMTOOLS_PATH mpileup -uf $INHOUSECONSENSUSFA $INHOUSESANITYBAM | $BCFTOOLS call -c --ploidy 1 - | $VCFTOOLS vcf2fq > $INHOUSESANITYCONSENSUSFQ
	$SEQTK_PATH seq -a $INHOUSESANITYCONSENSUSFQ > $INHOUSESANITYCONSENSUSFA
else
	echo "No InHouse hits identified"
fi

###Create new XML file with cigar strings for alignments to check
echo "Creating Alignment information and XML file"
echo "Rscript $XML_SHORTREAD_PATH $SITE $LIBRARY_PREP $DATE $VERSION $VERSIONDOC $CARDHITFASTA $CARDCONSENSUSFA $CARDSANITYCONSENSUSFA $PLASMIDHITFASTA $PLASMIDCONSENSUSFA $PLASMIDSANITYCONSENSUSFA $INHOUSEHITFASTA $INHOUSECONSENSUSFA $INHOUSESANITYCONSENSUSFA $CARDLIST $PLASMIDLIST $INHOUSELIST $OUTPUTXMLTEMP"
Rscript $XML_SHORTREAD_PATH $SITE $LIBRARY_PREP $DATE $VERSION $VERSIONDOC $CARDHITFASTA $CARDCONSENSUSFA $CARDSANITYCONSENSUSFA $PLASMIDHITFASTA $PLASMIDCONSENSUSFA $PLASMIDSANITYCONSENSUSFA $INHOUSEHITFASTA $INHOUSECONSENSUSFA $INHOUSESANITYCONSENSUSFA $CARDLIST $PLASMIDLIST $INHOUSELIST $OUTPUTXMLTEMP

###Takes flat xml file and restores tree structure
xmllint --format $OUTPUTXMLTEMP > $OUTPUTXML

###Clean up
echo "Cleaning up"
rm -r temp/
rm ${OUTPUTLOCATION}/*sam
gzip ${OUTPUTLOCATION}/*fasta

rm $OUTPUTXMLTEMP
