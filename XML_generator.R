##################
#USAGE - Rscript $XML_GENERATOR_PATH $SAMPLE_NAME $SITE $LIBRARY_PREP $DATE $VERSION $VERSIONDOC $RASPBERRY_RAW_R1 $RASPBERRY_RAW_R2 $RASPBERRY_TRIM_R1 $RASPBERRY_TRIM_R2 $QUAST_OUT $STRAINSEEKER_OUT $MLST_OUT $KLEBORATE_OUT $ABRICATE_CARD $ABRICATE_PLASMID $ABRICATE_INHOUSE $CARDLIST $PLASMIDLIST $INHOUSELIST $OUTPUT_XML_TEMP

#For use within the CRACKLEII_WGS_SingleIsolate.sh v1.0 pipeline 

#Author: Blake Hanson (blake.hanson [at] uth.tmc.edu)
#Version:1.1
#Date Last updated:2019/11/19
#Updates: Fixes bug where the name in the XML file is sometimes truncated
##################

###Arguments
args<-commandArgs(TRUE)

###Libraries
library(XML)

###Set variable names for files/data passed into script 
SAMPLE_NAME=args[1]
SITE=args[2]
LIBRARY_PREP=args[3]
DATE=args[4]
PIPELINEVER=args[5]
VERSIONDOC <- args[6]
RASP_RAW_R1 <- args[7]
RASP_RAW_R2 <- args[8]
RASP_TRIM_R1 <- args[9]
RASP_TRIM_R2 <- args[10]
QUASTDOC <- args[11]
STRAINDOC <- args[12]
MLSTDOC <- args[13]
KLEBORATEDOC <- args[14]
CARDDOC <- args[15]
PLASMIDDOC <- args[16]
INHOUSEDOC <- args[17]
CARDLIST <- args[18]
PLASMIDLIST <- args[19]
INHOUSELIST <- args[20]
OUTFILE=args[21]


###Import database specific mechanisms for the three  files (prepopulated)
#CARD
print("Reading in CARD gene list")
print(CARDLIST)
Card_gene_list <- read.table(CARDLIST, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
#Plasmidfinder
print("Reading in Plasmidfinder gene list")
print(PLASMIDLIST)
Plasmid_gene_list <- read.table(PLASMIDLIST, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
#In_house_db
print("Reading in Inhouse gene list")
print(INHOUSELIST)
Inhouse_gene_list <- read.table(INHOUSELIST, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")

###Import Program Versions
print("Reading in version file")
print(VERSIONDOC)
Version_doc <- read.csv(VERSIONDOC, header=FALSE, row.names = 1, col.names = c("","version"))

###Import Raspberry reports for raw data
#Raspberry raw read 1
print("Reading in raspberry raw R1 file")
print(RASP_RAW_R1)
Rasp_Raw_R1 <- read.table(RASP_RAW_R1, skip=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
Rasp_Raw_R1_formatted <- as.data.frame(do.call(rbind, strsplit(Rasp_Raw_R1$V1, ' : ')))
Rasp_Raw_R1_formatted$V1 <- gsub("%", "Percent", Rasp_Raw_R1_formatted$V1)
Rasp_Raw_R1_formatted$V1 <- gsub(" ", "_", Rasp_Raw_R1_formatted$V1)
Rasp_Raw_R1_formatted$V1 <- gsub(">=", "GTEQ_", Rasp_Raw_R1_formatted$V1)
rownames(Rasp_Raw_R1_formatted) <- Rasp_Raw_R1_formatted$V1
Rasp_Raw_R1_formatted$V1 <- NULL
names(Rasp_Raw_R1_formatted) <- "value"
#Raspberry raw read 2
print("Reading in raspberry raw R2 file")
print(RASP_RAW_R2)
Rasp_Raw_R2 <- read.table(RASP_RAW_R2, skip=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
Rasp_Raw_R2_formatted <- as.data.frame(do.call(rbind, strsplit(Rasp_Raw_R2$V1, ' : ')))
Rasp_Raw_R2_formatted$V1 <- gsub("%", "Percent", Rasp_Raw_R2_formatted$V1)
Rasp_Raw_R2_formatted$V1 <- gsub(" ", "_", Rasp_Raw_R2_formatted$V1)
Rasp_Raw_R2_formatted$V1 <- gsub(">=", "GTEQ_", Rasp_Raw_R2_formatted$V1)
rownames(Rasp_Raw_R2_formatted) <- Rasp_Raw_R2_formatted$V1
Rasp_Raw_R2_formatted$V1 <- NULL
names(Rasp_Raw_R2_formatted) <- "value"

###Import Raspberry reports for trimmed data
#Raspberry trimmed read 1
print("Reading in raspberry trimmed R1 file")
print(RASP_TRIM_R1)
Rasp_Trim_R1 <- read.table(RASP_TRIM_R1, skip=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
Rasp_Trim_R1_formatted <- as.data.frame(do.call(rbind, strsplit(Rasp_Trim_R1$V1, ' : ')))
Rasp_Trim_R1_formatted$V1 <- gsub("%", "Percent", Rasp_Trim_R1_formatted$V1)
Rasp_Trim_R1_formatted$V1 <- gsub(" ", "_", Rasp_Trim_R1_formatted$V1)
Rasp_Trim_R1_formatted$V1 <- gsub(">=", "GTEQ_", Rasp_Trim_R1_formatted$V1)
rownames(Rasp_Trim_R1_formatted) <- Rasp_Trim_R1_formatted$V1
Rasp_Trim_R1_formatted$V1 <- NULL
names(Rasp_Trim_R1_formatted) <- "value"
#Raspberry trimmed read 2
print("Reading in raspberry trimmed R2 file")
print(RASP_TRIM_R2)
Rasp_Trim_R2 <- read.table(RASP_TRIM_R2, skip=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
Rasp_Trim_R2_formatted <- as.data.frame(do.call(rbind, strsplit(Rasp_Trim_R2$V1, ' : ')))
Rasp_Trim_R2_formatted$V1 <- gsub("%", "Percent", Rasp_Trim_R2_formatted$V1)
Rasp_Trim_R2_formatted$V1 <- gsub(" ", "_", Rasp_Trim_R2_formatted$V1)
Rasp_Trim_R2_formatted$V1 <- gsub(">=", "GTEQ_", Rasp_Trim_R2_formatted$V1)
rownames(Rasp_Trim_R2_formatted) <- Rasp_Trim_R2_formatted$V1
Rasp_Trim_R2_formatted$V1 <- NULL
names(Rasp_Trim_R2_formatted) <- "value"

###Import Spades quality information
print("Reading in quast assembly quality file")
print(QUASTDOC)
Quast_doc <- read.table(QUASTDOC, sep="\t", header=TRUE, col.names = c("","value"))
Quast_doc$X <- gsub("%", "percent", Quast_doc$X)
Quast_doc$X <- gsub("\\(", "", Quast_doc$X)
Quast_doc$X <- gsub("\\)", "", Quast_doc$X)
Quast_doc$X <- gsub(" ", "_", Quast_doc$X)
Quast_doc$X <- gsub(">=", "GTEQ", Quast_doc$X)
rownames(Quast_doc) <- Quast_doc$X
Quast_doc$X <- NULL

###Import StrainSeeker Results
print("Reading in strainseeker file")
print(STRAINDOC)
err <- try(Strain_doc <- read.csv(STRAINDOC, header=FALSE, skip = 1))
if("try-error" %in% class(err))
{
  Strain_doc<-data.frame(V1="NA")
}

###Import MLST Results
print("Reading in MLST file")
print(MLSTDOC)
MLST_doc <- read.table(MLSTDOC, header = FALSE)
#Pull species out of MLST results (variable 2)
MLST_species <- MLST_doc$V2
#Pull ST type out of MLST results (variable 3)
MLST_ST <- MLST_doc$V3
#Create dataframe of only genes and gene ID numbers
drops <- c("V1","V2", "V3")
MLST_doc_sub <- MLST_doc[ , !(names(MLST_doc) %in% drops)]
tMLST_doc_sub <- as.data.frame(t(MLST_doc_sub))
tMLST_doc_sub$V2 <- gsub("\\(", " ", tMLST_doc_sub$V1)
tMLST_doc_sub$V3 <- gsub("\\)", "", tMLST_doc_sub$V2)
MLST_final_DF <- as.data.frame(do.call(rbind, strsplit(tMLST_doc_sub$V3, ' ')))
try(names(MLST_final_DF) <- c("gene", "allele_number"))

#When there is not an MLST scheme found, the dataframe is empty. This creates a null value DF to keep everything from breaking
if(nrow(MLST_final_DF) == 0){
  print("MLST does not contain data - Creating empty section for xml")
  gene <- c("NA")
  allele_number <- c("NA")
  MLST_final_DF <- data.frame(gene, allele_number)
}else{
  print("MLST contains data - Everything is good")
}

###Import Kleborate results 
#Note: Kleborate results are always in the same order and all genes are present in output
#so no need to merge with consistent list like done with Abricate output
print("Reading in Kleborate file")
print(KLEBORATEDOC)
Kleborate_doc <- read.table(KLEBORATEDOC, header=TRUE, quote = '', stringsAsFactors = FALSE, sep = "\t")

###Import card output from abricate for sample
print("Reading in CARD output file")
print(CARDDOC)
col_names <- c("file", "sequence", "start", "end", "strand", "gene", "coverage", "coverage_map",
              "gaps", "percent_coverage", "identity", "database", "accession", "product", "resistance")
Card_doc <- read.table(CARDDOC, sep="\t", header=FALSE, quote = "", stringsAsFactors = FALSE, col.names = col_names)
Card_doc_filter <- Card_doc[2:12]
#Merge card output from abricate with full gene list for card database
Card_final <- merge(Card_gene_list, Card_doc_filter, by.x="V1", by.y="gene", all = TRUE)
colnames(Card_final)[1] <- "gene"
colnames(Card_final)[2] <- "accession"
colnames(Card_final)[3] <- "resistance"
colnames(Card_final)[4] <- "product"

###Import plasmidfinder output from abricate for sample
print("Reading in plasmidfinder output file")
print(PLASMIDDOC)
col_names <- c("file", "sequence", "start", "end", "strand", "gene", "coverage", "coverage_map",
               "gaps", "percent_coverage", "identity", "database", "accession", "product", "resistance")
Plasmid_doc <- read.table(PLASMIDDOC, sep="\t", header=FALSE, quote = "", stringsAsFactors = FALSE, col.names = col_names)
Plasmid_doc_filter <- Plasmid_doc[2:12]
#Merge plasmidfinder output from abricate with full gene list for plasmidfinder database
Plasmid_final <- merge(Plasmid_gene_list, Plasmid_doc_filter, by.x="V1", by.y="gene", all = TRUE)
colnames(Plasmid_final)[1] <- "gene"
colnames(Plasmid_final)[2] <- "accession"
colnames(Plasmid_final)[2] <- "product"

###Import in_house_db output from abricate for sample
print("Reading in inhouse output file")
print(INHOUSEDOC)
col_names <- c("file", "sequence", "start", "end", "strand", "gene", "coverage", "coverage_map",
               "gaps", "percent_coverage", "identity", "database", "accession", "product", "resistance")
InHouse_doc <- read.table(INHOUSEDOC, sep="\t", header=FALSE, quote = "", stringsAsFactors = FALSE, col.names = col_names)
InHouse_doc_filter <- InHouse_doc[2:12]
#Merge inhouse db output from abricate with full gene list for inhouse database
InHouse_final <- merge(Inhouse_gene_list, InHouse_doc_filter, by.x="V1", by.y="gene", all = TRUE)
colnames(InHouse_final)[1] <- "gene"
colnames(InHouse_final)[2] <- "description"

###Create xml file
#Open empty xml tree
xml <- xmlTree()
#Create overarching node "data"
xml$addTag("data", close=FALSE)
#Create child node - lab
xml$addTag("sample_name", SAMPLE_NAME)
xml$addTag("lab", close=FALSE)
xml$addTag("site", SITE)
xml$addTag("library_prep", LIBRARY_PREP)
xml$addTag("date", DATE)
xml$closeTag("lab")
#create child node - pipeline
xml$addTag("pipeline", close=FALSE)
xml$addTag("version", PIPELINEVER)
xml$closeTag("pipeline")
#create child node - software_versions
xml$addTag("software_versions", close=FALSE)
for (i in 1:nrow(Version_doc)) {
  xml$addTag(rownames(Version_doc)[i], close=FALSE)
  for (j in names(Version_doc)) {
    xml$addTag(j, Version_doc[i, j])
  }
  xml$closeTag()
}
xml$closeTag("software_versions")
#create child node - raspberry_raw_read1 (results for R1 for untrimmed data)
xml$addTag("raspberry_raw_read1", close=FALSE)
for (i in 1:nrow(Rasp_Raw_R1_formatted)) {
  xml$addTag(rownames(Rasp_Raw_R1_formatted)[i], close=FALSE)
  for (j in names(Rasp_Raw_R1_formatted)) {
    xml$addTag(j, Rasp_Raw_R1_formatted[i, j])
  }
  xml$closeTag()
}
xml$closeTag("raspberry_raw_read1")
#create child node - raspberry_raw_read2 (results for R2 for untrimmed data)
xml$addTag("raspberry_raw_read2", close=FALSE)
for (i in 1:nrow(Rasp_Raw_R2_formatted)) {
  xml$addTag(rownames(Rasp_Raw_R2_formatted)[i], close=FALSE)
  for (j in names(Rasp_Raw_R2_formatted)) {
    xml$addTag(j, Rasp_Raw_R2_formatted[i, j])
  }
  xml$closeTag()
}
xml$closeTag("raspberry_raw_read2")
#create child node - raspberry_trim_read1 (results for R1 for trimmed data)
xml$addTag("raspberry_trim_read1", close=FALSE)
for (i in 1:nrow(Rasp_Trim_R1_formatted)) {
  xml$addTag(rownames(Rasp_Trim_R1_formatted)[i], close=FALSE)
  for (j in names(Rasp_Trim_R1_formatted)) {
    xml$addTag(j, Rasp_Trim_R1_formatted[i, j])
  }
  xml$closeTag()
}
xml$closeTag("raspberry_trim_read1")
#create child node - raspberry_trim_read2 (results for R2 for trimmed data)
xml$addTag("raspberry_trim_read2", close=FALSE)
for (i in 1:nrow(Rasp_Trim_R2_formatted)) {
  xml$addTag(rownames(Rasp_Trim_R2_formatted)[i], close=FALSE)
  for (j in names(Rasp_Trim_R2_formatted)) {
    xml$addTag(j, Rasp_Trim_R2_formatted[i, j])
  }
  xml$closeTag()
}
xml$closeTag("raspberry_trim_read2")
#create child node - assebly_stats (output from Quast assessment of Spades assembly)
xml$addTag("assembly_stats", close=FALSE)
for (i in 1:nrow(Quast_doc)) {
  xml$addTag(rownames(Quast_doc)[i], close=FALSE)
  for (j in names(Quast_doc)) {
    xml$addTag(j, Quast_doc[i, j])
  }
  xml$closeTag()
}
xml$closeTag("assembly_stats")
#create child node - strain_id (output from StrainSeeker)
xml$addTag("strain_id", close=FALSE)
for (i in 1:nrow(Strain_doc)) {
  xml$addTag("strain", close=FALSE)
  for (j in names(Strain_doc)) {
    xml$addTag(j, Strain_doc[i, j])
  }
  xml$closeTag()
}
xml$closeTag("strain_id")
#create child node - mlst (output from MLST)
xml$addTag("mlst", close = FALSE)
xml$addTag("species", MLST_species)
xml$addTag("ST_type", MLST_ST)
for (i in 1:nrow(MLST_final_DF)) {
  xml$addTag("loci", close=FALSE)
  for (j in names(MLST_final_DF)) {
    xml$addTag(j, MLST_final_DF[i, j])
  }
  xml$closeTag()
}
xml$closeTag("mlst")
#create child node - kleborate (output from kleborate)
xml$addTag("kleborate", close = FALSE)
for (i in 1:nrow(Kleborate_doc)) {
  xml$addTag("loci", close=FALSE)
  for (j in names(Kleborate_doc)) {
    xml$addTag(j, Kleborate_doc[i, j])
  }
  xml$closeTag()
}
xml$closeTag("kleborate")
#create child node - card db (output from abricate)
xml$addTag("card_resistance", close = FALSE)
for (i in 1:nrow(Card_final)) {
  xml$addTag("gene", close=FALSE)
  for (j in names(Card_final)) {
    xml$addTag(j, Card_final[i, j])
  }
  xml$closeTag()
}
xml$closeTag("card_resistance")
#create child node - card db (output from abricate)
xml$addTag("plasmidfinder_resistance", close = FALSE)
for (i in 1:nrow(Plasmid_final)) {
  xml$addTag("gene", close=FALSE)
  for (j in names(Plasmid_final)) {
    xml$addTag(j, Plasmid_final[i, j])
  }
  xml$closeTag()
}
xml$closeTag("plasmidfinder_resistance")
#create child node - card db (output from abricate)
xml$addTag("inhouse_resistance", close = FALSE)
for (i in 1:nrow(InHouse_final)) {
  xml$addTag("gene", close=FALSE)
  for (j in names(InHouse_final)) {
    xml$addTag(j, InHouse_final[i, j])
  }
  xml$closeTag()
}
xml$closeTag("inhouse_resistance")
xml$closeTag("data")

###Save XML file
saveXML(xml, file=OUTFILE, compression=0, indent=TRUE, prefix = '<?xml version="1.0"?>\n', doctype = NULL)


