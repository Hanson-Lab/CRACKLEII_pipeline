##################
#USAGE - Rscript $XML_ShortRead.R [XML File] [Resistance Database] [Output File]

#For use with the CRACKLEII_WGS_SingleIsolate.sh v1.0 pipeline 

#Author: Blake Hanson (blake.hanson [at] uth.tmc.edu)
#Version:1.0
#Date Last updated:2018/04/19
#Updates:
##################

###Arguments
args<-commandArgs(TRUE)

###Libraries
library(Biostrings)
library(ape)
library(XML)

###Variables
SAMPLE_NAME=args[1]
SITE=args[2]
LIBRARY_PREP=args[3]
DATE=args[4]
PIPELINEVER=args[5]
VERSIONDOC <- args[6]
CARD_DB <- args[7]
CARD_Consensus <- args[8]
CARD_Sanity <- args[9]
Plasmid_DB <- args[10]
Plasmid_Consensus <- args[11]
Plasmid_Sanity <- args[12]
InHouse_DB <- args[13]
InHouse_Consensus <- args[14]
InHouse_Sanity <- args[15]
CARDLIST <- args[16]
PLASMIDLIST <- args[17]
INHOUSELIST <- args[18]
OutputFile <- args[19]

###Code
#Get Sample Name
FileName_temp <- basename(OutputFile)
FileName_temp2 <- gsub(DATE, "", FileName_temp)
FileName_temp3 <- gsub(PIPELINEVER, "", FileName_temp2)
FileName <- gsub ("__ShortReadAligned_ResistanceGenes_.xml", "", FileName_temp3)

#Import version document
Version <- read.csv(VERSIONDOC, header=FALSE, row.names = 1, col.names = c("","version"))

###Import database specific mechanisms for the three  files (prepopulated)
#CARD
Card_gene_list <- read.table(CARDLIST, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
colnames(Card_gene_list) <- c("GeneName", "Description")
#Plasmidfinder
Plasmid_gene_list <- read.table(PLASMIDLIST, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
colnames(Plasmid_gene_list) <- c("GeneName", "Description")
#In_house_db
Inhouse_gene_list <- read.table(INHOUSELIST, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
colnames(Inhouse_gene_list) <- c("GeneName")

#Set mismatch parameters
SubMatrix <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
N <- c(0,0,0,0)
SubMatrix <- cbind(SubMatrix, N)
N <- c(0,0,0,0,0)
SubMatrix <- rbind(SubMatrix, N)

#Import CARD files
if(file.exists(CARD_DB)){
  CARD_RefGenes <- readDNAStringSet(CARD_DB)
  CARD_GeneName = names(CARD_RefGenes)
  CARD_Reference = paste(CARD_RefGenes)
  CARD_Reference_fasta <- data.frame(CARD_GeneName, CARD_Reference) 
  names(CARD_Reference_fasta) <- c("GeneName", "Reference")
  CARD_Reference_fasta$GeneName <- gsub("card~~~", "", CARD_Reference_fasta$GeneName)
  CARD_Reference_fasta$GeneName <- gsub("~~~.*", "", CARD_Reference_fasta$GeneName)
  
  CARD_ConsensusGenes <- readDNAStringSet(CARD_Consensus)
  CARD_ConsensusGeneName = names(CARD_ConsensusGenes)
  CARD_Consensus = paste(CARD_ConsensusGenes)
  CARD_Consensus_fasta <- data.frame(CARD_ConsensusGeneName, CARD_Consensus)  
  names(CARD_Consensus_fasta) <- c("GeneName", "Consensus")
  CARD_Consensus_fasta$GeneName <- gsub("card~~~", "", CARD_Consensus_fasta$GeneName)
  CARD_Consensus_fasta$GeneName <- gsub("~~~.*", "", CARD_Consensus_fasta$GeneName)
  
  CARD_SanityGenes <- readDNAStringSet(CARD_Sanity)
  CARD_SanityGeneName = names(CARD_SanityGenes)
  CARD_Sanity = paste(CARD_SanityGenes)
  CARD_Sanity_fasta <- data.frame(CARD_SanityGeneName, CARD_Sanity)  
  names(CARD_Sanity_fasta) <- c("GeneName", "Sanity")
  CARD_Sanity_fasta$GeneName <- gsub("card~~~", "", CARD_Sanity_fasta$GeneName)
  CARD_Sanity_fasta$GeneName <- gsub("~~~.*", "", CARD_Sanity_fasta$GeneName)
  
  CARD_Merge <- merge(CARD_Reference_fasta, CARD_Consensus_fasta, by="GeneName")
  CARD_Merge$Consensus <- as.character(CARD_Merge$Consensus)
  CARD_Merge$Reference <- as.character(CARD_Merge$Reference)
  
  CARD_Merge$Match <- CARD_Merge$Consensus==CARD_Merge$Reference
  CARD_Merge$ConsensusLength <- nchar(CARD_Merge$Consensus)
  CARD_Merge$ReferenceLength <- nchar(CARD_Merge$Reference)
  
  CARD_Merge$ConsensusAlign <- pattern(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                    gapOpening = -2, gapExtension = -4))
  CARD_Merge$ReferenceAlign <- subject(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                    gapOpening = -2, gapExtension = -4))
  CARD_Merge$ConsensusAlignLength <- nchar(CARD_Merge$ConsensusAlign)
  CARD_Merge$ReferenceAlignLength <- nchar(CARD_Merge$ReferenceAlign)
  CARD_Merge$PercentIdentity <- pid(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                 gapOpening = -2, gapExtension = -4))
  CARD_Merge$PercentCoverage <- (CARD_Merge$ReferenceAlignLength/CARD_Merge$ReferenceLength)*100
  CARD_Merge$Matches <- nmatch(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                            gapOpening = -2, gapExtension = -4))
  CARD_Merge$Insertions <- nindel(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                               gapOpening = -2, gapExtension = -4))@insertion[,1]
  CARD_Merge$Deletions <- nindel(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                              gapOpening = -2, gapExtension = -4))@deletion[,1]
  CARD_Merge$Mismatches <- nmismatch(pairwiseAlignment(pattern=CARD_Merge$Consensus, subject=CARD_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                  gapOpening = -2, gapExtension = -4))
  CARD_Merge_temp <- merge(CARD_Merge, CARD_Sanity_fasta, by="GeneName")
  CARD_Merge_temp$IssueFlag <- ifelse(CARD_Merge_temp$Consensus==CARD_Merge_temp$Sanity, "OK", "PotentialProblem")
  CARD_Merge_temp$ConsensusAlign <- as.character(CARD_Merge_temp$ConsensusAlign)
  CARD_Merge_temp$ReferenceAlign <- as.character(CARD_Merge_temp$ReferenceAlign)
  
  CARD_Merge_Final <- merge(Card_gene_list, CARD_Merge_temp, by="GeneName", all=TRUE)
} else {
  CARD_Merge_temp <- data.frame(matrix(ncol=18, nrow=1))
  ColNames <- c("GeneName","Reference","Consensus","Match","ConsensusLength","ReferenceLength" ,"ConsensusAlign",
                "ReferenceAlign","ConsensusAlignLength","ReferenceAlignLength","PercentIdentity","PercentCoverage",
                "Matches","Insertions","Deletions","Mismatches","Sanity","IssueFlag")
  colnames(CARD_Merge_temp) <- ColNames
  CARD_Merge_Final <- merge(Card_gene_list, CARD_Merge_temp, by="GeneName", all=TRUE)
}

#Import PlasmidFinder files
if(file.exists(Plasmid_DB)){
  Plasmid_RefGenes <- readDNAStringSet(Plasmid_DB)
  Plasmid_GeneName = names(Plasmid_RefGenes)
  Plasmid_Reference = paste(Plasmid_RefGenes)
  Plasmid_Reference_fasta <- data.frame(Plasmid_GeneName, Plasmid_Reference) 
  names(Plasmid_Reference_fasta) <- c("GeneName", "Reference")
  Plasmid_Reference_fasta$GeneName <- gsub("plasmidfinder~~~", "", Plasmid_Reference_fasta$GeneName)
  Plasmid_Reference_fasta$GeneName <- gsub("~~~.*", "", Plasmid_Reference_fasta$GeneName)
  
  Plasmid_ConsensusGenes <- readDNAStringSet(Plasmid_Consensus)
  Plasmid_ConsensusGeneName = names(Plasmid_ConsensusGenes)
  Plasmid_Consensus = paste(Plasmid_ConsensusGenes)
  Plasmid_Consensus_fasta <- data.frame(Plasmid_ConsensusGeneName, Plasmid_Consensus)  
  names(Plasmid_Consensus_fasta) <- c("GeneName", "Consensus")
  Plasmid_Consensus_fasta$GeneName <- gsub("plasmidfinder~~~", "", Plasmid_Consensus_fasta$GeneName)
  Plasmid_Consensus_fasta$GeneName <- gsub("~~~.*", "", Plasmid_Consensus_fasta$GeneName)
  
  Plasmid_SanityGenes <- readDNAStringSet(Plasmid_Sanity)
  Plasmid_SanityGeneName = names(Plasmid_SanityGenes)
  Plasmid_Sanity = paste(Plasmid_SanityGenes)
  Plasmid_Sanity_fasta <- data.frame(Plasmid_SanityGeneName, Plasmid_Sanity)  
  names(Plasmid_Sanity_fasta) <- c("GeneName", "Sanity")
  Plasmid_Sanity_fasta$GeneName <- gsub("plasmidfinder~~~", "", Plasmid_Sanity_fasta$GeneName)
  Plasmid_Sanity_fasta$GeneName <- gsub("~~~.*", "", Plasmid_Sanity_fasta$GeneName)
  
  Plasmid_Merge <- merge(Plasmid_Reference_fasta, Plasmid_Consensus_fasta, by="GeneName")
  Plasmid_Merge$Consensus <- as.character(Plasmid_Merge$Consensus)
  Plasmid_Merge$Reference <- as.character(Plasmid_Merge$Reference)
  
  Plasmid_Merge$Match <- Plasmid_Merge$Consensus==Plasmid_Merge$Reference
  Plasmid_Merge$ConsensusLength <- nchar(Plasmid_Merge$Consensus)
  Plasmid_Merge$ReferenceLength <- nchar(Plasmid_Merge$Reference)
  
  Plasmid_Merge$ConsensusAlign <- pattern(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                         gapOpening = -2, gapExtension = -4))
  Plasmid_Merge$ReferenceAlign <- subject(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                         gapOpening = -2, gapExtension = -4))
  Plasmid_Merge$ConsensusAlignLength <- nchar(Plasmid_Merge$ConsensusAlign)
  Plasmid_Merge$ReferenceAlignLength <- nchar(Plasmid_Merge$ReferenceAlign)
  Plasmid_Merge$PercentIdentity <- pid(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                      gapOpening = -2, gapExtension = -4))
  Plasmid_Merge$PercentCoverage <- (Plasmid_Merge$ReferenceAlignLength/Plasmid_Merge$ReferenceLength)*100
  Plasmid_Merge$Matches <- nmatch(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                 gapOpening = -2, gapExtension = -4))
  Plasmid_Merge$Insertions <- nindel(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                    gapOpening = -2, gapExtension = -4))@insertion[,1]
  Plasmid_Merge$Deletions <- nindel(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                   gapOpening = -2, gapExtension = -4))@deletion[,1]
  Plasmid_Merge$Mismatches <- nmismatch(pairwiseAlignment(pattern=Plasmid_Merge$Consensus, subject=Plasmid_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                       gapOpening = -2, gapExtension = -4))
  Plasmid_Merge_temp <- merge(Plasmid_Merge, Plasmid_Sanity_fasta, by="GeneName")
  Plasmid_Merge_temp$IssueFlag <- ifelse(Plasmid_Merge_temp$Consensus==Plasmid_Merge_temp$Sanity, "OK", "PotentialProblem")
  Plasmid_Merge_temp$ConsensusAlign <- as.character(Plasmid_Merge_temp$ConsensusAlign)
  Plasmid_Merge_temp$ReferenceAlign <- as.character(Plasmid_Merge_temp$ReferenceAlign)
  
  Plasmid_Merge_Final <- merge(Plasmid_gene_list, Plasmid_Merge_temp, by="GeneName", all=TRUE)
} else {
  Plasmid_Merge_temp <- data.frame(matrix(ncol=18, nrow=1))
  ColNames <- c("GeneName","Reference","Consensus","Match","ConsensusLength","ReferenceLength" ,"ConsensusAlign",
                "ReferenceAlign","ConsensusAlignLength","ReferenceAlignLength","PercentIdentity","PercentCoverage",
                "Matches","Insertions","Deletions","Mismatches","Sanity","IssueFlag")
  colnames(Plasmid_Merge_temp) <- ColNames
  Plasmid_Merge_Final <- merge(Plasmid_gene_list, Plasmid_Merge_temp, by="GeneName", all=TRUE)
}

#Import InHouse files
if(file.exists(InHouse_DB)){
  InHouse_RefGenes <- readDNAStringSet(InHouse_DB)
  InHouse_GeneName = names(InHouse_RefGenes)
  InHouse_Reference = paste(InHouse_RefGenes)
  InHouse_Reference_fasta <- data.frame(InHouse_GeneName, InHouse_Reference) 
  names(InHouse_Reference_fasta) <- c("GeneName", "Reference")
  
  InHouse_ConsensusGenes <- readDNAStringSet(InHouse_Consensus)
  InHouse_ConsensusGeneName = names(InHouse_ConsensusGenes)
  InHouse_Consensus = paste(InHouse_ConsensusGenes)
  InHouse_Consensus_fasta <- data.frame(InHouse_ConsensusGeneName, InHouse_Consensus)  
  names(InHouse_Consensus_fasta) <- c("GeneName", "Consensus")
  
  InHouse_SanityGenes <- readDNAStringSet(InHouse_Sanity)
  InHouse_SanityGeneName = names(InHouse_SanityGenes)
  InHouse_Sanity = paste(InHouse_SanityGenes)
  InHouse_Sanity_fasta <- data.frame(InHouse_SanityGeneName, InHouse_Sanity)  
  names(InHouse_Sanity_fasta) <- c("GeneName", "Sanity")
  
  InHouse_Merge <- merge(InHouse_Reference_fasta, InHouse_Consensus_fasta, by="GeneName")
  InHouse_Merge$Consensus <- as.character(InHouse_Merge$Consensus)
  InHouse_Merge$Reference <- as.character(InHouse_Merge$Reference)
  
  InHouse_Merge$Match <- InHouse_Merge$Consensus==InHouse_Merge$Reference
  InHouse_Merge$ConsensusLength <- nchar(InHouse_Merge$Consensus)
  InHouse_Merge$ReferenceLength <- nchar(InHouse_Merge$Reference)
  
  InHouse_Merge$ConsensusAlign <- pattern(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                         gapOpening = -2, gapExtension = -4))
  InHouse_Merge$ReferenceAlign <- subject(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                         gapOpening = -2, gapExtension = -4))
  InHouse_Merge$ConsensusAlignLength <- nchar(InHouse_Merge$ConsensusAlign)
  InHouse_Merge$ReferenceAlignLength <- nchar(InHouse_Merge$ReferenceAlign)
  InHouse_Merge$PercentIdentity <- pid(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                      gapOpening = -2, gapExtension = -4))
  InHouse_Merge$PercentCoverage <- (InHouse_Merge$ReferenceAlignLength/InHouse_Merge$ReferenceLength)*100
  InHouse_Merge$Matches <- nmatch(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                 gapOpening = -2, gapExtension = -4))
  InHouse_Merge$Insertions <- nindel(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                    gapOpening = -2, gapExtension = -4))@insertion[,1]
  InHouse_Merge$Deletions <- nindel(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                   gapOpening = -2, gapExtension = -4))@deletion[,1]
  InHouse_Merge$Mismatches <- nmismatch(pairwiseAlignment(pattern=InHouse_Merge$Consensus, subject=InHouse_Merge$Reference, substitutionMatrix=SubMatrix, type="local",
                                                       gapOpening = -2, gapExtension = -4))
  InHouse_Merge_temp <- merge(InHouse_Merge, InHouse_Sanity_fasta, by="GeneName")
  InHouse_Merge_temp$IssueFlag <- ifelse(InHouse_Merge_temp$Consensus==InHouse_Merge_temp$Sanity, "OK", "PotentialProblem")
  InHouse_Merge_temp$ConsensusAlign <- as.character(InHouse_Merge_temp$ConsensusAlign)
  InHouse_Merge_temp$ReferenceAlign <- as.character(InHouse_Merge_temp$ReferenceAlign)
  
  InHouse_Merge_Final <- merge(Inhouse_gene_list, InHouse_Merge_temp, by="GeneName", all=TRUE)
} else {
  InHouse_Merge_temp <- data.frame(matrix(ncol=18, nrow=1))
  ColNames <- c("GeneName","Reference","Consensus","Match","ConsensusLength","ReferenceLength" ,"ConsensusAlign",
                "ReferenceAlign","ConsensusAlignLength","ReferenceAlignLength","PercentIdentity","PercentCoverage",
                "Matches","Insertions","Deletions","Mismatches","Sanity","IssueFlag")
  colnames(InHouse_Merge_temp) <- ColNames
  InHouse_Merge_Final <- merge(Inhouse_gene_list, InHouse_Merge_temp, by="GeneName", all=TRUE)
}

###CREATE XML
#Open empty xml tree
xml <- xmlTree()
#Create overarching node "data"
xml$addTag("data", close=FALSE)
#Create child node - lab
xml$addTag("sample_name", FileName)
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
for (i in 1:nrow(Version)) {
  xml$addTag(rownames(Version)[i], close=FALSE)
  for (j in names(Version)) {
    xml$addTag(j, Version[i, j])
  }
  xml$closeTag()
}
xml$closeTag("software_versions")
#create child node - card_resistance_short_read (results for CARD genes)
xml$addTag("card_resistance_short_read", close = FALSE)
for (i in 1:nrow(CARD_Merge_Final)) {
  xml$addTag("gene", close=FALSE)
  for (j in names(CARD_Merge_Final)) {
    xml$addTag(j, CARD_Merge_Final[i, j])
  }
  xml$closeTag()
}
xml$closeTag("card_resistance_short_read")
#create child node - card_resistance_short_read (results for CARD genes)
xml$addTag("plasmid_finder_short_read", close = FALSE)
for (i in 1:nrow(Plasmid_Merge_Final)) {
  xml$addTag("gene", close=FALSE)
  for (j in names(Plasmid_Merge_Final)) {
    xml$addTag(j, Plasmid_Merge_Final[i, j])
  }
  xml$closeTag()
}
xml$closeTag("plasmid_finder_short_read")
#create child node - card_resistance_short_read (results for CARD genes)
xml$addTag("in_house_resistance_short_read", close = FALSE)
for (i in 1:nrow(InHouse_Merge_Final)) {
  xml$addTag("gene", close=FALSE)
  for (j in names(InHouse_Merge_Final)) {
    xml$addTag(j, InHouse_Merge_Final[i, j])
  }
  xml$closeTag()
}
xml$closeTag("in_house_resistance_short_read")
xml$closeTag("data")

###Save XML file
saveXML(xml, file=OutputFile, compression=0, indent=TRUE, prefix = '<?xml version="1.0"?>\n', doctype = NULL)
