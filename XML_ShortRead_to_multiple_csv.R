##################
#USAGE - Rscript $XML_ShortRead_to_multiple_csv.R [LIST Of XML Files] [Location of XML files] [% Identity Cutoff] [% Coverage Cutoff] [Output Directory]

#For use with the CRACKLEII_WGS_SingleIsolate.sh v1.0 pipeline 

#Make sure the xml files are all in the same folder (can also use symlinks)
#Example to get list of files: ls -1 *xml > ListFiles.txt

#Author: Blake Hanson (blake.hanson [at] uth.tmc.edu)
#Version:1.0
#Date Last updated:2018/04/30
#Updates:
##################

###Arguments
args<-commandArgs(TRUE)

###Libraries
library(XML)
library(dplyr)
library(plyr)

###Functions
ProcessXML <- function(file, FileLocation, OutputLocation, Coverage, Identity){
  print(paste("Starting file: ", file, sep=""))
  #Get pipeline version
  version <- gsub(".*ResistanceGenes_", "", file)
  version <- gsub(".xml", "", version)

  #Parse XML file 
  xmlfile <- xmlParse(paste(FileLocation, file, sep=""))
  
  #Extract sample name node
  sample_name <- as.character(xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/sample_name"), stringsAsFactors = FALSE))
  
  #Extract pipeline version node
  pipeline_version <- as.character(xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/pipeline"), stringsAsFactors = FALSE))
  
  #Extract lab node
  lab <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/lab"), stringsAsFactors = FALSE)
  lab <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, lab)
  
  #Set variable for data processing site
  site <- lab$site[1]
  
  #Extract software_versions node
  software_versions <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/software_versions"), stringsAsFactors = FALSE)
  software_versions <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, software_versions)
  
  #Extract card node
  card_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/card_resistance_short_read/*"), stringsAsFactors = FALSE)
  #Remove "NA" characters and set identity and percent_coverage as numeric values for filtering
  card_output$PercentIdentity[card_output$PercentIdentity == "NA"] <- NA
  
  card_output$PercentCoverage[card_output$PercentCoverage == "NA"] <- NA
  card_output$PercentIdentity <- as.numeric(card_output$PercentIdentity)
  card_output$PercentCoverage <- as.numeric(card_output$PercentCoverage)
  #Filter card results on Identity variable and if identity is less than value specified, set NA for match fields
  card_output_filter_i <- card_output %>% mutate_at(.vars = c("Description", "Reference",  "Consensus",
                                                              "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                                              "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                                              "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                                              "Deletions", "Mismatches", "Sanity", "IssueFlag"),
                                                    funs(ifelse(card_output$PercentCoverage < Coverage, NA, .)))
  card_output_filter <- card_output_filter_i %>% mutate_at(.vars = c("Description", "Reference",  "Consensus",
                                                                     "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                                                     "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                                                     "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                                                     "Deletions", "Mismatches", "Sanity", "IssueFlag"),
                                                           funs(ifelse(card_output_filter_i$PercentIdentity < Identity, NA, .)))
  #Check for duplicates, remove duplicates that are NA, and count number of gene duplications
  card_output_filter$duplicate <- duplicated(card_output_filter$GeneName)|duplicated(card_output_filter$GeneName, fromLast=TRUE)
  card_output_filter <-card_output_filter[!(card_output_filter$duplicate=="TRUE" & is.na(card_output_filter$Reference)),]
  card_output_count <- count(card_output_filter$GeneName)
  card_output_temp <- merge(card_output_filter, card_output_count, by.x="GeneName", by.y="x", all.x=TRUE)
  #Set identity for absent genes to NA
  card_output_temp$freq[is.na(card_output_temp$PercentIdentity)] <- "NA"
  #Keep only gene name and frequency for final table 
  drop_vars <- names(card_output_temp) %in% c("Description", "Reference",  "Consensus",
                                              "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                              "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                              "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                              "Deletions", "Mismatches", "Sanity", "IssueFlag", "duplicate")
  card_subset <- card_output_temp[!drop_vars]
  #Remove duplicate records (now that the count is maintained)
  card_subset <- unique(card_subset)
  #Transpose dataframe for future merging and set row name as sample ID
  card_final <- as.data.frame(t(card_subset))
  colnames(card_final) <- as.character(unlist(card_final[1,]))
  card_final <- card_final[-1,]
  card_final <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, card_final)
  names(card_final) <- sub("_$","",gsub('([_])\\1+', '\\1', gsub('([[:punct:]])|\\s+','_', names(card_final))))
  
  #Extract plasmidfinder node
  plasmidfinder_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/plasmid_finder_short_read/*"), stringsAsFactors = FALSE)
  #Remove "NA" characters and set identity and percent_coverage as numeric values for filtering
  plasmidfinder_output$PercentIdentity[plasmidfinder_output$PercentIdentity == "NA"] <- NA
  
  plasmidfinder_output$PercentCoverage[plasmidfinder_output$PercentCoverage == "NA"] <- NA
  plasmidfinder_output$PercentIdentity <- as.numeric(plasmidfinder_output$PercentIdentity)
  plasmidfinder_output$PercentCoverage <- as.numeric(plasmidfinder_output$PercentCoverage)
  #Filter plasmidfinder results on Identity variable and if identity is less than value specified, set NA for match fields
  plasmidfinder_output_filter_i <- plasmidfinder_output %>% mutate_at(.vars = c("Description", "Reference",  "Consensus",
                                                              "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                                              "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                                              "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                                              "Deletions", "Mismatches", "Sanity", "IssueFlag"),
                                                    funs(ifelse(plasmidfinder_output$PercentCoverage < Coverage, NA, .)))
  plasmidfinder_output_filter <- plasmidfinder_output_filter_i %>% mutate_at(.vars = c( "Description", "Reference",  "Consensus",
                                                                     "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                                                     "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                                                     "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                                                     "Deletions", "Mismatches", "Sanity", "IssueFlag"),
                                                           funs(ifelse(PercentIdentity < Identity, NA, .)))
  #Check for duplicates, remove duplicates that are NA, and count number of gene duplications
  plasmidfinder_output_filter$duplicate <- duplicated(plasmidfinder_output_filter$GeneName)|duplicated(plasmidfinder_output_filter$GeneName, fromLast=TRUE)
  plasmidfinder_output_filter <-plasmidfinder_output_filter[!(plasmidfinder_output_filter$duplicate=="TRUE" & is.na(plasmidfinder_output_filter$Reference)),]
  plasmidfinder_output_count <- count(plasmidfinder_output_filter$GeneName)
  plasmidfinder_output_temp <- merge(plasmidfinder_output_filter, plasmidfinder_output_count, by.x="GeneName", by.y="x", all.x=TRUE)
  #Set identity for absent genes to NA
  plasmidfinder_output_temp$freq[is.na(plasmidfinder_output_temp$PercentIdentity)] <- "NA"
  #Keep only gene name and frequency for final table 
  drop_vars <- names(plasmidfinder_output_temp) %in% c("Description", "Reference",  "Consensus",
                                              "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                              "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                              "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                              "Deletions", "Mismatches", "Sanity", "IssueFlag", "duplicate")
  plasmidfinder_subset <- plasmidfinder_output_temp[!drop_vars]
  #Remove duplicate records (now that the count is maintained)
  plasmidfinder_subset <- unique(plasmidfinder_subset)
  #Transpose dataframe for future merging and set row name as sample ID
  plasmidfinder_final <- as.data.frame(t(plasmidfinder_subset))
  colnames(plasmidfinder_final) <- as.character(unlist(plasmidfinder_final[1,]))
  plasmidfinder_final <- plasmidfinder_final[-1,]
  plasmidfinder_final <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, plasmidfinder_final)
  names(plasmidfinder_final) <- sub("_$","",gsub('([_])\\1+', '\\1', gsub('([[:punct:]])|\\s+','_', names(plasmidfinder_final))))
  
  #Extract in_house database node

  in_house_db_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/in_house_resistance_short_read/*"), stringsAsFactors = FALSE)
  #Remove "NA" characters and set identity and percent_coverage as numeric values for filtering
  in_house_db_output$PercentIdentity[in_house_db_output$PercentIdentity == "NA"] <- NA
  
  in_house_db_output$PercentCoverage[in_house_db_output$PercentCoverage == "NA"] <- NA
  in_house_db_output$PercentIdentity <- as.numeric(in_house_db_output$PercentIdentity)
  in_house_db_output$PercentCoverage <- as.numeric(in_house_db_output$PercentCoverage)
  #Filter in_house_db results on Identity variable and if identity is less than value specified, set NA for match fields
  in_house_db_output_filter_i <- in_house_db_output %>% mutate_at(.vars = c("Reference",  "Consensus",
                                                              "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                                              "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                                              "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                                              "Deletions", "Mismatches", "Sanity", "IssueFlag"),
                                                    funs(ifelse(in_house_db_output$PercentCoverage < Coverage, NA, .)))
  in_house_db_output_filter <- in_house_db_output_filter_i %>% mutate_at(.vars = c("Reference",  "Consensus",
                                                                     "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                                                     "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                                                     "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                                                     "Deletions", "Mismatches", "Sanity", "IssueFlag"),
                                                           funs(ifelse(PercentIdentity < Identity, NA, .)))
  #Check for duplicates, remove duplicates that are NA, and count number of gene duplications
  in_house_db_output_filter$duplicate <- duplicated(in_house_db_output_filter$GeneName)|duplicated(in_house_db_output_filter$GeneName, fromLast=TRUE)
  in_house_db_output_filter <-in_house_db_output_filter[!(in_house_db_output_filter$duplicate=="TRUE" & is.na(in_house_db_output_filter$Reference)),]
  in_house_db_output_count <- count(in_house_db_output_filter$GeneName)
  in_house_db_output_temp <- merge(in_house_db_output_filter, in_house_db_output_count, by.x="GeneName", by.y="x", all.x=TRUE)
  #Set identity for absent genes to NA
  in_house_db_output_temp$freq[is.na(in_house_db_output_temp$PercentIdentity)] <- "NA"
  #Keep only gene name and frequency for final table 
  drop_vars <- names(in_house_db_output_temp) %in% c("Description", "Reference",  "Consensus",
                                              "Match", "ConsensusLength", "ReferenceLength",  "ConsensusAlign",
                                              "ReferenceAlign", "ConsensusAlignLength", "ReferenceAlignLength",
                                              "PercentIdentity", "PercentCoverage", "Matches", "Insertions",
                                              "Deletions", "Mismatches", "Sanity", "IssueFlag", "duplicate")
  in_house_db_subset <- in_house_db_output_temp[!drop_vars]
  #Remove duplicate records (now that the count is maintained)
  in_house_db_subset <- unique(in_house_db_subset)
  #Transpose dataframe for future merging and set row name as sample ID
  in_house_db_final <- as.data.frame(t(in_house_db_subset))
  colnames(in_house_db_final) <- as.character(unlist(in_house_db_final[1,]))
  in_house_db_final <- in_house_db_final[-1,]
  in_house_db_final <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, in_house_db_final)
  names(in_house_db_final) <- sub("_$","",gsub('([_])\\1+', '\\1', gsub('([[:punct:]])|\\s+','_', names(in_house_db_final))))
  
  #Write csv files to temp directory
  print("Writing csv files")
  Lab_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "lab.csv", sep="")
  Software_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "software.csv", sep="")
  Card_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "card.csv", sep="")
  Plasmid_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "plasmid.csv", sep="")
  InHouse_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "inhouse.csv", sep="")
  
  write.csv(lab, Lab_Output, row.names = FALSE)
  write.csv(software_versions, Software_Output, row.names = FALSE)
  write.csv(card_final, Card_Output, row.names = FALSE)
  write.csv(plasmidfinder_final, Plasmid_Output, row.names = FALSE)
  write.csv(in_house_db_final, InHouse_Output, row.names = FALSE)
}

#Combine files and output tables for each node extracted above
MergeTables <- function(FileLocation, filesuffix){
  #Pattern to create a list of files to merge
  MergeFilePattern <- paste("^.*", filesuffix, ".csv", sep="")
  #Creates list of files to merge
  filenames <- list.files(path =FileLocation, pattern = MergeFilePattern, ignore.case = TRUE)
  PathFileNames = list()
  #Uses the filenames list to create a list of file names with absolute paths
  filenames <- for(filename in filenames){filename <- paste(FileLocation, "/", filename, sep="")
  PathFileNames[[filename]] <- filename}
  #Reads in each csv file and then merges them using the consistent format
  datalist <- lapply(PathFileNames, function(x){read.csv(file=x, header=TRUE, stringsAsFactors = FALSE)})
  Reduce(function(...) merge(..., all = T), datalist)
}

###Variables
XMLListFile <- args[1]
FileLocation <- args[2]
Identity <- as.numeric(args[3])
Coverage <- as.numeric(args[4])
OutputDir <- args[5]

###Code
#Create output directory (will output warning if directory exists)
dir.create(OutputDir)

#Create temporary directory 
temp_dir <- paste(OutputDir, "/temp", sep="")
dir.create(temp_dir)

#Setting time variable
time <- format(Sys.time(), "%Y-%m-%d-%H-%M")

#Read in list of XML files
XMLList <- read.table(XMLListFile)

#Parse over XML files and generate tables for each node in temp directory
for(i in 1:nrow(XMLList)) {
  datafr <- XMLList[i,]
  ProcessXML(datafr, FileLocation, temp_dir, Coverage, Identity)
}

#Create merged tables for each node 
print("Combining All Samples into tables by xml section")
LabMergedTable <- MergeTables(temp_dir, "lab")
SoftwareMergedTable <- MergeTables(temp_dir, "software")
CardMergedTable <- MergeTables(temp_dir, "card")
PlasmidMergedTable <- MergeTables(temp_dir, "plasmid")
InHouseMergedTable <- MergeTables(temp_dir, "inhouse")

###Create card, plasmid, and inhouse data frames that only include present genes
CardMergedTable_Filter <- CardMergedTable[,colSums(is.na(CardMergedTable))<nrow(CardMergedTable)]
PlasmidMergedTable_Filter <- PlasmidMergedTable[,colSums(is.na(PlasmidMergedTable))<nrow(PlasmidMergedTable)]
InHouseMergedTable_Filter <- InHouseMergedTable[,colSums(is.na(InHouseMergedTable))<nrow(InHouseMergedTable)]

###Write out merged dataframes 
print("Writing merged csv files")
Lab_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_lab_merged_V1.0_", time, ".csv", sep="")
Software_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_software_merged_V1.0_", time, ".csv", sep="")
Card_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_card_all_genes_merged_V1.0_", time, ".csv", sep="")
Plasmid_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_plasmid_all_genes_merged_V1.0_", time, ".csv", sep="")
InHouse_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_inhouse_all_genes_merged_V1.0_", time, ".csv", sep="")
Card_sub_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_card_only_detected_genes_merged_V1.0_", time, ".csv", sep="")
Plasmid_sub_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_plasmid_only_detected_genes_merged_V1.0_", time, ".csv", sep="")
InHouse_sub_Output <- paste(OutputDir, "/CRACKLE_II_ShortRead_inhouse_only_detected_genes_merged_V1.0_", time, ".csv", sep="")

write.csv(LabMergedTable, Lab_Output, row.names = FALSE)
write.csv(SoftwareMergedTable, Software_Output, row.names = FALSE)
write.csv(CardMergedTable, Card_Output, row.names = FALSE)
write.csv(PlasmidMergedTable, Plasmid_Output, row.names = FALSE)
write.csv(InHouseMergedTable, InHouse_Output, row.names = FALSE)
write.csv(CardMergedTable_Filter, Card_sub_Output, row.names = FALSE)
write.csv(PlasmidMergedTable_Filter, Plasmid_sub_Output, row.names = FALSE)
write.csv(InHouseMergedTable_Filter, InHouse_sub_Output, row.names = FALSE)


###Remove temp directory
setwd(OutputDir)
unlink("temp", recursive=TRUE, force=TRUE)

