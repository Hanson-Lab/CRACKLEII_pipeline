##################
#USAGE - Rscript $XML_to_multiple_csv.R [LIST Of XML Files] [Location of XML files] [% Identity Cutoff] [% Coverage Cutoff] [Output Directory]

#For use with the CRACKLEII_WGS_SingleIsolate.sh v1.0 pipeline 

#Make sure the xml files are all in the same folder (can also use symlinks)
#Example to get list of files: ls -1 *xml > ListFiles.txt

#Author: Blake Hanson (blake.hanson [at] uth.tmc.edu)
#Version:1.0
#Date Last updated:2018/02/14
#Updates:
##################

###Arguments
args<-commandArgs(TRUE)

###Libraries
library(XML)
library(plyr)
library(dplyr)

###Functions
ProcessXML <- function(file, FileLocation, OutputLocation, Coverage, Identity){
  print(paste("Starting file: ", file, sep=""))
  #Get pipeline version
  version <- gsub(".*CRACKLE_", "", file)
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
  
  #Extract raspberry result nodes
  raspberry_raw_read1 <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/raspberry_raw_read1"), stringsAsFactors = FALSE)
  raspberry_raw_read1 <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, raspberry_raw_read1)
  
  raspberry_raw_read2 <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/raspberry_raw_read2"), stringsAsFactors = FALSE)
  raspberry_raw_read2 <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, raspberry_raw_read2)
  
  raspberry_trim_read1 <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/raspberry_trim_read1"), stringsAsFactors = FALSE)
  raspberry_trim_read1 <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, raspberry_trim_read1)
  
  raspberry_trim_read2 <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/raspberry_trim_read2"), stringsAsFactors = FALSE)
  raspberry_trim_read2 <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, raspberry_trim_read2)
  
  #Extract Quast assembly stats node
  assembly_stats <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/assembly_stats"), stringsAsFactors = FALSE)
  assembly_stats <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, assembly_stats)
  
  #Extract strain seeker output node
  strain_seeker_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/strain_id/*"), stringsAsFactors = FALSE)
  strain_seeker_output <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, strain_seeker_output)
  
  #Extract mlst node
  mlst <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/mlst/*"), stringsAsFactors = FALSE)
  #If else statement to check if there is data in the mlst section, if not, deal with differently.
  if(nrow(mlst) <= 3){
    mlst_df <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, 
                     species = "NA", st_type = "NA")
  } else{
    #Transpose mlst data
    mlst_df <- as.data.frame(t(mlst))
    #Pull MLST species ID from the first row of the first column
    mlst_species <- mlst_df[1,1]
    #Pull MLST ST number from first row of second column
    mlst_ST <- mlst_df[1,2]
    #Remove first two columns that contain the species and ST number
    mlst_df$V1 <- mlst_df$V2 <- NULL
    #Remove empty first row from mlst loci information
    mlst_df <- mlst_df[-1,]
    #Set gene names to column headers
    colnames(mlst_df) <- as.character(unlist(mlst_df[1,]))
    #Remove row containing gene names
    mlst_df <- mlst_df[-1,]
    #Add sample name, pipeline version, species name, and ST
    mlst_df <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, 
                     species = mlst_species, st_type = mlst_ST, mlst_df)
    #Subsitute a pipe for any comma in the mlst_df dataframe
    mlst_df <- data.frame(lapply(mlst_df, function(x) {
      gsub(",", "|", x)
    }))
  }
  
  #Extract kleborate node
  kleborate_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/kleborate/*"), stringsAsFactors = FALSE)
  #Kleborate has the name of the file listed as strain, but that is not used in this pipeline so removed
  kleborate_output$strain <- NULL
  kleborate_output <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, kleborate_output)
  
  #Extract card node
  card_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/card_resistance/*"), stringsAsFactors = FALSE)
  #Remove "NA" characters and set identity and percent_coverage as numeric values for filtering
  card_output$percent_coverage[card_output$percent_coverage == "NA"] <- NA
  card_output$identity[card_output$identity == "NA"] <- NA
  card_output$percent_coverage <- as.numeric(card_output$percent_coverage)
  card_output$identity <- as.numeric(card_output$identity)
  #Filter card results on Identity variable and if identity is less than value specified, set NA for match fields
  card_output_filter_i <- card_output %>% mutate_at(.vars = c("sequence", "database", "accession", "start", "end", 
                                                              "coverage", "coverage_map", "gaps", "percent_coverage", 
                                                              "identity", "strand"),
                                                    funs(ifelse(card_output$percent_coverage < Coverage, NA, .)))
  card_output_filter <- card_output_filter_i %>% mutate_at(.vars = c("sequence", "database", "accession", "start", "end", 
                                                                     "coverage", "coverage_map", "gaps", "identity", 
								     "percent_coverage", "strand"),
                                                           funs(ifelse(identity < Identity, NA, .)))
  #Check for duplicates, remove duplicates that are NA, and count number of gene duplications
  card_output_filter$duplicate <- duplicated(card_output_filter$gene)|duplicated(card_output_filter$gene, fromLast=TRUE)
  card_output_filter <-card_output_filter[!(card_output_filter$duplicate=="TRUE" & is.na(card_output_filter$sequence)),]
  card_output_count <- count(card_output_filter, gene)
  card_output_temp <- merge(card_output_filter, card_output_count, by="gene", all.x=TRUE)
  #Set identity for absent genes to NA
  card_output_temp$freq[is.na(card_output_temp$identity)] <- "NA"
  #Keep only gene name and frequency for final table 
  drop_vars <- names(card_output_temp) %in% c("sequence", "database", "accession", "start", "end", 
                                                "coverage", "coverage_map", "gaps", "percent_coverage", 
                                                "description", "identity", "strand", "duplicate")
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
  plasmidfinder_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/plasmidfinder_resistance/*"), stringsAsFactors = FALSE)
  #Remove "NA" characters and set identity and percent_coverage as numeric values for filtering
  plasmidfinder_output$percent_coverage[plasmidfinder_output$percent_coverage == "NA"] <- NA
  plasmidfinder_output$identity[plasmidfinder_output$identity == "NA"] <- NA
  plasmidfinder_output$percent_coverage <- as.numeric(plasmidfinder_output$percent_coverage)
  plasmidfinder_output$identity <- as.numeric(plasmidfinder_output$identity)
  #Filter card results on Identity variable and if identity is less than value specified, set NA for match fields
  plasmidfinder_output_filter_i <- plasmidfinder_output %>% mutate_at(.vars = c("sequence", "database", "start", "end", 
                                                                                "coverage", "coverage_map", "gaps", "percent_coverage", 
                                                                                "identity", "strand"),
                                                                      funs(ifelse(percent_coverage < Coverage, NA, .)))
  plasmidfinder_output_filter <- plasmidfinder_output_filter_i %>% mutate_at(.vars = c("sequence", "database", "start", "end", 
                                                                                       "coverage", "coverage_map", "gaps", "identity", "strand", "percent_coverage"),
                                                                             funs(ifelse(identity < Identity, NA, .)))
  #Check for duplicates, remove duplicates that are NA, and count number of gene duplications
  plasmidfinder_output_filter$duplicate <- duplicated(plasmidfinder_output_filter$gene)|duplicated(plasmidfinder_output_filter$gene, fromLast=TRUE)
  plasmidfinder_output_filter <-plasmidfinder_output_filter[!(plasmidfinder_output_filter$duplicate=="TRUE" & is.na(plasmidfinder_output_filter$sequence)),]
  plasmidfinder_output_count <- count(plasmidfinder_output_filter, gene)
  plasmidfinder_output_temp <- merge(plasmidfinder_output_filter, plasmidfinder_output_count, by="gene", all.x=TRUE)
  #Set identity for absent genes to NA
  plasmidfinder_output_temp$freq[is.na(plasmidfinder_output_temp$identity)] <- "NA"
  #Keep only gene name and frequency for final table 
  drop_vars <- names(plasmidfinder_output_temp) %in% c("sequence", "database", "accession", "start", "end", 
                                              "coverage", "coverage_map", "gaps", "percent_coverage", 
                                              "description", "identity", "strand", "duplicate")
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
  in_house_db_output <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/inhouse_resistance/*"), stringsAsFactors = FALSE)
  #Remove "NA" characters and set identity and percent_coverage as numeric values for filtering
  in_house_db_output$percent_coverage[in_house_db_output$percent_coverage == "NA"] <- NA
  in_house_db_output$identity[in_house_db_output$identity == "NA"] <- NA
  in_house_db_output$percent_coverage <- as.numeric(in_house_db_output$percent_coverage)
  in_house_db_output$identity <- as.numeric(in_house_db_output$identity)
  #Filter card results on Identity variable and if identity is less than value specified, set NA for match fields
  in_house_db_output_filter_i <- in_house_db_output %>% mutate_at(.vars = c("description", "database", "start", "end", 
                                                                            "coverage", "coverage_map", "gaps", "percent_coverage", 
                                                                            "identity", "strand"),
                                                                  funs(ifelse(percent_coverage < Coverage, NA, .)))
  in_house_db_output_filter <- in_house_db_output_filter_i %>% mutate_at(.vars = c("description", "database", "start", "end", 
                                                                                   "coverage", "coverage_map", "gaps", "identity", "strand", "percent_coverage"),
                                                                         funs(ifelse(identity < Identity, NA, .)))
  #Check for duplicates, remove duplicates that are NA, and count number of gene duplications
  in_house_db_output_filter$duplicate <- duplicated(in_house_db_output_filter$gene)|duplicated(in_house_db_output_filter$gene, fromLast=TRUE)
  in_house_db_output_filter <- in_house_db_output_filter[!(in_house_db_output_filter$duplicate=="TRUE" & is.na(in_house_db_output_filter$start)),]
  in_house_db_output_count <- count(in_house_db_output_filter, gene)
  in_house_db_output_temp <- merge(in_house_db_output_filter, in_house_db_output_count, by="gene", all.x=TRUE)
  #Set identity for absent genes to NA
  in_house_db_output_temp$freq[is.na(in_house_db_output_temp$identity)] <- "NA"
  #Keep only gene name and frequency for final table 
  drop_vars <- names(in_house_db_output_temp) %in% c("database", "accession", "start", "end", 
                                              "coverage", "coverage_map", "gaps", "percent_coverage", 
                                              "description", "identity", "strand", "duplicate")
  in_house_subset <- in_house_db_output_temp[!drop_vars]
  #Remove duplicate records (now that the count is maintained)
  in_house_subset <- unique(in_house_subset)
  #Transpose dataframe for future merging and set row name as sample ID
  in_house_final <- as.data.frame(t(in_house_subset))
  colnames(in_house_final) <- as.character(unlist(in_house_final[1,]))
  in_house_final <- in_house_final[-1,]
  in_house_final <- cbind(sample_name = sample_name, pipeline_version = pipeline_version, site = site, in_house_final)
  names(in_house_final) <- sub("_$","",gsub('([_])\\1+', '\\1', gsub('([[:punct:]])|\\s+','_', names(in_house_final))))
  
  #Write csv files to temp directory
  print("Writing csv files")
  Lab_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "lab.csv", sep="")
  Software_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "software.csv", sep="")
  RaspberryR1_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "raspberry_R1.csv", sep="")
  RaspberryR2_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "raspberry_R2.csv", sep="")
  RaspberryT1_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "raspberry_T1.csv", sep="")
  RaspberryT2_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "raspberry_T2.csv", sep="")
  Assembly_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "assembly.csv", sep="")
  Strain_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "strain.csv", sep="")
  MLST_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "mlst.csv", sep="")
  Kleborate_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "kleborate.csv", sep="")
  Card_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "card.csv", sep="")
  Plasmid_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "plasmid.csv", sep="")
  InHouse_Output <- paste(OutputLocation, "/", sample_name, "_", version, "_", "inhouse.csv", sep="")
  
  write.csv(lab, Lab_Output, row.names = FALSE)
  write.csv(software_versions, Software_Output, row.names = FALSE)
  write.csv(raspberry_raw_read1, RaspberryR1_Output, row.names = FALSE)
  write.csv(raspberry_raw_read2, RaspberryR2_Output, row.names = FALSE)
  write.csv(raspberry_trim_read1, RaspberryT1_Output, row.names = FALSE)
  write.csv(raspberry_trim_read2, RaspberryT2_Output, row.names = FALSE)
  write.csv(assembly_stats, Assembly_Output, row.names = FALSE)
  write.csv(strain_seeker_output, Strain_Output, row.names = FALSE)
  write.csv(mlst_df, MLST_Output, row.names = FALSE)
  write.csv(kleborate_output, Kleborate_Output, row.names = FALSE)
  write.csv(card_final, Card_Output, row.names = FALSE)
  write.csv(plasmidfinder_final, Plasmid_Output, row.names = FALSE)
  write.csv(in_house_final, InHouse_Output, row.names = FALSE)
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
RaspberryR1MergedTable <- MergeTables(temp_dir, "raspberry_R1")
RaspberryR2MergedTable <- MergeTables(temp_dir, "raspberry_R2")
RaspberryT1MergedTable <- MergeTables(temp_dir, "raspberry_T1")
RaspberryT2MergedTable <- MergeTables(temp_dir, "raspberry_T2")
AssemblyMergedTable <- MergeTables(temp_dir, "assembly")
StrainMergedTable <- MergeTables(temp_dir, "strain")
MLSTMergedTable <- MergeTables(temp_dir, "mlst")
KleborateMergedTable <- MergeTables(temp_dir, "kleborate")
CardMergedTable <- MergeTables(temp_dir, "card")
PlasmidMergedTable <- MergeTables(temp_dir, "plasmid")
InHouseMergedTable <- MergeTables(temp_dir, "inhouse")

###Create card, plasmid, and inhouse data frames that only include present genes
CardMergedTable_Filter <- CardMergedTable[,colSums(is.na(CardMergedTable))<nrow(CardMergedTable)]
PlasmidMergedTable_Filter <- PlasmidMergedTable[,colSums(is.na(PlasmidMergedTable))<nrow(PlasmidMergedTable)]
InHouseMergedTable_Filter <- InHouseMergedTable[,colSums(is.na(InHouseMergedTable))<nrow(InHouseMergedTable)]

###Write out merged dataframes 
print("Writing merged csv files")
Lab_Output <- paste(OutputDir, "/CRACKLE_II_lab_merged_V1.0_", time, ".csv", sep="")
Software_Output <- paste(OutputDir, "/CRACKLE_II_software_merged_V1.0_", time, ".csv", sep="")
RaspberryR1_Output <- paste(OutputDir, "/CRACKLE_II_raspberry_raw_R1_merged_V1.0_", time, ".csv", sep="")
RaspberryR2_Output <- paste(OutputDir, "/CRACKLE_II_raspberry_raw_R2_merged_V1.0_", time, ".csv", sep="")
RaspberryT1_Output <- paste(OutputDir, "/CRACKLE_II_raspberry_trim_R1_merged_V1.0_", time, ".csv", sep="")
RaspberryT2_Output <- paste(OutputDir, "/CRACKLE_II_raspberry_trim_R2_merged_V1.0_", time, ".csv", sep="")
Assembly_Output <- paste(OutputDir, "/CRACKLE_II_assembly_merged_V1.0_", time, ".csv", sep="")
Strain_Output <- paste(OutputDir, "/CRACKLE_II_strain_merged_V1.0_", time, ".csv", sep="")
MLST_Output <- paste(OutputDir, "/CRACKLE_II_mlst_merged_V1.0_", time, ".csv", sep="")
Kleborate_Output <- paste(OutputDir, "/CRACKLE_II_kleborate_merged_V1.0_", time, ".csv", sep="")
Card_Output <- paste(OutputDir, "/CRACKLE_II_card_all_genes_merged_V1.0_", time, ".csv", sep="")
Plasmid_Output <- paste(OutputDir, "/CRACKLE_II_plasmid_all_genes_merged_V1.0_", time, ".csv", sep="")
InHouse_Output <- paste(OutputDir, "/CRACKLE_II_inhouse_all_genes_merged_V1.0_", time, ".csv", sep="")
Card_sub_Output <- paste(OutputDir, "/CRACKLE_II_card_only_detected_genes_merged_V1.0_", time, ".csv", sep="")
Plasmid_sub_Output <- paste(OutputDir, "/CRACKLE_II_plasmid_only_detected_genes_merged_V1.0_", time, ".csv", sep="")
InHouse_sub_Output <- paste(OutputDir, "/CRACKLE_II_inhouse_only_detected_genes_merged_V1.0_", time, ".csv", sep="")

write.csv(LabMergedTable, Lab_Output, row.names = FALSE)
write.csv(SoftwareMergedTable, Software_Output, row.names = FALSE)
write.csv(RaspberryR1MergedTable, RaspberryR1_Output, row.names = FALSE)
write.csv(RaspberryR2MergedTable, RaspberryR2_Output, row.names = FALSE)
write.csv(RaspberryT1MergedTable, RaspberryT1_Output, row.names = FALSE)
write.csv(RaspberryT2MergedTable, RaspberryT2_Output, row.names = FALSE)
write.csv(AssemblyMergedTable, Assembly_Output, row.names = FALSE)
write.csv(StrainMergedTable, Strain_Output, row.names = FALSE)
write.csv(MLSTMergedTable, MLST_Output, row.names = FALSE)
write.csv(KleborateMergedTable, Kleborate_Output, row.names = FALSE)
write.csv(CardMergedTable, Card_Output, row.names = FALSE)
write.csv(PlasmidMergedTable, Plasmid_Output, row.names = FALSE)
write.csv(InHouseMergedTable, InHouse_Output, row.names = FALSE)
write.csv(CardMergedTable_Filter, Card_sub_Output, row.names = FALSE)
write.csv(PlasmidMergedTable_Filter, Plasmid_sub_Output, row.names = FALSE)
write.csv(InHouseMergedTable_Filter, InHouse_sub_Output, row.names = FALSE)


###Remove temp directory
setwd(OutputDir)
unlink("temp", recursive=TRUE, force=TRUE)


