##################
#USAGE - Rscript $XML_to_gene_hit_list.R [XML File] [Resistance Database] [Output File]

#For use with the CRACKLEII_WGS_SingleIsolate.sh v1.0 pipeline 

#Author: Blake Hanson (blake.hanson [at] uth.tmc.edu)
#Version:1.0
#Date Last updated:2018/04/01
#Updates:
##################

###Arguments
args<-commandArgs(TRUE)

###Libraries
library(XML)

###Functions

###Variables
XMLFile <- args[1]
Database <- args[2]
OutputFile <- args[3]

###Code
#Parse XML file 
xmlfile <- xmlParse(XMLFile)

#Read resistance gene section based upon specified Database
ifelse(Database == "CARD", resistance <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/card_resistance/*"), stringsAsFactors = FALSE),
       ifelse(Database == "PLASMIDFINDER", resistance <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/plasmidfinder_resistance/*"), stringsAsFactors = FALSE),
              ifelse(Database == "INHOUSE", resistance <- xmlToDataFrame(nodes = getNodeSet(xmlfile, "//*/inhouse_resistance/*"), stringsAsFactors = FALSE), 
                     print("Error: specified database not supported"))))

#Subset only genes with a hit
resistance_sub <- resistance[which(resistance$start != "NA"),]

#Write list of hits to csv file
write.csv(resistance_sub, OutputFile, row.names = FALSE)


