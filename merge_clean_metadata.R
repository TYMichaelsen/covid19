#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)


option_list = list(
  make_option(c("-l", "--local_meta"), type="character", default=NULL, 
              help="path to local metafile", metavar="character"),
  make_option(c("-g", "--global_meta"), type="character", default=NULL, 
              help="path to gisaid metafile", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="path to output directory", metavar="character")
) 

if(FALSE){ 
# Specify data directories for interactive testing of script
opt <- list("g" = "/srv/rbd/covid19/global_data/metadata_2020-07-06_22-15.tsv",
            "l" = "/srv/rbd/covid19/genomes/2020-07-01-18-35_export/metadata.tsv",
            "o" = ".")
}

opt <- parse_args(OptionParser(option_list=option_list))

###################
### GISAID data ###
###################

## Load Global data
gisaid_meta  <- read_tsv(opt$g, na = c("?", ""))
gisaid_meta$date <- as.Date(gisaid_meta$date) 

## Remove Danish samples from global data
gisaid_meta  <- gisaid_meta[grep("Denmark", gisaid_meta$strain,invert = T), ]  

##################
### Local data ###
##################

# Load local data
local_meta <- read_tsv(opt$l)
local_meta$location <- local_meta$zipcode_name

## Extract data from folder name
local_date <- as.Date(strsplit(strsplit(opt$l , "\\/")[[1]][6], "_")[[1]][1])
local_meta$date_submission <- local_date

####################
### Combine data ###
####################

## Columns needed for nextStrain
ns_cols <- c("strain",	"virus",	"gisaid_epi_isl",	"genbank_accession",	"date",	"country",
             "region_exposure",	"country_exposure",	"division_exposure", "segment",	"length",
             "host",	"age",	"sex",	"originating_lab",	"submitting_lab",	"authors",
             "url",	"title", "paper_url",	"date_submitted",	"region",	"division",	"location")

## Columns to select for SSI
ssi_cols <- c("SampleDate","Sex","Address","ZipCodeCity","Region","SymptomsStartDate",
              "Symptoms","Travel","PlaceOfInfection_EN","ContactWithCase","Occupation",
              "ReportAge","ReportAgeGrp","CodR_Death30Days","CPR_Death30Days","DateOfDeath_final",
              "Death60Days_final","Death30Days_final","CitizenshipText","RegHospital","Doctor",
              "Nurse","HealthAssist","Plejehjemsnavn","branche1","branche2","branche3")

comb_meta_full <- bind_rows(gisaid_meta[,ns_cols] , local_meta) %>% 
  filter(is.na(date) == F)

#################
### Dump data ###
#################

write_tsv(comb_meta_full, paste0(opt$o ,"/metadata_full_", local_date,".tsv"))


