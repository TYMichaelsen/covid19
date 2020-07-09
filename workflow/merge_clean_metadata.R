#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(stringr)

swapDKlet <- function(x){stringr::str_replace_all(x,c("ø" = "oe","Ø" = "Oe","å" = "aa","Å" = "Aa","æ" = "ae","Æ" = "Ae"))}

option_list = list(
  make_option(c("-l", "--local_meta"), type="character", default=NULL, 
              help="path to local metafile", metavar="character"),
  make_option(c("-g", "--global_meta"), type="character", default=NULL, 
              help="path to gisaid metafile", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".", 
              help="path to output directory", metavar="character")
) 

opt <- parse_args(OptionParser(option_list=option_list))

if(FALSE){ 
# Specify data directories for interactive testing of script
opt <- list("g" = "/srv/rbd/covid19/global_data/metadata_2020-07-06_22-15.tsv",
            "l" = "/srv/rbd/covid19/genomes/2020-07-01-18-35_export/metadata.tsv",
            "o" = ".")
}

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

## Extract data from folder name
local_date <- as.Date(strsplit(strsplit(opt$l , "\\/")[[1]][6], "_")[[1]][1])

## Rename local variables to match gisaid data
local_meta <- rename(local_meta,
                     strain           = gisaid_id,
                     date             = SampleDate,
                     location         = zipcode_name,
                     division         = Region,
                     age              = SampleAge,
                     country_exposure = CountryOfTravel,
                     sex              = Sex) %>%
  mutate(
    virus   = "ncov",
    region  = "Europe",
    country = "Denmark",
    division = swapDKlet(division),
    sex     = ifelse(sex == "M","Male","Female"),
    date_submitted = local_date)

############################
### Combine & Clean data ###
############################

## Columns to keep in gisaid data for nextstrain
ns_cols <- c("strain",	"virus",	"gisaid_epi_isl",	"genbank_accession",	"date",	"country",
             "region_exposure",	"country_exposure",	"division_exposure", "segment",	"length",
             "host",	"age",	"sex",	"originating_lab",	"submitting_lab",	"authors",
             "url",	"title", "paper_url",	"date_submitted",	"region",	"division",	"location")

## Combine & Remove samples with missing date
comb_meta_full <- bind_rows(gisaid_meta[,ns_cols], local_meta) %>% 
  filter(is.na(date) == F)

## Clean Sex
comb_meta_full$sex <- recode(comb_meta_full$sex,
                             Woman = "Female",
                             FEmale = "Female")
comb_meta_full$sex <- if_else(comb_meta_full$sex %in% c("Male", "Female"), comb_meta_full$sex, NA_character_)

## Clean Host
comb_meta_full$host <- recode(comb_meta_full$host,
                              Female = "Human",
                              human = "Human",
                              unknown = NA_character_)

## Remove samples with very old dates
comb_meta_full <- comb_meta_full[comb_meta_full$date > as.Date("2010-01-01"),] 

###################################
### Extract NextStrain metadata ###
###################################

## Columns to select for SSI
ssi_cols <- c("Address","SymptomsStartDate","Symptoms","Travel","PlaceOfInfection_EN","ContactWithCase","Occupation",
              "ReportAgeGrp","CodR_Death30Days","CPR_Death30Days","DateOfDeath_final",
              "Death60Days_final","Death30Days_final","CitizenshipText","RegHospital","Doctor",
              "Nurse","HealthAssist","Plejehjemsnavn","branche1","branche2","branche3")

comb_meta_nextstrain <- select(comb_meta_full, any_of(c(ns_cols, ssi_cols)))

#################
### Dump data ###
#################

write_tsv(comb_meta_full, paste0(opt$o ,"/metadata_full_", local_date,".tsv"))
write_tsv(comb_meta_nextstrain, paste0(opt$o ,"/metadata_nextstrain_", local_date,".tsv"))

