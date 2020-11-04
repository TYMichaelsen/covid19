#!/usr/bin/env Rscript

## Fix library resoultion problem
## Exclude user's library
my_home = Sys.getenv("HOME")
libvec = .libPaths()
.libPaths(libvec[!grepl(my_home, libvec)])

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
opt <- list("g" = "/srv/rbd/covid19/global_data/metadata_2020-11-04_09-37.tsv",
            "l" = "/srv/rbd/covid19/metadata/2020-10-31-22-54_metadata_nextstrain.tsv",
            "o" = ".")
}

###################
### GISAID data ###
###################

## Load Global data
gisaid_meta      <- read_tsv(opt$g, na = c("?", ""))
gisaid_meta$date <- as.Date(gisaid_meta$date) 

## Remove Danish samples from global data
gisaid_meta  <- gisaid_meta[grep("Denmark", gisaid_meta$strain,invert = T), ]

# Correct sample names.
gisaid_meta$strain <- recode(gisaid_meta$strain,
                             `Benin/197/03.2020` = "Benin/197/2020")
##################
### Local data ###
##################

# Load local data
local_meta <- read_tsv(opt$l,guess_max = 10000)

## Extract date from file name

local_date <- as.Date(str_match(opt$l, "metadata/(.*?)_metadata")[2])

## Rename local variables to match gisaid data
local_meta <- dplyr::rename(local_meta,
                     strain           = ssi_id,
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
    location = swapDKlet(location),
    sex     = ifelse(sex == "M","Male","Female"),
    date_submitted = local_date)

############################
### Combine & Clean data ###
############################

## Columns to keep in gisaid data for nextstrain
ns_cols <- c("strain",	"virus",	"date",	"country",
             "region_exposure",	"country_exposure",	"division_exposure", "length",
             "host",	"age",	"sex",	"date_submitted",	"region",	"division",	"location")

# "originating_lab",	"submitting_lab",	"authors", "url",	"title", "paper_url",
# "gisaid_epi_isl",	"genbank_accession", "segment",

## Combine & Remove samples with missing date
comb_meta_full <- bind_rows(gisaid_meta[,ns_cols], local_meta) %>% 
  filter(!is.na(date))

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

## Clean locations.
comb_meta_full$location <- recode(comb_meta_full$location,
                                  Ukendt = NA_character_)

#################
### Dump data ###
#################

write_tsv(comb_meta_full, paste0(opt$o ,"/metadata_nextstrain_DKglobal.tsv"),na = "")
