library(tidyverse)
library(seqinr)

###################
### GISAID data ###
###################

## Load Global data
gisaid_meta  <- read_tsv("/srv/rbd/covid19/global_data/metadata_2020-06-29_12-15.tsv", na = c("?", ""))
gisaid_meta$date <- as.Date(gisaid_meta$date) 
gisaid_fasta <- read.fasta("/srv/rbd/covid19/global_data/sequences_2020-06-29_06-18.fasta")

## Remove Danish samples from global data
gisaid_meta  <- gisaid_meta[grep("Denmark", gisaid_meta$strain,invert = T), ]  
gisaid_fasta <- gisaid_fasta[grep("Denmark", names(gisaid_fasta), invert = T)]  

## Filter fasta and metadata to have matching names
gisaid_meta <- gisaid_meta[gisaid_meta$strain %in% names(gisaid_fasta), ] 
gisaid_fasta <- gisaid_fasta2[names(gisaid_fasta) %in% gisaid_meta$strain] 



##################
### Local data ###
##################
## Columns needed for nextStrain
nscols <- c("strain",	"virus",	"gisaid_epi_isl",	"genbank_accession",	"date",	"country",
            "region_exposure",	"country_exposure",	"division_exposure", "segment",	"length",
            "host",	"age",	"sex",	"originating_lab",	"submitting_lab",	"authors",
            "url",	"title", "paper_url",	"date_submitted",	"region",	"division",	"location")

# Load local data
#local_meta <- read_tsv("/srv/rbd/covid19/genomes/2020-06-29-14-47_export/metadata.tsv")
local_meta <- read_tsv("/srv/rbd/covid19/metadata/2020-06-29-14-47_metadata_nextstrain.tsv")
local_meta$location <- local_meta$zipcode_name   

local_fasta <- read.fasta("/srv/rbd/covid19/genomes/2020-06-29-14-47_export/sequences.fasta")

selected_local  <- select(local_meta, any_of(nscols))
selected_local$date_submission <- as.Date("2020-06-29") 

####################
### Combine data ###
####################

comb_meta <- bind_rows(gisaid_meta[,nscols] , selected_local) %>% 
  filter(is.na(date) == F)
comb_meta <- comb_meta[, nscols]


if(TRUE){
  ## Subsample
  sampleRows <- c("Wuhan/Hu-1/2019", 'Wuhan/WH01/2019', sample(comb_meta$strain, 100))
  comb_meta  <- comb_meta[comb_meta$strain %in% sampleRows, ]
}  

comb_fasta <- c(gisaid_fasta, local_fasta)
comb_fasta <- comb_fasta[comb_meta$strain] 

outdir <- "/srv/rbd/covid19/testing/gisaid/ncov/data/"
write_tsv(comb_meta, paste0(outdir,"metadata_", "2020-06-29",".tsv"))
write.fasta(comb_fasta, names = names(comb_fasta), file = paste0(outdir,"sequences_", "2020-06-29",".fasta"))


