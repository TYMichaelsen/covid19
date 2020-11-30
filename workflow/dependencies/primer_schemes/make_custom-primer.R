# This R-script formats the scheme definition file into a format that can be used for nextclade.

library(tidyverse)

# Set wd according to location of script.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### aau_long_v3.1 ###

aau_long_v3.1 <- read_tsv("nCoV-2019/aau_long_v3.1/nCoV-2019.tsv") %>%
  mutate('Country (Institute)' = "AAU (DK)",
         Target = "N",
         Oligonucleotide = name,
         Sequence = seq) %>% 
  select('Country (Institute)', Target, Oligonucleotide, Sequence)

write_csv(aau_long_v3.1,file = "nCoV-2019/aau_long_v3.1/custom_primer2.csv")

### v3 ###

v3 <- read_tsv("nCoV-2019/v3/nCoV-2019.tsv") %>%
  mutate('Country (Institute)' = "UK (Artic Network)",
         Target = "N",
         Oligonucleotide = name,
         Sequence = seq) %>% 
  select('Country (Institute)', Target, Oligonucleotide, Sequence)

write_csv(v3,file = "nCoV-2019/v3/custom_primer.csv")