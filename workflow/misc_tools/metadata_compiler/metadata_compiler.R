# AAU COVID-19 metadata compiler

cat("\n### COVID-19 metadata compiler ###\n \n")

# R-packages
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(readxl)
library(utils)

# Define helper functions
step_stamp <- function(step){
  paste(
    "\n[",
    format(Sys.time(),"%H:%M:%S"),
    "] ",
    step,
    "\n\n",
    sep = ""
  )
}

check_extractionid <- function(df){
  cd <- row.names(df)[is.na(df$ExtractionID) | df$ExtractionID == ""]
  if (length(cd) > 0){
    cat("Empty 'ExtractionID' found in following data rows:\n")
    print(cd)
    cat("Curate and rerun. Exiting....")
    Sys.sleep(10)
    q()
  }
}

# Define journal settings

cat(step_stamp("Importing registration metadata"))

CJ <- readline(prompt = "Enter journal ID (CJXXX): ") %>%
  ifelse(
    grepl("CJ[0-9]{3}", .),
    .,
    stop(
      paste("'", ., "' is not correct journal ID format. Should be 'CJXXX'", sep =""),
      call. = F)
  )


COV <- readline(prompt = "Enter project ID (COVXXX): ") %>%
  ifelse(
    grepl("COV[0-9]{3}", .),
    .,
    stop(
      paste("'", ., "' is not correct project ID format. Should be 'COVXXX'", sep =""),
      call. = F)
  )

LibraryMethod <- menu(
  c("artic_std", "aau_long"),
  title = "Choose library method:"
) %>%
  c("artic_std", "aau_long")[.]

Primer <- menu(
  c("V3", "V3.1"),
  title = "Choose primer pool version:"
) %>%
  c("V3", "V3.1")[.]

# Import registration files
cat(step_stamp("Import registration metadata"))
      
reg_files <- {
  cat("Choose registration metadata file(s):\n")
  choose.files()
}
print(reg_files)

reg <- 
  lapply(
    reg_files,
    function(file){
      read_csv(
        file,
        col_types = cols(.default = "c")
      )
    }
  ) %>%
  bind_rows() %>%
  select(
    -AmpliconConc,       
    -LibraryPoolAmount,
    -LibraryMethod,
    -Primer,
    -Barcode,
    -SampleContent
  ) %>%
  mutate(ExtractionID = sub("-A$", "", ExtractionID))

check_extractionid(reg)

# Import journal files
cat(step_stamp("Import journal metadata"))

sample <- {
  cat("Choose journal sample overview file:\n")
  sample_file <- file.choose()
  print(sample_file)
  sample_file %>%
    read_excel() %>%
    select(ExtractionID, SampleContent)  %>%
    filter(!grepl("empty*", tolower(ExtractionID))) %>%
    mutate(ExtractionID = sub("-A$", "", ExtractionID)) %>%
    mutate(
      sample_group =  case_when(
        grepl(".*-NEG[0-9]$", ExtractionID, perl = T) ~ "NC",
        grepl("EXT-CJ002-4", ExtractionID, perl = T) ~ "PC",
        TRUE ~ "RS"
      )
    ) %>%
    group_by(sample_group) %>%
    mutate(sample_group_nr = 1:n())
}
check_extractionid(sample)

qubit <- {
  cat("Choose journal qubit and barcode file:\n")
  qubit_file <- file.choose()
  print(qubit_file)
  qubit_file %>%
    read_excel() %>%
    select(ExtractionID, AmpliconConc = matches(".*[cC]onc.*"), Barcode) %>%
    filter(!grepl("empty*", tolower(ExtractionID))) %>%
    mutate(ExtractionID = sub("-A$", "", ExtractionID))  %>%
    mutate(
      sample_group =  case_when(
        grepl(".*-NEG[0-9]$", ExtractionID, perl = T) ~ "NC",
        grepl("EXT-CJ002-4", ExtractionID, perl = T) ~ "PC",
        TRUE ~ "RS"
      )
    ) %>%
    group_by(sample_group) %>%
    mutate(sample_group_nr = 1:n())
}
check_extractionid(qubit)

norm <- {
  cat("Choose journal normalization file:\n")
  norm_file <- file.choose()
  print(norm_file)
  norm_file  %>%
    read_excel()%>%
    select(ExtractionID, LibraryPoolAmount = `Pool amount (ng)`)   %>%
    filter(!grepl("empty*", tolower(ExtractionID))) %>%
    mutate(ExtractionID = sub("-A$", "", ExtractionID))  %>%
    mutate(
      sample_group =  case_when(
        grepl(".*-NEG[0-9]$", ExtractionID, perl = T) ~ "NC",
        grepl("EXT-CJ002-4", ExtractionID, perl = T) ~ "PC",
        TRUE ~ "RS"
      )
    ) %>%
    group_by(sample_group) %>%
    mutate(sample_group_nr = 1:n())
}
check_extractionid(norm)

# Compile draft metadata sheet
cat(step_stamp("Merging metadata"))

ms <- 
  sample %>%
  left_join(reg, by = c("ExtractionID")) %>%
  left_join(qubit, by = c("ExtractionID", "sample_group", "sample_group_nr")) %>%
  left_join(norm, by = c("ExtractionID", "sample_group", "sample_group_nr"))


# Fill missing data for controls and generate libraru ids
cat(step_stamp("Adding metadata for controls and generating LibraryID's"))

msf <- 
  ms %>%
  # define sample groups for counting purposes
  # Fill out data
  mutate(
    ProjectName = if_else(
      grepl(".*-NEG[0-9]$|EXT-CJ002-4", ExtractionID, perl = T),
      COV,
      ProjectName 
    ),
    InternalContact = if_else(
      grepl(".*-NEG[0-9]$", ExtractionID, perl = T),
      "SMK",
      InternalContact
    ),
    ExternalContact = if_else(
      grepl(".*-NEG[0-9]$", ExtractionID, perl = T),
      "AAU",
      ExternalContact
    ),
    LibraryID = case_when(
      grepl(".*-NEG[0-9]$", ExtractionID, perl = T) ~ sub("^EXT-", "LIB-", ExtractionID),
      grepl("EXT-CJ002-4", ExtractionID, perl = T) ~ paste("LIB-", CJ, "-POS", 1:n(), sep =""),
      TRUE ~ paste("LIB-", CJ, "-", 1:n(), sep = "")
    ),
    SampleName = case_when(
      grepl(".*-NEG[0-9]$", ExtractionID, perl = T) ~ LibraryID,
      grepl("EXT-CJ002-4", ExtractionID, perl = T) ~ LibraryID,
      TRUE ~ SampleName
    ),
    SampleSite = if_else(
      grepl(".*-NEG[0-9]$", ExtractionID, perl = T),
      "AAU",
      SampleSite
    )
  ) %>%
  mutate(
    LibraryMethod = LibraryMethod,
    Primer = Primer
  )


# Trim dataframe
mo <- msf %>%
  ungroup() %>%
  select(
    ProjectName,
    InternalContact,
    ExternalContact,
    SampleName,
    SampleID,
    ExtractionID,
    LibraryID,
    SampleContent,
    SampleSite,
    RecieveDate,
    qPCR_Ct,
    AmpliconConc,
    LibraryPoolAmount,
    SamplingMethod,
    ExtractionMethod,
    LibraryMethod,
    Primer,
    Barcode,
    SampleStorageID,
    ExtractionStorageID,
    LibraryStorageID,
    Comments
  )

# Run checks
cat(step_stamp("Checking for missing values"))

check1_list <-
  c(
    "ProjectName",
    "InternalContact",
    "ExternalContact",
    "SampleName",
    "ExtractionID",
    "LibraryID",
    "SampleContent",
    "SampleSite",
    "AmpliconConc",
    "LibraryPoolAmount",
    "LibraryMethod",
    "Primer",
    "Barcode"
  )

mo_check1 <- 
  mo %>%
  select(check1_list) %>%
  group_by(ExtractionID) %>%
  filter_all(any_vars(is.na(.))) %>%
  select_if(function(x) any(is.na(x)))

if (nrow(mo_check1) > 0){
  cat("Following ExtractionID's have missing values in one or more columns:\n")
  print(mo_check1)
} else {
  cat("No missing values found.\n")
}
  

# Write output
out_dir <- paste(dirname(reg_files[1]), "/output", sep = "")
dir.create(out_dir)
cat(step_stamp(paste("Writing output files to ", out_dir, sep = "")))

write_csv(
  mo,
  paste(out_dir, "/", CJ, "_metadata_sequencing.csv", sep="")
)

write_csv(
  mo %>% select(LibraryID, Barcode, AmpliconConc),
  paste(out_dir, "/sample_sheet.csv", sep=""),
  col_names = F
)

cat(step_stamp("Done..."))
q()