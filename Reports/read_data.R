## Misc functions.---------------------------------------------------------------
median_range <- function(x) paste0(median(x, na.rm = T), " (", min(x, na.rm = T), "-", max(x, na.rm = T), ")")
median_iqr   <- function(x) paste0(median(x, na.rm = T), " (", IQR(x, na.rm = T), ")")

## Custom colors
branded_colors <- list(
  "blue"   = "#00798c",
  "red"    = "#d1495b",
  "yellow" = "#edae49",
  "green"  = "#66a182",
  "navy"   = "#2e4057", 
  "grey"   = "#8d96a3"
)

clade_colors <- list(
  "19A"   = "#00798c",
  "19B"    = "#d1495b",
  "20A" = "#edae49",
  "20B"  = "#66a182",
  "20C"   = "#2e4057")

## Danish month names
dk_month <- factor(c("Januar", "Februar", "Marts", "April", "Maj", "Juni",
                     "Juli", "August", "September", "Oktober", "November", "December"),
                   levels = c("Januar", "Februar", "Marts", "April", "Maj", "Juni",
                              "Juli", "August", "September", "Oktober", "November", "December"))


# Select the most recent nextstrain.
timestmp_nxt <- list.files("/srv/rbd/covid19/nextstrain",pattern = "_nextstrain") %>%
  sub("_nextstrain","",x = .) %>%
  gsub("_","-",x = .) %>%
  {.[!grepl("[A-z]",x = .)]} %>%
  strptime(format = "%Y-%m-%d-%H-%M") %>%
  max(na.rm = T) %>%
  format("%Y-%m-%d-%H-%M")

# Select the most recent linelist data.
tmp <- list.files("/srv/rbd/covid19/metadata",pattern = "linelist") %>%  
  str_split(string = .,pattern = "_linelist") %>%  
  .[[1]] 

timestmp      <- tmp[1]
timestmp_ll   <- tmp[2] %>%   sub(".tsv","",x = .) %>%  as.Date()

### Load the data.--------------------------------------------------------------

# Linelist
meta_ll <- read_tsv(file = list.files("/srv/rbd/covid19/metadata",
                                      pattern = "linelist",
                                      full.names = T),
                    guess_max = 100000) %>% 
  #filter(host == "Human") %>% 
  filter(host == "Human"| is.na(host)) %>% 
  mutate(
    genome_qc    = factor(genome_qc,levels = c("HQ","MQ","Fail")),
    firstDayWeek = {floor_date(date_linelist,"week", week_start = 1)} %>% as.Date(),
    clade        = sub("/.*","",x = clade))

# Tree.
tree     <- read.tree(file = paste0("/srv/rbd/covid19/nextstrain/",timestmp_nxt,"_nextstrain/results/Denmark_Full/tree_raw.nwk"))
timetree <- read.tree(file = paste0("/srv/rbd/covid19/nextstrain/",timestmp_nxt,"_nextstrain/results/Denmark_Full/tree.nwk"))

# intersect tree and metadata.
wh1 <- match(tree$tip.label,meta_ll$ssi_id) %>%
  `names<-`(tree$tip.label) %>%
  na.omit()

wh2 <- match(timetree$tip.label,meta_ll$ssi_id) %>%
  `names<-`(timetree$tip.label) %>%
  na.omit()

wh <- intersect(names(wh1),names(wh2))

tree     <- keep.tip(tree,tip = wh)
timetree <- keep.tip(timetree,tip = wh)


# Misc data.--------------------------------------------------------------------
# For plotting denmark.
dk_nuts2 <- read_delim(file = "/srv/rbd/ma/test/maps/DK_NUTS2.txt", delim ="\t")


