
library(httr)

# Swap DK wierd letters function.
swapDKlet <- function(x){stringr::str_replace_all(x,c("ø" = "oe","Ø" = "Oe","å" = "aa","Å" = "Aa","æ" = "ae","Æ" = "Ae"))}


path <- "https://dawa.aws.dk/postnumre?landpostnumre"

request <- GET(url = path)

response <- content(request, as = "text", encoding = "UTF-8")

kommune <- fromJSON(response, flatten = TRUE) %>% 
  data.frame() %>%
  separate(visueltcenter, c("long", "lat"), sep = ",") %>%
  mutate(long = gsub('[c()]', '',long) %>% as.numeric()) %>%
  mutate(lat  = gsub('[c()]', '',lat) %>% as.numeric()) %>%
  mutate(muncipal = swapDKlet(navn)) %>%
  select(nr,muncipal,navn,long,lat)

# Average by muncipal. 
latlong_muncipal <- group_by(kommune,muncipal) %>%
  summarise(lat = mean(lat),long = mean(long)) %>%
  mutate(type = "location",loc = muncipal) %>%
  select(type,loc,lat,long)

# Add regions.
latlong_regions <- rbind(
  c("Syddanmark",55.378456,9.131928),
  c("Nordjylland",57.047218,9.920100),
  c("Midtjylland",56.453121,9.402010),
  c("Sjaelland",55.503750,11.836965),
  c("Hovedstaden",55.861710,12.313743)) %>%
  data.frame() %>%
  `colnames<-`(c("loc","lat","long")) %>%
  mutate(type = "division") %>%
  select(type,loc,lat,long)

# Add country.
latlong_country <- data.frame(type = "country",loc = "Denmark",lat = 56.2639,long = 9.5018)

latlong_nxt <- rbind(latlong_muncipal,latlong_regions,latlong_country)

fwrite(latlong_nxt,"/srv/rbd/covid19/git/covid19/workflow/dependencies/nextstrain/latlong_nextstrain.tsv",col.names = F,sep = "\t")
