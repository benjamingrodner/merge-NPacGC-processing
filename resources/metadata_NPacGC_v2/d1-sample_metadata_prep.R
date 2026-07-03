#aim, prepare standardized metadata files for RNASpades assemblies
#format: ID	LookupID	sequencingID	Status	extractionID	SampleID	Cruise	Gradients	Type	project	Dataset	Station	Cast	Replicate	Datetime.local	Datetime.utc	Latitude.dec	Longitude.dec	Depth.m	Filter.um	Volume.L	Treatment	Incubation.time	Notes

library(dplyr)
library(lubridate)

setwd("/Users/sachacoesel/Documents/NPac_v2/NPAc_v2_metadata/sample_metadata")

#load naming convention
convention <- read.csv("../NPac_v2_sample_nomenclature.csv")

###############################################################################
#d1-pa
###############################################################################

#Set dataset specifics
Gradients = "d1"
SampleType = "st-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/Lookup_diel1_from_poolA_poolB_barcode_lookup_SC.csv")
head(lookup)
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status       = 4,
    extractionID = 5,
    ID           = 1
  ) %>%
  mutate(LookupID = "Lookup_diel1_from_poolA_poolB_barcode") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta <- read.csv("/Users/sachacoesel/Documents/Gradients/NPAc_metadata/D1_diel_sample_metadata.csv")
head(meta)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse as UTC; your original string is "yy-mm-dd HH:MM"
    Datetime.utc = as.POSIXct(Datetime_GMT,
                              format = "%y-%m-%d %H:%M",
                              tz = "UTC"),
    
    # Convert to local Hawaii Standard Time (HST, UTC-10)
    Datetime.local = with_tz(Datetime.utc, "Pacific/Honolulu"),
    
    # Extract month-day for convenience
    MonthDay = format(Datetime.local, "%m-%d"),
    
    # (if you really want to reset Cast, but be careful)
    Cast = 1
  )
colnames(meta)
#get meta in right format: month.day.filter.rep
Meta <- meta %>%
  transmute(
    ID             = Alias1,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".",Alias2, "h.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "d1",                      
    Type           = "diel",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Alias2,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth,
    Filter.um      = Filter,
    Volume.L       = Volume_L,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup, Meta, by.x = "extractionID", by.y = "ID", all = T)
#colnames in same order as others
df <- df[,c(5,2:4,1,6:24)]
#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

