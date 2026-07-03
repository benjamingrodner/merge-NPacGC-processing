#aim, prepare standardized metadata files for RNASpades assemblies
#format: ID	LookupID	sequencingID	Status	extractionID	SampleID	Cruise	Gradients	Type	project	Dataset	Station	Cast	Replicate	Datetime.local	Datetime.utc	Latitude.dec	Longitude.dec	Depth.m	Filter.um	Volume.L	Treatment	Incubation.time	Notes

library(dplyr)
library(lubridate)

setwd("/Users/sachacoesel/Documents/NPac_v2/NPAc_v2_metadata/sample_metadata")

#load naming convention
convention <- read.csv("../NPac_v2_sample_nomenclature.csv")

###############################################################################
#g2-st-am-pa
###############################################################################

#Set dataset specifics
Gradients = "g1"
SampleType = "st-am-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/morales_gradient_lookup.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = Sample.ids,
    Status       = NA,
    extractionID = NA,
    ID           = Investigator.Sample.Name
  ) %>%
  mutate(LookupID = "morales_gradient") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("../../gradients1/g1_station_pa_metat/sample_metadata.csv")
head(meta1)
meta1$Alias2 <- gsub("_0_2", "_0.2", meta1$Alias2)

meta2 <- read.csv("../../Count_standards_workflow/Normalization_factors/G1PA_SUM_norm_factors.csv")
colnames(meta2)
colnames(meta1)
head(meta2)
meta2$Alias2 <- gsub("_02", "_0.2", meta2$sample_name_short)


meta <- merge(meta1[,-11], meta2[,c(9,4)], by = "Alias2")

head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.utc = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "UTC"), 
    # Convert to local (HST)
    Datetime.local   = with_tz(Datetime.utc, "HST"),
    MonthDay = format(Datetime.local, "%m-%d"),
    Longitude = -1* Longitude,
    Cast = 1
  )

#get meta in right format: month.day.filter.rep
Meta <- meta %>%
  transmute(
    ID             = Alias1,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station,".C", Cast, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G1",                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = "AM",
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth,
    Filter.um      = Filter,
    Volume.L       = volume_filtered,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )

#merge and save metadata file
head(lookup)
colnames(Meta)

df <- merge(lookup, Meta, by.x = "ID", by.y = "ID")

#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g1-st-am-ns
###############################################################################

#Set dataset specifics
Gradients = "g1"
SampleType = "st-am-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/2256_Morales_barcode_lookup.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = Sample.Name,
    Status       = NA,
    extractionID = NA,
    ID           = Investigator.Sample.Name
  ) %>%
  mutate(LookupID = "2256_Morales") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
#load metadata
meta1 <- read.csv("../../gradients1/g1_station_ns_metat/sample_metadata.csv")
head(meta1)
meta2 <- read.csv("../../Count_standards_workflow/Normalization_factors/G1NS_SUM_norm_factors.csv")
colnames(meta2)
colnames(meta1)
head(meta2)
#meta2$Alias2 <- gsub("_02", "_0.2", meta2$sample_name_short)



meta1$Alias2
meta2$sample_name_short
colnames(meta1)
colnames(meta2)
meta <- merge(meta1[,-11], meta2[,c(2,4)], by.x = "Alias2", by.y = "sample_name_short")

head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.utc = as.POSIXct(Datetime,
                              format = "%m/%d/%y %H:%M",
                              tz = "UTC"),  
    # Convert to local (HST)
    Datetime.local   = with_tz(Datetime.utc, "HST"),
    MonthDay = format(Datetime.local, "%m-%d")
  )
colnames(meta)
#get meta in right format: month.day.filter.rep
Meta <- meta %>%
  transmute(
    ID             = Alias1,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station,".C", Cast, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G1",                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = "AM",
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth,
    Filter.um      = Filter,
    Volume.L       = volume_filtered,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)

df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)
