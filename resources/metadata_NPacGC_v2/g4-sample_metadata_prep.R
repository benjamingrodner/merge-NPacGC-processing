#aim, prepare standardized metadata files for RNASpades assemblies
#format: ID	LookupID	sequencingID	Status	extractionID	SampleID	Cruise	Gradients	Type	project	Dataset	Station	Cast	Replicate	Datetime.local	Datetime.utc	Latitude.dec	Longitude.dec	Depth.m	Filter.um	Volume.L	Treatment	Incubation.time	Notes

library(dplyr)
setwd("/Users/sachacoesel/Documents/NPac_v2/NPac_v2_metadata/")

#load naming convention
convention <- read.csv("NPac_v2_sample_nomenclature.csv")

###############################################################################
#g4-uw-am-pa
###############################################################################
#Set dataset specifics
Gradients = "g4"
SampleType = "uw-am-pa"

#load lookup table
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_11.csv")
head(lookup)
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_11.csv") %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_11") %>%
  select(LookupID, sequencingID, Status, extractionID, ID) %>%
  filter(!(is.na(ID) | ID == "" | Status == "Failed Low volume"))
head(lookup)

meta <- read.csv("sample_info_google_drive/g4.uw.am.google.drive.csv")
head(meta)

#prepare Datetime.local column
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2021"),
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system's local tz, prevents unwanted UTC tag
    )|> format("%Y-%m-%dT%H:%M:%S")
  )


#get meta in right format
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample.Id,
    SampleID       = paste0(Gradients, "-", SampleType, ".", sub("G4.UW-AM\\.", "", Sample.Name)),
    Cruise         = "TN397",                   
    Gradients      = Gradients,                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = "UW",                   
    Cast           = NA,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime..UTC.,            
    Latitude.dec   = Latitude..Decimal.,
    Longitude.dec  = Longitude..Decimal.,
    Depth.m        = Depth..m.,
    Filter.um      = Filter.Pore.Size..um.,
    Volume.L       = Total.Filter.Volume..L.,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )

#merge and save metadata file
df <- merge(lookup, Meta, by = "ID") 

#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0("sample_metadata/", Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g4-uw-am-ns
###############################################################################

#Set dataset specifics
Gradients = "g4"
SampleType = "uw-am-ns"

#load lookup table
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_18.csv")
head(lookup)
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_18.csv") %>%
  rename_with(~ c("sequencingID", "extractionID", "ID"), .cols = 1:3) %>%
  mutate(Status = NA) %>%   # add a proper Status column
  mutate(LookupID = "armbrust_grc_rnaseq_18") %>%
  select(LookupID, sequencingID, Status, extractionID, ID) %>%
  filter(!(is.na(ID) | ID == ""))
head(lookup)

meta <- read.csv("sample_info_google_drive/g4.uw.am.google.drive.csv")
head(meta)

#prepare Datetime.local column
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2021"),
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system's local tz, prevents unwanted UTC tag
    )|> format("%Y-%m-%dT%H:%M:%S")
  )

#get meta in right format
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample.Id,
    SampleID       = paste0(Gradients, "-", SampleType, ".", sub("G4.UW-AM\\.", "", Sample.Name)),
    Cruise         = "TN397",                   
    Gradients      = Gradients,                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = "UW",                   
    Cast           = NA,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime..UTC.,            
    Latitude.dec   = Latitude..Decimal.,
    Longitude.dec  = Longitude..Decimal.,
    Depth.m        = Depth..m.,
    Filter.um      = Filter.Pore.Size..um.,
    Volume.L       = Total.Filter.Volume..L.,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )

#merge and save metadata file
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0("sample_metadata/", Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g4-uw-pm-ns
###############################################################################

#Set dataset specifics
Gradients = "g4"
SampleType = "uw-pm-ns"

#load lookup table
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_23.csv")
head(lookup)
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_23.csv") %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_23") %>%
  select(LookupID, sequencingID, Status, extractionID, ID) %>%
  filter(!(is.na(ID) | ID == "")) %>%
  filter(Status == "Complete; Pass")
head(lookup)

meta <- read.csv("sample_info_google_drive/g4.inc.google.drive.csv")
head(meta)

#prepare Datetime.local column
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2021"),
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system's local tz, prevents unwanted UTC tag
    )|> format("%Y-%m-%dT%H:%M:%S")
  )

#get meta in right format : station.cast.treatment.time.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample.ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station,".", CAST,".", TREATMENT,".", TIME,".", Replicate),
    Cruise         = "TN412",                   
    Gradients      = Gradients,                      
    Type           = "Incubation",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth..m.,
    Filter.um      = Filter.Pore.Size..um.,
    Volume.L       = Total.Filter.Volume..L.,
    Treatment      = Treatment,                        
    Incubation.time= Incubation.Time,                        
    Notes          = Notes
  )

#merge and save metadata file
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "UW")) %>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0("sample_metadata/", Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g4-ctd-ns
###############################################################################

#Set dataset specifics
Gradients = "g4"
SampleType = "ctd-ns"

#load lookup table
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_23.csv")
head(lookup)
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_23.csv") %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_23") %>%
  select(LookupID, sequencingID, Status, extractionID, ID) %>%
  filter(!(is.na(ID) | ID == "")) %>%
  filter(Status == "Complete; Pass")
head(lookup)

meta <- read.csv("sample_info_google_drive/g4.ctd.google.drive.csv")
head(meta)

#prepare Datetime.local column
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2021"),
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system's local tz, prevents unwanted UTC tag
    )|> format("%Y-%m-%dT%H:%M:%S")
  )
colnames(meta)
#get meta in right format : type.station.cast.depth.rep
Meta <- meta %>%
  transmute(
    ID             = Sample.Id,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".", Depth..m., "m.", Replicate),
    Cruise         = "TN412",                   
    Gradients      = Gradients,                      
    Type           = "CTD",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime..UTC.,            
    Latitude.dec   = Latitude..Decimal.,
    Longitude.dec  = Longitude..Decimal.,
    Depth.m        = Depth..m.,
    Filter.um      = Filter.Pore.Size..um.,
    Volume.L       = Total.Filter.Volume..L.,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )
#merge and save metadata file
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "CTD")) %>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0("sample_metadata/", Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g4-inc-ns
###############################################################################

#Set dataset specifics
Gradients = "g4"
SampleType = "inc-ns"

#load lookup table
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_23.csv")
head(lookup)
lookup <- read.csv("lookup_armbrust_files/lookup_armbrust_grc_rnaseq_23.csv") %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_23") %>%
  select(LookupID, sequencingID, Status, extractionID, ID) %>%
  filter(!(is.na(ID) | ID == "")) %>%
  filter(Status == "Complete; Pass")
head(lookup)

meta <- read.csv("sample_info_google_drive/g4.inc.google.drive.csv")
head(meta)

#prepare Datetime.local column
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2021"),
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system's local tz, prevents unwanted UTC tag
    )|> format("%Y-%m-%dT%H:%M:%S")
  )

#get meta in right format : station.cast.treatment.time.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample.ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station,".", CAST,".", TREATMENT,".", Incubation.Time,".", Replicate),
    Cruise         = "TN412",                   
    Gradients      = Gradients,                      
    Type           = "Incubation",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime..UTC.,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth..m.,
    Filter.um      = Filter.Pore.Size..um.,
    Volume.L       = Total.Filter.Volume..L.,
    Treatment      = Treatment,                        
    Incubation.time= Incubation.Time,                       
    Notes          = Notes
  )

#merge and save metadata file
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "SDR")) %>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0("sample_metadata/", Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)
