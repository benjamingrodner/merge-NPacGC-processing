#aim, prepare standardized metadata files for RNASpades assemblies
#format: ID	LookupID	sequencingID	Status	extractionID	SampleID	Cruise	Gradients	Type	project	Dataset	Station	Cast	Replicate	Datetime.local	Datetime.utc	Latitude.dec	Longitude.dec	Depth.m	Filter.um	Volume.L	Treatment	Incubation.time	Notes

library(dplyr)
library(lubridate)

setwd("/Users/sachacoesel/Documents/NPac_v2/NPAc_v2_metadata/sample_metadata")

#load naming convention
convention <- read.csv("../NPac_v2_sample_nomenclature.csv")

###############################################################################
#g5-uw-am-pa
###############################################################################

#Set dataset specifics
Gradients = "g5"
SampleType = "uw-am-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_20.csv") %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  filter(Status != "Failed; Low mass") %>%
  mutate(LookupID = "armbrust_grc_rnaseq_20") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta <- read.csv("../sample_info_google_drive/g5.uw.am.google.drive.csv")
head(meta$Datetime..UTC.)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2025"),
    # Local time (system tz)
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system tz
    ),
    # UTC time (parse m/d/yy HH:MM explicitly)
    Datetime.utc = as.POSIXct(
      Datetime..UTC.,
      format = "%m/%d/%y %H:%M",
      tz = "UTC"
    )
  )

#get meta in right format
Meta <- meta %>%
  transmute(
    ID             = Sample.Id,
    SampleID       = paste0(Gradients, "-", SampleType, ".", sub("G5.UW-AM\\.", "", Sample.Name)),
    Cruise         = "TN412",                   
    Gradients      = Gradients,                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = NA,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
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
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g5-uw-am-ns
###############################################################################
#Set dataset specifics
Gradients = "g5"
SampleType = "uw-am-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_21.csv") %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  filter(Status != "Failed; Low mass") %>%
  mutate(LookupID = "armbrust_grc_rnaseq_21") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#get meta in right format
Meta <- meta %>%
  transmute(
    ID             = Sample.Id,
    SampleID       = paste0(Gradients, "-", SampleType, ".", sub("G5.UW-AM\\.", "", Sample.Name)),
    Cruise         = "TN412",                   
    Gradients      = Gradients,                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = NA,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
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
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g5-ctd-pa
###############################################################################
#Set dataset specifics
Gradients = "g5"
SampleType = "ctd-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_22.csv")
head(lookup)
lookup <- lookup %>%
  rename(
    sequencingID = 1,
    Status       = 2,
    extractionID = 3,
    ID           = 4
  ) %>%
  filter(Status != "Failed; Low mass", ID != "") %>%
  mutate(LookupID = "armbrust_grc_rnaseq_22") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta <- read.csv("../sample_info_google_drive/g5.ctd.google.drive.csv")
head(meta)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2025"),
    # Local time (system tz)
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system tz
    ),
    # UTC time (parse m/d/yy HH:MM explicitly)
    Datetime.utc = as.POSIXct(
      Datetime..UTC.,
      format = "%m/%d/%y %H:%M",
      tz = "UTC"
    )
  )

#get meta in right format : type.station.cast.depth.rep
Meta <- meta %>%
  transmute(
    ID             = Sample.Id,
    SampleID       = paste0(Gradients, "-", SampleType, ".", sub("G5.CTD\\.", "", Sample.Name)),
    Cruise         = "TN412",                   
    Gradients      = Gradients,                      
    Type           = "CTD",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
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
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g5-inc-PA
###############################################################################
#Set dataset specifics
Gradients = "g5"
SampleType = "inc-pa"

#lookup is "armbrust_grc_rnaseq_22"

#load metadata
meta <- read.csv("../sample_info_google_drive/g5.inc.google.drive.csv")
head(meta)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2025"),
    # Local time (system tz)
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system tz
    ),
    # UTC time (parse m/d/yy HH:MM explicitly)
    Datetime.utc = as.POSIXct(
      Datetime..UTC.,
      format = "%m/%d/%y %H:%M",
      tz = "UTC"
    )
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
  filter(startsWith(extractionID, "INC")) %>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g5-inc-ns
###############################################################################
#Set dataset specifics
Gradients = "g5"
SampleType = "inc-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_15.csv")
head(lookup)
lookup <- lookup %>%
  rename(
    sequencingID = 1,
    Status       = 4,
    extractionID = 2,
    ID           = 3
  ) %>%
  filter(Status != "Failed; Low mass", ID != "N/A") %>%
  mutate(LookupID = "armbrust_grc_rnaseq_15") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta <- read.csv("../sample_info_google_drive/g5.inc.google.drive.csv")
head(meta)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    Date_full = paste(Date..Local., "2025"),
    # Local time (system tz)
    Datetime.local = as.POSIXct(
      paste(Date_full, Time..Local.),
      format = "%d-%b %Y %H:%M",
      tz = ""   # system tz
    ),
    # UTC time (parse m/d/yy HH:MM explicitly)
    Datetime.utc = as.POSIXct(
      Datetime..UTC.,
      format = "%m/%d/%y %H:%M",
      tz = "UTC"
    )
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

nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)
