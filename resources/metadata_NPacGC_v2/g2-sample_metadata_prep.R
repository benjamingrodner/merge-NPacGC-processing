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
Gradients = "g2"
SampleType = "st-am-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/Morales_polyA lookup.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = NWGC.ID,
    Status       = NA,
    extractionID = NA,
    ID           = Investigator.ID
  ) %>%
  mutate(LookupID = "Morales_polyA") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("/Users/sachacoesel/Documents/Gradients/NPAc_metadata/G2NS_sample_metadata.csv")
head(meta1)
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G2NS_15m_standards/G2NS_15m_norm_factors.csv")
colnames(meta2)
head(meta2)
meta2$Alias2 <- gsub("G2NS.", "", meta2$sample_name)

meta1$Alias2
meta2$Alias2
meta <- merge(meta1[,-11], meta2[,c(19,25)], by = "Alias2")

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
    Gradients      = "G2",                      
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
    Volume.L       = volume,
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

###############################################################################
#g2-st-am-ns
###############################################################################

#Set dataset specifics
Gradients = "g2"
SampleType = "st-am-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/3690_Morales_60plex_depletion_lookup.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = NWGC.ID,
    Status       = NA,
    extractionID = NA,
    ID           = Investigator.Sample.Id
  ) %>%
  mutate(LookupID = "3690_Morales_60plex_depletion") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("/Users/sachacoesel/Documents/Gradients/NPAc_metadata/G2NS_sample_metadata.csv")
head(meta1)
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G2NS_15m_standards/G2NS_15m_norm_factors.csv")
head(meta2)
meta2$Alias2 <- gsub("G2NS.", "", meta2$sample_name)

meta1$Alias2
meta2$Alias2
colnames(meta1)
colnames(meta2)
meta <- merge(meta1[,-11], meta2[,c(19,25)], by = "Alias2")

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
    Gradients      = "G2",                      
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
    Volume.L       = volume,
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
###############################################################################
#g2-ctd-pa
###############################################################################

#Set dataset specifics
Gradients = "g2"
SampleType = "ctd-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/morales_grc_rnaseq_6_lookup.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = NWGC.Sample.ID,
    Status       = NA,
    extractionID = NA,
    ID           = Investigator.Sample.ID
  ) %>%
  mutate(LookupID = "morales_grc_rnaseq_6") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)
lookup$extractionID <- gsub("\\ ", "", lookup$ID)

#load metadata
meta1 <- read.csv("../../gradients2/g2_dcm_rr_ns_metat/sample_metadata.csv") #load ns data because it has datetime column
head(meta1)
meta2 <- read.csv("../../Count_standards_workflow/Normalization_factors/G2PA_DCM_norm_factors.csv")
head(meta2)
meta1$Alias1 <- gsub("\\ ", "", meta1$Alias1)

colnames(meta1)
colnames(meta2)

meta <- merge(meta2[,c(3,4,21)],meta1, by.x = "Sample_ID", by.y = "Alias1" )
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
    Longitude = Longitude * -1,
    dataset = "DCM"
  )

#get meta in right format:
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample,
    extractionID   = Sample_ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".", dataset, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G2",                      
    Type           = "CTD",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth,
    Filter.um      = Filter,
    Volume.L       = volume,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = NA
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup[,-5], Meta, by = "extractionID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "BD"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

df <- df[,c(5,2:4, 1, 6:24)]

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g2-ctd-ns
###############################################################################

#Set dataset specifics
Gradients = "g2"
SampleType = "ctd-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/G2NS.DCM_and_resource_ratios-lookup_armbrust_grc_rnaseq_2.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = Sample.Name,
    Status       = Notes,
    extractionID = Investigator.Sample.Name,
    ID           = NA
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_2") %>%
  filter(Status != "Seq failed") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)
lookup$extractionID <- gsub("\\ ", "", lookup$extractionID)

#load metadata
meta1 <- read.csv("../../gradients2/g2_dcm_rr_ns_metat/sample_metadata.csv") #load ns data because it has datetime column
head(meta1)
meta2 <- read.csv("../../Count_standards_workflow/Normalization_factors/G2PA_DCM_norm_factors.csv")
head(meta2)
meta1$Alias1 <- gsub("\\ ", "", meta1$Alias1)

colnames(meta1)
colnames(meta2)

meta <- merge(meta2[,c(3,4,21)],meta1, by.x = "Sample_ID", by.y = "Alias1" )
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
    Longitude = Longitude * -1,
    dataset = "DCM"
  )

#get meta in right format:
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample,
    extractionID   = Sample_ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".", dataset, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G2",                      
    Type           = "CTD",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth,
    Filter.um      = Filter,
    Volume.L       = volume,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = NA
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup[,-5], Meta, by = "extractionID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "BD"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

df <- df[,c(5,2:4, 1, 6:24)]

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g2-inc-pa
###############################################################################

#Set dataset specifics
Gradients = "g2"
SampleType = "inc-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/morales_grc_rnaseq_6_lookup.csv")
head(lookup)
lookup <- lookup %>%
  transmute(
    sequencingID = NWGC.Sample.ID,
    Status       = NA,
    extractionID = NA,
    ID           = Investigator.Sample.ID
  ) %>%
  mutate(LookupID = "morales_grc_rnaseq_6") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)
lookup <- lookup %>%
  mutate(extractionID = gsub(" ", "", ID)) %>%       # remove spaces from ID
  filter(!grepl("^BD", ID))     

#load metadata
meta1 <- read.csv("/Users/sachacoesel/Documents/Gradients/Gradients2/g2_dcm_rr_pa_metat/sample_metadata_Sept2025.csv") #updated with Angie/Randie/Bryn input on sampling datetime
head(meta1)
meta1 <- meta1 %>%
  filter(!grepl("^BD", Alias1)) #mta1 is missing some data! 3 lines.
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G2PA_RR_DCM_standards/G2NS_DCM_RR_sampling_volumes.csv")

head(meta2)
meta1$SampleID <- gsub("\\ ", "", meta1$SampleID)

colnames(meta1)
colnames(meta2)

meta <- merge(meta1,meta2, by.x = "SampleID", by.y = "Alias1" )
head(meta$Datetime)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse original local time (e.g. HST)
    Datetime.local = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "HST"),
    
    # Convert to UTC
    Datetime.utc = with_tz(Datetime.local, "UTC"),
    
    # Other derived columns
    MonthDay  = format(Datetime.local, "%m-%d"),
    Longitude = Longitude,
    dataset   = "inc"
  )

#get meta in right format:
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Alias2,
    extractionID   = SampleID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".", EXP, ".", time, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G2",                      
    Type           = "RRexp",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = dataset,
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
    Treatment      = EXP,                        
    Incubation.time= time,                        
    Notes          = Notes
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup[,-5], Meta, by = "extractionID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(!startsWith(extractionID, "BD"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

df <- df[,c(5,2:4, 1, 6:24)]

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g2-inc-ns
###############################################################################

#Set dataset specifics
Gradients = "g2"
SampleType = "inc-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/G2NS.DCM_and_resource_ratios-lookup_armbrust_grc_rnaseq_2.csv")
head(lookup)
lookup <- lookup %>%
  filter(!Notes == "Seq failed") %>%
  transmute(
    sequencingID = Sample.Name,
    Status       = Notes,
    extractionID = Investigator.Sample.Name,
    ID           = Investigator.Sample.Name
  ) %>%
  mutate(LookupID = "G2NS.DCM_and_resource_ratios-lookup_armbrust_grc_rnaseq_2") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)
lookup <- lookup %>%
  mutate(extractionID = gsub(" ", "", ID)) %>%       # remove spaces from ID
  filter(!grepl("^BD", ID))     

#load metadata
meta1 <- read.csv("/Users/sachacoesel/Documents/Gradients/Gradients2/g2_dcm_rr_pa_metat/sample_metadata_Sept2025.csv") #updated with Angie/Randie/Bryn input on sampling datetime
head(meta1)
meta1 <- meta1 %>%
  filter(!grepl("^BD", Alias1)) #mta1 is missing some data! 3 lines.
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G2PA_RR_DCM_standards/G2NS_DCM_RR_sampling_volumes.csv")

head(meta2)
meta1$SampleID <- gsub("\\ ", "", meta1$SampleID)

colnames(meta1)
colnames(meta2)

meta <- merge(meta1,meta2, by.x = "SampleID", by.y = "Alias1" )
head(meta$Datetime)

#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse original local time (e.g. HST)
    Datetime.local = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "HST"),
    
    # Convert to UTC
    Datetime.utc = with_tz(Datetime.local, "UTC"),
    
    # Other derived columns
    MonthDay  = format(Datetime.local, "%m-%d"),
    Longitude = Longitude,
    dataset   = "inc",
    Alias2 = gsub("\\.PA\\.", ".ns.", Alias2)
  )

#get meta in right format:
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Alias2,
    extractionID   = SampleID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".", EXP, ".", time, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G2",                      
    Type           = "RRexp",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = dataset,
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
    Treatment      = EXP,                        
    Incubation.time= time,                        
    Notes          = Notes
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup[,-5], Meta, by = "extractionID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(!startsWith(extractionID, "BD"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

df <- df[,c(5,2:4, 1, 6:24)]

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)




