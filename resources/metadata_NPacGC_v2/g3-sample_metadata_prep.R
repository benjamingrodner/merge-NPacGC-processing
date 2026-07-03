#aim, prepare standardized metadata files for RNASpades assemblies
#format: ID	LookupID	sequencingID	Status	extractionID	SampleID	Cruise	Gradients	Type	project	Dataset	Station	Cast	Replicate	Datetime.local	Datetime.utc	Latitude.dec	Longitude.dec	Depth.m	Filter.um	Volume.L	Treatment	Incubation.time	Notes

library(dplyr)
library(lubridate)

setwd("/Users/sachacoesel/Documents/NPac_v2/NPAc_v2_metadata/sample_metadata")

#load naming convention
convention <- read.csv("../NPac_v2_sample_nomenclature.csv")

###############################################################################
#g3-uw-am-pa
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "uw-am-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_7.csv")
head(lookup)
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status       = 1,
    extractionID = 3,
    ID           = 4
  ) %>%
  filter(Status != "Failed library") %>%
  mutate(LookupID = "armbrust_grc_rnaseq_7") %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("/Users/sachacoesel/Documents/Gradients/Gradients3/g3_uw_pa_metat/sample_metadata.csv")
head(meta1)
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G3PA_underway_standards/G3PA_underway_dawn_norm_factors.csv")
colnames(meta2)

meta <- merge(meta1, meta2[,c(17,24)], by = "Alias2")
head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    MonthDay = format(Datetime.local, "%m-%d"),
    Station = "UW"
  )
colnames(meta)
#get meta in right format: month.day.filter.rep
Meta <- meta %>%
  transmute(
    ID             = Alias2,
    SampleID       = paste0(Gradients, "-", SampleType, ".", MonthDay, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = "AM",
    Station        = Station,                   
    Cast           = NA,                        
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

###############################################################################
#g3-uw-pm-pa
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "uw-pm-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_12.csv")
head(lookup)
lookup$Inv_Sample_ID == lookup$Alias2 #use Alias2 for merging

lookup <- lookup %>%
  rename(
    sequencingID = 1,
    extractionID = 4,
    ID           = 3
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_12", Status = NA) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("/Users/sachacoesel/Documents/Gradients/Gradients3/g3_uw_pa_pm_metat/sample_metadata.csv")
head(meta1)
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G3PA_PM_underway_standards/G3PA_underway_dusk_norm_factors.csv")
head(meta2)
colnames(meta2)

meta <- merge(meta1, meta2[,c(21:22)], by.x = "SampleID", by.y = "Sample.ID")
head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    MonthDay = format(Datetime.local, "%m-%d")
  )
colnames(meta)
#get meta in right format: month.day.filter.rep
Meta <- meta %>%
  transmute(
    ID             = Alias2,
    SampleID       = paste0(Gradients, "-", SampleType, ".", MonthDay, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = "PM",
    Station        = Station,                   
    Cast           = NA,                        
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
head(Meta)
df <- merge(lookup, Meta, by.x = "extractionID", by.y = "ID", all = T)

#colnames in same order as others
df <- df[,c(5,2:4,1,6:24)]
#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g3-uw-am-ns
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "uw-am-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_10.csv")
head(lookup)
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status = 1,
    extractionID = 3
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_10", Status = NA, ID = extractionID) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("../../gradients3/g3_uw_ns_metat/sample_metadata.csv")
head(meta1)
meta2 <- read.csv("../../Count_standards_workflow/Normalization_factors/G3NS_underway_norm_factors.csv")
head(meta2)
colnames(meta2)

meta <- merge(meta1, meta2[,c(2,17)], by.x = "SampleID", by.y = "sample_name")
head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    MonthDay = format(Datetime.local, "%m-%d"),
    Station = "UW"
  )
colnames(meta)

#get meta in right format: month.day.filter.rep
Meta <- meta %>%
  transmute(
    ID             = Alias2,
    SampleID       = paste0(Gradients, "-", SampleType, ".", MonthDay, ".", Filter, "um.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
    Type           = "transect",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = "AM",
    Station        = Station,                   
    Cast           = NA,                        
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
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
nrow(lookup) == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g3-ctd-pa
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "ctd-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_6.csv")
head(lookup[,1:4])
colnames(lookup[,1:4])
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status = 1,
    extractionID = 3,
    ID = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_6", Status = NA) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("../../gradients3/g3_depth_pa_metat/G3Depth_polyA_sample_metadata.csv")
head(meta1)
meta2 <- read.csv("../../gradients3/g3_depth_pa_metat/G3Depth_polyA_norm_factors.csv")
head(meta2)
meta2 <- merge(lookup[,4:5],meta2, by.x = "ID"  ,by.y = "Add.l.ID", all.x = T)
colnames(meta1)
colnames(meta2)


meta <- merge(meta2[,c(1,2,21)],meta1, by.x = "extractionID", by.y = "Alias1", all.y = T )
head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(DatetimeHST,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    dataset = sapply(strsplit(SampleID, "\\."), `[`, 4)
  )


#get meta in right format: month.day.filter.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".", dataset, ".", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
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
    Volume.L       = V_total..L.,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = NA
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "DP"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g3-ctd-pa
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "ctd-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_6.csv")
head(lookup[,1:4])
colnames(lookup[,1:4])
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status = 1,
    extractionID = 3,
    ID = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_6", Status = NA) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("../../gradients3/g3_depth_pa_metat/G3Depth_polyA_sample_metadata.csv")
head(meta1)
meta2 <- read.csv("../../gradients3/g3_depth_pa_metat/G3Depth_polyA_norm_factors.csv")
head(meta2)
meta2 <- merge(lookup[,4:5],meta2, by.x = "ID"  ,by.y = "Add.l.ID", all.x = T)
colnames(meta1)
colnames(meta2)


meta <- merge(meta2[,c(1,2,21)],meta1, by.x = "extractionID", by.y = "Alias1", all.y = T )
head(meta$Datetime)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(DatetimeHST,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    dataset = sapply(strsplit(SampleID, "\\."), `[`, 4)
  )

#get meta in right format: month.day.filter.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".", dataset, ".", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
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
    Volume.L       = V_total..L.,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = NA
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "DP"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g3-inc-pa
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "inc-pa"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_6.csv")
head(lookup[,1:4])
colnames(lookup[,1:4])
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status = 1,
    extractionID = 3,
    ID = 4
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_6", Status = NA) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("../../gradients3/g3_lightdark_pa_metat/G3LightDark_polyA_sample_metadata.csv")
head(meta1)
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G3PA_LD_standards/LD_Gradients3_discrete_sampling_googlesheet.csv")
head(meta2)

colnames(meta1)
colnames(meta2)
colnames(meta2[,c(9,15)])

meta <- merge(meta2[,c(9,15)],meta1, by.x = "Sample.ID", by.y = "Alias1", all.y = T )
head(meta$DatetimeHST_TimeOfFiltering.)
colnames(meta)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(DatetimeHST_TimeOfFiltering.,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    Datetime.local.rounded = round_date(Datetime.local, unit = "hour"),
    Hour = format(Datetime.local.rounded, "%H"),  # 2-digit hour
    dataset = sapply(strsplit(SampleID, "\\."), `[`, 4),
    Treatment = case_when(
      Experimental.treatment == "Dark"  ~ "d",
      Experimental.treatment == "Light" ~ "l",
      TRUE ~ Experimental.treatment   # keep anything else unchanged
    )
  )

#get meta in right format: month.day.filter.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Sample.ID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".ld.", Treatment , ".",Hour, "h.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
    Type           = "Incubation",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = dataset,
    Station        = Station,                   
    Cast           = NA,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude,
    Longitude.dec  = Longitude,
    Depth.m        = Depth,
    Filter.um      = Filter,
    Volume.L       = V_total..L.,
    Treatment      = Experimental.treatment ,                        
    Incubation.time= Experimental.timepoint,                        
    Notes          = NA
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup, Meta, by = "ID")

#check dimentions df against lookup and against unique
lookup %>%
  filter(startsWith(ID, "LD"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g3-diel-ns
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "diel-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_9.csv")
head(lookup[,1:4])
colnames(lookup[,1:4])
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status = 1,
    extractionID = 3
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_9", ID = NA) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("../../gradients3/g3_diel_pa_metat/sample_metadata.csv")
head(meta1)
meta2 <- read.csv("/Users/sachacoesel/Documents/Gradients/Standards/G3PA_diel_standards/G3PA.diel.V.filtered.csv")
head(meta2)
meta2$Sample.ID <- gsub("\\.", " ", meta2$Sample.ID)

colnames(meta1)
colnames(meta2)
colnames(meta2[,c(3:4)])

meta <- merge(meta2[,c(3:4)],meta1, by.x = "Sample.ID", by.y = "Alias1", all.y = T )
head(meta$Datetime)
colnames(meta)
#prepare Datetime.local and Datetime.utc columns in ISO 8601 format
meta <- meta %>%
  mutate(
    # Parse your original column (currently "m/d/y H:M")
    Datetime.local = as.POSIXct(Datetime,
                                format = "%m/%d/%y %H:%M",
                                tz = "Pacific/Honolulu"),  # HST = Honolulu tz
    # Convert to UTC
    Datetime.utc   = with_tz(Datetime.local, "UTC"),
    Datetime.local.rounded = round_date(Datetime.local, unit = "hour"),
    Hour = format(Datetime.local.rounded, "%H")  # 2-digit hour
    )

#get meta in right format: month.day.filter.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = Alias2,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast , ".",Hour, "h.", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
    Type           = "diel",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Hour,
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
    Treatment      = NA ,                        
    Incubation.time= NA,                        
    Notes          = Notes
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup, Meta, by.x = "extractionID", by.y = "ID")
#re-order
colnames(df)
df <- df[,c(5,2:4,1,6:24)]

#check dimentions df against lookup and against unique
lookup %>%
  nrow() == nrow(df) #44 entrees start with D, followed by nr
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

###############################################################################
#g3-ctd-ns
###############################################################################

#Set dataset specifics
Gradients = "g3"
SampleType = "ctd-ns"

#load the lookup table
lookup <- read.csv("../lookup_armbrust_files/lookup_armbrust_grc_rnaseq_9.csv")
head(lookup[,1:4])
colnames(lookup[,1:4])
lookup <- lookup %>%
  rename(
    sequencingID = 2,
    Status = 1,
    extractionID = 3
  ) %>%
  mutate(LookupID = "armbrust_grc_rnaseq_9", ID = NA) %>%
  select(LookupID, sequencingID, Status, extractionID, ID)
head(lookup)

#load metadata
meta1 <- read.csv("g3-ctd-pa.sample.metadata.csv")
head(meta1)
meta2 <- read.csv("../../gradients3/g3_depth_pa_metat/G3Depth_polyA_norm_factors.csv")
head(meta2)
colnames(meta2)
colnames(meta1)
meta <- merge(meta1[,-c(3,6,10)],meta2[,c(18,19:20)], by.x="ID", by.y = "Add.l.ID")

#get meta in right format: month.day.filter.rep
colnames(meta)
Meta <- meta %>%
  transmute(
    ID             = ID,
    extractionID   = extractionID,
    SampleID       = paste0(Gradients, "-", SampleType, ".S", Station, ".C", Cast, ".", Dataset, ".", Replicate),
    Cruise         = Cruise,                   
    Gradients      = "G3",                      
    Type           = "CTD",
    project        = paste0(Gradients, "-", SampleType),
    Dataset        = Dataset,
    Station        = Station,                   
    Cast           = Cast,                        
    Replicate      = Replicate,
    Datetime.local = Datetime.local,            
    Datetime.utc   = Datetime.utc,            
    Latitude.dec   = Latitude.dec,
    Longitude.dec  = Longitude.dec,
    Depth.m        = Depth.m,
    Filter.um      = Filter.um,
    Volume.L       = Volume.L,
    Treatment      = NA,                        
    Incubation.time= NA,                        
    Notes          = NA
  )

#merge and save metadata file
colnames(lookup)
colnames(Meta)
df <- merge(lookup[,-5], Meta, by = "extractionID")
#re-order
colnames(df)
df <- df[,c(5,2:4,1,6:24)]

#check dimensions df against lookup and against unique
lookup %>%
  filter(startsWith(extractionID, "DP"))%>%
  nrow() == nrow(df)
length(unique(df$SampleID)) == nrow(df)

#save file
write.csv(df, file = paste0(Gradients, "-", SampleType, ".sample.metadata.csv"), row.names=FALSE)

