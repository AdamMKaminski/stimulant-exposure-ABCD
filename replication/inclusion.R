

setwd("/Users/adamkaminski/Desktop/stimproj/data/")
ksads <- read.csv('abcd_ksad01.csv')

# Remove "Symptom..." columns
cols <- colnames(ksads)
symptom_cols <- grepl("Symptom...", cols)
ksads <- ksads[, !symptom_cols]

# Remove "Past" columns
cols <- colnames(ksads)
Past_cols <- grepl("Past", cols)
ksads <- ksads[, !past_cols]

# Remove "PAST" columns
cols <- colnames(ksads)
PAST_cols <- grepl("PAST", cols)
ksads <- ksads[, !PAST_cols]

# Remove "In.Partial.Remission" columns
cols <- colnames(ksads)
remission_cols <- grepl("In.Partial.Remission", cols)
ksads <- ksads[, !remission_cols]

# write
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
write.csv(ksads,file="ksads_filtered.csv")

# read
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
ksads <- read.csv('ksads_filtered.csv')

ksads <- ksads[ksads$any==1,]

# Recommended use
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
imgincl <- read.csv('abcd_imgincl01.csv')
imgincl <- imgincl[imgincl$nBack_SST==1,]
use_IDs <- imgincl$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject

IDs <- ksads$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject[ksads$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject%in%use_IDs] 



