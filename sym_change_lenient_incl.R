


# Get IDs with ADHD at baseline based on KSADS
setwd("/Users/adamkaminski/Desktop/stimproj/data/ABCD_data/")
ksads <- read.csv('abcd_ksad01_adhd_present.csv')
ksads <- ksads[,-2:-3]
ksads <- reshape(ksads, idvar = "ID", timevar = "time", direction = "wide")
origlen <- dim(ksads)[1]
ksads <- ksads[,-3:-4]
ksads <- ksads[ksads$adhd_present.0==1,]
ksads <- ksads[complete.cases(ksads), ]

baseline_adhd_IDs <- ksads$ID
writeLines(paste(dim(ksads)[1],"/",origlen," (", (dim(ksads)[1]/origlen)*100,"%) of people have ADHD at baseline according to KSADS",sep=""))

# Get CBCL data -- is there CBCL data for these subjects at baseline and Y2? Exclude if not
setwd("/Users/adamkaminski/Desktop/stimproj/data/ABCD_data/")
cbcl <- read.csv('abcd_cbcls01.csv')
cbcl <- cbcl[cbcl$A.ID%in%baseline_adhd_IDs,]

cbcl_baseline <- cbcl[cbcl$time=="0_baseline_year_1_arm_1",]
cbcl_1 <- cbcl[cbcl$time=="1_year_follow_up_y_arm_1",]
cbcl_2 <- cbcl[cbcl$time=="2_year_follow_up_y_arm_1",]

full_IDs <- cbcl_baseline$A.ID[cbcl_baseline$A.ID%in%cbcl_2$A.ID]
writeLines(paste(length(full_IDs),"/",dim(cbcl_baseline)[1]," (", (length(full_IDs)/dim(cbcl_baseline)[1])*100,"%) of people who have CBCL baseline data also have CBCL Y2 data",sep=""))

# Exclude subjects whose med data was previously decoded
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
orig_dat <- read.csv('stimulant_master_df.csv')
orig_IDs <- orig_dat$ID

full_IDs_new <- full_IDs[!(full_IDs%in%orig_IDs)]

# Get Medication data, save csv to be decoded
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
med <- read.csv('medsy01.csv')

med_subs <- med[med$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject%in%full_IDs_new,]
write.csv(med_subs,"medsy01_moreIDs.csv",na = "",row.names = FALSE)


