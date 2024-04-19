


#
# No KSADS diagnosis
# Match based on: Sex, Gender, IQ, head motion
# CBCL less than t-score of 60 across the board for DSM-oriented scales
# 

# ksads
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
ksads <- read.csv('ksads_filtered.csv')
ksads_filtered <- ksads[ksads$any==0,]
filtered_ids <- ksads_filtered$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject

# cbcl
setwd("/Users/adamkaminski/Desktop/stimproj/data/ABCD_data/")
cbcl <- read.csv('abcd_cbcls01.csv')
cbcl <- cbcl[cbcl$A.ID%in%filtered_ids,]

# has data at baseline and Y2
cbcl_baseline <- cbcl[cbcl$time=="0_baseline_year_1_arm_1",]
cbcl_1 <- cbcl[cbcl$time=="1_year_follow_up_y_arm_1",]
cbcl_2 <- cbcl[cbcl$time=="2_year_follow_up_y_arm_1",]

full_IDs <- cbcl_baseline$A.ID[cbcl_baseline$A.ID%in%cbcl_2$A.ID]
cbcl <- cbcl[cbcl$A.ID%in%full_IDs,]

# gas below 60 on all DSM-oriented scales
cbcl$times_60_above <- 0
for(col in colnames(cbcl)) {
  if(any(!(col %in% c("A.ID","age","sex","time")))) {
    cbcl_integer <- as.integer(cbcl[[col]]>=60)
    cbcl_integer <- replace(cbcl_integer, is.na(cbcl_integer), 0)
    cbcl$times_60_above <- cbcl$times_60_above + cbcl_integer
  }
}

cbcl_filtered <- cbcl[cbcl$times_60_above==0,]

# Remove Y1 time
cbcl_filtered2 <- cbcl_filtered[cbcl_filtered$time!="1_year_follow_up_y_arm_1",]

