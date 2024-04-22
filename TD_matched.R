


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

# clean up and reshape df
cbcl$time[cbcl$time=="0_baseline_year_1_arm_1"] <- "0"
cbcl$time[cbcl$time=="1_year_follow_up_y_arm_1"] <- "1"
cbcl$time[cbcl$time=="2_year_follow_up_y_arm_1"] <- "2"
cbcl <- cbcl[cbcl$time!="1",]
cbcl_wide <- reshape(cbcl, idvar = "A.ID", timevar = "time", direction = "wide")

# has below 60 on all DSM-oriented scales
cbcl_wide$times_60_above <- 0

for(col in colnames(cbcl_wide)) {
  if(any(!(col %in% c("A.ID","age.2","sex.2","age.0","sex.0")))) {
    cbcl_integer <- as.integer(cbcl_wide[[col]]>=60)
    cbcl_integer <- replace(cbcl_integer, is.na(cbcl_integer), 0)
    cbcl_wide$times_60_above <- cbcl_wide$times_60_above + cbcl_integer
  }
}

cbcl_wide_filtered <- cbcl_wide[cbcl_wide$times_60_above==0,]

# Remove overlapping IDs
origin_IDs <- read.csv('original_analysis_IDs.csv')
cbcl_wide_filtered2 <- cbcl_wide_filtered[!(cbcl_wide_filtered$A.ID%in%origin_IDs$x),]

# Remove people with NOT recommended resting state data
setwd("/Users/adamkaminski/Desktop/stimproj/data/ABCD_data/")
imgincl <- read.csv('abcd_imgincl01_rsfmri.csv')
imgincl <- imgincl[,-2:-3]
imgincl <- reshape(imgincl, idvar = "ID", timevar = "time", direction = "wide")
imgincl <- imgincl[imgincl$imgincl_rsfmri_include.0==1,]
imgincl <- imgincl[imgincl$imgincl_rsfmri_include.2==1,]
imgincl <- imgincl[complete.cases(imgincl), ]
cbcl_wide_filtered3 <- cbcl_wide_filtered2[cbcl_wide_filtered2$A.ID%in%imgincl$ID,]

# Get imaging data
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
netdat <- read.csv('mrirscor02_v1.csv')

