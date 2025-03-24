
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(brms)

#
#
# SELECTING TD GROUP
#
#

# Criteria:
#  - No KSADS diagnosis
#  - CBCL less than t-score of 60 across the board for DSM-oriented scales
#  - Compare on: Sex, Gender, IQ, head motion

# ksads
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
ksads <- read.csv('ksads_filtered.csv')
ksads_filtered <- ksads[ksads$any==0,]
filtered_ids <- ksads_filtered$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject

# cbcl
setwd("/Users/adamkaminski/Desktop/stimulant-exposure-ABCD")
cbcl <- read.csv('abcd_cbcls01_2.csv')
cbcl <- cbcl[cbcl$ID%in%filtered_ids,]

# remove unwanted cols
cbcl <- cbcl %>%
  select(-contains("..raw.score"))

cbcl <- cbcl %>%
  select(-contains("..missing.values"))

cbcl <- cbcl %>%
  select(-contains("..number.of.missing.values"))

cbcl <- cbcl %>%
  select(-contains("collection_title"))

cbcl <- cbcl %>%
  select(-contains("study_cohort_name"))

cbcl <- cbcl %>%
  select(-contains("Syndrome"))

# has data at baseline and Y2
cbcl_baseline <- cbcl[cbcl$time=="baseline_year_1_arm_1",]
cbcl_1 <- cbcl[cbcl$time=="1_year_follow_up_y_arm_1",]
cbcl_2 <- cbcl[cbcl$time=="2_year_follow_up_y_arm_1",]

full_IDs <- cbcl_baseline$ID[cbcl_baseline$ID%in%cbcl_2$ID]
cbcl <- cbcl[cbcl$ID%in%full_IDs,]

# clean up and reshape df
cbcl$time[cbcl$time=="baseline_year_1_arm_1"] <- "0"
cbcl$time[cbcl$time=="1_year_follow_up_y_arm_1"] <- "1"
cbcl$time[cbcl$time=="2_year_follow_up_y_arm_1"] <- "2"
cbcl <- cbcl[cbcl$time!="1",]
cbcl_wide <- reshape(cbcl, idvar = "ID", timevar = "time", direction = "wide")

# has below 60 on all DSM-oriented scales
cbcl_wide$times_60_above <- 0

for(col in colnames(cbcl_wide)) {
  if(any(!(col %in% c("ID","age.2","sex.2","age.0","sex.0")))) {
    cbcl_integer <- as.integer(cbcl_wide[[col]]>=60)
    cbcl_integer <- replace(cbcl_integer, is.na(cbcl_integer), 0)
    cbcl_wide$times_60_above <- cbcl_wide$times_60_above + cbcl_integer
  }
}

cbcl_wide_filtered <- cbcl_wide[cbcl_wide$times_60_above==0,]

# Remove overlapping IDs
origin_IDs <- read.csv('original_analysis_IDs.csv')
cbcl_wide_filtered2 <- cbcl_wide_filtered[!(cbcl_wide_filtered$ID%in%origin_IDs$x),]

# Remove people with NOT recommended resting state data
imgincl <- read.csv('abcd_imgincl01_rsfmri.csv')
imgincl <- imgincl[,-2:-3]
imgincl <- reshape(imgincl, idvar = "ID", timevar = "time", direction = "wide")
imgincl <- imgincl[imgincl$imgincl_rsfmri_include.0==1,]
imgincl <- imgincl[imgincl$imgincl_rsfmri_include.2==1,]
imgincl <- imgincl[complete.cases(imgincl), ]
cbcl_wide_filtered3 <- cbcl_wide_filtered2[cbcl_wide_filtered2$ID%in%imgincl$ID,]

#
#
# RS-FC DATA
#
#


# Get imaging data
setwd("/Users/adamkaminski/Downloads")
netdat <- read.csv('mrirscor02_v1.csv')

# filter based on included subjects up to this point
netdat_filt <- netdat[netdat$ID%in%cbcl_wide_filtered3$ID,]

# neaten time column
netdat_filt$time[netdat_filt$time=="baseline_year_1_arm_1"] <- 0
netdat_filt$time[netdat_filt$time=="2_year_follow_up_y_arm_1"] <- 2

# reshape
netdat_filt_wide <- reshape(netdat_filt, idvar = "ID", timevar = "time", direction = "wide") 

# remove people if FD is not < 0.5mm 
netdat_filt_wide$baseline_good_motion <- as.integer(netdat_filt_wide$avg_FD.0<0.5)
netdat_filt_wide$Y2_good_motion <- as.integer(netdat_filt_wide$avg_FD.2<0.5)

netdat_filt2_wide <- netdat_filt_wide[(netdat_filt_wide$baseline_good_motion==1 & netdat_filt_wide$Y2_good_motion==1),]

# rename df
netdat <- netdat_filt2_wide
netdat <- netdat[,-c(ncol(netdat)-1,ncol(netdat))]

setwd("/Users/adamkaminski/Downloads/")
write.csv(netdat, "td_IDs.csv", row.names = FALSE)


# combine networks
combine <- 1
if (combine == 1) {
  CON <- as.data.frame(netdat[,grepl("_cerc_",colnames(netdat),fixed=TRUE)])
  SAL <- as.data.frame(netdat[,grepl("_sa_",colnames(netdat),fixed=TRUE)])
  SMH <- netdat[,grepl("_smh_",colnames(netdat),fixed=TRUE)]
  SMM <- netdat[,grepl("_smm_",colnames(netdat),fixed=TRUE)]
  
  avgconsal <- matrix(0,dim(netdat)[1],dim(CON)[2])
  avgsmhsmm <- matrix(0,dim(netdat)[1],dim(SMH)[2])
  for (i in 1:dim(CON)[2]) {
    avgconsal[,i] <- (as.numeric(CON[,i]) + as.numeric(SAL[,i]))/2
    avgsmhsmm[,i] <- (as.numeric(SMH[,i]) + as.numeric(SMM[,i]))/2
  }
  colnames(avgconsal) <-  gsub('_cerc_','_cercsa_',colnames(CON))
  colnames(avgsmhsmm) <- gsub('_smh_','_smhsmm_',colnames(SMH))
  
  netdat <- netdat[,!grepl("_cerc_",colnames(netdat),fixed=TRUE)]
  netdat <- netdat[,!grepl("_sa_",colnames(netdat),fixed=TRUE)]
  netdat <- netdat[,!grepl("_smh_",colnames(netdat),fixed=TRUE)]
  netdat <- netdat[,!grepl("_smm_",colnames(netdat),fixed=TRUE)]
  
  netdat <- cbind(netdat,avgconsal)
  netdat <- cbind(netdat,avgsmhsmm)
}

# calculate FC difference:
# define net codes
if (combine == 1) {
  net_codes <- c("_au_","_cercsa_","_copa_","_df_","_dsa_","_fopa_","_rst_",
                 "_smhsmm_","_vta_","vs")
} else {
  net_codes <- c("_au_","_cerc_","_copa_","_df_","_dsa_","_fopa_","_rst_",
                 "_sa_","_smh_","_smm_","_vta_","_vs_")
}

if (combine == 1) {n<-(10*6)} else {n<-(12*6)}
diff <- as.data.frame(matrix(0,dim(netdat)[1],n))
count <- 1
for (i in 1:length(net_codes)) {
  temp_net <- netdat[,grepl(net_codes[i],colnames(netdat),fixed=TRUE)]
  for (se in 1:6) {
    str1 <- colnames(temp_net)[se]
    str2 <- colnames(temp_net)[se+6]
    len <- nchar(str1)
    if (substr(str1,1,len-2) == substr(str2,1,len-2)) {
      new_str <- gsub('rsfmri_cor_ngd_','',substr(str1,1,len-2))
      new_str <- gsub('_scs','',new_str)
      colnames(diff)[count] <- new_str
      diff[,count] <- as.numeric(temp_net[,se]) - as.numeric(temp_net[,se+6])
      count <- count + 1
    } else {
      stop(paste("These colnames do not match:",se,"and",se+6))
    }
  }
}

# remove unnecessary cols (0 + 2 FC) from netdat and add change in FC
# if interested in difference scores
# if (combine == 1) {netdat <- netdat[,c(1:10,59:67)]}
# netdat <- cbind(netdat,diff)

# remove unnecessary cols (all FC besides sign. from original analysis for comparison, also desired covariates)
netdat_TD_orig_cxs <- netdat[,c(1:10,18,31,33,35,36,51,57,59:67,75,88,90,92,93,108,114,116,122,129,135)]

# if you don't want to remove
netdat_TD_orig_cxs <- netdat

# make df long
netdat_TD_orig_cxs_long <- netdat_TD_orig_cxs %>%
  pivot_longer(cols = -ID, 
               names_to = c(".value", "time"), 
               names_pattern = "(.*?)\\.(\\d+)$")

#
#
# COVARIATES
#
#

# add other covariates
# site
setwd("/Users/adamkaminski/Desktop/stimulant-exposure-ABCD/")
sites <- read.csv("sites.csv")
sites_filtered <- sites[sites$ID %in% netdat_TD_orig_cxs_long$ID,]
sites_filtered <- sites_filtered[order(sites_filtered$ID),]
netdat_TD_orig_cxs_long <- netdat_TD_orig_cxs_long[order(netdat_TD_orig_cxs_long$ID),]

if (all(sites_filtered$ID == netdat_TD_orig_cxs_long$ID)) {
  netdat_TD_orig_cxs_long$site <- sites_filtered$Site
}

# cbcl
cbcl_filtered <- cbcl[cbcl$ID %in% netdat_TD_orig_cxs_long$ID,]
cbcl_filtered <- cbcl_filtered[order(cbcl_filtered$ID),]
netdat_TD_orig_cxs_long <- netdat_TD_orig_cxs_long[order(netdat_TD_orig_cxs_long$ID),]

if (all(cbcl_filtered$ID == netdat_TD_orig_cxs_long$ID)) {
  netdat_TD_orig_cxs_long$ADHD.Problems <- cbcl_filtered$ADHD.CBCL.DSM5.Scale..t.score.
}

# IQ
IQ <- read.csv('nih_toolbox.csv')

IQ_TD <- IQ[IQ$src_subject_id %in% cbcl_filtered$ID,]
IQ_TD <- IQ_TD[,c(1,4,5)]
colnames(IQ_TD) <- c("ID","time","iq")
IQ_TD$time[IQ_TD$time=='baseline_year_1_arm_1'] <- 0
IQ_TD$time[IQ_TD$time=='2_year_follow_up_y_arm_1'] <- 2
IQ_TD_wide <- reshape(IQ_TD, idvar = "ID", timevar = "time", direction = "wide") 

IQ_adhd <- IQ[IQ$src_subject_id %in% orig_dat$ID,]  # from below
IQ_adhd <- IQ_adhd[,c(1,4,5)]
colnames(IQ_adhd) <- c("ID","time","iq")
IQ_adhd$time[IQ_adhd$time=='baseline_year_1_arm_1'] <- 0
IQ_adhd$time[IQ_adhd$time=='2_year_follow_up_y_arm_1'] <- 2
IQ_adhd_wide <- reshape(IQ_adhd, idvar = "ID", timevar = "time", direction = "wide") 

# SES
#setwd("/Users/adamkaminski/Desktop/stimproj/data/nda_downloads/")
demos <- read.csv("demographics.csv",header=TRUE)
demos <- demos[demos$src_subject_id %in% cbcl_filtered$ID,]
demos$parent_ed[demos$parent_ed==777] <- NA
demos$parent_ed_z <- (demos$parent_ed-mean(demos$parent_ed,na.rm=TRUE))/sd(demos$parent_ed,na.rm=TRUE)
demos$partner_ed[demos$partner_ed==999] <- NA
demos$partner_ed[demos$partner_ed==777] <- NA
demos$partner_ed_z <- (demos$partner_ed-mean(demos$partner_ed,na.rm=TRUE))/sd(demos$partner_ed,na.rm=TRUE) 
demos$comb_inc[demos$comb_inc==999] <- NA
demos$comb_inc[demos$comb_inc==777] <- NA
demos$comb_inc_z <- (demos$comb_inc-mean(demos$comb_inc,na.rm=TRUE))/sd(demos$comb_inc,na.rm=TRUE)
demos$ses <- rowMeans(demos[,c(8:10)], na.rm = TRUE)

demos <- demos[order(demos$src_subject_id),]
netdat_TD_orig_cxs_long <- netdat_TD_orig_cxs_long[order(netdat_TD_orig_cxs_long$time),]
netdat_TD_orig_cxs_long$ses <- 0
if (all(demos$src_subject_id == netdat_TD_orig_cxs_long$ID[1:1420])) {
  netdat_TD_orig_cxs_long$ses[1:1420] <- demos$ses
}
if (all(demos$src_subject_id == netdat_TD_orig_cxs_long$ID[1421:2840])) {
  netdat_TD_orig_cxs_long$ses[1421:2840] <- demos$ses
}

# Comorbidities



# plot FC change
# Left Caudate - Frontoparietal Network
ggplot(netdat_TD_orig_cxs_long, aes(x=time, y=rsfmri_cor_ngd_fopa_scs_cdelh, color=time)) +
  geom_line(aes(x = time, group = ID), 
            size = 0.5, color = 'gray') + 
  geom_point(aes(x = time), size = 2, position = position_dodge(width = 0.75)) + 
  geom_boxplot() +
  scale_color_manual(values=c("gray", "black"), labels = c("Baseline", "Year 2")) +
  xlab("Time") + ylab("rs-FC") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  guides(color=guide_legend(title="Time")) + 
  ggtitle("Left Caudate - Frontoparietal Network") +
  theme_minimal_grid(12)

# Left Putamen - Frontoparietal Network
ggplot(netdat_TD_orig_cxs_long, aes(x=time, y=rsfmri_cor_ngd_fopa_scs_ptlh, color=time)) +
  geom_line(aes(x = time, group = ID), 
            size = 0.5, color = 'gray') + 
  geom_point(aes(x = time), size = 2, position = position_dodge(width = 0.75)) + 
  geom_boxplot() +
  scale_color_manual(values=c("gray", "black"), labels = c("Baseline", "Year 2")) +
  xlab("Time") + ylab("rs-FC") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  guides(color=guide_legend(title="Time")) + 
  ggtitle("Left Putamen - Frontoparietal Network") +
  theme_minimal_grid(12)

# Right Putamen - Visual Network
ggplot(netdat_TD_orig_cxs_long, aes(x=time, y=rsfmri_cor_ngd_vs_scs_ptrh, color=time)) +
  geom_line(aes(x = time, group = ID), 
            size = 0.5, color = 'gray') + 
  geom_point(aes(x = time), size = 2, position = position_dodge(width = 0.75)) + 
  geom_boxplot() +
  scale_color_manual(values=c("gray", "black"), labels = c("Baseline", "Year 2")) +
  xlab("Time") + ylab("rs-FC") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  guides(color=guide_legend(title="Time")) + 
  ggtitle("Right Putamen - Visual Network") +
  theme_minimal_grid(12)

# Paired t tests for rs-FC change
netdat_TD_orig_cxs_long_0 <- netdat_TD_orig_cxs_long[netdat_TD_orig_cxs_long$time=="0",]
netdat_TD_orig_cxs_long_2 <- netdat_TD_orig_cxs_long[netdat_TD_orig_cxs_long$time=="2",]

t.test(netdat_TD_orig_cxs_long_0$rsfmri_cor_ngd_fopa_scs_cdelh,netdat_TD_orig_cxs_long_2$rsfmri_cor_ngd_fopa_scs_cdelh,paired=TRUE)
t.test(netdat_TD_orig_cxs_long_0$rsfmri_cor_ngd_fopa_scs_ptlh,netdat_TD_orig_cxs_long_2$rsfmri_cor_ngd_fopa_scs_ptlh,paired=TRUE)
t.test(netdat_TD_orig_cxs_long_0$rsfmri_cor_ngd_vs_scs_ptrh,netdat_TD_orig_cxs_long_2$rsfmri_cor_ngd_vs_scs_ptrh,paired=TRUE)

#
#
# ADHD data
#
#

# Get original covariate data to compare (make long)
setwd("/Users/adamkaminski/Desktop/stimulant-exposure-ABCD")
orig_dat <- read.csv('stimulant_master_df.csv')
orig_dat_filitered <- orig_dat[,c(1,2,3,6,16,21,22,25,40,56,59,66,69,132,135)]
colnames(orig_dat_filitered)[5] <- c("site")
colnames(orig_dat_filitered)[c(10:13,15)] <- c("FD.under.2","meanFD.2","FD.under.0","meanFD.0","group")
orig_dat_filitered$group[orig_dat_filitered$group==0] <- 'ADHD_naive'
orig_dat_filitered$group[orig_dat_filitered$group==1] <- 'ADHD_stimulant'

# make long
orig_dat_filitered_long <- orig_dat_filitered %>%
  pivot_longer(cols = -c(ID,group,ses,comorbs,site), 
               names_to = c(".value", "time"), 
               names_pattern = "(.*?)\\.(\\d+)$")

# Get long resting state data

setwd("/Users/adamkaminski/Desktop/stimulant-exposure-ABCD")
rsfc <- read.csv('orig_netdat_lt.csv')
colnames(rsfc)[87] <- c("smhsmm_ptlh")
rsfc_filtered <- rsfc[,c(2,27,40,42,44,45,78,84,87,92)]

rsfc_filtered <- rsfc_filtered[order(rsfc_filtered$src_subject_id),]
orig_dat_filitered_long <- orig_dat_filitered_long[order(orig_dat_filitered_long$ID),]
if (all(rsfc_filtered$src_subject_id==orig_dat_filitered_long$ID)) {
  orig_dat_filitered_long <- cbind(orig_dat_filitered_long,rsfc_filtered[2:10])
}

colnames(orig_dat_filitered_long)[c(12:18)] <- c("copa_ptlh","dsa_aalh","dsa_ptrh","fopa_cdelh","fopa_ptlh","vta_ptrh","vs_ptrh")

#
#
# compare TD, ADHD stimulant-exposed, and ADHD naive
#
#

# get desired variables from TD dataset
# just significant conections:
netdat_TD_orig_cxs_long_filtered <- netdat_TD_orig_cxs_long[,c(1,2,3,4,8,11,12:23)]
# all connections:
netdat_TD_orig_cxs_long_filtered <- netdat_TD_orig_cxs_long[,c(1,2,3,4,8,11,12:86)]
colnames(netdat_TD_orig_cxs_long_filtered)[c(5,6)] <- c("FD.under","meanFD")
colnames(netdat_TD_orig_cxs_long_filtered)[c(7:15)] <- c("copa_ptlh","dsa_aalh","dsa_ptrh","fopa_cdelh","fopa_ptlh","vta_ptrh","vs_ptrh","cercsa_cdelh","smhsmm_ptlh")
netdat_TD_orig_cxs_long_filtered$group <- 'TD'
netdat_TD_orig_cxs_long_filtered$comorbs <- 0

# rearrange columns
netdat_TD_orig_cxs_long_filtered <- netdat_TD_orig_cxs_long_filtered %>%
  select(ID, site, comorbs, ses, group, time, age, sex, ADHD.Problems, FD.under, meanFD, copa_ptlh, dsa_aalh, dsa_ptrh, fopa_cdelh, fopa_ptlh, vta_ptrh, vs_ptrh, smhsmm_ptlh, cercsa_cdelh)

# simplify df names
TD_dat <- netdat_TD_orig_cxs_long_filtered
ADHD_dat <- orig_dat_filitered_long

# add covariate columns 
IQ_TD <- IQ_TD[order(IQ_TD$time),]
IQ_TD <- IQ_TD[order(IQ_TD$ID),]
TD_dat <- TD_dat[order(TD_dat$time),]
TD_dat <- TD_dat[order(TD_dat$ID),]
if (all(TD_dat$ID == IQ_TD$ID)) {
  if (all(TD_dat$time == IQ_TD$time)) {
    TD_dat$IQ <- IQ_TD$iq
  }}

IQ_adhd <- IQ_adhd[order(IQ_adhd$ID),]
if (all(ADHD_dat$ID == IQ_adhd$ID)) {
  ADHD_dat$IQ <- IQ_adhd$iq
}

# make FD.under sqrt if not
if (mean(ADHD_dat$FD.under)>1000) {
  ADHD_dat$FD.under <- sqrt(ADHD_dat$FD.under)
}
if (mean(TD_dat$FD.under)>1000) {
  TD_dat$FD.under.sqrt <- sqrt(TD_dat$FD.under)
}

combined_df <- rbind(ADHD_dat,TD_dat)

# test for group differences:

# length of good data (sqrt)
# combined_df$FD.under <- sqrt(combined_df$FD.under)

t.test(combined_df$FD.under[combined_df$group=='TD'],combined_df$FD.under[combined_df$group=='ADHD_stimulant'],paired=FALSE)
t.test(combined_df$FD.under[combined_df$group=='TD'],combined_df$FD.under[combined_df$group=='ADHD_naive'],paired=FALSE)

# mean FD
t.test(combined_df$meanFD[combined_df$group=='TD'],combined_df$meanFD[combined_df$group=='ADHD_stimulant'],paired=FALSE)
t.test(combined_df$meanFD[combined_df$group=='TD'],combined_df$meanFD[combined_df$group=='ADHD_naive'],paired=FALSE)

# sex
chisq.test(combined_df$sex[combined_df$group!='ADHD_naive'],combined_df$group[combined_df$group!='ADHD_naive'])
chisq.test(combined_df$sex[combined_df$group!='ADHD_stimulant'],combined_df$group[combined_df$group!='ADHD_stimulant'])

# IQ
t.test(combined_df$IQ[combined_df$group=='TD'],combined_df$IQ[combined_df$group=='ADHD_stimulant'],paired=FALSE)
t.test(combined_df$IQ[combined_df$group=='TD'],combined_df$IQ[combined_df$group=='ADHD_naive'],paired=FALSE)

#
#
# Final
# Bayesian hierarchical logistic regressions
#
#

# No naive
combined_df_no_naive <- combined_df[combined_df$group!="ADHD_naive",]

combined_df_no_naive_lt_lc <- combined_df_no_naive %>%
  pivot_longer(cols = c(copa_ptlh, dsa_aalh, dsa_ptrh, fopa_cdelh, fopa_ptlh, vta_ptrh, vs_ptrh, smhsmm_ptlh, cercsa_cdelh), 
               names_to = "Connection", values_to = "rsFC")

write.csv(combined_df_no_naive_lt_lc, "combined_df_no_naive_lt_lc.csv", row.names = FALSE)

# No exposed
combined_df_no_stim <- combined_df[combined_df$group!="ADHD_stimulant",]

combined_df_no_stim_lt_lc <- combined_df_no_stim %>%
  pivot_longer(cols = c(copa_ptlh, dsa_aalh, dsa_ptrh, fopa_cdelh, fopa_ptlh, vta_ptrh, vs_ptrh, smhsmm_ptlh, cercsa_cdelh), 
               names_to = "Connection", values_to = "rsFC")

write.csv(combined_df_no_stim_lt_lc, "combined_df_no_stim_lt_lc.csv", row.names = FALSE)


slope.prior.all <- c(
  # Covariates
  prior(normal(0,10), class=b, coef=meanFD),
  prior(normal(0,10), class=b, coef=FD.under),
  prior(normal(0,10), class=b, coef=sexM),
  prior(normal(0,10), class=b, coef=IQ),
  prior(normal(0,10), class=b, coef=ADHD.Problems),
  prior(normal(0,10), class=b, coef=groupTD),
  prior(normal(0,10), class=b, coef=time2),
  prior(normal(0,10), class=b, coef=Connectioncopa_ptlh),
  prior(normal(0,10), class=b, coef=Connectiondsa_aalh),
  prior(normal(0,10), class=b, coef=Connectiondsa_ptrh),
  prior(normal(0,10), class=b, coef=Connectionfopa_cdelh),
  prior(normal(0,10), class=b, coef=Connectionfopa_ptlh),
  prior(normal(0,10), class=b, coef=Connectionsmhsmm_ptlh),
  prior(normal(0,10), class=b, coef=Connectionvs_ptrh),
  prior(normal(0,10), class=b, coef=Connectionvta_ptrh),
  # 2-way interactions
  prior(normal(0,10), class=b, coef=groupTD:Connectioncopa_ptlh),
  prior(normal(0,10), class=b, coef=groupTD:Connectiondsa_aalh),
  prior(normal(0,10), class=b, coef=groupTD:Connectiondsa_ptrh),
  prior(normal(0,10), class=b, coef=groupTD:Connectionfopa_cdelh),
  prior(normal(0,10), class=b, coef=groupTD:Connectionfopa_ptlh),
  prior(normal(0,10), class=b, coef=groupTD:Connectionsmhsmm_ptlh),
  prior(normal(0,10), class=b, coef=groupTD:Connectionvs_ptrh),
  prior(normal(0,10), class=b, coef=groupTD:Connectionvta_ptrh),
  prior(normal(0,10), class=b, coef=groupTD:time2),
  prior(normal(0,10), class=b, coef=time2:Connectioncopa_ptlh),
  prior(normal(0,10), class=b, coef=time2:Connectiondsa_aalh),
  prior(normal(0,10), class=b, coef=time2:Connectiondsa_ptrh),
  prior(normal(0,10), class=b, coef=time2:Connectionfopa_cdelh),
  prior(normal(0,10), class=b, coef=time2:Connectionfopa_ptlh),
  prior(normal(0,10), class=b, coef=time2:Connectionsmhsmm_ptlh),
  prior(normal(0,10), class=b, coef=time2:Connectionvs_ptrh),
  prior(normal(0,10), class=b, coef=time2:Connectionvta_ptrh),
  # 3-way interactions
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectioncopa_ptlh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectiondsa_aalh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectiondsa_ptrh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectionfopa_cdelh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectionfopa_ptlh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectionsmhsmm_ptlh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectionvs_ptrh),
  prior(normal(0,10), class=b, coef=groupTD:time2:Connectionvta_ptrh)
)

# covariates to add
# - SES
# - non-stim med bi -- not plausible 
# - comorbidities -- all 0s for TD

three_way_interaction_TD_stim_no_naive <- brm(rsFC ~ 
                             meanFD + FD.under + sex +
                             IQ + ADHD.Problems +
                             group + time + Connection +
                             group*time +
                             group*Connection +
                             time*Connection +
                             group*time*Connection +
                             (1|site),
                           data=combined_df_no_naive_lt_lc,      
                           prior = slope.prior.all,
                           silent = FALSE,
                           warmup = 2000, iter = 10000)

writeLines(capture.output(summary(three_way_interaction_TD_stim_no_naive)),"step_2_exp_TD.txt")

save.image("env_5_10_2024.RData")

# undersampling 

reps <- 100
results <- array(rep(0, 8*4*reps), dim=c(8,4,reps))

temp_df <- combined_df_no_naive_lt_lc[combined_df_no_naive_lt_lc$group=="TD",]
TD_ids <- unique(temp_df$ID)
STIM_ids <- combined_df_no_naive_lt_lc$ID[combined_df_no_naive_lt_lc$group=="ADHD_stimulant"]
STIM_ids <- unique(STIM_ids)

for(i in 1:reps) {
  
  print(paste('Working on rep ',i,'/100...',sep=''))
  
  temp_TD_ids <- sample(TD_ids,81)
  ids <- c(temp_TD_ids,STIM_ids)
  
  temp_df <- combined_df_no_naive_lt_lc[combined_df_no_naive_lt_lc$ID %in% ids,]
  
  temp <- brm(rsFC ~ meanFD + FD.under + sex +
              IQ + ADHD.Problems +
              group + time + Connection +
              group*time +
              group*Connection +
              time*Connection +
              group*time*Connection + (1|site),
              data=temp_df,      
              prior = slope.prior.all,
              silent = FALSE,
              warmup = 2000, iter = 10000)

  results[,,i] <- posterior_summary(temp)[34:41,]
  
  print('Finished')
}


# paired t tests

td <- combined_df_no_stim[combined_df_no_stim$group=="TD",]
td_0 <- td[td$time==0,]
td_2 <- td[td$time==2,]

# Left caudate-Frontoparietal network
cat("Left Caudate - FPN FC\ TD Group, effect of Time")
t.test(td_0$fopa_cdelh,td_2$fopa_cdelh,paired=TRUE)
mean(td_0$fopa_cdelh,na.rm=TRUE)
mean(td_2$fopa_cdelh,na.rm=TRUE)

# Left putamen-Frontoparietal network
cat("Left Putamen - FPN FC\ TD Group, effect of Time")
t.test(td_0$fopa_ptlh,td_2$fopa_ptlh,paired=TRUE)
mean(td_0$fopa_ptlh)0
mean(td_2$fopa_ptlh)

# Right putamen-Visual network
cat("Right Putamen - FPN FC\ TD Group, effect of Time")
t.test(td_0$vs_ptrh,td_2$vs_ptrh,paired=TRUE)
mean(td_0$vs_ptrh)
mean(td_2$vs_ptrh)

# plot TD change
x_time <- as.factor(as.numeric(td$time))

y_fc <- as.numeric(td$fopa_cdelh)
label <- "Left Caudate - FPN FC"
p <- "p=0.286"

#y_fc <- as.numeric(td$fopa_ptlh)
#label <- "Left Putamen - FPN FC"
#p <- ""

#y_fc <- as.numeric(td$vs_ptrh)
#label <- "Right Putamen - VN FC"
#p <- "p<0.001"

fig6 <- ggplot(td, aes(x=x_time, y=y_fc, color=x_time)) +
  geom_line(aes(x = x_time, group = ID), 
            size = 0.5, color = 'gray') + 
  geom_violin() +
  geom_point(aes(x = x_time), size = 2, position = position_dodge(width = 0.75)) + 
  scale_color_manual(values=c("khaki", "goldenrod"), labels = c("TD — Baseline", "TD — Year 2")) +
  xlab("Group") + ylab("Functional Connectivity") +
  scale_x_discrete(labels = c("TD","TD")) +
  theme(axis.text.x = element_blank(),axis.ticks.margin=unit(0,'cm')) +
  guides(color=guide_legend(title="Time")) + 
  ggtitle(label) + 
  geom_signif(y_position=0.5, xmin=0.9, xmax=2.1,
              annotation=p, tip_length = 0, vjust = -1, color = "black") +
  theme_minimal_grid(12) +
  ylim(-0.5,0.75)


# ANOVAs (3x2)
anova_result <- aov(fopa_cdelh ~ time * group + sex + ADHD.Problems + IQ + meanFD + FD.under, data = combined_df)
summary(anova_result)

anova_result <- aov(fopa_ptlh ~ time * group + sex + ADHD.Problems + IQ + meanFD + FD.under, data = combined_df)
summary(anova_result)

anova_result <- aov(vs_ptrh ~ time * group + sex + ADHD.Problems + IQ + meanFD + FD.under, data = combined_df)
summary(anova_result)

# ANOVAs (2x2)
# no naive
combined_df_no_naive <- combined_df[combined_df$group!="ADHD_naive",]

anova_result <- aov(fopa_cdelh ~ time * group + sex + ADHD.Problems + FD.under, data = combined_df_no_naive)
summary(anova_result)

anova_result <- aov(fopa_ptlh ~ time * group + sex + ADHD.Problems + FD.under, data = combined_df_no_naive)
summary(anova_result)

anova_result <- aov(vs_ptrh ~ time * group + sex + ADHD.Problems + FD.under, data = combined_df_no_naive)
summary(anova_result)

# plot results
ggplot(combined_df_no_naive, aes(x = time, y = fopa_cdelh, color = group)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.5) +  
  geom_line() +
  geom_line(aes(group = group, y = fopa_ptlh), stat = "summary", size = 0.5, linetype="dashed") +
  scale_color_manual(values=c("red", "darkblue"), labels = c("ADHD\nStimulant Exposed","Healthy Control")) +
  labs(x = "Time", y = "Resting State Functional Connectivity", color = "Group") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  ggtitle("Left Caudate - Frontoparietal Network") +
  theme_minimal() 

ggplot(combined_df_no_naive, aes(x = time, y = fopa_ptlh, color = group)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.5) +  
  geom_line() +
  geom_line(aes(group = group, y = fopa_ptlh), stat = "summary", size = 0.5, linetype="dashed") +
  scale_color_manual(values=c("red", "darkblue"), labels = c("ADHD\nStimulant Exposed","Healthy Control")) +
  labs(x = "Time", y = "Resting State Functional Connectivity", color = "Group") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  ggtitle("Left Putamen - Frontoparietal Network") +
  theme_minimal() 

ggplot(combined_df_no_naive, aes(x = time, y = vs_ptrh, color = group)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.5) +  
  geom_line() +
  geom_line(aes(group = group, y = fopa_ptlh), stat = "summary", size = 0.5, linetype="dashed") +
  scale_color_manual(values=c("red", "darkblue"), labels = c("ADHD\nStimulant Exposed","Healthy Control")) +
  labs(x = "Time", y = "Resting State Functional Connectivity", color = "Group") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  ggtitle("Right Putamen - Visual Network") +
  theme_minimal() 

# no stim
combined_df_no_stim <- combined_df[combined_df$group!="ADHD_stimulant",]

anova_result <- aov(fopa_cdelh ~ time * group + sex + ADHD.Problems + FD.under + meanFD + IQ, data = combined_df_no_stim)
summary(anova_result)

anova_result <- aov(fopa_ptlh ~ time * group + sex + ADHD.Problems + FD.under + meanFD + IQ, data = combined_df_no_stim)
summary(anova_result)

anova_result <- aov(vs_ptrh ~ time * group + sex + ADHD.Problems + FD.under + meanFD + IQ, data = combined_df_no_stim)
summary(anova_result)

# plot results
ggplot(combined_df_no_stim, aes(x = time, y = fopa_cdelh, color = group)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.5) +  
  geom_line() +
  geom_line(aes(group = group, y = fopa_ptlh), stat = "summary", size = 0.5, linetype="dashed") +
  scale_color_manual(values=c("orange", "darkblue"), labels = c("ADHD\nStimulant Naive","Healthy Control")) +
  labs(x = "Time", y = "Resting State Functional Connectivity", color = "Group") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  ggtitle("Left Caudate - Frontoparietal Network") +
  theme_minimal() 

ggplot(combined_df_no_stim, aes(x = time, y = fopa_ptlh, color = group)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.5) +  
  geom_line() +
  geom_line(aes(group = group, y = fopa_ptlh), stat = "summary", size = 0.5, linetype="dashed") +
  scale_color_manual(values=c("orange", "darkblue"), labels = c("ADHD\nStimulant Naive","Healthy Control")) +
  labs(x = "Time", y = "Resting State Functional Connectivity", color = "Group") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  ggtitle("Left Putamen - Frontoparietal Network") +
  theme_minimal() 

ggplot(combined_df_no_stim, aes(x = time, y = vs_ptrh, color = group)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.5) +  
  geom_line() +
  geom_line(aes(group = group, y = fopa_ptlh), stat = "summary", size = 0.5, linetype="dashed") +
  scale_color_manual(values=c("orange", "darkblue"), labels = c("ADHD\nStimulant Naive","Healthy Control")) +
  labs(x = "Time", y = "Resting State Functional Connectivity", color = "Group") +
  scale_x_discrete(labels = c("Baseline","Year 2")) +
  ggtitle("Right Putamen - Visual Network") +
  theme_minimal() 

# t-tests
combined_df_0 <- combined_df[combined_df$time=="0",]
combined_df_2 <- combined_df[combined_df$time=="2",]

# for baseline
t.test(combined_df_0$fopa_ptlh[combined_df_0$group=='TD'],combined_df_0$fopa_ptlh[combined_df_0$group=='ADHD_stimulant'],paired=FALSE)
t.test(combined_df_0$fopa_ptlh[combined_df_0$group=='TD'],combined_df_0$fopa_ptlh[combined_df_0$group=='ADHD_naive'],paired=FALSE)
t.test(combined_df_0$fopa_ptlh[combined_df_0$group=='ADHD_stimulant'],combined_df_0$fopa_ptlh[combined_df_0$group=='ADHD_naive'],paired=FALSE)

# for Y2
t.test(combined_df_2$fopa_ptlh[combined_df_2$group=='TD'],combined_df_2$fopa_ptlh[combined_df_2$group=='ADHD_stimulant'],paired=FALSE)
t.test(combined_df_2$fopa_ptlh[combined_df_2$group=='TD'],combined_df_2$fopa_ptlh[combined_df_2$group=='ADHD_naive'],paired=FALSE)
t.test(combined_df_2$fopa_ptlh[combined_df_2$group=='ADHD_stimulant'],combined_df_2$fopa_ptlh[combined_df_2$group=='ADHD_naive'],paired=FALSE)

# bayesian analyses

require(brms)

combined_df_wide <- reshape(as.data.frame(combined_df), idvar = "ID", timevar = "time", direction = "wide", drop = FALSE)
combined_df_wide$group.0 <- as.factor(combined_df_wide$group.0)
combined_df_wide$group.2 <- as.factor(combined_df_wide$group.2)

combined_df_wide$fopa_ptlh_change <- combined_df_wide$fopa_ptlh.2 - combined_df_wide$fopa_ptlh.0
combined_df_wide$vs_ptrh_change <- combined_df_wide$vs_ptrh.2 - combined_df_wide$vs_ptrh.0

combined_df_wide_no_naive <- combined_df_wide[combined_df_wide$group.0!="ADHD_naive",]

slope.prior.all <- c(
  prior(normal(0,10), class=b, coef=meanFD.0),
  prior(normal(0,10), class=b, coef=meanFD.2),
  prior(normal(0,10), class=b, coef=FD.under.0),
  prior(normal(0,10), class=b, coef=FD.under.2),
  prior(normal(0,10), class=b, coef=sex.0M),
  prior(normal(0,10), class=b, coef=IQ.0),
  prior(normal(0,10), class=b, coef=IQ.2),
  prior(normal(0,10), class=b, coef=ADHD.Problems.0),
  prior(normal(0,10), class=b, coef=ADHD.Problems.2),
  prior(normal(0,10), class=b, coef=fopa_ptlh_change),
  prior(normal(0,10), class=b, coef=vs_ptrh_change)
)

model_no_naive<-brm(as.numeric(group.0) ~ 
               meanFD.0 + meanFD.2 +
               FD.under.0 + FD.under.2 + sex.0 +
               IQ.0 + IQ.2 +
               ADHD.Problems.0 + ADHD.Problems.2 +
               fopa_ptlh_change + vs_ptrh_change,
             data=combined_df_wide_no_naive, family="bernoulli",           
             prior = slope.prior.all,
             silent = FALSE,
             warmup = 2000, iter = 10000)



