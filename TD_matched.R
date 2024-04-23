


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
setwd("/Users/adamkaminski/Downloads")
netdat <- read.csv('mrirscor02_v1.csv')

# filter based on included subjects up to this point
netdat_filt <- netdat[netdat$ID%in%cbcl_wide_filtered3$A.ID,]

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

# define net codes
if (combine == 1) {
  net_codes <- c("_au_","_cercsa_","_copa_","_df_","_dsa_","_fopa_","_rst_",
                 "_smhsmm_","_vta_","vs")
} else {
  net_codes <- c("_au_","_cerc_","_copa_","_df_","_dsa_","_fopa_","_rst_",
                 "_sa_","_smh_","_smm_","_vta_","_vs_")
}

# calculate FC difference
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
if (combine == 1) {netdat <- netdat[,c(1:10,59:67)]}
netdat <- cbind(netdat,diff)

# remove unnecessary cols (all FC besides sign. from original analysis for comparison)
netdat_TD_orig_cxs <- netdat[,c(1:10,36,57,59:67,93,114)]

# make df long
netdat_TD_orig_cxs_long <- netdat_TD_orig_cxs %>%
  pivot_longer(cols = -ID, 
               names_to = c(".value", "time"), 
               names_pattern = "(.*?)\\.(\\d+)$")

# plot FC change
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

# t tests for rs-FC change
netdat_TD_orig_cxs_long_0 <- netdat_TD_orig_cxs_long[netdat_TD_orig_cxs_long$time=="0",]
netdat_TD_orig_cxs_long_2 <- netdat_TD_orig_cxs_long[netdat_TD_orig_cxs_long$time=="2",]

t.test(netdat_TD_orig_cxs_long_0$rsfmri_cor_ngd_fopa_scs_ptlh,netdat_TD_orig_cxs_long_2$rsfmri_cor_ngd_fopa_scs_ptlh,paired=TRUE)
t.test(netdat_TD_orig_cxs_long_0$rsfmri_cor_ngd_vs_scs_ptrh,netdat_TD_orig_cxs_long_2$rsfmri_cor_ngd_vs_scs_ptrh,paired=TRUE)






