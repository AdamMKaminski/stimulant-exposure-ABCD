
require(tidyr)
require(ggplot2)
require(cowplot)
require(dplyr)

#
#
# Obtain IDs for 1. baseline ADHD (KSADS), and 2. has CBCL data at baseline and Y2
#
#

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
setwd("/Users/adamkaminski/Desktop/stimulant-exposure-ABCD")
cbcl <- read.csv('abcd_cbcls01.csv')
cbcl <- cbcl[cbcl$ID%in%baseline_adhd_IDs,]

cbcl_baseline <- cbcl[cbcl$time=="baseline_year_1_arm_1",]
cbcl_1 <- cbcl[cbcl$time=="1_year_follow_up_y_arm_1",]
cbcl_2 <- cbcl[cbcl$time=="2_year_follow_up_y_arm_1",]

full_IDs <- cbcl_baseline$ID[cbcl_baseline$ID%in%cbcl_2$ID]
writeLines(paste(length(full_IDs),"/",dim(cbcl_baseline)[1]," (", (length(full_IDs)/dim(cbcl_baseline)[1])*100,"%) of people who have CBCL baseline data also have CBCL Y2 data",sep=""))

# Exclude subjects whose med data was previously decoded
orig_dat <- read.csv('stimulant_master_df.csv')
orig_IDs <- orig_dat$ID

full_IDs_new <- full_IDs[!(full_IDs%in%orig_IDs)]

# Get Medication data, save csv to be decoded
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
med <- read.csv('medsy01.csv')

med_subs <- med[med$The.NDAR.Global.Unique.Identifier..GUID..for.research.subject%in%full_IDs_new,]
write.csv(med_subs,"medsy01_moreIDs.csv",na = "",row.names = FALSE)

#
#
# Following med decoding...
#
#

# get stimulant exposure for new IDs
setwd("/Users/adamkaminski/Desktop/stimproj/data/")
new_exp <- read.csv('new_subs.csv')
new_exp <- new_exp[,-4]

# get stimulant exposure for old IDs
old_dat <- read.csv('stimulant_master_df.csv')
old_exp <- as.data.frame(matrix(0,179,3))
colnames(old_exp) <- c('sub_id','stim_exposed','combined_other')
old_exp$sub_id <- old_dat$ID
old_exp$stim_exposed <- old_dat$stim_group_bi
old_exp$combined_other <- old_dat$nonstim_med_bi

# combine old and new IDs
total_exp <- rbind(new_exp,old_exp)


#
#
# For obtained IDs, look at symptom change (DSM-oriented CBCL scales)
#
#

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

cbcl <- cbcl[cbcl$ID%in%full_IDs,]
cbcl$time[cbcl$time=='baseline_year_1_arm_1'] <- 0
cbcl$time[cbcl$time=='1_year_follow_up_y_arm_1'] <- 1
cbcl$time[cbcl$time=='2_year_follow_up_y_arm_1'] <- 2

cbcl <- cbcl[,c(1,2,3,4,16:24)]
colnames(cbcl) <- c('ID','age','sex','time','Depress','AnxDisord','SomaticPr','ADHD','Opposit','Conduct','Sluggish','Obsessive','Stress')

# for wide
cbcl_wide <- reshape(cbcl, idvar = "ID", timevar = "time", direction = "wide")

# add stimulant exposure
cbcl_wide_sorted <- cbcl_wide[order(cbcl_wide$ID),]
total_exp_sorted <- total_exp[order(total_exp$sub_id),]

if(all(cbcl_wide_sorted$ID==total_exp_sorted$sub_id)){
  cbcl_wide_sorted_exp <- cbind(cbcl_wide_sorted,total_exp_sorted)
}

# function for returning long
make_long <- function(wide_dat) {
  long_dat <- as.data.frame(pivot_longer(wide_dat, 
                            cols = matches("\\.\\d+$"),
                            names_to = c(".value", "Time"),  
                            names_pattern = "(\\w+)\\.(\\d+)", 
                            values_to = "Value"))
  return(long_dat)
}
  
cbcl_long_sorted_exp <- make_long(cbcl_wide_sorted_exp)

# remove Y1 rows
#cbcl_long_no1 <- cbcl_long[cbcl_long$Time!="1",]
#cbcl_long_sorted_exp <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time!="1",]

# Function for visualizing symptom change from baseline to Y2
plot_symptom_change <- function(df,col,title,grouped) { 
  
  if (grouped==0) {
    fig <- ggplot(df, aes(x=Time, y=.data[[col]], color=Time)) +
      geom_line(aes(x = Time, group = df$ID), 
                size = 0.5, color = 'gray') + 
      geom_point(aes(x = Time), size = 2, position = position_dodge(width = 0.75)) + 
      geom_boxplot() +
      scale_color_manual(values=c("gray", "black"), labels = c("Baseline", "Year 2")) +
      xlab("Group") + ylab("T-Score") +
      scale_x_discrete(labels = c("Baseline","Year 2")) +
      theme(axis.text.x = element_blank(),axis.ticks.margin=unit(0,'cm')) +
      guides(color=guide_legend(title="Time")) + 
      ggtitle(title) +
      theme_minimal_grid(12)
  } else if (grouped==1) {
    fig <- ggplot(df, aes(x=interaction(Time,stim_exposed), y=.data[[col]], color=Time)) +
      geom_line(aes(x = interaction(Time,stim_exposed), group = df$ID), 
                size = 0.5, color = 'gray') + 
      geom_point(aes(x = interaction(Time,stim_exposed)), size = 2, position = position_dodge(width = 0.75)) + 
      geom_boxplot() +
      scale_color_manual(values=c("gray", "black"), labels = c("Baseline", "Year 2")) +
      xlab("Group") + ylab("T-score") +
      scale_x_discrete(labels = c("Naive","Naive","Exposed","Exposed")) +
      theme(axis.text.x = element_blank(),axis.ticks.margin=unit(0,'cm')) +
      guides(color=guide_legend(title="Time")) + 
      ggtitle(title) + 
      theme_minimal_grid(12)
  }
  
  print(fig)
}

# Loop through DSM-oriented scales and plot symptom change
dat <- cbcl_long_sorted_exp
cols <- colnames(dat)[8:16]
par(mfrow = c(1, 8))
for(col in cols){
  dev.new()
  plot_symptom_change(dat,col,col,1)
}

# Just plot ADHD symptom change
ggplot(cbcl_long_sorted_exp, aes(x=interaction(as.factor(Time),as.factor(stim_exposed)), y=ADHD, color=as.factor(stim_exposed))) +
  geom_line(aes(x = interaction(as.factor(Time),as.factor(stim_exposed)), group = ID), 
            size = 0.5, color = 'gray') + 
  geom_point(aes(x = interaction(as.factor(Time),as.factor(stim_exposed))), size = 2, position = position_dodge(width = 0.75)) + 
  geom_boxplot() +
  scale_color_manual(values=c("gray", "black"), labels = c("Stimulant Naive", "Stimulant Exposed")) +
  xlab("Time") + ylab("ADHD Problems (t-score)") +
  scale_x_discrete(labels = c("Baseline","Y1","Y2","Baseline","Y1","Y2")) +
  theme(axis.text.x = element_blank(),axis.ticks.margin=unit(0,'cm')) +
  guides(color=guide_legend(title="Group")) + 
  ggtitle("Change in ADHD Symptoms with More Lenient Inclusion") + 
  theme_minimal_grid(12)

# supp material anova
#cbcl_long_sorted_exp_no1 <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time!="1",]
result <- aov(ADHD ~ stim_exposed*Time + Error(sub_id/Time), data = cbcl_long_sorted_exp)
summary(result)

cbcl_long_sorted_exp_0 <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time==0,]
cbcl_long_sorted_exp_1 <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time==1,]
cbcl_long_sorted_exp_2 <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time==2,]

t.test(cbcl_long_sorted_exp_0$ADHD~cbcl_long_sorted_exp_0$stim_exposed,paired=FALSE)
t.test(cbcl_long_sorted_exp_1$ADHD~cbcl_long_sorted_exp_1$stim_exposed,paired=FALSE)
t.test(cbcl_long_sorted_exp_2$ADHD~cbcl_long_sorted_exp_2$stim_exposed,paired=FALSE)

cbcl_long_sorted_exp_stim <- cbcl_long_sorted_exp_no1[cbcl_long_sorted_exp_no1$stim_exposed==1,]
cbcl_long_sorted_exp_naive <- cbcl_long_sorted_exp_no1[cbcl_long_sorted_exp_no1$stim_exposed==0,]

t.test(cbcl_long_sorted_exp_stim$ADHD~cbcl_long_sorted_exp_stim$Time,paired=TRUE)
t.test(cbcl_long_sorted_exp_naive$ADHD~cbcl_long_sorted_exp_naive$Time,paired=TRUE)

# ANOVAs
result <- aov(Depress ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(AnxDisord ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(SomaticPr ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(ADHD ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(Opposit ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(Conduct ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(Sluggish ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(Obsessive ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

result <- aov(Stress ~ stim_exposed*Time, data = cbcl_long_sorted_exp)
summary(result)

# t tests for ADHD symptoms

cbcl_long_sorted_exp_0 <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time=="0",]
cbcl_long_sorted_exp_2 <- cbcl_long_sorted_exp[cbcl_long_sorted_exp$Time=="2",]

t.test(cbcl_long_sorted_exp_0$ADHD~cbcl_long_sorted_exp_0$stim_exposed)
t.test(cbcl_long_sorted_exp_2$ADHD~cbcl_long_sorted_exp_2$stim_exposed)

# Cohen's d

#
# Baseline
#

mean1 <- mean(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_0$stim_exposed==1])
mean2 <- mean(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_0$stim_exposed==0])
sd1 <- sd(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_0$stim_exposed==1])
sd2 <- sd(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_0$stim_exposed==0])
n1 <- length(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_0$stim_exposed==1])
n2 <- length(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_0$stim_exposed==0])

sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

cohens_d_baseline <- (mean1 - mean2) / sp

#
# Y2
#

mean1 <- mean(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
mean2 <- mean(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
sd1 <- sd(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
sd2 <- sd(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
n1 <- length(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
n2 <- length(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])

sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

cohens_d_y2 <- (mean1 - mean2) / sp

#
# Paired
#

#
# Stim
#

mean1 <- mean(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
mean2 <- mean(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
sd1 <- sd(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
sd2 <- sd(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
n1 <- length(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])
n2 <- length(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_2$stim_exposed==1])

sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

cohens_d_stim <- (mean1 - mean2) / sp

#
# No stim
#

mean1 <- mean(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
mean2 <- mean(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
sd1 <- sd(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
sd2 <- sd(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
n1 <- length(cbcl_long_sorted_exp_2$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])
n2 <- length(cbcl_long_sorted_exp_0$ADHD[cbcl_long_sorted_exp_2$stim_exposed==0])

sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

cohens_d_nostim <- (mean1 - mean2) / sp

