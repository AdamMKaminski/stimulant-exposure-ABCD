
require(tidyr)
require(ggplot2)

#
# 
#
#
#
# Obtain IDs for 1. baseline ADHD (KSADS), and 2. has CBCL data at baseline and Y2
#
#
#
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

#
# 
#
#
#
# For obtained IDs, look at symptom change (DSM-oriented CBCL scales)
#
#
#
#
#

cbcl <- cbcl[cbcl$ID%in%full_IDs,]
cbcl$time[cbcl$time=='0_baseline_year_1_arm_1'] <- '0'
cbcl$time[cbcl$time=='1_year_follow_up_y_arm_1'] <- '1'
cbcl$time[cbcl$time=='2_year_follow_up_y_arm_1'] <- '2'
colnames(cbcl) <- c('ID','age','sex','time','Depress','AnxDisord','SomaticPr','ADHD','Opposit','Conduct','Sluggish','Obsessive','Stress')

# for wide
cbcl_wide <- reshape(cbcl, idvar = "ID", timevar = "time", direction = "wide")
# for long
cbcl_long <- as.data.frame(pivot_longer(cbcl_wide, 
                        cols = matches("\\.\\d+$"),
                        names_to = c(".value", "Time"),  
                        names_pattern = "(\\w+)\\.(\\d+)", 
                        values_to = "Value"))

cbcl_long_no1 <- cbcl_long[(cbcl_long$Time!="1"),]


# Function for visualizing symptom change from baseline to Y2
plot_symptom_change <- function(df,col,title) { 
  
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
  
  print(fig)
}

# Loop through DSM-oriented scales and plot symptom change
cols <- colnames(cbcl_long_no1)[5:13]
par(mfrow = c(1, 8))
for(col in cols){
  dev.new()
  plot_symptom_change(cbcl_long_no1,col,col)
}

