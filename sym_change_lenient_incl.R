


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

writeLines(paste(dim(ksads)[1],"/",origlen," (", round(origlen/dim(ksads)[1]),"%) of people have ADHD at baseline according to KSADS",sep=""))

