library("tidyverse")
library("ggplot2")
library("readxl")
library("gridExtra")
library("dpseg")


df <- read_xlsx("SEAP_raw_data.xlsx") # read excel file

r_names <- df[as.numeric(which(df == "Well")+1):nrow(df),2] # extract row names
c_names <- df[as.numeric(which(df == "Well")),3:ncol(df)] # extract column names

df2 <- df[as.numeric(which(df == "Well")+1):nrow(df),3:as.numeric(ncol(df))] # index out the dataframe with data only

rownames(df2) <- t(r_names) # rename row & column
colnames(df2) <- c(1:96)

df2 <- as.data.frame(lapply(df2,as.numeric)) # make it all numeric


res = list()
for (i in 1:96) {
  segs <- dpseg(x = c(1:60), y = df2[,i], minl = 10,P = 0.0001); segs
  s <- max(subset(segs$segments, var < 2000000)$slope)
  if (s < 0){
    res[[i]] <- NA
  } else {
    res[[i]] <- (s / (18600 * 0.5)) * (200 / 40) * 10^6
  }
  
}

res

######
#Try dpseg in the dpseg package. We restrict the minimum length to 50 to avoid short linear stretches which may occur by chance. There are other tuning parameters available. See ?dpseg and the vignette that comes with the package for more information.
#To make the input reproducible we need to use set.seed and have done this in the Note at the end.
##########

?dpseg

i = 61
segs <- dpseg(x = c(1:60), y = df2[,i], minl = 10,P = 0.0001); segs
## ... this output is shown just before the image ...
subset(segs$segments, var < 2000000)
##    x1  x2 start end intercept    slope        r2      var
## 3 116 203   116 203  1458.242 10.15865 0.8613225 10844.35
plot(segs)
max(res$slope)
