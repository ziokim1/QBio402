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
  segs <- dpseg(x = c(1:60), y = df2[,i], minl = 10,P = 0.0001); segs # P: break-point penalty, increase to get longer segments with lower scores; minl: minimal segment length
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

#####



############

f_path <- file.path("SEAP_experimental_layout.xlsx")

if (!file.exists(f_path)) {
  f_path <- file.choose()
}

desired_sheet <- 1
temp_read <- readxl::read_excel(f_path,sheet = desired_sheet)

skip_rows1 <- NULL
col_skip <- 0
search_string <- "Untreated"

max_cols_to_search <- ncol(temp_read)
max_rows_to_search <- nrow(temp_read)

# Note, for the - 0, you may need to add/subtract a row if you end up skipping too far later.
while (length(skip_rows1) == 0) {
  col_skip <- col_skip + 1
  if (col_skip == max_cols_to_search) break
  skip_rows1 <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string)) - 0
  
}

skip_rows2 <- NULL
col_skip <- 0
search_string2 <- "Treated"

max_cols_to_search <- ncol(temp_read)
max_rows_to_search <- nrow(temp_read)

# Note, for the - 0, you may need to add/subtract a row if you end up skipping too far later.
while (length(skip_rows2) == 0) {
  col_skip <- col_skip + 1
  if (col_skip == max_cols_to_search) break
  skip_rows2 <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string2)) - 0
  
}

skip_rows <- min(c(skip_rows1[1], skip_rows2[1]))

max_rows <- NULL
col_skip <- 0
search_string3 <- "Steps for automated analysis"

while (length(max_rows) == 0) {
  col_skip <- col_skip + 1
  if (col_skip == max_cols_to_search) break
  max_rows <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string3)) - 0
  
}


# ... now we re-read from the known good starting point.
real_data <- readxl::read_excel(
  f_path,
  sheet = desired_sheet,
  skip = skip_rows-1,
  n_max = max_rows-skip_rows-2,
  row_names = TRUE
)

df_t <- real_data[,-1]
cv <- as.list(c(t(df_t)))

for (i in )res

