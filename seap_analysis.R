# Load packages
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("readxl")
library("dpseg")

###############################################################################----
# Function for reading the machine data, then converting it into SEAP conc. in U/L

seap_calc <- function(filename,
                     search_string1="A01", # starting position define by sample labelling
                     search_string2="59 min", # ending position defined by time
                     min_l, # minimal segment length
                     break_pt_pen, # increase to get longer segments with lower scores
                     ext_coeff, # extinction coefficient
                     light_path, # light path in 96-well plate
                     cell_vol, # cell culture supernatant volume added to the well
                     inact_sup_vol, # inactivated supernatant / dilution
                     desired_sheet = 1){
  
  f_path <- file.path(filename)
  temp_read <- readxl::read_excel(f_path,sheet = desired_sheet)
  
  # Find the starting position by searching for the search_string1
  skip_rows <- NULL
  col_skip <- 0
  max_cols_to_search <- ncol(temp_read)
  max_rows_to_search <- nrow(temp_read)
  while (length(skip_rows) == 0) {
    col_skip <- col_skip + 1
    if (col_skip == max_cols_to_search) break
    skip_rows <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string1))
  }
  
  # Find the ending position by searching for the search_string2
  max_rows <- NULL
  col_skip <- 0
  while (length(max_rows) == 0) {
    col_skip <- col_skip + 1
    if (col_skip == max_cols_to_search) break
    max_rows <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string2))
  }
  
  # Now re-read from the known good starting point
  real_data <- readxl::read_excel(
    f_path,
    sheet = desired_sheet,
    skip = skip_rows,
    n_max = max_rows-skip_rows
  )
  df <- as.data.frame(lapply(real_data[,-c(1,2)],as.numeric)) # Dataframe read-out
  
  # Calculate the concentration from the raw data
  res = list()
  for (i in 1:ncol(df)) {
    segs <- dpseg(x = c(1:nrow(df)), y = df[,i], minl = min_l,P = 0.00001) 
    abs <- max(subset(segs$segments, var < 0.02)$slope) # setting variance threshold
    if (abs < 0){
      res[[i]] <- NA # for empty wells
    } else {
      res[[i]] <- (abs / (ext_coeff * light_path)) * (cell_vol / inact_sup_vol) * 10^6 # Beer Lambert Law
    }
  }
  return(as.numeric(c(t(res))))
}

conc <- seap_calc(filename = "SEAP_raw_data.xlsx",
                  min_l = 10,
                  break_pt_pen = 0.00001,
                  ext_coeff = 18600, # M^-1/cm^-1
                  light_path = 0.5, # cm
                  cell_vol = 200, # ul
                  inact_sup_vol = 40 # ul
                  )

###############################################################################----
# dpseg Tutorial: 

df_tut <- as.data.frame(lapply(readxl::read_excel(path = "SEAP_raw_data.xlsx",skip = 9,n_max = 60)[,-c(1,2)],as.numeric))
?dpseg # minl & P

# Restrict the minimum length to 10 to avoid short linear stretches which may occur by chance
# Add break point penalty to see one continuous line rather than many segmented

i = 96 # select the sample
segs <- dpseg(x = c(1:60), y = df_tut[,i], minl = 10, P=0.00001)
subset(segs$segments)
plot(segs) # Graphical visualisation
max(segs$segments$slope) # selecting the slope == dA_405/dt

###############################################################################----
# Function for reading the design matrix

design_matrix <- function(filename,
                      search_string1="A",
                      search_string2="H",
                      desired_sheet = 1){
  
  f_path <- file.path(filename)
  temp_read <- readxl::read_excel(f_path,sheet = desired_sheet)
  
  # Find the starting position by searching for the search_string1
  skip_rows <- NULL
  col_skip <- 0
  search_string1 <- "A"
  max_cols_to_search <- ncol(temp_read)
  max_rows_to_search <- nrow(temp_read)
  
  while (length(skip_rows) == 0) {
    col_skip <- col_skip + 1
    if (col_skip == max_cols_to_search) break
    skip_rows <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string1))
  }
  skip_row1 <- skip_rows[1]
  skip_row2 <- skip_rows[2]
  
  max_rows <- NULL
  col_skip <- 0
  search_string2 <- "H"
  
  while (length(max_rows) == 0) {
    col_skip <- col_skip + 1
    if (col_skip == max_cols_to_search) break
    max_rows <- which(stringr::str_detect(temp_read[1:max_rows_to_search,col_skip][[1]],search_string2))
  }
  max_row1 <- max_rows[1]
  max_row2 <- max_rows[2]
  
  # Re-read from the known good starting point & extract sample names and treatments
  samples <- readxl::read_excel(
    f_path,
    sheet = desired_sheet,
    skip = skip_row1-1,
    n_max = max_row1-skip_row1+1
  )
  samples <- samples[,-1]
  samples <- as.factor(c(t(samples)))
  
  trt <- readxl::read_excel(
    f_path,
    sheet = desired_sheet,
    skip = skip_row2-1,
    n_max = max_row2-skip_row2+1
  )
  trt <- trt[,-1]
  trt <- as.factor(c(t(trt)))
  
  return(data.frame(samples, trt))
}

###############################################################################----
# Building the final dataframe
foo <- design_matrix(filename = "SEAP_experimental_layout.xlsx")
df <- cbind(foo, conc)
df <- na.omit(df)

# Now graph with mean value, with error bar, dot points of individual values for each sample
ggbarplot(
  df, x = "samples", y = "conc", 
  add = c("mean_sd", "jitter"), 
  add.params = list(shape = "trt"),
  fill= "trt", palette = c("#807F7F", "#BF504D"),
  position = position_dodge(0.8),
  ylab = "[SEAP] (U/L)",
  xlab = "Samples",
  order = unique(df$samples)) + labs(fill = "Treatment", shape = "Treatment")

# Some samples show greater contrast between treatment (e.g. 5, 6, 8, 10)
# Others show only a small difference (e.g. 4, 9, 11)
