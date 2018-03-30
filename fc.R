library(readr)
library(rstudioapi)

rm(list=ls())

setwd(dirname(getSourceEditorContext()$path)) # Set working directory

df_exp  <- as.data.frame(read_csv("data/GSE23561.csv", col_names = TRUE)) # gene expressions
df_states  <- as.data.frame(read_csv("data/states.csv", col_names = TRUE)) # disease states
df_platform  <- as.data.frame(read_csv("data/GPL10775.csv", col_names = TRUE)) # platform

df_platform <- df_platform[, c("ID", "Symbol v12")] # remove unnecessary columns
df_exp <- df_exp[, -1] # remove id_ref column
df_exp <- cbind("Symbol v12" = df_platform[, "Symbol v12"], df_exp) # add the symbols column

# CONVERT EACH INVALID VALUE TO NA:
for(row in 1:nrow(df_exp)) {
  val <- df_exp[row, "Symbol v12"]
  if(!is.null(val) && !is.na(val))
    if(val == "." || val == "")
      df_exp[row, "Symbol v12"] <- NA
}

# REMOVE NA VALUES:
df_exp <- na.omit(df_exp)

# GROUP ROWS BY GENE SYMBOLS:
df_exp <- aggregate(df_exp[, -1], list(df_exp[, 1]), mean)
symbols <- unlist(df_exp[, 1])
df_exp <- round(df_exp[, -1], 2)
row.names(df_exp) <- symbols

# GROUP COLUMNS BY DISEASE STATES:
df_group <- NULL
df_group <- data.frame(Control = rowMeans(df_exp[ ,1:9]))                 # 9 Control
df_group <- cbind(df_group, data.frame(RA = rowMeans(df_exp[ ,10:15])))   # 6 RA
df_group <- cbind(df_group, data.frame(MetS = rowMeans(df_exp[ ,16:21]))) # 6 MetS
df_group <- cbind(df_group, data.frame(CAD = rowMeans(df_exp[ ,22:27])))  # 6 CAD
df_group <- cbind(df_group, data.frame(T2D = rowMeans(df_exp[ ,28:35])))  # 8 T2D
df_group <- round(df_group, 2)
row.names(df_group) <- symbols

################### FOLD CHANGE ###########################

# If there is any attribute which is differentially 
# expressed as against the target attribute, keep that row.
fc_on_target <- function(df, target_index, cutoff) {

  df_temp <- df
  
  if(target_index == 1)
    range <- 2:ncol(df_temp)
  else if(target_index == ncol(df_temp)) 
      range <- 1:(ncol(df_temp)-1)
  else
    range <- c(1:(target_index-1),(target_index+1):ncol(df_temp))
  
  for(row in 1:nrow(df_temp)) {
    flag <- FALSE
    
    for(pair_index in range) {
      fc <- abs(log(df_temp[row, target_index]/df_temp[row, pair_index], base = 2))

      if(fc >= cutoff) {
        flag <- TRUE
        break
      }
    }
    
    # Remove the gene unless it is a DEG for any disease group
    if(flag)
      next
    else
      df_temp[row, target_index] <- NA
  }
  
  df_temp <- na.omit(df_temp)
  return(df_temp)
}

df_deg_tg <- fc_on_target(df_group, 1, 2) # TEST: DEGs as against CONTROL group

# If there are any pair of attributes in the target range which are differentially 
# expressed as against each other, keep that row. 
fc_pairwise <- function(df, target_range, cutoff) {

  df_temp <- df

  for(row in 1:nrow(df_temp)) {
    flag <- FALSE
    last_index = length(target_range)
    
    for(curr_index in 1:last_index) {

      last_col = target_range[last_index] # last column for the target_range
      curr_col = target_range[curr_index]
      
      if(curr_col == last_col)
        break
      
      for(next_index in (curr_index + 1):last_index) {
        next_col = target_range[next_index]
        fc <- abs(log(df_temp[row, curr_col]/df_temp[row, next_col], base = 2))
        
        if(fc >= cutoff) {
          flag <- TRUE
          break
        }
      }
      if(flag)
        break
    }
    
    # Remove the gene unless it is a DEG for any disease group
    if(flag)
      next
    else
      df_temp[row, 1] <- NA
      # row.names(df_temp)[row] <- paste(row.names(df_temp)[row],"***") # TO TEST
  }
  
  df_temp <- na.omit(df_temp)
  return(df_temp)
}

df_deg_pw <- fc_pairwise(df_group, c(2:5), 2) # TEST: DEGs as against each other


