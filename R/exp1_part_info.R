all_files <- list.files("/Users/simon/congruency/data/responses/exp_1_info")

# loop over participants
msi_scores <- vector("list",length(all_files))
n <- 1
for (i in all_files) {
  df <- read.csv(paste("/Users/simon/congruency/data/responses/exp_1_info/",i,
                       sep = ""),
                 skip = 3)
  start_row <- which(df$condition1 == "msiPerception5", arr.ind = TRUE)
  end_row <- start_row + 15
  
  msi_scores[[n]] <- df$responseCode[start_row:end_row]
  
  n <- n + 1
}

msi_scores <- do.call(rbind, msi_scores)
msi_scores[,c(3, 5, 8, 11)] <- 7 - msi_scores[,c(3, 5, 8, 11)] + 1
msi_perception <- rowSums(msi_scores[,1:9])
msi_training <- rowSums(msi_scores[,10:16])
df <- data.frame(msi_perception,msi_training)