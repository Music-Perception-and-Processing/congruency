library(lme4)
library(lmerTest)
library(ggplot2)
library(MuMIn)
library(boot)
library(xtable)
library(performance)

# define directories
script_dir <- dirname(sys.frame(1)$ofile)
main_folder <- file.path(script_dir, "..")

# load functions
if(!exists("estimates_table")) {
  
  source(file.path(main_folder, "R", "estimates_table.R"))
}

# load pleasantness data
pleasantness <- read.csv(file.path(main_folder, "data", "compiled_responses", "exp2_pleasantness.csv"))

# retrieve instrument names
instruments <- as.list(unique(pleasantness$instrument))

# retrieve condition names
conditions <- as.list(unique(pleasantness$condition))

# loop over instruments
for (i in instruments) {
  # create subset for instrument
  df <- subset(pleasantness, instrument == i)
  
  # find existing pitch values
  idx <- which(is.na(df$response))
  missing_pitches <- unique(df$noteIdx[idx])
  min_pitch <- min(missing_pitches) - 1
  df1 <- subset(df, ! noteIdx %in% missing_pitches)
  df1$noteIdx <- df1$noteIdx - min_pitch
  
  # aggregate df for plotting
  agg_df <- aggregate(response ~ noteIdx + condition,
                      data = df,
                      FUN = mean)
  print(ggplot(agg_df, aes(x=factor(noteIdx),
                     y=response,
                     group=condition,
                     color=condition)) + geom_line() + ylim(1, 6))
  
  # select only data within range of instrument
  #df <- subset(df, withinRange == 1)
  
  # reorder levels and set references
  # df$condition <- factor(df$condition,
                           # levels = c("LRSE", "MRSE", "HRSE", "congr."))
  # df$condition <- relevel(df$condition, ref = "congr.")
  
  # deviation (effects) coding
  df$condition <- factor(df$condition,
                         levels = c("congr.", "LRSE", "MRSE", "HRSE"))
  # df$condition <- relevel(df$condition, ref = "congr.")
  contrasts(df$condition) <- contr.sum(levels(df$condition))
  
  df$withinRange <- factor(df$withinRange,
                           levels = c(1, 0),
                           labels = c("inRange", "outRange"))
  df$withinRange <- relevel(df$withinRange, ref = "inRange")
  contrasts(df$withinRange) <- contr.treatment(levels(df$withinRange))
  
  # F0 effects
  # F0 <- I((df$noteIdx-1)/3)
  F0 <- I((df$noteIdx-10)/3)
  F02 <- I(F0^2)
  
  # create LME
  model <- lmer(response ~ condition * F0 + F02 + (1 | participants),
                data = df)
  
  row_names <- c("Intercept", "LRSE-congr.", "MRSE-congr.", "HRSE-congr.",
                 "F0", "F02", "F0:LRSE-congr.", "F0:MRSE-congr.", "F0:HRSE-congr.")
  
  # # compute bootstrapped means
  # statfun <- function(data, num_obs, i){
  #   d <- data[i, ]
  #   c(
  #     mean(d, na.rm=TRUE),
  #     sd(d, na.rm=TRUE)
  #   )
  # }
  # 
  # for (j in conditions) {
  #   cond_df <- subset(df, condition == j)
  #   num_obs <- length(unique(cond_df$noteIdx[which(!is.na(cond_df$response))]))
  #   
  #   boot_stat <- boot(cond_df[, "response", drop = FALSE], statfun, R = 1000)
  #   print(boot_stat)
  # }
  # boot_mean <- boot(df[, "response", drop = FALSE], meanfun, R = 1000)
  
  print(i)
  mus <- aggregate(df$response, list(df$condition), FUN=mean, na.rm=TRUE)
  print(mus)
  print(mean(mus[,2]))
  print(sd(mus[,2]))
  estimates_table(model, row_names)
  print(model_performance(model))
  # print(r.squaredGLMM(model))
  # print(anova(model))
  # print(summary(model))
  
  # ######## df1
  # 
  # # conding sceme df1
  # df1$condition <- factor(df1$condition,
  #                         levels = c("LRSE", "MRSE", "HRSE", "congr."),
  #                         labels = c("LRSE", "MRSE", "HRSE", "congr."))
  # df1$condition <- relevel(df1$condition, ref = "congr.")
  # contrasts(df1$condition) <- contr.treatment(levels(df1$condition))
  # 
  # F0 <- I((df1$noteIdx-1)/3)
  # F02 <- I(F0^2)
  # 
  # # create LME
  # model <- lmer(response ~ condition * F0 + F02 + (1 | participants),
  #               data = df1)
  # 
  # row_names <- c("Intercept", "LRSE-congr.", "MRSE-congr.", "HRSE-congr.",
  #                "F0", "F02", "F0:LRSE-congr.", "F0:MRSE-congr.", "F0:HRSE-congr.")
  # 
  # print(i)
  # print(aggregate(df1$response, list(df1$condition), FUN=mean, na.rm=TRUE))
  # estimates_table(model, row_names)
  # print(anova(model))
  
}
