library(lme4)
library(lmerTest)
library(ggplot2)

# load plausibility data
plausibility <- read.csv("/Users/simon/congruency/data/exp2_plausibility.csv")

# retrieve instrument names
instruments <- as.list(unique(plausibility$instrument))

# retrieve condition names
conditions <- as.list(unique(plausibility$condition))

# retrieve register names
registers <- as.list(unique(plausibility$register))

# loop over instruments
for (i in instruments) {
  # create subset for instrument
  df <- subset(plausibility, instrument == i)
  
  # aggregate df for plotting
  agg_df <- aggregate(response ~ register + condition,
                      data = df,
                      FUN = mean)
  print(ggplot(agg_df, aes(x=factor(register),
                           y=response,
                           group=condition,
                           color=condition)) + geom_line() + ylim(1, 6))
  
  # plot <- ggplot(df, aes(x = msi_score, y = response, group = participants)) + 
  #   geom_line(aes(color = participants), alpha = 0.5) + 
  #   stat_summary(aes(group = interaction(participants, instrument)), fun = mean,
  #                geom = "point", color = "black", size = 3) + 
  #   facet_wrap(~ instrument, scales = "free_y") +
  #   labs(x = "MSI score", y = "Mean Response")
  # print(plot)
  
  # select only data within range of instrument
  #df <- subset(df, withinRange == 1)
  
  # reorder levels and set references
  # df$condition <- factor(df$condition,
  # levels = c("LRSE", "MRSE", "HRSE", "congr."))
  # df$condition <- relevel(df$condition, ref = "congr.")
  
  # deviation (effects) coding
  df$condition <- factor(df$condition,
                         levels = c("LRSE", "MRSE", "HRSE", "congr."),
                         labels = c("LRSE", "MRSE", "HRSE", "congr."))
  df$condition <- relevel(df$condition,
                          ref = "LRSE")
  
  df$register <- factor(df$register,
                        levels = c("low", "mid", "high"))
  # df$register <- relevel(df$register,
  #                        ref = "high")
  
  contrasts(df$condition) <- contr.sum(levels(df$condition))
  contrasts(df$register) <- contr.sum(levels(df$register))
  
  # create LME
  model <- lmer(response ~ register * condition + (1 | participants),
                data = df)
  
  # define row names
  row_names <- c("Intercept", "low-mid", "high-mid", "LRSE-congr.", "MRSE-congr.",
                 "HRSE-congr.", "low-mid:LRSE-congr.", "high-mid:LRSE-congr.",
                 "low-mid:MRSE-congr.", "high-mid:MRSE-congr.", "low-mid:HRSE-congr.",
                 "high-mid:HRSE-congr.")
  # row_names <- c("Intercept", "low-mid", "high-mid", "LRSE-congr.", "MRSE-congr.",
  #                "HRSE-congr.")
  
  print(i)
  print(mean(df$response))
  print(sd(df$response))
  estimates_table(model, row_names)
  print(model_performance(model))
  
}
