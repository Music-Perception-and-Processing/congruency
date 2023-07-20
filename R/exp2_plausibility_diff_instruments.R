library(lme4)
library(lmerTest)

# read plausibility data for Exp. 2
df <- read.csv("/Users/simon/congruency/data/exp2_plausibility.csv")

# set reference values for treatment coding
df$instrument <- factor(df$instrument,
                        levels = c("Violin", "Alto voice", "Clarinet", "Tuba"))
df$instrument <- relevel(df$instrument,
                         ref = "Violin")

contrasts(df$instrument) <- contr.sum(levels(df$instrument))

model <- lmer(response ~ instrument + (1 | participants),
              data = df)

row_names <- c("Intercept", "Violin", "Alto voice", "Clarinet")

estimates_table(model, row_names)
