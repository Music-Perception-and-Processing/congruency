library(lme4)
library(lmerTest)
library(xtable)

# define directories
script_dir <- dirname(sys.frame(1)$ofile)
main_folder <- file.path(script_dir, "..")

# load functions
if(!exists("estimates_table")) {
  
  source(file.path(main_folder, "R", "estimates_table.R"))
}

# load pleasantness data
df = read.csv(file.path(main_folder, "data", "compiled_responses", "exp1_pleasantness.csv"))

# specify custom order of levels and set reference level
df$space <- factor(df$space,
                   levels = c("inner", "outer"))
df$space <- relevel(df$space, ref = "outer")
df$condition <- factor(df$condition,
                       levels = c("congruent", "incongruent", "fixed"))
df$condition <- relevel(df$condition, ref = "congruent")

contrasts(df$condition) <- contr.sum(levels(df$condition))
contrasts(df$space) <- contr.sum(levels(df$space))

# F0 effects
F0 <- I(log2(df$f0))
F02 <- I(F0^2)

# setup linear mixed-effects model
model <- lmer(score ~ space * condition + F0 + F02 +  (1 + space + condition | participant), 
              data = df, REML = TRUE)

row_names <- c("Intercept", "space", "congr.", "incongr.", "F0", "F02",
               "space:congr.", "space:incongr.")
estimates_table(model, row_names)
