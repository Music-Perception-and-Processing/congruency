library(lme4)
library(lmerTest)

# define directories
script_dir <- dirname(sys.frame(1)$ofile)
main_folder <- file.path(script_dir, "..")

# load functions
if(!exists("estimates_table")) {
  
  source(file.path(main_folder, "R", "estimates_table.R"))
}

# read plausibility data for Exp. 2
df <- read.csv(file.path(main_folder, "data", "compiled_responses", "exp2_pleasantness.csv"))

# set reference values for treatment coding
df$instrument <- factor(df$instrument,
                        levels = c("violin", "vocal", "clarinet", "tuba"))
df$instrument <- relevel(df$instrument,
                         ref = "violin")

contrasts(df$instrument) <- contr.sum(levels(df$instrument))

model <- lmer(response ~ instrument + (1 | participants),
              data = df)

row_names <- c("Intercept", "Violin", "Alto voice", "Clarinet")

estimates_table(model, row_names)
