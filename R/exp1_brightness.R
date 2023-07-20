library(lme4)
library(lmerTest)
library(xtable)

# load brightness data
brightness <- read.csv("/Users/simon/congruency/data/exp1_brightness.csv")

## fix pc1
# search for fixed F0 data and pc1
ind_fix_pc1 <- which(brightness$condition == "fixed" & brightness$pc2 == 0)

# setup linear mixed-effects model
model_fix_pc1 <- lmer(response ~ pc1 + (1 | participant),
                  data=brightness[ind_fix_pc1, ])
row_names = c("Intercept", "PC2")

# print estimates
estimates_table(model_fix_pc1, row_names)

## fix pc2
# search for fixed F0 data and pc2
ind_fix_pc2 <- which(brightness$condition == "fixed" & brightness$pc1 == 0)

# setup linear mixed-effects model
model_fix_pc2 <- lmer(response ~ pc2 + (1 | participant),
                      data=brightness[ind_fix_pc2, ])

row_names = c("Intercept", "PC2")

# print estimates
estimates_table(model_fix_pc2, row_names)





## congr/incongr pc1
# search for congruent/incongruent data and pc1
ind_con_pc1 <- which(brightness$condition != "fixed" & brightness$pc2 == 0)
brightness_pc1 <- brightness[ind_con_pc1, ]
# specify custom order of levels and set reference level
brightness_pc1$condition <- factor(brightness_pc1$condition,
                               levels = c("congruent", "incongruent"))
brightness_pc1$condition <- relevel(brightness_pc1$condition, ref = "incongruent")
contrasts(brightness_pc1$condition) <- contr.sum(levels(brightness_pc1$condition))

# setup linear mixed-effects model
model_con_pc1 <- lmer(response ~ condition*pc1 + log2(f0) + (1 | participant),
                      data=brightness_pc1)

row_names = c("Intercept", "congr.", "PC1", "F0", "congr:PC1")
# print estimates
estimates_table(model_con_pc1, row_names)

## congr/incongr pc2
# search for congruent/incongruent data and pc1
ind_con_pc2 <- which(brightness$condition != "fixed" & brightness$pc1 == 0)
brightness_pc2 <- brightness[ind_con_pc2, ]
# specify custom order of levels and set reference level
brightness_pc2$condition <- factor(brightness_pc2$condition,
                                   levels = c("congruent", "incongruent"))
brightness_pc2$condition <- relevel(brightness_pc2$condition, ref = "incongruent")
contrasts(brightness_pc2$condition) <- contr.sum(levels(brightness_pc2$condition))

# setup linear mixed-effects model
model_con_pc2 <- lmer(response ~ condition*pc2 + log2(f0) + (1 | participant),
                      data=brightness_pc2)

row_names = c("Intercept", "congr.", "PC2", "F0", "congr:PC2")

# print estimates
estimates_table(model_con_pc2, row_names)