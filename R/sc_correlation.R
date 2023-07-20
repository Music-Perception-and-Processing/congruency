# Define data frame
df <- read.csv("/Users/simon/congruency/data/sc_correlation.csv")

# Calculate Spearman rank correlation
cor.test(df$PC1, df$SC, method = "spearman")