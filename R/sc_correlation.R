# define directories
script_dir <- dirname(sys.frame(1)$ofile)
main_folder <- file.path(script_dir, "..")

# load functions
if(!exists("estimates_table")) {
  
  source(file.path(main_folder, "R", "estimates_table.R"))
}

# Define data frame
df <- read.csv(file.path(main_folder, "data", "sc_correlation.csv"))

# Calculate Spearman rank correlation
cor.test(df$PC1, df$SC, method = "spearman")