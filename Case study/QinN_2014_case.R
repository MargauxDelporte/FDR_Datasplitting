# Install if needed
# install.packages("remotes")
# remotes::install_github("MRCIEU/ieugwasr")

library(ieugwasr)

# Check if the dataset exists
datasets <- gwasinfo("QinN_2014")
datasets

# Download full summary statistics to a file
outfile <- "QinN_2014.txt.gz"
write_gwas("QinN_2014", outfile = outfile)

cat("Downloaded file:", outfile, "\n")