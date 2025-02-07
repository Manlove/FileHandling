# Load necessary libraries
library(dplyr)

# Define the directory where your files are located
directory <- ""

# Define the pattern to match filenames
pattern <- ".txt"

# List all files that match the pattern
files <- list.files(directory, pattern = pattern, full.names = TRUE)
f <- sub("_annotated","",files)
f <- sub(".txt","",f)
f <- unique(f)

for (file in f) {
  print(file)
  data <- read.table(paste0(file,".txt"),sep="\t",header=TRUE)
  data$ID <- paste(data$Chromosome, data$Start, data$End, sep="_")
  
  data_annot <- read.table(paste0(file,"_annotated.txt"),sep="\t",header=TRUE)
  data_annot$ID <- paste(data_annot$Chromosome, data_annot$Start, data_annot$End, sep="_")
  data_annot <- data_annot[,(ncol(data_annot)-7):ncol(data_annot)]
  
  merged = merge(data,data_annot,by="ID")
  merged$ID <- NULL
  
  write.table(merged, file = paste0(file,"_updated.txt"), row.names = FALSE, sep = "\t",quote = FALSE)
}

