assembly_fasta = "C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\ddCarHirs1.hap1.1.fasta"
repeats <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\ddCarHirs1.hap1.1_repeats_filtered.csv")
arrays <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\ddCarHirs1.hap1.1_arrays_filtered.csv")
genome_classes_data <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\genome_classes_ddCarHirs1.hap1.1.csv")
metadata_data <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\metadata_ddCarHirs1.hap1.1.fasta.csv")
edta <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\ddCarHirs1.hap1.1_edta_filtered.csv")
genes <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\ddCarHirs1.hap1.1_helixer_filtered.csv")
genes <- genes[genes$V3 == "gene",]


write.csv(chr_classes, file = "C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\ddCarHirs1.hap1.1.fasta\\ddCarHirs1.hap1.1_centromeric_scores.csv", row.names = FALSE)  # Assuming scores_data is created in your code

# sync from windows to linux: f
rsync -rltDvh --delete --exclude='.git/' --exclude='.Rproj.user/' "/mnt/c/Users/Piotr Włodzimierz/Documents/GitHub/CAP/" ~/CAP/
  
# back from linux to windows:
rsync -rltDvh --delete \
--exclude='.git/' \
--exclude='.Rproj.user/' \
~/CAP/ "/mnt/c/Users/Piotr Włodzimierz/Documents/GitHub/CAP/" 


install.packages(c("stringr", 
                   "stringdist", 
                   "seqinr",
                   "doParallel",
                   "getopt",
                   "ape",
                   "BiocManager"))

install.packages(c("doParallel",
                   "getopt",
                   "ape",
                   "BiocManager"))

BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("msa")




install.packages("renv")
install.packages("yaml")
