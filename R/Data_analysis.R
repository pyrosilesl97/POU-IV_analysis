
# Title:
  #"Cnidarian hair cell development illuminates an ancient role for the class IV POU transcription factor in defining mechanoreceptor identity"
# Code author:
  # Pablo Yamild Rosiles Loeza
  # For information please go to https://pyrosilesl97.github.io/
# Date: "9/28/2021"

# NOTES: Within the code, there are square brackets in square brackets those paths that must be modified by the user for its use.
################################################################################
#################### DOWNLOAD DATA #############################################
################################################################################

# To initiate the analysis, we required to process the RNA-seq data from Tournière, et al. (2020).

## To achieve this, reference genomic data must be downloaded.
## To download the complete fasta file of the sequence genome
## Run in terminal
#### wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-51/fasta/nematostella_vectensis/dna/Nematostella_vectensis.ASM20922v1.dna.toplevel.fa.gz

# To all annotation files, please go to https://doi.org/https://doi.org/10.6084/m9.figshare.807696.v2
# Or run in terminal
### wget https://figshare.com/ndownloader/articles/807696/versions/2
### unzip 2

# In order to run this script, you can clone the repository from Github. The annotation data downloaded
# in the previous step should be located in the folder ref_genomes

########################################
# Secondly, RNA-seq data were downloaded
#To download all files interactively, go to https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8658/
# Or run in terminal
#### cd data
# POU mutants data
# Run in terminal
####  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/003/ERR3809533/ERR3809533.fastq.gz
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/002/ERR3809532/ERR3809532.fastq.gz
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/004/ERR3809534/ERR3809534.fastq.gz
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/005/ERR3809535/ERR3809535.fastq.gz

# POU wild type data
# Run in terminal
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/006/ERR3809536/ERR3809536.fastq.gz
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/007/ERR3809537/ERR3809537.fastq.gz
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/008/ERR3809538/ERR3809538.fastq.gz
#### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/009/ERR3809539/ERR3809539.fastq.gz

# Run in terminal to uncompress files
#### gzip *.fastq.gz
#### mv *.fastq /before_trim
#These data downloaded should be located in the folder data, in before_trim folder

########################################
# Finally, we need t download ChIP-seq data
# Go to [link]
# Or run in terminal
##### wget [link]

################################################################################
#################### DATA QUALITY ##############################################
################################################################################

# For RNA-seq data, we use fastqc,
# In terminal run
#### fastqc -O quality_files/before_trim data/*fastq

# The html files, will be located in quality_files/before_trim

# Since k-mer content, sequence duplication levels and overrepresented sequences may fail due to
# adaptor content or poor sequence quality, we used the Trimmomatic tool to process the datasets.
#For this, we download illumina adapters from https://gist.github.com/photocyte/3edd9401d0b13476e60f8b104c2575f8

#### cd data/
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809539.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep4.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809538.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep3.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809537.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep2.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809536.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep1.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809533.fastq Nematostella_vectensis_CH2_12d_POU--_rep1 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809532.fastq Nematostella_vectensis_CH2_12d_POU--_rep2 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809534.fastq Nematostella_vectensis_CH2_12d_POU--_rep3 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
#### java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809535.fastq Nematostella_vectensis_CH2_12d_POU--_rep4 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

# We included the code for all replicates in order to ensure that the file names are correct.


# We can check data quelity again, using this code in terminal

#### fastqc -O quality_files/after_trim data/*fastq

# The html files, will be located in quality_files/after_trim


# For ChIP-seq data, we use ChipQC in R
# You can install it, using BiocManager::install("ChIPQC")
library(ChIPQC)
sample1 = ChIPQCsample('../data/[]')
sample2 = ChIPQCsample('../data/[]')
sample3 = ChIPQCsample('../data/[]')

# Report from sample 1
ChIPQCreport(sample1)
# Report from sample 2
ChIPQCreport(sample2)
# Report from sample 3
ChIPQCreport(sample3)

################################################################################
#################### GUIDED SEQUENCE ALIGNMENT #################################
################################################################################



################################################################################
#################### SYSTEM INFORMATION #################################
################################################################################

# The system information can be found in R/system.txt and was generated with the following commands.

options(width = 120)
pkgs <- loadedNamespaces()
pkgs <- installed.packages()[, "Package"]
sessioninfo::session_info(pkgs)

