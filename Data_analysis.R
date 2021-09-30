#MIT License

#Copyright (c) 2021 Pablo Yamild Rosiles Loeza

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

 # The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.



## wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-51/fasta/nematostella_vectensis/dna/Nematostella_vectensis.ASM20922v1.dna.toplevel.fa.gz


## wget https://figshare.com/ndownloader/articles/807696/versions/2

## unzip 2


## cd data

## # POU mutants data

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/003/ERR3809533/ERR3809533.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/002/ERR3809532/ERR3809532.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/004/ERR3809534/ERR3809534.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/005/ERR3809535/ERR3809535.fastq.gz

##
## # POU wild type data

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/006/ERR3809536/ERR3809536.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/007/ERR3809537/ERR3809537.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/008/ERR3809538/ERR3809538.fastq.gz

## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR380/009/ERR3809539/ERR3809539.fastq.gz

##

## gzip *.fastq.gz

## #Then we move the files

## mv *.fastq /before_trim


##  wget [link]


## fastqc -O quality_files/before_trim data/*fastq

##

##
## cd data/

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809539.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep4.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809538.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep3.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809537.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep2.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809536.fastq Nematostella_vectensis_CH2_12d_POU-+_wild_rep1.fastq ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

##
## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809533.fastq Nematostella_vectensis_CH2_12d_POU--_rep1 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809532.fastq Nematostella_vectensis_CH2_12d_POU--_rep2 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809534.fastq Nematostella_vectensis_CH2_12d_POU--_rep3 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10

## java -jar [#PATH_TO_TRIMMOMATIC_JAVA_FILE]/trimmomatic-0.32.jar SE -threads 4 before_trim/ERR3809535.fastq Nematostella_vectensis_CH2_12d_POU--_rep4 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10


## fastqc -O quality_files/after_trim data/*fastq

##

## #$

## Rscript [/path_to_spp_script]run_spp.R -c=[ChIPseq_rep1] -savp -out=[ChIPseq_rep1.outphantom]

## #$


## ----eval=FALSE---------------------------------------------------------------------------------
## BiocManager::install("ChIPQC")


## ---- eval=FALSE--------------------------------------------------------------------------------
## library(ChIPQC)
## sample1 = ChIPQCsample('../data/[ChIPseq_rep1]')
## sample2 = ChIPQCsample('../data/[ChIPseq_rep2]')
## sample3 = ChIPQCsample('../data/[ChIPseq_rep3]')
##
## # Report from sample 1
## ChIPQCreport(sample1)
## # Report from sample 2
## ChIPQCreport(sample2)
## # Report from sample 3
## ChIPQCreport(sample3)


## sed 's/NEMVE//g' [fasta_file_with_sequences] > Nvectensis.fa


## STAR --runThreadN 6 --runMode genomeGenerate --genomeDir STAR_genome_index --genomeFastaFiles ref_genomes/[Nvectensis.fa] --sjdbGTFfile ref_genomes/[nveGenes.vienna130208.nemVec1.gtf] --sjdbOverhang 99


## #@ run for each replicate

## STAR --runThreadN 6 --genomeDir  ref_genomes/ --readFilesIn data/Nematostella_vectensis_CH2_12d_POU-+_wild_rep1.fastq --outFileNamePrefix alignment/Nematostella_vectensis_CH2_12d_POU-+_wild_rep1

## #@ run for each replicate


## cd alignment

## #@ run for each replicate

## # To compress file to bam

## samtools view -S -b alignment/Nematostella_vectensis_CH2_12d_POU-+_wild_rep1.out.sam > Nvect_POU-+_rep1_Aligned.bam

## # Sort

## samtools sort Nvect_POU-+_rep1_Aligned.bam -o Nvect_POU-+_rep1_sorted.bam

## # Counts

## htseq-count -f bam -s yes -r pos Nvect_POU-+_rep1_sorted.bam ../ref_genomes/[nveGenes.vienna130208.nemVec1.gtf] > POU-+_rep1.counts

## #@ run for each replicate


## cd alignment

## paste POU--_rep*.counts POU-+_rep*.counts > counts_table_DGE.txt


## ---- message=FALSE, warning=FALSE--------------------------------------------------------------
library(DESeq2)
library(apeglm)
#Read data
data.dge <- read.delim(file = '../alignment/counts_table_DGE.txt',sep = '\t', header = F, row.names = 1)
#Assign names
condition <- factor(c("POU_ko","POU_ko","POU_ko","POU_ko","POU_wt","POU_wt","POU_wt","POU_wt"), levels = c('POU_wt','POU_ko'))
colData <- data.frame(row.names=colnames(data.dge), condition)

dds <- DESeqDataSetFromMatrix(countData = data.dge,
                              colData = colData,
                              design = ~ condition )

#Erase data with no counts
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)

#You can preprocces all again, and you should get this figure
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
plotMA(resLFC, ylim=c(-4,4))


## ---- message=FALSE-----------------------------------------------------------------------------
library(reshape)
library(SummarizedExperiment)
vst <- varianceStabilizingTransformation(dds)
vst@colData@listData[["condition"]] <- c("POU_ko","POU_ko","POU_ko","POU_ko","POU_wt","POU_wt","POU_wt","POU_wt")
plotPCA(vst, intgroup = "condition")


## -----------------------------------------------------------------------------------------------
#Para 0.05
nrow(data.frame(resLFC@rownames[resLFC@listData[["pvalue"]] <0.05]))
#Para 0.01
nrow(data.frame(resLFC@rownames[resLFC@listData[["pvalue"]] <0.01]))


## -----------------------------------------------------------------------------------------------
#Keep only p-value<0.01
differential_genes <- data.frame(resLFC@rownames[resLFC@listData[["pvalue"]] <0.01])
#Add column name
colnames(differential_genes) <- c('gene')
#Add Fold change
differential_genes$fold_change <- resLFC@listData[["log2FoldChange"]][resLFC@listData[["pvalue"]] <0.01]
#Add base mean
differential_genes$mean_counts <- resLFC@listData[["baseMean"]][resLFC@listData[["pvalue"]] <0.01]
#Add p value
differential_genes$pvalue <- resLFC@listData[["pvalue"]][resLFC@listData[["pvalue"]] <0.01]
# Eliminate empty rows
differential_genes <- differential_genes[complete.cases(differential_genes),]


## ----eval=FALSE---------------------------------------------------------------------------------
## #Read file
## GO_terms <- read.delim(file = '../data/nveGenes.vienna130208.GO_annotation_141017.txt',sep = '\t', header = F)
## # Extract those that are differentially expressed
## GO_terms_differential_expressed <- GO_terms[GO_terms$V1%in%differential_genes$gene, ]
## #Quiero tener los valores de p, guardados en differential_genes$pvalue en el df de GO
## for(j in 1:nrow(GO_terms_differential_expressed)) {
##   for(k in 1:nrow(differential_genes)) {
##      if (GO_terms_differential_expressed$V1[j] == differential_genes$gene[k]) {
##          GO_terms_differential_expressed$pvalue[j] <- differential_genes$pvalue[k]
##      }
##     k <- k+1
##   }
##   j<- j+1
##   }


## -----------------------------------------------------------------------------------------------
#reading table
equival <- read.delim(file = '../data/equivalence_table.tsv',sep = '\t', header = F)
###Loop for extracting data
differential_genes$JGI_ID <- 'No data'
j <-1
k<-1

for(j in 1:nrow(differential_genes)) {
  for(k in 1:nrow(equival)) {
     if (differential_genes$gene[j] == equival$V1[k]) {
         differential_genes$JGI_ID[j] <- equival$V3[k]
     }
    k <- k+1
  }
  j<- j+1
}


upregulated_genes <- differential_genes[differential_genes$fold_change>0,]
downregulated_genes <- differential_genes[differential_genes$fold_change<0,]


## ----eval=FALSE---------------------------------------------------------------------------------
## library(rrvgo)
## library(org.Ce.eg.db) #C. elegans data as reference
## library(clusterProfiler)
## #We need some dfs to do the job
##
## #Getting terms to gene dataframe
## GOterms_to_genes <- GOterms_to_genes[,c(3,1)]
## colnames(GOterms_to_genes) <- c('term','gene')
## #Getting term to name df
## GOterms_to_names <- GOterms_to_names[,c(3,4)]
## colnames(GOterms_to_names) <- c('term','name')
##
## #GO enrichment for downregulated genes
## ggo <- enricher(
##   downregulated_genes$gene,
##   pvalueCutoff = 0.05,
##   pAdjustMethod = "BH",
##   TERM2GENE=GOterms_to_genes,
##   TERM2NAME = GOterms_to_names
## )
##
## #GO enrichment result
## ggo_result <- ggo@result
## #Filtering GO enrichment result
## ggo_result <- ggo_result[ggo_result$p.adjust<0.05,]
##
## #Generating matrix to plot
## simMatrix <- calculateSimMatrix(ggo_result$ID,
##                                 orgdb="org.Ce.eg.db",
##                                 ont="BP",
##                                 method="Rel")
##
## GO_terms_names <- ggo_result[ggo_result$ID%in%row.names(simMatrix),]
## scores <- setNames(-log10(GO_terms_names$p.adjust), GO_terms_names$ID)
## #Reducing dimensions to plot
## reducedTerms <- reduceSimMatrix(simMatrix, scores,
##                                 orgdb="org.Ce.eg.db", threshold = 0.8)
## treemapPlot(reducedTerms)


## macs2 callpeak -t data/[ChIPseq_rep1] -c data/[ChIPseq_rep1_input] -f BAMPE -g 361523379 -n cp -B


## # Code used to obtain consensus peaks.

## [path_to_mspc]/mspc -i data/*.narrowPeak -r tec -w 1e-4 -s 1e-8 -o ./MSPC_outdir


## bedtools getfasta -fi ../ref_genomes/Nvectensis.fa -bed MSPC_outdir/ConsensusPeaks.bed -fo consensus_peaks_sequence


## ---- eval=FALSE--------------------------------------------------------------------------------
## library('GenomicFeatures')
## library("GenomicRanges")
## library("AnnotationHub")
## library("rtracklayer")
##
## #Reading annotation files
## anota_data <-  makeTxDbFromGFF(file = '../data/[annotation_file.gtf]', format = "gff3")
##
## #getting promoters
## promoters_genome <- data.frame(promoters(genes(anota_data), upstream = 350, downstream = 100))
## promoters_genome$start[promoters_genome$start<0] <- 0
##
##
##
## #To get the genes
## genes_anotation <- data.frame(genes(anota_data))
## #Write file to bedtools
## write.table(x = genes_anotation,file = '../data/genes.bed',quote = F,sep = '\t',row.names = F,col.names = F)
## write.table(x = promoters_total,file = '../data/promoter.bed',quote = F,sep = '\t',row.names = F,col.names = F)
##


## bedtools intersect -wb -a  MSPC_outdir_1/ConsensusPeaks.bed -b promoter.bed  > peaks_on_promoters.tsv

## bedtools intersect -wb -a  MSPC_outdir_1/ConsensusPeaks.bed -b genes.bed  > peaks_on_genes.tsv


## -----------------------------------------------------------------------------------------------
#Change the name of the file to peaks_on_promoters.tsv and peaks_on_genes.tsv
peaks_on_promoters <- read.table(file = '../data/peaks_on_prom_example.tsv', sep = '\t')
peaks_on_genes <- read.table(file = '../data/peaks_on_gene_example.tsv', sep = '\t')


## -----------------------------------------------------------------------------------------------
#Promoter regions of genes with peaks
#There should not be repeated data
peaks_on_promoters <- peaks_on_promoters[,c(6,7,8,9,10,11)]
peaks_on_promoters <- peaks_on_promoters[!duplicated(peaks_on_promoters),]

print('We have 1271 genes with peaks on the promoter region')

#Genes with peaks
#Should not be duplicated data
genes_with_peaks <- peaks_on_genes[,c(6,7,8,9,10,11)]
genes_with_peaks <- genes_with_peaks[!duplicated(genes_with_peaks),]
print('We have 3867 genes with peaks on the promoter region')


## ----warning=FALSE------------------------------------------------------------------------------
library(dplyr)
###########Peaks on differentially expressed genes ########
differential_genes_with_peaks_on_genes <-differential_genes[differential_genes$gene%in%genes_with_peaks$V11,]

###########Peaks on promoter ########
differential_genes_with_peaks_on_promoter <-differential_genes[differential_genes$gene%in%peaks_on_promoters$V11,]


#List of all genes
differential_genes_with_peaks_all <- unique.data.frame(data.frame(rbind(differential_genes_with_peaks_on_promoter,differential_genes_with_peaks_on_genes)))

differential_genes_with_peaks_upregulated <- differential_genes_with_peaks_all[differential_genes_with_peaks_all$fold_change>0,]

differential_genes_with_peaks_downregulated <- differential_genes_with_peaks_all[differential_genes_with_peaks_all$fold_change<0,]

#Differences in the data
genes_with_peaks_in_gene_body <- differential_genes_with_peaks_on_genes[!differential_genes_with_peaks_on_genes$gene%in%differential_genes_with_peaks_on_promoter$gene,]

genes_with_peaks_in_promoter_region <- differential_genes_with_peaks_on_promoter[!differential_genes_with_peaks_on_promoter$gene%in%differential_genes_with_peaks_on_genes$gene,]

genes_with_peaks_on_gene_body_and_promoter_region <- differential_genes_with_peaks_on_genes[differential_genes_with_peaks_on_genes$gene%in%differential_genes_with_peaks_on_promoter$gene,]


## ---- eval=FALSE--------------------------------------------------------------------------------
## library(ggvenn)
##
## x <- list(Gene_body = differential_genes_with_peaks_on_genes$gene,Promoters = differential_genes_with_peaks_on_promoter$gene)
## ggvenn(x,
##          fill_color = c("#0073C2FF", "#CD534CFF"),
##   stroke_size = 0.5, set_name_size = 4
##   )+
## ggtitle('Venn diagram of ChIP-seq peaks in regions of \n differentially expressed genes.')


## ---- eval=F------------------------------------------------------------------------------------
## #Graphic
## library(rrvgo)
## library(org.Ce.eg.db)
## ddgt <- differential_genes_with_peaks_on_promoter[differential_genes_with_peaks_on_promoter$fold_change<0,]
## gene_enricher <- ddgt$gene
## ggo <- enricher(
##   gene_enricher,
##   pvalueCutoff = 0.05,
##   pAdjustMethod = "BH",
##   minGSSize = 10,
##   maxGSSize = 500,
##   qvalueCutoff = 0.2,
##   TERM2GENE=GOterms_to_genes,
##   TERM2NAME = GOterms_to_names
## )
##
## ggo_result <- ggo@result
## ggo_result <- ggo_result[ggo_result$p.adjust<0.05,]
## simMatrix <- calculateSimMatrix(ggo_result$ID,
##                                 orgdb="org.Ce.eg.db",
##                                 ont="BP",
##                                 method="Rel")
##
## GO_terms_names <- ggo_result[ggo_result$ID%in%row.names(simMatrix),]
## scores <- setNames(-log10(GO_terms_names$p.adjust), GO_terms_names$ID)
## reducedTerms <- reduceSimMatrix(simMatrix, scores,
##                                 orgdb="org.Ce.eg.db", threshold = 0.7)
## treemapPlot(reducedTerms)


## -----------------------------------------------------------------------------------------------
options(width = 120)
pkgs <- loadedNamespaces()
pkgs <- installed.packages()[, "Package"]
sessioninfo::session_info(pkgs)

