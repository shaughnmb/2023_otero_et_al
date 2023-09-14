#!/bin/env Rscript

#### License, Title, Authors, and Info  ####
  #### License Notice ####

##
# Copyright (c) 2023 Cedars-Sinai Medical Center
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

  #### Title and Authors ####

##
# Title: Organ-Chips Enhance the Maturation of Human iPSC-Derived Dopamine Neurons
#
# Authors:  MAria G. Otero 1, Shaughn Bell 1, Alexander H. Laperle 1, George Lawless 1, 
#           Zachary NMyers 1, Marian A. Castro 1, Jaquelyn M. Villalba, Clive N. Svendsen*1
#
# Affiliations: 1 Cedars-Sinai Board of Governors Regenerative Medicine 
#                 Institute, Cedars-Sinai Medical Center, Los Angeles CA
# *Corresponding author email:  Clive.Svendsen@cshs.org 
##

  #### Script Information ####

##
# R version 4.2.1
# R Script Title:  2023_otero_et_al_v3.R
# R Script Author:  Shaughn Bell
# R Script Corresponding Email:  Shaughn.Bell@cshs.org
#
# Notes: 
#   A) Script makes use of the variables set up under "project information" as 
#      well as additional "prefixes" throughout the script for ease of saving
#      files with a similar path and naming structure.  When reloading data 
#      (see note "B"), you must either use the full path or reload the prefixes.
#   B) Script saves intermediate steps at each major manipulation of the seurat
#      object via the "saveRDS" function.  If needed, these RDS objects can then 
#      be reloaded to easily restart at one of these save points without needing 
#      to start from scratch.  However, these are not required for analysis, and
#      they can be skipped to save time and disk space.
##

#### Samples and Setup ####
  #### Sample Information #### 

##
# GEO Accession Number:  GSE239951
#
# aligned to homo sapiens hg38 via CellRanger v6.1.2
# used the 10x Genomics supplied reference file "refdata-gex-GRCh38-2020-A.tar.gz"
#
# fastq files used:
#
# EDi044_A_2DA_S5_L001_I1_001.fastq.gz
# EDi044_A_2DA_S5_L001_R1_001.fastq.gz
# EDi044_A_2DA_S5_L001_R2_001.fastq.gz
# EDi044_A_2DA_S5_L002_I1_001.fastq.gz
# EDi044_A_2DA_S5_L002_R1_001.fastq.gz
# EDi044_A_2DA_S5_L002_R2_001.fastq.gz
# EDi044_A_chip_S1_L001_I1_001.fastq.gz
# EDi044_A_chip_S1_L001_R1_001.fastq.gz
# EDi044_A_chip_S1_L001_R2_001.fastq.gz
# EDi044_A_chip_S1_L002_I1_001.fastq.gz
# EDi044_A_chip_S1_L002_R1_001.fastq.gz
# EDi044_A_chip_S1_L002_R2_001.fastq.gz
#
# Cellranger outs used:  
#
# hg38_EDi044_A_2DA_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# hg38_EDi044_A_2DA_outs/filtered_feature_bc_matrix/features.tsv.gz
# hg38_EDi044_A_2DA_outs/filtered_feature_bc_matrix/matrix.mtx.gz 
#
# hg38_EDi044_A_chip_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# hg38_EDi044_A_chip_outs/filtered_feature_bc_matrix/features.tsv.gz
# hg38_EDi044_A_chip_outs/filtered_feature_bc_matrix/matrix.mtx.gz
#
# Cellranger outs are not included in the GEO record.
#
# Genesummed matricies of counts:  
# 
# 1_EDi044_A_DA_filtered_genesummed_mtx.csv
# 2_EDi044_A_chip_filtered_genesummed_mtx.csv
#
# Genesummed matrix files are included in the GEO record. 
#
##

  #### Project Information ####

date <- "20221109"
project <- "chip_vs_2Dastro"
datadir <- "/chip_vs_2Dastro" # directory to save files generated
sourcedir <- "/20221026_hg38_outs_copy" # save Cellranger outs here

  #### Load required packages ####

library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(irlba)
library(gridExtra)
library(scCustomize)
library(nichenetr)
library(EnhancedVolcano)
library(Vennerable)
library(paletteer)

  #### Set working dir and save session info ####
setwd(datadir)

sessionInfo()

sink(paste0(date,"_",project,"_devtools_sessionInfo.txt"))
devtools::session_info()
sink()

sink(paste0(date,"_",project,"_sessionInfo.txt"))
sessionInfo()
sink()

  #### Functions for saving images ####
hirestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white"
  )
}

lowrestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 100,
    limitsize = TRUE,
    bg = "white"
  )
}

hirestiffsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white",
    width = 9,
    height = 9,
    units = c("in")
  )
}

lowrestiffsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 100,
    limitsize = TRUE,
    bg = "white",
    width = 9,
    height = 9,
    units = c("in")
  )
}
  #### set color scheme ####

volc_pallette <- c("#DAA520","#088DA5","#ABABAB")
scales::show_col(volc_pallette)

color_palette <- c("#00ff00","#0099ff")
scales::show_col(color_palette)

color_palette2 <- paletteer_d("ggthemes::Classic_Color_Blind")
scales::show_col(color_palette2)

color_palette3 <- c("#DAA520","#088da5")
scales::show_col(color_palette3)

#### Gene Sum and Seurat Object ####
  #### Gene Sum Setup ####

# Some gene symbols refer to more than one ENSG ID.  
# In order to ensure all ENSG IDs are represented while using the much more
# user-friendly gene symbol, the cellranger outpus were processed as follows:
#    1) Cellranger output was loaded into a matrix object
#    2) All genes with more than one ENSG ID were put into a new sub-matrix
#    3) Counts from all instances of a given gene symbol were summed
#    4) The genesummed sub-matrix was rejoined to the main matrix
#    5) The ENSG ID column was removed
# Genesummed matrix files are included in the GEO record.  
# Cellranger outs are not included in the GEO record.

# Create matrix files from the cellranger outs

filelist <- list.files(path = sourcedir)
filelist <- as.data.frame(filelist, files = list.files(path = sourcedir))
path_list <- paste0(sourcedir,"/",filelist[,1],"/filtered_feature_bc_matrix")

### IMPORTANT ###
# Files with numbers in their name may be sorted in a 
# different order depending on the way files are listed.  
# Some systems sort by number instead of numeric value.
# Make sure to use the same sort order as in the path_list object.

path_list

# check that all the files exist
all(file.exists(path_list)) 

  #### Combine Cellranger output files ####
#Read the Cellranger output files and save expression matrix csv
for (k in 1:length(path_list)) {
  
  # read in and transpose the barcode file
  bc <- as.data.frame(
    read.delim(
      file = paste0(
        path_list[k], "/barcodes.tsv.gz"
      ),
      header = F
    )
  ) %>% t(.)
  
  # read in the features file and select just the first two columns
  feat <- as.data.frame(
    read.delim(
      file = paste0(
        path_list[k], "/features.tsv.gz"
      ),
      header = F
    )
  ) %>% .[, 1:2]
  
  # read in the counts
  matx <- as.data.frame(
    readMM(
      file = paste0(
        path_list[k], "/matrix.mtx.gz"
      )
    )
  )
  
  # assign barcodes as colnames for the matrix
  colnames(matx) <- bc[1, ]
  
  # rename the columns of features table
  colnames(feat) <- c("ensmbl.id", "symbol")
  
  # cbind the features table and the matrix
  matx <- cbind(feat, matx)
  
  write.csv(matx,
            file = paste0(path_list[k], "/filtered_combined_mtx.csv"),
            row.names = FALSE
  )
}

rm(matx, bc, feat, filelist, k)

# check the order of the files in path_list
path_list
  #### Load expression matrix ####
# Load the data
# Make sure the index number after path_list matches the correct file

data.2DA <- read.csv(
  file = paste0(
    path_list[1],
    "/filtered_combined_mtx.csv"
  ),
  header = TRUE
)
data.CH <- read.csv(
  file = paste0(
    path_list[2],
    "/filtered_combined_mtx.csv"
  ),
  header = TRUE
)
  #### Gene Sum ####
# Put all the loaded matrices into one list
obj.list <- list(data.2DA, data.CH)

# Iterate the summing over all samples in obj.list
for (j in 1:length(obj.list)) {
  
  # Create a list of duplicated genes
  gene.duplicates <- obj.list[[j]]$symbol[
    duplicated(obj.list[[j]]$symbol, incomparables = NA)
  ] %>%
    as.list(.) %>%
    unique(.)
  
  # Create empty dataframe to rbind duplicated rows to
  dup.rows <- data.frame()
  
  # Create identical expression matrix to subtract duplicated rows from
  sub.rows <- obj.list[[j]]
  
  # Extract rows containing the duplicated genes.
  for (i in 1:length(gene.duplicates)) {
    # subsets duplicated rows
    dups <- subset(obj.list[[j]], subset = symbol == gene.duplicates[[i]])
    # binds duplicated rows together
    dup.rows <- rbind(dup.rows, dups)
    # deletes duplicated rows from expression matrix
    sub.rows <- subset(sub.rows, subset = symbol != gene.duplicates[[i]])
  }
  
  # Remove ensembl ID columns
  dup.rows <- dup.rows[, 2:length(colnames(dup.rows))]
  sub.rows <- sub.rows[, 2:length(colnames(sub.rows))]
  
  # Find the gsum of duplicated rows
  dup.rows <- aggregate(. ~ symbol, data = dup.rows, sum)
  
  # rbind the gsum of duplicated rows with the 
  # expression matrix that had these genes removed
  new.mat <- rbind(sub.rows, dup.rows)
  
  if (nrow(new.mat) == length(unique(obj.list[[j]]$symbol))) {
    new.mat <- as.data.frame(new.mat)
    rownames(new.mat) <- new.mat[, 1]
    new.mat <- subset(new.mat, select = -c(symbol))
    obj.list[[j]] <- new.mat
    message("Successfully merged!")
  } else {
    message("There's an error!")
  }
}

# Put genesummed objects back into the data matrix objects
data.2DA <- obj.list[[1]]
data.CH <- obj.list[[2]]

  #### Save genesummed matrices ####

write.csv(data.2DA,
          file = paste0(path_list[1], "/11852DA_filtered_genesummed_mtx.csv"),
          row.names = FALSE
)
write.csv(data.CH,
          file = paste0(path_list[2], "/1185chip_filtered_genesummed_mtx.csv"),
          row.names = FALSE
)

# save R objects for easy reloading
save(data.2DA, data.CH,
     file = paste0(datadir, "/", date, "_", project, "_genesummed.rdata")
)

# load(file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

  #### Create Seurat objects ####

# Initialize the Seurat object with the genesummed data.
# Keep all genes expressed in >= 1 cell

s2DA <- CreateSeuratObject(
  counts = obj.list[[1]],
  project = "EDi044_A_2DA",
  min.cells = 1,
  min.features = 0
)

sCH <- CreateSeuratObject(
  counts = obj.list[[2]],
  project = "EDi044_A_CH",
  min.cells = 1,
  min.features = 0
)

rm(
  data.2DA, data.CH,
  dup.rows, dups, gene.duplicates, new.mat,
  obj.list, sub.rows, i, j, path_list
)

save(s2DA, sCH,
     file = paste0(datadir, "/", date, "_", project, "_seuratobjs.rdata")
)

# load(file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

  #### Add metadata ####

s2DA[["Format"]] <- "2D"
sCH[["Format"]] <- "Chip"

s2DA[["orig.ident"]] <- "EDi044_A_2DA"
sCH[["orig.ident"]] <- "EDi044_A_CH"

  #### Merge all Seurat objects into a single Seurat object ####

tenx1 <- merge(s2DA, y = c(sCH),
               add.cell.ids = c("s2DA","sCH"), 
               project = project)

tenx1

# save merged seurat object 
saveRDS(tenx1, 
        file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

# tenx1 <- read_rds(file = 
# paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

#### QC and Filtering ####
  #### QC setup ####

# Temporarily copy "tenx1" in "data" just in case need to revert to prefiltered
data <- tenx1

# Add mito info to metadata
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")

  #### Prefiltered plots ####
##Plots of UMI#, Gene#, and %mito
vp1 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nCount_RNA prefiltered")
vp1
pdf(paste0("./",date,"_",project,"_prefilter_","nCount_RNA",".pdf"))
print(vp1)
dev.off()

vp2 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nFeature_RNA prefiltered")
vp2
pdf(paste0("./",date,"_",project,"_prefilter_","nFeature_RNA",".pdf"))
print(vp2)
dev.off()

vp3 <- VlnPlot(data, features = "percent.mt", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt prefiltered")
vp3
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt",".pdf"))
print(vp3)
dev.off()

vp4 <- VlnPlot(data, features = "percent.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.ribo prefiltered")
vp4
pdf(paste0("./",date,"_",project,"_prefilter_","percent.ribo",".pdf"))
print(vp4)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_prefilter_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()
  #### Add Z-Score ####
data[["nUMI.z"]] <- scale(data$nCount_RNA)
data[["nGene.z"]] <- scale(data$nFeature_RNA)
data[["percent.mt.z"]] <- scale(data$percent.mt)
data[["percent.ribo.z"]] <- scale(data$percent.ribo)

  #### Prefiltered QC data ####
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

length(WhichCells(data, expression = Format == "2D"))
length(WhichCells(data, expression = Format == "Chip"))

# save prefilter QC data
sink(paste0("./",date,"_",project,"_prefilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

  #### Filter cells based on Z-score ####
data <- subset(data, subset = percent.mt.z < 3)
data <- subset(data, subset = percent.ribo.z < 3)
data <- subset(data, subset = nGene.z < 3)
data <- subset(data, subset = nUMI.z < 3)

  #### Postfilter QC data ####
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

length(WhichCells(data, expression = Format == "2D"))
length(WhichCells(data, expression = Format == "Chip"))

# save postfilter QC data
sink(paste0("./",date,"_",project,"_postfilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

  #### Postfiltered plots ####
##Plots of UMI#, Gene#, and %mito
vp5 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")  + 
  NoLegend() + 
  ggtitle("nCount_RNA filtered")
vp5
pdf(paste0("./",date,"_",project,"_filtered_","nCount_RNA",".pdf"))
print(vp5)
dev.off()

vp6 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
  NoLegend() + 
  ggtitle("nFeature_RNA filtered")
vp6
pdf(paste0("./",date,"_",project,"_filtered_","nFeature_RNA",".pdf"))
print(vp6)
dev.off()

vp7 <- VlnPlot(data, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  NoLegend() + 
  ggtitle("percent.mt filtered")
vp7
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt",".pdf"))
print(vp7)
dev.off()

vp8 <- VlnPlot(data, features = "percent.ribo", pt.size = 0, group.by = "orig.ident") +
  NoLegend() + 
  ggtitle("percent.ribo filtered")
vp8
pdf(paste0("./",date,"_",project,"_filtered_","percent.ribo",".pdf"))
print(vp8)
dev.off()

p3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p3 + p4

pdf(paste0("./",date,"_",project,"_filtered_","scatterplots",".pdf"))
print(p3 + p4)
dev.off()

  #### Additional QC plots ####
# Additional QC metrics
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo")]],
                        rownames="Cell.Barcode")

p5 <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics")
p5
pdf(paste0("./",date,"_",project,"QC_Metrics",".pdf"))
print(p5)
dev.off()

p6 <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mitochondria")
p6
pdf(paste0("./",date,"_",project,"Dist_mito",".pdf"))
print(p6)
dev.off()

p7 <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Ribosomes")
p7
pdf(paste0("./",date,"_",project,"Dist_ribo",".pdf"))
print(p7)
dev.off()

# store filtered data back in original object
tenx1 <- data

rm(data)

  #### Set factor levels ####
seurat_sct <- tenx1
seurat_sct$orig.ident <- factor(x = seurat_sct$orig.ident, levels = c("EDi044_A_2DA","EDi044_A_CH"))
seurat_sct$Format <- factor(x = seurat_sct$Format, levels = c("2D","Chip"))

rm(tenx1)

  #### Save post-QC object ####

saveRDS(seurat_sct, file = paste0("./",date,"_",project,"_seurat_post_qc"))

# seurat_sct <- read_rds(file = paste0("./",date,"_",project,"_seurat_post_qc"))

rm(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,
   p1,p2,p3,p4,p5,p6,p7,qc.metrics)

#### SCTransform, PCA, UMAP, Cluster, and Normalize ####
  #### SCTransform ####
seurat_sct <- SCTransform(seurat_sct, vars.to.regress = c("percent.mt"), verbose = TRUE)

saveRDS(seurat_sct, file = paste0("./",date,"_",project,"_post_sct_pre_pca.rds"))

# seurat_sct <- readRDS(file = paste0("./",date,"_",project,"_post_sct_pre_pca.rds"))

  #### Run PCA ####

# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)
# write.table(xx.coord, file = paste0("./",date,"_",project,"_PCAcoord_R.txt"), col.names = NA, sep = "\t", row.names = T)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)
# write.table(xx.gload, file = paste0("./",date,"_",project,"_PCAgload_R.txt"), col.names = NA, sep = "\t", row.names = T)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), Format = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), Format = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

# meaningful.PCs <- 13

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,eig,
   eig.percent,expected.contribution,i,scree,sum.eig)

  #### Run UMAP ####

seurat_sct <- RunUMAP(seurat_sct, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident",
        cols = color_palette3)

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "orig.ident", 
        group.by = "orig.ident", 
        cols = color_palette3,
        combine = TRUE) + 
  NoLegend()

# by Format 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Format",
        cols = color_palette3)

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Format", 
        group.by = "Format",
        cols = color_palette3) + 
  NoLegend()

  #### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              
seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.1))

res <- "_res.0.1"
ident_res <- paste0("SCT_snn",res)

Idents(seurat_sct) <- ident_res
DimPlot(seurat_sct, reduction = "umap", label = TRUE, label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# set cluster colors
clust_cols <- scales::hue_pal()(length(levels(seurat_sct)))
names(clust_cols) <- levels(seurat_sct)

# Plot the UMAP
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 1,
        cols = color_palette2) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

# UMAP of cells in each cluster by Format without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE,
        cols = color_palette2, 
        split.by = "Format",
        pt.size = 1) + 
  NoLegend() +
  ggtitle(paste0(ident_res))

# UMAP of cells in each cluster by Format with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        cols = color_palette2,
        split.by = "Format",
        label.size = 5,
        pt.size = 1) +
  NoLegend() +
  ggtitle(paste0(ident_res))

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

# save RDS containing reduction and cluster idents
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### Normalize RNA, find variable features, scale data for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct, verbose = FALSE)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

# save RDS containing RNA normalized data
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

#### Cluster Quality Check####
  #### Check cluster quality ####      
VlnPlot(seurat_sct,
        features = "nCount_RNA") + NoLegend() #cluster 2 low
hirestiff(paste(prefixPCres,"VLN","nCount_RNA","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","nCount_RNA","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = "nFeature_RNA") + NoLegend() #cluster 2 low
hirestiff(paste(prefixPCres,"VLN","nFeatures_RNA","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","nFeatures_RNA","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = "percent.mt") + NoLegend() #cluster 2 high
hirestiff(paste(prefixPCres,"VLN","percent.mt","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","percent.mt","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = "percent.ribo") + NoLegend() # cluster 2 high
hirestiff(paste(prefixPCres,"VLN","percent.ribo","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","percent.ribo","lowres.tiff", sep = "_"))

view(seurat_sct@meta.data)

  #### Remove poor quality cluster(s) ####
seurat_sct
seurat_sct <- subset(seurat_sct, idents = c(0,1,3,4,5))
seurat_sct

length(WhichCells(seurat_sct, expression = Format == "2D"))
length(WhichCells(seurat_sct, expression = Format == "Chip"))

#### SCTranscform, PCA, UMAP, Cluster and Normalize remaining cells ####
  #### SCTransform ####
seurat_sct <- SCTransform(seurat_sct, vars.to.regress = c("percent.mt"), verbose = TRUE)

  #### Run PCA ####

# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)
# write.table(xx.coord, file = paste0("./",date,"_",project,"_PCAcoord_R.txt"), col.names = NA, sep = "\t", row.names = T)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)
# write.table(xx.gload, file = paste0("./",date,"_",project,"_PCAgload_R.txt"), col.names = NA, sep = "\t", row.names = T)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

# meaningful.PCs <- 13

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,eig,
   eig.percent,expected.contribution,i,scree,sum.eig)

  #### Run UMAP ####

seurat_sct <- RunUMAP(seurat_sct, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = 1, 
        group.by = "orig.ident",
        cols = color_palette3)

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = 1, 
        split.by = "orig.ident", 
        group.by = "orig.ident",
        cols = color_palette3) + 
  NoLegend()

# by Format 
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = 1, 
        group.by = "Format",
        cols = color_palette3)

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = 1, 
        split.by = "Format", 
        group.by = "Format",
        cols = color_palette3) + 
  NoLegend()

  #### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              
seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.1))
res <- "_res.0.1"
ident_res <- paste0("SCT_snn",res)

Idents(seurat_sct) <- ident_res

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 1,
        cols = color_palette2) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

  #### Counts after last filtering #### 

length(seurat_sct@meta.data$orig.ident)
mean(seurat_sct@meta.data$percent.mt)
mean(seurat_sct@meta.data$percent.ribo)
mean(seurat_sct@meta.data$nCount_RNA)
median(seurat_sct@meta.data$nCount_RNA)
mean(seurat_sct@meta.data$nFeature_RNA)
median(seurat_sct@meta.data$nFeature_RNA)
max(seurat_sct@meta.data$nCount_RNA)

length(WhichCells(seurat_sct, expression = Format == "2D"))
length(WhichCells(seurat_sct, expression = Format == "Chip"))

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### Normalize RNA, find variable features, scale data for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct, verbose = FALSE)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

# save RDS containing RNA normalized data
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

  #### Rename clusters ####
seurat_sct <- RenameIdents(object = seurat_sct, 
                           '0' = "4",
                           '1' = "1", 
                           '2' = "2", 
                           '3' = "3",
                           '4' = "5",
                           '5' = "6"
                           )

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 6,
        pt.size = 1,
        cols = color_palette2) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

# save new cluster numbers in metadata

seurat_sct[["renumbered.clusters"]] <- seurat_sct@active.ident

# reorder the cluster numbers numerically

seurat_sct$renumbered.clusters <- factor(
  x = seurat_sct$renumbered.clusters, 
  levels = c("1","2","3","4","5","6"))

Idents(seurat_sct) <- "renumbered.clusters"

#### Figures ####
  #### Figure 2A ####

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = F,  # numbers added manually with illustrator
        pt.size = 1,
        cols = color_palette2) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

hirestiff(paste(prefixPCres,"Fig2A","DimPlot","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2A","DimPlot","lowres.tiff", sep = "_"))

  #### Figure 2B ####

DoHeatmap(seurat_sct,
          features = c("POU5F1","NANOG","KLF4","MKI67","TOP2A",
                       "CENPF","AQP4","S100B","GFAP","VIM",
                       "HES1","NFIA","MAPT","SNAP25","TUBB3")) +
  scale_fill_gradientn(colors = c("grey","black","yellow"))
NoLegend()

hirestiff(paste(prefixPCres,"Fig2B","heatmap","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2B","heatmap","lowres.tiff", sep = "_"))

  #### Figure 2C ####
    #### vln plots ####

VlnPlot(seurat_sct,
        features = c("NFIA"),
        pt.size = 0,
        cols = color_palette2) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig2C","vln","NFIA","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","vln","NFIA","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("HES1"),
        pt.size = 0,
        cols = color_palette2) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig2C","vln","HES1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","vln","HES1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("SOX6"),
        pt.size = 0,
        cols = color_palette2) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig2C","vln","SOX6","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","vln","SOX6","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("SNAP25"),
        pt.size = 0,
        cols = color_palette2) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig2C","vln","SNAP25","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","vln","SNAP25","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("RBFOX3"),
        pt.size = 0,
        cols = color_palette2) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig2C","vln","RBFOX3","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","vln","RBFOX3","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("MAPT"),
        pt.size = 0,
        cols = color_palette2) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig2C","vln","MAPT","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","vln","MAPT","lowres.tiff", sep = "_"))

    #### Feature Plots ####

FeaturePlot(seurat_sct,
            reduction = "umap",
            features = c("NFIA"),
            order = TRUE,
            label = FALSE,
            pt.size = 0.5,
            cols = color_palette2[3:4])

hirestiff(paste(prefixPCres,"Fig2C","featureplot","NFIA","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","featureplot","NFIA","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct,
            reduction = "umap",
            features = c("HES1"),
            order = TRUE,
            label = FALSE,
            pt.size = 0.5,
            cols = color_palette2[3:4])

hirestiff(paste(prefixPCres,"Fig2C","featureplot","HES1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","featureplot","HES1","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct,
            reduction = "umap",
            features = c("SOX6"),
            order = TRUE,
            label = FALSE,
            pt.size = 0.5,
            cols = color_palette2[3:4])

hirestiff(paste(prefixPCres,"Fig2C","featureplot","SOX6","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","featureplot","SOX6","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct,
            reduction = "umap",
            features = c("SNAP25"),
            order = TRUE,
            label = FALSE,
            pt.size = 0.5,
            cols = color_palette2[3:4])

hirestiff(paste(prefixPCres,"Fig2C","featureplot","SNAP25","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","featureplot","SNAP25","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct,
            reduction = "umap",
            features = c("RBFOX3"),
            order = TRUE,
            label = FALSE,
            pt.size = 0.5,
            cols = color_palette2[3:4])

hirestiff(paste(prefixPCres,"Fig2C","featureplot","RBFOX3","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","featureplot","RBFOX3","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct,
            reduction = "umap",
            features = c("MAPT"),
            order = TRUE,
            label = FALSE,
            pt.size = 0.5,
            cols = color_palette2[3:4])

hirestiff(paste(prefixPCres,"Fig2C","featureplot","MAPT","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2C","featureplot","MAPT","lowres.tiff", sep = "_"))

  #### Figure 2D ####
# UMAP of cells in each cluster by type 

for (i in 1:length(levels(seurat_sct$Format))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = F,  # clusters added manually
                pt.size = 1,
                cols = color_palette2,
                cells = c(
                  WhichCells(seurat_sct, 
                             expression = Format == levels(seurat_sct$Format)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Format)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"Fig2D","UMAP","split_by","Format",
                  levels(seurat_sct$Format)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"Fig2D","UMAP","split_by","Format",
                   levels(seurat_sct$Format)[i],"lowres.tiff",sep = "_"))
}

rm(p)

  #### Figure 2E ####

DotPlot(seurat_sct,
        features = c("TH"),
        cols = color_palette2[c(10,6)],
        dot.scale = 10,
        scale = T,
        scale.by = "size",
        col.min = 0) + 
  scale_y_discrete(limits=rev) 

hirestiff(paste(prefixPCres,"Fig2E","TH","dotplot","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig2E","TH","dotplot","lowres.tiff", sep = "_"))

  #### Figure 3A ####

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = 1, 
        group.by = "Format",
        cols = color_palette3)

hirestiff(paste(prefixPCres,"Fig3A","DimPlot","by","Format","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3A","DimPlot","by","Format","lowres.tiff", sep = "_"))

  #### Figure 3B ####

VlnPlot(seurat_sct,
        features = c("TH"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","TH","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","TH","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("CHAT"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","CHAT","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","CHAT","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("SST"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","SST","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","SST","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("GAD1"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","GAD1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","GAD1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("TPH1"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","TPH1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","TPH1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("TPH2"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","TPH2","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","TPH2","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("SLC17A7"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","SLC17A7","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","SLC17A7","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("SLC17A8"),
        split.by = "Format",
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3B","vln","SLC17A8","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3B","vln","SLC17A8","lowres.tiff", sep = "_"))

  #### Figure 3C ####

VlnPlot(seurat_sct,
        features = c("NES"),
        split.by = "Format",
        idents = c(5,6),
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3C","vln","NES","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3C","vln","NES","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("SYP"),
        split.by = "Format",
        idents = c(5,6),
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3C","vln","SYP","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3C","vln","SYP","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("MAP2"),
        split.by = "Format",
        idents = c(5,6),
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3C","vln","MAP2","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3C","vln","MAP2","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("DCX"),
        split.by = "Format",
        idents = c(5,6),
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3C","vln","DCX","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3C","vln","DCX","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("DLG4"),
        split.by = "Format",
        idents = c(5,6),
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3C","vln","DLG4","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3C","vln","DLG4","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("RBFOX3"),
        split.by = "Format",
        idents = c(5,6),
        cols = color_palette3) +
  NoLegend()

hirestiff(paste(prefixPCres,"Fig3C","vln","RBFOX3","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3C","vln","RBFOX3","lowres.tiff", sep = "_"))

  #### Figure 3D & 3E ####

expression_data <- FetchData(object = seurat_sct, vars = 
                               c("orig.ident","SCT_snn_res.0.1",
                                 "TH","KCNJ6","SOX6", "LMO3","NR4A2"))

write.csv(expression_data, 
          file = paste(prefixPCres,"Fig3E","Expression_selected_genes","by","genotype.csv", sep = "_"))

  #### Figure 3F ####

# cluster markers 
cluster_markers <- FindAllMarkers(seurat_sct, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.4)

write.csv(cluster_markers, file = paste(prefixPCres,"clustermarkers_pos_neg_minpct0.10_logfcthresh0.4_renumbered.csv", sep = "_"))

# load cluster markers into a dataframe

volc <- read.csv(file = paste(prefixPCres,"clustermarkers_pos_neg_minpct0.10_logfcthresh0.4_renumbered.csv", sep = "_"))

# make df with just cluster 6 data
volc6 <- volc[volc$cluster == "6",]

# set all p_val_adj = 0 to the lowest nonzero value
volc6[volc6$p_val_adj == 0,]$p_val_adj <- min(volc6[volc6$p_val_adj > 0,]$p_val_adj)

# volcano plot
EnhancedVolcano(volc6,
                lab = volc6$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Cluster 6",
                subtitle = NULL,
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ p_val_adj),
                xlim = c(min(volc6$avg_log2FC, na.rm = TRUE) - 0.5, 
                         max(volc6$avg_log2FC, na.rm = TRUE) + 0.5),
                ylim = c(0, max(-log10(volc6$p_val_adj), na.rm = TRUE) + 5),
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "p_val_adj", 
                                 expression(p_val_adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 7,
                legendIconSize = 2,
                drawConnectors = TRUE,
                labSize = 2.75,
                max.overlaps = 5,
                col = c("grey75", "grey65", "grey55", "#DAA520"),
                colAlpha = 2/2,
                pCutoff = 0.05,
                FCcutoff = 1.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 3.5)

hirestiff(paste(prefixPCres,"Fig3F","volcano","cluster","6","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3F","volcano","cluster","6","lowres.tiff", sep = "_"))

  #### Figure 3G ####

# make df with just cluster 5 data 
volc5 <- volc[volc$cluster == "5",]

# set all p_val_adj = 0 to the lowest nonzero value
volc5[volc5$p_val_adj == 0,]$p_val_adj <- min(volc5[volc5$p_val_adj > 0,]$p_val_adj)

# volcano plot
EnhancedVolcano(volc5,
                lab = volc5$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Cluster 5",
                subtitle = NULL,
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ p_val_adj),
                xlim = c(min(volc5$avg_log2FC, na.rm = TRUE) - 0.5, 
                         max(volc5$avg_log2FC, na.rm = TRUE) + 0.5),
                ylim = c(0, max(-log10(volc5$p_val_adj), na.rm = TRUE) + 5),
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "p_val_adj", 
                                 expression(p_val_adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 7,
                legendIconSize = 2,
                drawConnectors = TRUE,
                labSize = 2.75,
                max.overlaps = 10,
                col = c("grey75", "grey65", "grey55", "#088DA5"),
                colAlpha = 2/2,
                pCutoff = 0.05,
                FCcutoff = 1.2,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 3.5)

hirestiff(paste(prefixPCres,"Fig3G","volcano","cluster","5","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Fig3G","volcano","cluster","5","lowres.tiff", sep = "_"))


  #### Supplementary Figure 2C ####

VlnPlot(seurat_sct,
        features = c("CORIN")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2C","vln","CORIN","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2C","vln","CORIN","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("LMX1A")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2C","vln","LMX1A","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2C","vln","LMX1A","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("FOXA2")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2C","vln","FOXA2","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2C","vln","FOXA2","lowres.tiff", sep = "_"))

  #### Supplementary Figure 2D ####

VlnPlot(seurat_sct,
        features = c("BARHL1")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2D","vln","BARHL1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2D","vln","BARHL1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("DBX1")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2D","vln","DBX1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2D","vln","DBX1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("WNT8B")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2D","vln","WNT8B","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2D","vln","WNT8B","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("PITX2")) +
  NoLegend()

hirestiff(paste(prefixPCres,"SupFig2D","vln","PITX2","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"SupFig2D","vln","PITX2","lowres.tiff", sep = "_"))

