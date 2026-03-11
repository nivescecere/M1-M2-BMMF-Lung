#Macrophages and Monocytes

rm(list=ls())
Sys.setenv(LANG = "en")
library(Seurat)
library(readxl)
library(umap)
library(edgeR)

#Subset Macrophages and Monocytes
load("IM01_Immune_Seurat_object_nodups-013.RData")

tiss_immune
mf_monocytes <- subset(tiss_immune, subset = immune_subtype_annotation == "MF-Monocytes")

#View(mf_monocytes@meta.data)

#New column in data frame

exceldata=read_excel("List_for_metadata_only_lung_bmmf1.xlsx")
macrophages_list <- data.frame(exceldata)

#View(mf_monocytes@meta.data[["bmmf"]])
#mf_monocytes
#View(mf_monocytes@meta.data)

pdat=mf_monocytes@meta.data
#View(pdat)
#pdat[["bmmf"]] <- macrophages_list$bmmf
table(pdat$bmmf)


# Subset MF in lung
mf_monocytes <- subset(mf_monocytes, subset = biopsy_site == "Lung")

mf_monocytes@meta.data[["bmmf"]] <- macrophages_list$bmmf
pdat=mf_monocytes@meta.data
table(pdat$bmmf)
#mf_monocytes=NormalizeData(mf_monocytes) already done
mf_monocytes <- SetIdent(mf_monocytes, value='bmmf')

#most variable features
mf_monocytes <- FindVariableFeatures(mf_monocytes, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mf_monocytes), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mf_monocytes)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

-----
  library(ggplot2)

ggplot(df_combined, aes(x = tSNE_1, y = tSNE_2, color = bmmf)) +
  # Plot points with the color aesthetic mapped to bmmf
  geom_point(data = subset(df_combined, bmmf == "Negative"), alpha = 0.4, size = 1.0) +
  geom_point(data = subset(df_combined, bmmf == "Positive"), size = 1.0) +
  theme_minimal() +
  theme_classic() +
  theme(panel.grid.major = element_blank(),   # Remove major grid lines
        panel.grid.minor = element_blank(),   # Remove minor grid lines
        legend.title = element_blank()) +  # Remove the legend title
  scale_color_manual(values = c("Negative" = "grey", "Positive" = "red")) +  # Manually specify the colors
  labs(color = NULL) +  # Remove the "Group" label from the legend
  guides(color = guide_legend(override.aes = list(size = 3)))  # Increase the size of the dots in the legend



-------
DimPlot(mf_monocytes, group.by = "bmmf", cols = c("grey", "red"), )
DimPlot(mf_monocytes, group.by = "sample_name")
DimPlot(mf_monocytes,  group.by = "patient_id")

#detection rate covariate
L=GetAssayData(object = mf_monocytes, slot = "counts") 
cdr <- scale(colMeans(L> 0))
mf_monocytes[["cdr"]] <- cdr

# Different expression test
levels(mf_monocytes)
markers <- FindMarkers(mf_monocytes, ident.1 = "Positive",test.use="MAST",latent.vars =c("sample_name","cdr"))
View(markers)

library("writexl")

markers <- cbind(" "=rownames(markers), markers)
write_xlsx(markers,"Markers_BMMF1_only_lung.xlsx")


# M1 and M2
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

#genes
M1 = c("CCR7", "IL2RA", "IL15RA", "IL7R", "CXCL11", "CCL19", "CXCL10", "CXCL9", "TNF", "CCL5", "CCL15", "IL12B", "IL15RA", "TRAIL", "IL6", "CCL20", "PBEF1", "ECGF1", "BCL2A1", "FAS", "BIRC3", "GADD45G", "HSXIAPAF1", "SLC7A5", "SLC21A15", "SLC2A6", "SLC31A2", "INDO", "PLA1A", "OASL", "CHI3L2", "HSD11B1", "AK3", "SPHK1", "PFKFB3", "PSME2", "PFKP", "PSMB9", "PSMA2", "OAS2", "PTX3", "CSPG2", "APOL3", "IGFBP4", "APOL1", "PDGFA", "EDN1", "APOL2", "INHBA", "APOL6", "HESX1", "IRF1", "ATF3", "IRF7")
M2 = c("GPR86", "P2RY5", "TGFBR2", "HRH1", "TLR5", "DCL-1", "MSR1", "CXCR4", "DECTIN1", "P2RY14", "DCSIGN", "CLECSF13", "MS4A4A", "MRC1", "IGF1", "CCL23", "CCL18", "CCL13", "SLC21A9", "SLC4A7", "SLC38A6", "CTSC", "HEXB", "LIPA", "ADK", "HNMT", "TPST2", "CERK", "HS3ST2", "LTA4H", "CA2", "ALOX15", "HS3ST1", "TGFBI", "SEPP1", "CHN2", "FN1", "FGL2", "GAS7", "EGR2", "MAF")


mf_monocytes <- AddModuleScore(mf_monocytes,
                               features = list(M1,M2),
                               name="M")
#FeaturePlot(mf_monocytes,
#features = c("M1","M2"), label = TRUE, repel = TRUE) +
#scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Compare M1 and M2 for BMMF
mf_monocytes$seurat_clusters <- NULL
FeaturePlot(mf_monocytes,
            features = "M2", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(mf_monocytes,
            features = "M1", label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

par(mfrow=c(1,2))

datM1M2=data.frame('Score'=c(mf_monocytes@meta.data$M1,mf_monocytes@meta.data$M2),
                   'Type'=rep(c('M1','M2'),each=nrow(mf_monocytes@meta.data)),
                   'BMMF'=rep(mf_monocytes@meta.data$bmmf,2))

p1 <- ggplot(datM1M2, aes(x=BMMF, y=Score)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + facet_wrap(as.factor(datM1M2$Type))
p1

library(robustbase)

#robust linear regression
lmM1=lmrob(mf_monocytes@meta.data$M1~mf_monocytes@meta.data$bmmf+mf_monocytes@meta.data$sample_name,setting = "KS2014",fast.s.large.n = Inf)
summary(lmM1)$coef[2,]
lmM2=lmrob(mf_monocytes@meta.data$M2~mf_monocytes@meta.data$bmmf+mf_monocytes@meta.data$sample_name,setting = "KS2014",fast.s.large.n = Inf)
summary(lmM2)$coef[2,]
