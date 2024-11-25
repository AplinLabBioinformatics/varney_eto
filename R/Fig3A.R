


library(Seurat)
library(GSEABase)

# seu <- readRDS("/GSE139829_seuratObj.Rds")

gmtIn <- getGmt("PSscore.2.gmt")

gsetList <- ls()
for( i in 1:length(gmtIn@.Data)){
  gsetList[[gmtIn@.Data[[i]]@setName]] <- gmtIn@.Data[[i]]@geneIds
}

seu <- AddModuleScore(seu, features = gsetList, name = names(gsetList))

signatures_list <- list()
signatures_list[["Scores_df"]] <- seu$meta.data



# Fig 3A

library(ggplot2)
library(dplyr)

# Scored Data
signatures_list = readRDS("GSE139829_Varney_ETO.Rds")


Scores_df <- signatures_list$Scores_df

# Name of metadata column with grouping values (usually treatment timepoint or treated/untreated) - Factor level
treatmentGroupLbl <- "BAP1_Status"

# Name of metadata column with sample identifier (a violin will be generated for each unique value) - NOT a factor
cellSampleLbl <- "Tumor"

# Name of new metadata annotation that will have the samples ordered for plotting
factorLabelName <- "Tumor.by.BAP1_Status"

# The metadata column name with the values to use for sorting the samples
getMedVarLbl <- "CPT1AinhibDN3"

# Create object with each sample's median value calculated
checkCellannots <- Scores_df %>% 
  group_by(!!sym(cellSampleLbl)) %>%
  summarise(median_Val = median(!!sym(getMedVarLbl)))

# Clean up object to only have one sample per row by removing duplicates
uCA <- unique(Scores_df[,c(treatmentGroupLbl, cellSampleLbl)])
rownames(uCA) <- uCA[,cellSampleLbl]

# Add grouping variable to object
checkCellannots$Treatment <- uCA[unname(unlist(checkCellannots[,cellSampleLbl])), treatmentGroupLbl]

# Re-order rows by group and sample median values
tumorOrderViolin <- checkCellannots[with(checkCellannots, order(Treatment, -median_Val)), cellSampleLbl]
tumorOrderViolin <- unname(unlist(tumorOrderViolin[,cellSampleLbl]))

# Add sample order to annotation
Scores_df[,factorLabelName] <- factor(Scores_df[,cellSampleLbl], levels = tumorOrderViolin)


# sample median
sampleMedians <- as.data.frame(checkCellannots)
rownames(sampleMedians) <- sampleMedians[,cellSampleLbl]
colnames(sampleMedians)[colnames(sampleMedians) %in% "median_Val"] <- paste0(getMedVarLbl, "_median")
colnames(sampleMedians)[colnames(sampleMedians) %in% "Treatment"] <- treatmentGroupLbl

rm(tumorOrderViolin, checkCellannots, uCA)

tmpFNout1 <- file.path("/Varney et. al. ETO", 
                       paste(format(Sys.time() ,"%Y%m%d.%H%M"), "VlnPlot", "groupBy", factorLabelName, getMedVarLbl, "png", sep = "."))

ggp2ordrd <- ggplot(Scores_df)

ggp2ordrd +  
  geom_violin(aes_string(x = factorLabelName, y = getMedVarLbl, fill = treatmentGroupLbl), 
              scale = 'width', adjust = 1, trim = TRUE) +
  geom_boxplot(aes_string(x = factorLabelName, y = getMedVarLbl), width = .2 ) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # +
  # labs(y="CPT1A-inhib. down signature score")

ggsave(tmpFNout1, width = 7, height = 5, dpi = "retina")

