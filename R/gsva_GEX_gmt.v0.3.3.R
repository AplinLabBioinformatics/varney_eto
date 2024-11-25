#!/usr/bin/env Rscript

# EXPRFILE="~/RNA-CancerCell-MORRISON1-combat_batch_corrected-logcpm-all_samples.tsv"
# GMTFILE="~/gmts/PSscore.gmt"
# OUTPUTDIR="~/MORRISON_CancerCell"
# GSVATOFILE="PSscore.MORRISON_CancerCell.tsv"

# Rscript ./gsva_GEX_gmt.v0.3.3.R --ExprMat $EXPRFILE --gmt $GMTFILE --outDir $OUTPUTDIR --gsvaToFile $GSVATOFILE

pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                            "--ExprMat /path/to/geneBySampleMatrix.tsv ",
                                            "--gmt /path/to/geneSets.gmt ",
                                            "--outDir /path/to/results",
                                            "--gsvaToFile filenameOut "
))

pargs <- optparse::add_option(pargs, c("--ExprMat"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="exprMat",
                              metavar="exprMat",
                              help=paste(".txt or .tsv file with genes as rows and samples as columns."))

pargs <- optparse::add_option(pargs, c("--gmt"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="gmtFile",
                              metavar="gmtFile",
                              help=paste(".gmt file with gene sets."))

pargs <- optparse::add_option(pargs, c("--outDir"),
                              type="character",
                              action="store",
                              default=NULL,
                              dest="baseOutputDirectory",
                              metavar="baseOutputDirectory",
                              help=paste("Output directory for generating heatmaps."))

pargs <- optparse::add_option(pargs, c("--gsvaToFile"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="gsvaToFile",
                              metavar="gsvaToFile",
                              help=paste("name of file to save gsva object"))

args <- optparse::parse_args(pargs)

tempTimestampTime <- format(Sys.time(), "%Y%m%d.%H%M")

outputDir <- file.path(args$baseOutputDirectory)
if(!dir.exists(outputDir)){ dir.create(outputDir, recursive = T) }

library("GSEABase")
library("GSVA")

ntcts2 <- read.table(args$exprMat, sep = "\t", header = T, row.names = 1)
ntcts2check <- read.table(args$exprMat, sep = "\t", header = T, row.names = 1, check.names = F)

storeNames <- data.frame(inNames = colnames(ntcts2), origNames = colnames(ntcts2check), row.names = colnames(ntcts2))

gmtIn <- getGmt(args$gmtFile)

gsvaOutT <- gsva(as.matrix(ntcts2), gmtIn)

colnames(gsvaOutT) <- storeNames[colnames(gsvaOutT), "origNames"]

save(list = gsvaOutT, file = file.path(outputDir, paste0(tempTimestampTime, ".", args$gsvaToFile,".Rdata")))

writeLines(paste0("Complete. Files are in ", outputDir))
