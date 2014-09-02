library(Biobase)

#'#' Load data
getDatasetDpulex246daysAllpg <- function(basedir=getwd()) {
    
    data_file <- paste0(basedir, '/data/Dpulex-MXQ-all-A-S/proteinGroups.txt')
    
    ## Load data from MQ iTRAQ
    temp_proteins <- read.delim(data_file, stringsAsFactors=FALSE, 
                                row.names=NULL, header=TRUE)
    
    ## Select intensities
    reporterIntHeaders2dDMSOM4 <- c(
        "LFQ.intensity.DMSO.2d.G", "LFQ.intensity.DMSO.2d.H", 
        "LFQ.intensity.M4.2d.A", "LFQ.intensity.M4.2d.B"
    )
    exprs2dDMSOM4 <- temp_proteins[, reporterIntHeaders2dDMSOM4]
    reporterIntHeaders2dDMSOC1 <- c(
        "LFQ.intensity.DMSO.2d.G", "LFQ.intensity.DMSO.2d.H", 
        "LFQ.intensity.Tamox.C1.2d.M", "LFQ.intensity.Tamox.C1.2d.N"
    )
    exprs2dDMSOC1 <- temp_proteins[, reporterIntHeaders2dDMSOC1]
    reporterIntHeaders2dDMSOC2 <- c(
        "LFQ.intensity.DMSO.2d.G", "LFQ.intensity.DMSO.2d.H", 
        "LFQ.intensity.Tamox.C2.2d.Pp", "LFQ.intensity.Tamox.C2.2d.Q"
    )
    exprs2dDMSOC2 <- temp_proteins[, reporterIntHeaders2dDMSOC2]
    
    reporterIntHeaders4dDMSOM4 <- c(
        "LFQ.intensity.DMSO.4d.I", "LFQ.intensity.DMSO.4d.J",
        "LFQ.intensity.M4.4d.C", "LFQ.intensity.M4.4d.D"
    )
    exprs4dDMSOM4 <- temp_proteins[, reporterIntHeaders4dDMSOM4]

    reporterIntHeaders6dDMSOM4 <- c(
        "LFQ.intensity.DMSO.6d.K", "LFQ.intensity.DMSO.7d.L", 
        "LFQ.intensity.M4.7d.E", "LFQ.intensity.M4.7d.F"
    )
    exprs6dDMSOM4 <- temp_proteins[, reporterIntHeaders6dDMSOM4]
    reporterIntHeaders6dDMSOC1 <- c(
        "LFQ.intensity.DMSO.6d.K", "LFQ.intensity.DMSO.7d.L", 
        "LFQ.intensity.Tamox.C1.7d.O", "LFQ.intensity.Tamox.C1.7d.P"
    )
    exprs6dDMSOC1 <- temp_proteins[, reporterIntHeaders6dDMSOC1]
    reporterIntHeaders6dDMSOC2 <- c(
        "LFQ.intensity.DMSO.6d.K", "LFQ.intensity.DMSO.7d.L", 
        "LFQ.intensity.Tamox.C2.8d.R", "LFQ.intensity.Tamox.C2.8d.S"
    )
    exprs6dDMSOC2 <- temp_proteins[, reporterIntHeaders6dDMSOC2]
        
    ## Phenotypic data
    pData2dDMSOM4 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData2dDMSOM4) <- reporterIntHeaders2dDMSOM4
    pData2dDMSOC1 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData2dDMSOC1) <- reporterIntHeaders2dDMSOC1
    pData2dDMSOC2 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData2dDMSOC2) <- reporterIntHeaders2dDMSOC2

    pData4dDMSOM4 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData4dDMSOM4) <- reporterIntHeaders4dDMSOM4

    pData6dDMSOM4 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData6dDMSOM4) <- reporterIntHeaders6dDMSOM4
    pData6dDMSOC1 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData6dDMSOC1) <- reporterIntHeaders6dDMSOC1
    pData6dDMSOC2 <- data.frame(Design=c(rep("Control", 2), rep("Experiment", 2)))
    rownames(pData6dDMSOC2) <- reporterIntHeaders6dDMSOC2
    
    ## Feature data
    infoHeaders <- c("Majority.protein.IDs", "Fasta.headers", "PEP", 
                     "Peptides", "Only.identified.by.site", "Reverse", 
                     "Contaminant")
    featuredata <- temp_proteins[, infoHeaders]
    names(featuredata) <- c("Majority.protein.IDs", "Fasta.headers", "PEP", 
                            "Peptides", "Only.identified.by.site", "Reverse", 
                            "Contaminant")
    
    ## ExpressionSet construction
    data_proteins2dDMSOM4 <- new("ExpressionSet", 
                           exprs = exprs2dDMSOM4, 
                           phenoData = new("AnnotatedDataFrame", pData2dDMSOM4), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    data_proteins2dDMSOC1 <- new("ExpressionSet", 
                           exprs = exprs2dDMSOC1, 
                           phenoData = new("AnnotatedDataFrame", pData2dDMSOC1), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    data_proteins2dDMSOC2 <- new("ExpressionSet", 
                           exprs = exprs2dDMSOC2, 
                           phenoData = new("AnnotatedDataFrame", pData2dDMSOC2), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )

    data_proteins4dDMSOM4 <- new("ExpressionSet", 
                           exprs = exprs4dDMSOM4, 
                           phenoData = new("AnnotatedDataFrame", pData4dDMSOM4), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )

    data_proteins6dDMSOM4 <- new("ExpressionSet", 
                           exprs = exprs6dDMSOM4, 
                           phenoData = new("AnnotatedDataFrame", pData6dDMSOM4), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    data_proteins6dDMSOC1 <- new("ExpressionSet", 
                           exprs = exprs6dDMSOC1, 
                           phenoData = new("AnnotatedDataFrame", pData6dDMSOC1), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    data_proteins6dDMSOC2 <- new("ExpressionSet", 
                           exprs = exprs6dDMSOC2, 
                           phenoData = new("AnnotatedDataFrame", pData6dDMSOC2), 
                           featureData = new("AnnotatedDataFrame", featuredata)
    )
    
    return(
        list(
            dpulex2dDMSOM4 = list(dataset=data_proteins2dDMSOM4, type="prot"),
            dpulex2dDMSOC1 = list(dataset=data_proteins2dDMSOC1, type="prot"),
            dpulex2dDMSOC2 = list(dataset=data_proteins2dDMSOC2, type="prot"),
            dpulex4dDMSOM4 = list(dataset=data_proteins4dDMSOM4, type="prot"),
            dpulex6dDMSOM4 = list(dataset=data_proteins6dDMSOM4, type="prot"),
            dpulex6dDMSOC1 = list(dataset=data_proteins6dDMSOC1, type="prot"),
            dpulex6dDMSOC2 = list(dataset=data_proteins6dDMSOC2, type="prot")
        )
    )
    
}



filterDataset <- function(expSet, threshold=0) {
    
    ## Remove contaminants
    whichOkProteins <- 
        pData(featureData(expSet))[["Only.identified.by.site"]] != '+' & 
        pData(featureData(expSet))[["Reverse"]] != '+' &
        pData(featureData(expSet))[["Contaminant"]] != '+'
    
    ## Remove non-complete observations
    whichOkIntensitites <- 
        apply(exprs(expSet), 1, FUN=function(x){all(!is.na(x))}) &
        apply(exprs(expSet), 1, FUN=function(x){all(x>threshold)})   
    
    return(expSet[(whichOkProteins & whichOkIntensitites), ])
    
}

## Save as tsv
exportDataset <- function(expSet, normExpSet, computedVals, outFolder) {
    
    results <- cbind(
        pData(featureData(expSet))[, c("Majority.protein.IDs", "Fasta.headers", 
                                      "PEP", "Peptides")],
        exprs(expSet),
        exprs(normExpSet),
        computedVals
    )
    
    #Modify format of UniProt IDs for further displays
    results[,"Majority.protein.IDs"] <- gsub(";", "; ", 
                                             results[,"Majority.protein.IDs"])
    
    #Improve headers
    renamedHeaders <- c(
        "UniProt_IDs", "Fasta_headers", "MaxQuant_PEP", "Peptides", 
        NA, NA, NA, NA, NA, NA, NA, NA, 
        "p-value", "Log2(fold_change)", 
        "Outlier_Ctrl", "Outlier_Exp", 
        "p-value_BH_corrected"
    )
    names(results)[1:4] <- renamedHeaders[1:4]
    names(results)[13:17] <- renamedHeaders[13:17]
    
    tempOutputFull <- paste(
        c(outFolder, '/MaxQuant_LFQ_non-linked_vsnlpe_full.csv'), 
        collapse='')
    write.table(results, tempOutputFull, sep="\t") 

    ## For clients
    clientHeaders <- c(
        "UniProt_IDs", "Fasta_headers", "MaxQuant_PEP", "Peptides", 
        "Log2(fold_change)", "p-value_BH_corrected",
        "Outlier_Ctrl", "Outlier_Exp"
    )
    tempOutputLight <- paste(
        c(outFolder, '/MaxQuant_LFQ_non-linked_vsnlpe_light.csv'), 
        collapse='')
    ## round p-values
    results[, "Log2(fold_change)"] <- 
        round(results[, "Log2(fold_change)"], digits=4)
    results[, "p-value_BH_corrected"] <- 
        round(results[, "p-value_BH_corrected"], digits=4)
    write.table(results[, clientHeaders], tempOutputLight, sep="\t") 
    
}
