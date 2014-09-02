## Full pipeline of the analysis


#########################################################################
## Clear workspace 
rm(list = ls(all.names=TRUE))

#########################################################################
## Run analysis
setwd('/home/chernan/Workspace/DataAnalysis/2013_03_DaphniaPulex')
source("./R/statistical_analysis.R")

#########################################################################
## Generate reports
setwd('/home/chernan/Workspace/DataAnalysis/2013_03_DaphniaPulex')
inputReports <- paste(c(getwd(), '/inst/reports'), collapse='')
outputReports <- paste(c(getwd(), '/outputs/reports'), collapse='')
dir.create(outputReports)

library("knitr")
library("markdown")
oldBasedir <- opts_knit$get("base.dir")
opts_knit$set(base.dir = outputReports)

report <- paste(c(inputReports, '/workflow2dDMSOM4.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow2dDMSOM4.Md'), collapse='')
report2dDMSOM4Html <- paste(c(outputReports, '/workflow2dDMSOM4.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report2dDMSOM4Html)

report <- paste(c(inputReports, '/workflow2dDMSOC1.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow2dDMSOC1.Md'), collapse='')
report2dDMSOC1Html <- paste(c(outputReports, '/workflow2dDMSOC1.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report2dDMSOC1Html)

report <- paste(c(inputReports, '/workflow2dDMSOC2.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow2dDMSOC2.Md'), collapse='')
report2dDMSOC2Html <- paste(c(outputReports, '/workflow2dDMSOC2.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report2dDMSOC2Html)


report <- paste(c(inputReports, '/workflow4dDMSOM4.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow4dDMSOM4.Md'), collapse='')
report4dDMSOM4Html <- paste(c(outputReports, '/workflow4dDMSOM4.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report4dDMSOM4Html)


report <- paste(c(inputReports, '/workflow6dDMSOM4.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow6dDMSOM4.Md'), collapse='')
report6dDMSOM4Html <- paste(c(outputReports, '/workflow6dDMSOM4.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report6dDMSOM4Html)

report <- paste(c(inputReports, '/workflow6dDMSOC1.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow6dDMSOC1.Md'), collapse='')
report6dDMSOC1Html <- paste(c(outputReports, '/workflow6dDMSOC1.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report6dDMSOC1Html)

report <- paste(c(inputReports, '/workflow6dDMSOC2.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/workflow6dDMSOC2.Md'), collapse='')
report6dDMSOC2Html <- paste(c(outputReports, '/workflow6dDMSOC2.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = report6dDMSOC2Html)


report <- paste(c(inputReports, '/comparisonDMSOM4.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/comparisonDMSOM4.Md'), collapse='')
reportCompDMSOM4Html <- paste(c(outputReports, '/comparisonDMSOM4.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = reportCompDMSOM4Html)

report <- paste(c(inputReports, '/comparison2d-6d.Rmd'), collapse='')
reportMd <- paste(c(outputReports, '/comparison2d-6d.Md'), collapse='')
reportCompAllHtml <- paste(c(outputReports, '/comparison2d-6d.html'), collapse='')
knit(report, output = reportMd)
markdownToHTML(file = reportMd, output = reportCompAllHtml)

opts_knit$set(base.dir = oldBasedir)

##########################################################################
## Save session status

setwd('/home/chernan/Workspace/DataAnalysis/2013_03_DaphniaPulex')
outputsFolder <- paste(c(getwd(), '/outputs'), collapse='')

instPackages <- installed.packages()[.packages(),]
write.table(instPackages, paste0(outputsFolder, "/loadedPackages.txt"), quote=TRUE)

sessionInfoFile <- paste0(outputsFolder, "/sessionInfo.txt")
capture.output(sessionInfo(), file=sessionInfoFile)

#########################################################################
## Move/copy files to send to clients

setwd('/home/chernan/Workspace/DataAnalysis/2013_03_DaphniaPulex')
folderToSend <- paste(c(getwd(), '/toSend'), collapse='')
dir.create(folderToSend)

## Copy session info
file.copy(sessionInfoFile, folderToSend)

## Copy analysis results
analysisFolders <- paste0(
    folderToSend, 
    c('/dpulex2dDMSOM4', '/dpulex4dDMSOM4', '/dpulex6dDMSOM4',
      '/dpulex2dDMSOC1', '/dpulex2dDMSOC2',
      '/dpulex6dDMSOC1', '/dpulex6dDMSOC2'))
sapply(analysisFolders, dir.create)

workflowReports <- c(report2dDMSOM4Html, report4dDMSOM4Html, report6dDMSOM4Html,
                     report2dDMSOC1Html, report2dDMSOC2Html,
                     report6dDMSOC1Html, report6dDMSOC2Html)
file.copy(report2dDMSOM4Html, paste(c(analysisFolders[1], 
                                '/summary_analysis_day2_DMSOM4_READFIRST.html'), 
                              collapse=''))
file.copy(report4dDMSOM4Html, paste(c(analysisFolders[2], 
                                '/summary_analysis_day4_DMSOM4_READFIRST.html'), 
                              collapse=''))
file.copy(report6dDMSOM4Html, paste(c(analysisFolders[3], 
                                '/summary_analysis_day6_DMSOM4_READFIRST.html'), 
                              collapse=''))
file.copy(report2dDMSOC1Html, paste(c(analysisFolders[4], 
                                      '/summary_analysis_day2_DMSOC1_READFIRST.html'), 
                                    collapse=''))
file.copy(report2dDMSOC2Html, paste(c(analysisFolders[5], 
                                      '/summary_analysis_day2_DMSOC2_READFIRST.html'), 
                                    collapse=''))
file.copy(report6dDMSOC1Html, paste(c(analysisFolders[6], 
                                      '/summary_analysis_day6_DMSOC1_READFIRST.html'), 
                                    collapse=''))
file.copy(report6dDMSOC2Html, paste(c(analysisFolders[7], 
                                      '/summary_analysis_day6_DMSOC2_READFIRST.html'), 
                                    collapse=''))

analysisFiles <- c(
    "/dpulex2dDMSOM4/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex2dDMSOM4/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex2dDMSOM4/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[1],
                 c("/full_output_day2_DMSOM4_data.csv", 
                   "/summary_analysis_day2_DMSOM4_data.csv", 
                   "/full_output_day2_DMSOM4.html"))
)

analysisFiles <- c(
    "/dpulex4dDMSOM4/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex4dDMSOM4/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex4dDMSOM4/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[2],
                 c("/full_output_day4_DMSOM4_data.csv", 
                   "/summary_analysis_day4_DMSOM4_data.csv", 
                   "/full_output_day4_DMSOM4.html"))
)

analysisFiles <- c(
    "/dpulex6dDMSOM4/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex6dDMSOM4/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex6dDMSOM4/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[3],
                 c("/full_output_day6_DMSOM4_data.csv", 
                   "/summary_analysis_day6_DMSOM4_data.csv", 
                   "/full_output_day6_DMSOM4.html"))
)

analysisFiles <- c(
    "/dpulex2dDMSOC1/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex2dDMSOC1/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex2dDMSOC1/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[4],
                 c("/full_output_day2_DMSOC1_data.csv", 
                   "/summary_analysis_day2_DMSOC1_data.csv", 
                   "/full_output_day2_DMSOC1.html"))
)
analysisFiles <- c(
    "/dpulex2dDMSOC2/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex2dDMSOC2/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex2dDMSOC2/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[5],
                 c("/full_output_day2_DMSOC2_data.csv", 
                   "/summary_analysis_day2_DMSOC2_data.csv", 
                   "/full_output_day2_DMSOC2.html"))
)

analysisFiles <- c(
    "/dpulex6dDMSOC1/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex6dDMSOC1/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex6dDMSOC1/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[6],
                 c("/full_output_day6_DMSOC1_data.csv", 
                   "/summary_analysis_day6_DMSOC1_data.csv", 
                   "/full_output_day6_DMSOC1.html"))
)
analysisFiles <- c(
    "/dpulex6dDMSOC2/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_full.csv",
    "/dpulex6dDMSOC2/vsn_lpe/MaxQuant_LFQ_non-linked_vsnlpe_light.csv",
    "/dpulex6dDMSOC2/vsn_lpe/output_vsn_lpe.html")
outputReports <- paste0(getwd(), '/outputs/0.05', analysisFiles)
file.copy(outputReports, 
          paste0(analysisFolders[7],
                 c("/full_output_day6_DMSOC2_data.csv", 
                   "/summary_analysis_day6_DMSOC2_data.csv", 
                   "/full_output_day6_DMSOC2.html"))
)

## Copy general summary reports
file.copy(reportCompDMSOM4Html, 
          paste0(folderToSend, '/comparison_DMSO-M4_2-4-6days.html'))
file.copy(paste0(outputsFolder, '/significantResultsDMSOM4.csv'), 
          paste0(folderToSend, '/comparison_DMSO-M4_2-4-6days_data.csv'))

file.copy(reportCompAllHtml, 
          paste0(folderToSend, '/comparison_All_2-4-6days.html'))
file.copy(paste0(outputsFolder, '/significantResults2-4-6dClustered.csv'), 
          paste0(folderToSend, '/comparison_All_2-4-6days_data.csv'))

## Manual step to be done:
## Add UniProt annotations from 'Retrieve' service.
## Tried using biomaRt but this THING is working 1/1000 times... Not reliable for a pipeline.