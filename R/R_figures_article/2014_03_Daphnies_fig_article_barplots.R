setwd("/home/chernan/Workspace/DataAnalysis/2014_03_Daphnies_fig_article_2")

################################################################################

## For the barplots

library(ggplot2)

# plot_title <- "GOBP - All"
# data_file <- data_file_all_GOBP
# pval_threshold <- pval_threshold
# col_terms <- "Biological.Process"
# col_exp <- "Dpulex.189.all..expected."
# col_obs <- "Dpulex.189.all..180."
# col_pval <- "Dpulex.189.all..P.value."


create_bar_plot <- function (plot_title, data_file, pval_threshold, col_terms, col_exp, col_obs, col_pval) {
    
    ## Load data
    data_all <- read.delim(data_file, stringsAsFactors=FALSE, sep="\t" ,
                           row.names=NULL, header=TRUE, comment.char="", 
                           na.strings='')
    
    which_significant <- data_all[ ,col_pval] < pval_threshold & !(data_all[, col_terms]=='Unclassified')
    sign_all <- data_all[which_significant , c(col_pval, col_terms, col_exp, col_obs)]
    
    ggplot2_data <- data.frame(Terms=c(sign_all[, col_terms], sign_all[, col_terms]),
                               Values=c(sign_all[, col_exp], sign_all[, col_obs]),
                               Categ=c(rep(col_exp, sum(which_significant)), rep(col_obs, sum(which_significant))),
                               stringsAsFactors=FALSE
    )
    ## Order Terms by P-Value
    ggplot2_data[, "Terms"] <- factor(ggplot2_data[, "Terms"], levels=sign_all[ order(sign_all[, col_pval]), col_terms] )
    ## Order bars by first Expected then Observed
    ggplot2_data[, "Categ"] <- factor(ggplot2_data[, "Categ"], levels=c(col_exp, col_obs))
    
    ggplot(ggplot2_data, aes(x=factor(Terms), y=Values, fill=Categ)) + 
        geom_bar(stat="identity", position="dodge") + 
        xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1)) +
        ggtitle(plot_title) +
        scale_fill_discrete(name='', breaks=c(col_exp, col_obs), labels=c("Expected", "Observed"))
}

pval_threshold <- 0.05

# data_file <- data_file_all_GOBP
# pval_threshold <- pval_threshold
# col_terms <- "Biological.Process"
# col_exp <- "Dpulex.189.all..expected."
# col_obs <- "Dpulex.189.all..180."
# col_pval <- "Dpulex.189.all..P.value."


pdf("Rplots_barplots.pdf")

data_file_all_GOBP <- paste0(getwd(), '/data/Panther-189-hits-all-GOBP-new.txt')
print(create_bar_plot("GOBP - All", data_file_all_GOBP, pval_threshold, 
                      "Biological.Process", "Dpulex.189.all..expected.", "Dpulex.189.all..180.", "Dpulex.189.all..P.value."))

data_file_all_GOMF <- paste0(getwd(), '/data/Panther-189-hits-all-GOMF-new.txt')
print(create_bar_plot("GOMF - All", data_file_all_GOMF, pval_threshold, 
                      "Molecular.Function", "IDs.189.hits.txt..expected.", "IDs.189.hits.txt..180.", "IDs.189.hits.txt..P.value."))

dev.off()

