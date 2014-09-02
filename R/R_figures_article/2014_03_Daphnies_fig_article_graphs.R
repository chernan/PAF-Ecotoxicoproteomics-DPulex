setwd("/home/chernan/Workspace/DataAnalysis/2014_03_Daphnies_fig_article_2")

do_GOMF <- FALSE

################################################################################

## Load side annotations
data_file_annot <- paste0(getwd(), '/data/Dpulex-189-sign-quant-annotated.txt')
data_annot <- read.delim(data_file_annot, stringsAsFactors=FALSE, sep="\t" ,
                         row.names=NULL, header=TRUE, comment.char="", 
                         na.strings='')
## Labels

col_ac <- "UniProt.ID"
col_names <- "Protein.name"
col_panther_names <- "PANTHER.Family.Subfamily"

## Ratios

col_ratios <-  c("Log.ratio.C1.DMSO.2d", "Log.ratio.C2.DMSO.2d", "Log.ratio.C1.DMSO.7d", "Log.ratio.C2.DMSO.7d") ## "Log.ratio.C2.DMSO.7d" ## 

################################################################################


## For the Cytoscape graphs
## GO

## Full MF:
col_MF <- "GO.Molecular.Function"
interesting_MF_terms <- c("structural constituent of ribosome" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/structural_constituent_of_ribosome.txt'), 
                          "translation regulator activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/translation_regulator_activity.txt'),
                          "translation factor activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/translation_factor_activity.txt'), 
                          "aminoacyl-tRNA ligase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/aminoacyl_tRNA_ligase_activity.txt'), 
                          "nucleic acid binding" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/nucleic_acid_binding.txt'), 
                          "translation initiation factor activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/translation_initiation_factor_activity.txt'), 
                          "RNA binding" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/RNA_binding.txt'), 
                          "oxidoreductase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/oxidoreductase_activity.txt'), 
                          "lyase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/lyase_activity.txt'), 
                          "isomerase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/isomerase_activity.txt'), 
                          "ligase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/ligase_activity.txt'), 
                          "actin binding" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/actin_binding.txt'), 
                          "structural constituent of cytoskeleton" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/structural_constituent_of_cytoskeleton.txt'), 
                          "acetyltransferase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/acetyltransferase_activity.txt'), 
                          "transferase activity" = paste0(getwd(), '/data/panther_sign_go/GOMF-significant/transferase_activity.txt'))

## Full BP:
col_BP <- "GO.Biological.Process"
interesting_BP_terms <- c("translation" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/translation.txt'), 
                          "carbohydrate metabolic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/carbohydrate_metabolic_process.txt'), 
                          "cellular amino acid metabolic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/cellular_amino_acid_metabolic_process.txt'), 
                          "regulation of translation" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/regulation_of_translation.txt'), 
                          "generation of precursor metabolites and energy" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/generation_of_precursor_metabolites_and_energy.txt'), 
                          "cellular amino acid biosynthetic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/cellular_amino_acid_biosynthetic_process.txt'), 
                          "nitrogen compound metabolic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/nitrogen_compound_metabolic_process.txt'), 
                          "oxidative phosphorylation" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/oxidative_phosphorylation.txt'), 
                          "respiratory electron transport chain" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/respiratory_electron_transport_chain.txt'), 
                          "lipid metabolic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/lipid_metabolic_process.txt'), 
                          "fatty acid beta-oxidation" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/fatty_acid_beta-oxidation.txt'), 
                          "coenzyme metabolic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/coenzyme_metabolic_process.txt'), 
                          "glycolysis" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/glycolysis.txt'), 
                          "nucleobase-containing compound transport" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/nucleobase-containing_compound_transport.txt'), 
                          "monosaccharide metabolic process" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/monosaccharide_metabolic_process.txt'),
                          "cellular component organization or biogenesis" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/cellular_component_organization_or_biogenesis.txt'),
                          "tricarboxylic acid cycle" = paste0(getwd(), '/data/panther_sign_go/GOBP-significant/tricarboxylic_acid_cycle.txt'))


## Depending on terms to look for (change do_GOMF at top of the script)
if(do_GOMF) {
    col_terms <- col_MF
    interesting_terms <- interesting_MF_terms
} else {
    col_terms <- col_BP
    interesting_terms <- interesting_BP_terms
}


# file <- interesting_terms[1]
terms_and_proteins <- sapply(interesting_terms, USE.NAMES=TRUE,
       FUN=function(file) {
           data <- read.delim(file, stringsAsFactors=FALSE, sep="\t" ,
                              row.names=NULL, header=FALSE, comment.char="", 
                              na.strings='')
           # DAPPU|Gene=DAPPUDRAFT_307788|UniProtKB=E9G121
           uniprot <- sapply(data[, 1], FUN=function(x){unlist(strsplit(x, split='='))[3]})
           return(uniprot)
       })

unique_protein_ac <- unique(unlist(terms_and_proteins))
nb_unique_prots <- length(unique_protein_ac)

## Table containing which protein is related to which term
# protein <- "E9G927"
prot_in_which_term <- sapply(
    unique_protein_ac,
    FUN=function(protein){
        unlist(lapply(lapply(terms_and_proteins, grepl, pattern=protein, fixed=TRUE), any))
    })


################################################################################
## Compute edge weight

## Distance matrix : percentage of shared terms
percent_shared_terms <- matrix(rep(0, nb_unique_prots**2), 
                               nrow=nb_unique_prots,
                               dimnames=list(unique_protein_ac, unique_protein_ac))
for(iIndex in 1:nb_unique_prots) {
    for(jIndex in 1:nb_unique_prots) {
        percent_shared_terms[iIndex, jIndex] <- sum(as.numeric(prot_in_which_term[, iIndex] & prot_in_which_term[, jIndex])) / sum(as.numeric(prot_in_which_term[, iIndex]))
    }
}

################################################################################
##

## Filter annotations for unique_protein_ac proteins
is_matching_any_term <- data_annot[, col_ac] %in% unique_protein_ac
subset_prots <- data.frame(
    data_annot[is_matching_any_term, c(col_ac, col_names, col_ratios, col_panther_names)],
    row.names=data_annot[is_matching_any_term, col_ac])


################################################################################
## Build graph

# source("http://bioconductor.org/biocLite.R")
# biocLite("RCytoscape")

library(RCytoscape)


# g_simple <- makeSimpleGraph()
# noa.names(g_simple)


## Default attributes
g_daphnies <- new ("graphNEL", edgemode = "directed")

g_daphnies <- initNodeAttribute (g_daphnies, "label", "char", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "nodeType", "char", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "ratio2d1", "numeric", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "ratio2d2", "numeric", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "ratio7d1", "numeric", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "ratio7d2", "numeric", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "terms", "char", "undefined")
g_daphnies <- initNodeAttribute (g_daphnies, "panther", "char", "undefined")

g_daphnies <- initEdgeAttribute (g_daphnies, "edgeType", "char", "undefined")
g_daphnies <- initEdgeAttribute (g_daphnies, "weight", "numeric", "undefined")

## Add nodes : proteins
# protein_name <- "E9GYD9"
for(protein_name in unique_protein_ac) {
    g_daphnies <- addNode(protein_name, g_daphnies)
    nodeData(g_daphnies, protein_name, 'label') <- subset_prots[protein_name, col_names]
    nodeData(g_daphnies, protein_name, 'nodeType') <- "protein"
    nodeData(g_daphnies, protein_name, 'terms') <- paste(names(interesting_terms)[unlist(prot_in_which_term[, protein_name])], collapse='; ')
    nodeData(g_daphnies, protein_name, 'ratio2d1') <- ifelse(is.na(subset_prots[protein_name, col_ratios[1]]), -0, subset_prots[protein_name, col_ratios[1]]) 
    nodeData(g_daphnies, protein_name, 'ratio2d2') <- ifelse(is.na(subset_prots[protein_name, col_ratios[2]]), -0, subset_prots[protein_name, col_ratios[2]]) 
    nodeData(g_daphnies, protein_name, 'ratio7d1') <- ifelse(is.na(subset_prots[protein_name, col_ratios[3]]), -0, subset_prots[protein_name, col_ratios[3]]) 
    nodeData(g_daphnies, protein_name, 'ratio7d2') <- ifelse(is.na(subset_prots[protein_name, col_ratios[4]]), -0, subset_prots[protein_name, col_ratios[4]]) 
    nodeData(g_daphnies, protein_name, 'panther') <- ifelse(is.na(subset_prots[protein_name, col_panther_names]), -0, subset_prots[protein_name, col_panther_names]) 
}

## Add edges between proteins if they share terms
for(iIndex in 1:nb_unique_prots) {
    
    from_prot <- unique_protein_ac[iIndex]
    
    for(jIndex in 1:nb_unique_prots) {
        
        if(iIndex==jIndex) { next }
        
        to_prot <- unique_protein_ac[jIndex]
        
        if(percent_shared_terms[iIndex, jIndex] > 0 & percent_shared_terms[iIndex, jIndex] <= percent_shared_terms[jIndex, iIndex]) {
            if(percent_shared_terms[iIndex, jIndex] < percent_shared_terms[jIndex, iIndex]){
                g_daphnies <- addEdge(from=from_prot, to=to_prot, g_daphnies)
                edgeData(g_daphnies, from=from_prot, to=to_prot, "edgeType") <- "share_terms_with"
                edgeData(g_daphnies, from=from_prot, to=to_prot, "weight") <- percent_shared_terms[iIndex, jIndex]
            }
            else {
                if(iIndex<jIndex) {
                    g_daphnies <- addEdge(from=from_prot, to=to_prot, g_daphnies)
                    edgeData(g_daphnies, from=from_prot, to=to_prot, "edgeType") <- "share_terms_with"
                    
                    if(percent_shared_terms[iIndex, jIndex] == 1) {
                        edgeData(g_daphnies, from=from_prot, to=to_prot, "weight") <- 2
                    }
                    else {
                        edgeData(g_daphnies, from=from_prot, to=to_prot, "weight") <- percent_shared_terms[iIndex, jIndex]
                    }
                }
            }
        }
        
    }
}

g_daphnies
noa.names(g_daphnies)
eda.names(g_daphnies)


################################################################################
## Cytoscape with R

# Start Cytoscape
# Start CytoscapeRCP plugin in Cytoscape 2.8.3


## Create a connection with Cytoscape (CytoscapeRCP must be already started)
cy <- CytoscapeConnection()

################################################################################
## Display new graph in Cytoscape


## Default attributes
setDefaultBackgroundColor(cy, "#FFFFFF")
setDefaultNodeFontSize(cy, 10)
setDefaultEdgeColor(cy, "lightgrey")
setDefaultEdgeLineWidth(cy, 1)


## Create a window

## Delete existing window if any
window.title <- gsub('.', '_', col_terms, fixed=TRUE)
if (window.title %in% as.character (getWindowList(cy)))
    deleteWindow (cy, window.title)
## New window
cw <- new.CytoscapeWindow(window.title, g_daphnies)
# if (window.title %in% as.character (getWindowList(cy)))
#     deleteWindow (cy, window.title)
# cw <- new.CytoscapeWindow(window.title, g_daphnies) # Leave it! In case of creation of another window in an existing session, plugin generates an error the first time... 
setWindowSize(cw, width=800, height=500)
## Link graph display to window
displayGraph(cw)


## Node attributes

## Assign displayed text for nodes

setNodeLabelRule(cw, 'canonicalName')
redraw(cw)

## Set node color

color_attribute <- 'ratio7d2'
data.values <- c(-2, 0, 2)
node.colors <- c("#00EE00", "#BBBBBB","#EE0000")
# setNodeColorRule(cw, 'nodeType', 'protein', c("#9797A5"), mode='lookup')
setNodeColorRule(cw, 
                 node.attribute.name = color_attribute,
                 control.points=data.values,
                 colors=node.colors,
                 mode="interpolate",
                 default="gray")
redraw(cw)

## Set node border color

setNodeBorderColorRule(cw, 'nodeType', 'protein', c("#000000"), mode='lookup')
redraw(cw)

## Edges

## Set bigger edge size for proteins sharing a lot of terms

size_attribute <- 'weight'
data.values <- as.character(sort(unique(eda(getGraph(cw), size_attribute))))
edge.sizes <- seq.int(from=1, to=10, length.out=length(data.values))
setEdgeLineWidthRule(
    cw,
    edge.attribute.name = size_attribute,
    attribute.values=data.values,
    line.widths=edge.sizes,
    default=1
)
redraw(cw)

## Set warmer color for edge among terms sharing a lot of proteins

color_attribute <- 'weight'
data.values <- c(0, max(eda(getGraph(cw), color_attribute)))
edge.colors <- c("#EEEEEE", "#CC0000")
setEdgeColorRule(cw, 
                 edge.attribute.name = color_attribute,
                 control.points=data.values,
                 colors=edge.colors,
                 mode="interpolate",
                 default="gray")

redraw(cw)


## Layout... 
# List of all mappings : unlist(getLayoutNameMapping(cw))

# For a start
# getLayoutNameMapping(cw)[['Spring Embedded Layout']]
# layoutNetwork(cw, getLayoutNameMapping(cw)[['Spring Embedded Layout']])

# Real interesting layout
# getLayoutPropertyNames(cw, getLayoutNameMapping(cw)[['Edge-Weighted Spring Embedded']])
setLayoutProperties(cw, getLayoutNameMapping(cw)[['Edge-Weighted Spring Embedded']], 
                    list(edge_attribute = 'weight'))
layoutNetwork(cw, getLayoutNameMapping(cw)[['Edge-Weighted Spring Embedded']])

## Set window and network sizes
fitContent(cw)
# setZoom(cw, 1.4 * getZoom (cw))

redraw(cw)



## Save visualisation style

new.unique.style.name <- paste(window.title, format(Sys.time(), "%m%d%H%M%S"), round(proc.time()[['elapsed']]), sep='_')
copyVisualStyle(cy, 'default', new.unique.style.name) 
# getVisualStyleNames(cy)
setVisualStyle(cy, new.unique.style.name)
