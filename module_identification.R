library(Seurat)
library(SeuratDisk)
library(scater)
library(SCINET)
library(ACTIONet)
library(SingleCellExperiment)
library(stringr)
options(stringsAsFactors=FALSE)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(repr)
library(igraph)
library(ComplexHeatmap)
source("./network_construction.R")
library("gprofiler2")

perform_pathway_enrichment <- function(gene_list, pathway_gene_list, pop_size){
    # Calculate the intersection of gene lists
    intersection = intersect(gene_list, pathway_gene_list)
    len_geneList = length(gene_list)
    len_pathway_geneList = length(pathway_gene_list)
    len_intersection = length(intersection)
    
    # Perform hypergeometric test
    p_value = phyper(q= len_intersection -1, m= len_pathway_geneList, 
                     n= pop_size- len_pathway_geneList, k= len_geneList, lower.tail=FALSE)
    returnObj <- list(p_value, intersection)
    return (returnObj)
}

get_GEM_enrichment <- function(cluster_num){
    db_size <-  # Number of features after scTransform normalisation
    pvalue_lst <- c()
    pathway_lst <- c()
    common_gene_lst <- list()
    pcnet.module.genes =V(pcnet_malignant_network)$name[pcnet_clusters == cluster_num]
    for (idx in (1:length(names(gem_pathways_lst)))){
        enrichment <- perform_pathway_enrichment(pcnet.module.genes, gem_pathways_lst[[idx]], db_size)
        pvalue <- enrichment[[1]]
        pvalue_lst <- c(pvalue_lst, pvalue)
        pathway_lst <- c(pathway_lst, names(gem_pathways_lst)[idx])
        if (length(enrichment[[2]]) > 2){
            common_gene_lst[[ names(gem_pathways_lst)[idx]]] <-  enrichment[[2]]
        }
        
    }
    adjp_values <- p.adjust(pvalue_lst, method = p.adjust.methods[5], n = length(pvalue_lst))
    enrichment.df<- data.frame(pathway_lst, pvalue_lst, adjp_values)
    enrichment.df <- filter(enrichment.df, pvalue_lst < 0.05)
    returnObj <- list(enrichment.df, common_gene_lst)
    return (returnObj)
    
}

get_kegg_enrichment <- function(cluster_num){
    pcnet.module.genes = V(pcnet_malignant_network)$name[pcnet_clusters == cluster_num]

    # Perform functional enrichment on the largest gene module
    pcnet.Mcell.module.enrichment = gost(pcnet.module.genes, sources = "KEGG")
    # Here , p-value represents  hypergeometric p-value after correction for multiple testing
    pcnet.terms.tbl = pcnet.Mcell.module.enrichment$result[, c("term_name", "p_value")]
    pcnet.terms.tbl
}

print_genes <- function(genes){
    for (g in genes){
        cat(g, "\n")
    }
}

outDir <- file.path("EEC_results")

## Preprocessing using ACTIONet framework
sce = readRDS(file = file.path(outDir, "reduced_sce.RDS"))
ACTIONet.out = readRDS(file = file.path(outDir, "ACTIONet_output.RDS"))
ACTIONet.out <- clusterFeatureSpecificity(ACTIONet.out, ACTIONet.out$cell_type_annotation, 
                                          "celltype_specificity_scores")

# Network Construction
pcnet_Celltype.specific.networks = run.SCINET.clusters(ACTIONet.out,
                                        specificity.slot.name = "celltype_specificity_scores_feature_specificity")
pcnet_malignant_network <- pcnet_Celltype.specific.networks$Malignant
# Save the network in Labeled edgelist
# write_graph(pcnet_malignant_network, file = file.path(outDir, "pcnet_malignant_ppi.txt"), format = "ncol") 

# Module Identification                      
pcnet_clusters = unsigned_cluster(as(get.adjacency(pcnet_malignant_network), "sparseMatrix"), 
                                  resolution_parameter = 2)
print(table(pcnet_clusters))



