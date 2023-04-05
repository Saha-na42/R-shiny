install.packages("Seurat")
install.packages("ggplot")
install.packages("shinybusy")
install.packages("tidyr")
install.packages("glue")
install.packages("ggthemes")
install.packages("DT")
install.packages("shinydashboardPlus")


library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
 library(shinythemes)
 library(markdown)
 library(Seurat)
library(ggplot2)
library(shinybusy)
library(tools)
library(tidyr)
library(ggthemes)
library(DT)
library(glue)


 load_seurat_obj <- function(path){
    errors <- c()
    # check file extension
    if (!tolower(tools::file_ext(path)) == "rds") { # ignores case
        errors <- c(errors, "Invalid rds file.")
        return(errors)
    }
 

 # to read in file
    tryCatch(
        {
        obj <- readRDS(path)
        },
        error = function(e) {
            errors <- c(errors, "Invalid rds file.")
            return(errors)
        }
    )
# Validate obj is a seurat object
    if (!inherits(obj, "Seurat")){
        errors <- c(errors, "File is not a seurat object")
        return(errors)
    }

    return(obj)
}

 

 create_metadata_UMAP <- function(obj, col){
    if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")){
    col_df <- data.frame(obj@reductions$umap@cell.embeddings, data = obj@meta.data[,col])
    umap <- ggplot(data = col_df) +
        geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(data)), size = 0.01) + 
        scale_colour_gradientn(colours = rainbow(7))
    } else if (col %in% colnames(obj@meta.data)) {
    umap <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "umap")
    } else {
    umap <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    return(umap)
}


create_metadata_tSNE <- function(obj, col){
    if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")){
        col_df <- data.frame(obj@reductions$tsne@cell.embeddings, data = obj@meta.data[,col])
        tsne <- ggplot(data = col_df) +
            geom_point(mapping = aes(TSNE_1, TSNE_2, color = log10(data)), size = 0.01) + 
            scale_colour_gradientn(colours = rainbow(7))
    } else if (col %in% colnames(obj@meta.data)) {
        tsne <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "tsne")
    } else {
        tsne <- ggplot() +
            theme_void() +
            geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    return(tsne)
}
