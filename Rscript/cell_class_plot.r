#!/usr/bin/env/ Rscript

#Author : yangfei
#data : 2022.6.24


library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(reshape2)
library(SPARK)

echo <- function(...){
    cat(...,"\n")
}


print.usage <- function(){
    echo("Plot the cell class")
    echo()
    echo("Usage:")
    echo("   $ ./spark_analysis.r counts.txt  cell_location.txt  gene.json")
    echo()
}

read.exptable <- function(path){
    exptable <- read.table(path)
    return(exptable)
}


read.cell_location <- function(path){
    cell_location <- read.table(path)
    return(cell_location)
}

read.config <- function(path){
    library("rjson")
    json.data <- fromJSON(file=path)
    return(json.data)
}

color <- function(class_num){
    i = 1
    color_list = c()
    while(i <= class_num){
        i = i+1
        code_list = c('0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F')
        frist_code = sample(code_list,1)
        second_code = sample(code_list,1)
        thrid_code = sample(code_list,1)
        forth_code = sample(code_list,1)
        fifth_code = sample(code_list,1)
        sixth_code = sample(code_list,1)
        seventh_code = sample(code_list,1)
        eighth_code = sample(code_list,1)
        A = paste0("#",frist_code,second_code,thrid_code,forth_code,fifth_code,sixth_code,seventh_code,eighth_code)
        color_list <- append(color_list,A)
    }
    return(color_list)
}


cell_class <- function(){
    
    args <- commandArgs(trailingOnly = TRUE)
    print.usage()
    cell_location.path <- args[2]

    #read cell location table
    cell_location <- read.cell_location(cell_location.path)
    pd <- cbind.data.frame(x = cell_location$x,y = cell_location$y,cell_label = cell_location$cell_label)
    cellname <- as.character(unique(pd$cell_label))

    tmp <- LETTERS[1:length(cellname)]
    names(tmp) <- cellname
    pd$label_idx <- tmp[as.character(pd$cell_label)]
    print(pd[1:4,])

    #set figure size
    min.pand = 0.99
    max.pand = 1.01
    pointsize = 1
    titlesize = 1

    library(ggplot2)
    library(RColorBrewer)
    
    pal <- colorRampPalette(c("mediumseagreen","lightyellow2","deeppink"))(length(cellname))

    class_dist <- ggplot(pd, aes(x = x, y = y, color = label_idx)) + geom_point(size = pointsize) + 
    scale_colour_manual(name = "", values = pal, labels = cellname) + 
    scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0, 
    1)) + expand_limits(x = c(min(pd$x) * min.pand, max(pd$x) * max.pand), 
    y = c(min(pd$y) * min.pand, max(pd$y) * max.pand)) + labs(title = "Bregma 0.11 (mm)", 
    x = NULL, y = NULL) + theme_bw() + theme(legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5, size = rel(titlesize)))

    ggsave(class_dist,file="./cell_class.png")
}

cell_class()


cell_individual_class <- function(pltdat,iclass,xy=TRUE,main=FALSE,titlesize=2,
    pointsize=3,min.pand=0.99,max.pand=1.01,title=NULL){

    pd <- pltdat
    cellname <- as.character(unique(pd$cell_label))

    pd$label_idx <- tmp[as.character(pd$cell_label)]
    pal <- color(length(cellname))

    gpt <- ggplot(pd, aes(x = x, y = y, color = (cell_label == cellname[iclass]))) + 
        geom_point(size = pointsize) + scale_colour_manual(name = "", values = c("gray90", 
        pal[iclass])) + scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0, 
        1)) + expand_limits(x = c(min(pd$x) * min.pand, max(pd$x) * max.pand), 
        y = c(min(pd$y) * min.pand, max(pd$y) * max.pand)) + theme_bw()

    if (main) {
        if (is.null(title)) {
            title = colnames(pd)[igene + 2]
        }
        out = gpt + labs(title = title, x = NULL, y = NULL) + theme(legend.position = "", 
            plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
    } else {
        out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = "")
    }
    return(out)       
}

args <- commandArgs(trailingOnly = TRUE)
cell_location.path <- args[2]
#read cell location table
cell_location <- read.cell_location(cell_location.path)
pd <- cbind.data.frame(x = cell_location$x,y = cell_location$y,cell_label = cell_location$cell_label)
cellname <- as.character(unique(pd$cell_label))

tmp <- LETTERS[1:length(cellname)]
names(tmp) <- cellname
pd$label_idx <- tmp[as.character(pd$cell_label)]

min.pand = 0.99
max.pand = 1.01
pointsize = 1
titlesize = 1

pp <- lapply(1:length(cellname),function(x){
    cell_individual_class(pltdat = pd,iclass = x,main=T,title = cellname[x],pointsize = 1,titlesize = 1.5)
})

library(gridExtra)
p <- grid.arrange(grobs=pp,ncol=3)
ggsave(p,file="./test.png")


# Spatial Distribution of Representative Genes

var_stabilize <- function(x, sv = 1) {
    varx = apply(x, 1, var)
    meanx = apply(x, 1, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = sv)))
    return(log(x + 1/(2 * phi)))
}

relative_func <- function(expres) {
    maxd = max(expres) - min(expres)
    rexpr = (expres - min(expres))/maxd
    return(rexpr)
}

pattern_plot <- function(pltdat, igene, xy = T, main = F, titlesize = 2, 
    pointsize = 3, min.pand = 0.99, max.pand = 1.01, title = NULL) {
    if (!xy) {
        xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 
            1]), split = "x"))), ncol = 2)
        rownames(xy) <- as.character(pltdat[, 1])
        colnames(xy) <- c("x", "y")
        pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
    } else {
        pd <- pltdat
    }
    
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    gpt <- ggplot(pd, aes(x = x, y = y, color = pd[, igene + 2])) + geom_point(size = pointsize) + 
        scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(0, 
        1)) + scale_y_discrete(expand = c(0, 1)) + expand_limits(x = c(min(pd$x) * 
        min.pand, max(pd$x) * max.pand), y = c(min(pd$y) * min.pand, max(pd$y) * 
        max.pand)) + theme_bw()
    if (main) {
        if (is.null(title)) {
            title = colnames(pd)[igene + 2]
        }
        out = gpt + labs(title = title, x = NULL, y = NULL) + theme(legend.position = "none", 
            plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
    } else {
        out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
    }
    return(out)
}

individual_gene <- function(){
        
    args <- commandArgs(trailingOnly=TRUE)
    exptable.path <- args[1]
    cell_location.path <- args[2] 
    gene.config <- args[3]

    exptable <- read.exptable(exptable.path)
    cell_location <- read.cell_location(cell_location.path)
    gene_config <- read.config(gene.config)

    colnames(exptable) <- paste0("Cell","_",1:ncol(exptable))
    rownames(cell_location) <- colnames(exptable)
    
    spark <- CreateSPARKObject(counts = exptable, location = cell_location[, 1:2], 
    percentage = 0.1, min_total_counts = 10)
    
    rm(exptable,cell_location)

    cell_location <- spark@location
    cell_location$total_counts <- apply(spark@counts,2,sum)
    exptable <- spark@counts

    gene_plot <- gene_config$genes
    sig_ct <- exptable[gene_plot,]
    sig_vst_ct <- var_stabilize(sig_ct)
    rel_vst_ct <- apply(sig_vst_ct,1,relative_func)
    pltdat <- cbind.data.frame(cell_location[,1:2],rel_vst_ct) 

    pp <- lapply(1:(ncol(pltdat) - 2), function(x) {
        pattern_plot(pltdat, x, main = T, titlesize = 1.5, pointsize = 1)
	    })
    grid.arrange(grobs = pp, ncol = 2)
    p <- grid.arrange(grobs=pp,ncol=3)
    ggsave(p,file="./gene.png")
}

individual_gene()

