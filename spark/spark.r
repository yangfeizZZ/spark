#!/usr/bin/env/ Rscript

#Author : yangfei
#data : 2022.6.24




echo <- function(...){
    cat(...,"\n")
}


print.usage <- function(){
    echo("Perform the SPARK analyze.")
    echo()
    echo("Usage:")
    echo("   $ ./spark_analysis.R  Experssion_table cell_location")
    echo()
}

read.exptable <- function(path){
    exptable <- read.table(path)
    return(exptable)
}


read.cell_location <- function(path){
    cell_location <- read.table(path)
    #if()
    return(cell_location)
}


main <- function(){
    
    args <- commandArgs(trailingOnly=TRUE)
    print.usage()
    exptable.path <- args[1]
    cell_location.path <- args[2]
    
    #read spatial expression table
    exptable <- read.exptable(exptable.path)
    echo("Shape of full spatial expression table:")
    echo(dim(exptable))
    echo()

    #read cell loaction table
    cell_location <- read.cell_location(cell_location.path)

    if(length(colnames(exptable)) != length(rownames(cell_location))){
        stop("The file is not valid,the colname of exptable and the rowname of cell_location must have same size")
    }   

    colnames(exptable) <- paste0("Cell","_",1:ncol(exptable))
    rownames(cell_location) <- colnames(exptable)


    #run spark to analysis
    library(SPARK)
    SPARK <- CreateSPARKObject(counts = exptable,loaction = cell_location[,1:2],percentage = 0,min_total_counts = 10)
    spark@lib_size <- apply(spark@counts,2,sum)
    spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
        num_core = 10, verbose = FALSE, fit.maxiter = 500)
    spark <- spark.test(spark, check_positive = T, verbose = FALSE)
    SVG <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
    SVG <- SVG[which(SVG$adjusted_pvalue <= 0.05),]
    write.table("./SVG.txt")
    print(SVG[1:4,])

}

main()


