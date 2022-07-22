# Compute SRT svg

The aim of this section is to compute the SVG of SRT which is use a R package of [spark](https://www.nature.com/articles/s41592-019-0701-7)  



## Work Flow
![](https://github.com/yangfeizZZ/spark/blob/master/image/pipeline.png)

## Requirements
This script was run on R,and following R package are required:
- spark(the latest)

## Use the script 
This script can compute SVG of spatial transcriptome.If you want to generate the SVG of your own spatoal transcriptome data,you can use this script.

If you want tu use this script,you should provide two txt of your own spatial transcriptome data.

1. The spatial transcriptome gene expression matrix. see [expression-matrix](https://github.com/yangfeizZZ/spark/blob/master/example/count.txt)
2. The cell location and celltype matrix. see  [location-matrix](https://github.com/yangfeizZZ/spark/blob/master/example/info.txt)

If you prepare the tow matrix,you can compute the SVG list by follow command

```R
$ Rscript spark.r count.txt info.txt
```

The result of comand is [SVG.txt](https://github.com/yangfeizZZ/spark/blob/master/example/SVG.txt).THis txt contain two colnume: "combined_pvalue" and "adjusted_pvalue" . The "adjusted_pvalue" is we need.

# Plot

The aim of this section is to plot the cell

## Requirement
This script was run on the R,and following R package are required:
- spare
- ggplot2
- gridExtra
- reshape2
- grid
- RColorBrewer

## Use the script
This script can plot the cell by celltype and plot gene of SVG

If you want tu use this script,you should provide two txt of your own spatial transcriptome data.

1. The spatial transcriptome gene expression matrix. see [expression-matrix](https://github.com/yangfeizZZ/spark/blob/master/example/count.txt)
2. The cell location and celltype matrix. see  [location-matrix](https://github.com/yangfeizZZ/spark/blob/master/example/info.txt)
3. The SVG json file.see [gene_config.json](https://github.com/yangfeizZZ/spark/blob/master/example/gene_config.json)

If you prepare the tow matrix,you can compute the SVG list by follow command

```R
$ Rscript spark.r count.txt info.txt gene_config.json
```

