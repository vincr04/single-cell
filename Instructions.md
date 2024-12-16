# Single-cell RNA-seq workshop
January 2025 version

In this workshop, we will learn how to analyse single-cell RNA-seq data using the Seurat software. This simplified exercise is inspired by the [original Seurat tutorial by Rahul Satija et al.](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
>Tasks or questions that you can try and answer, either on your own or by groups of 2-3.  

There is no report to submit, so you can just focus on **understanding** how the analysis works. There is no need to speed through the instructions. Instead, it is important that you make sure you understand the importance of each step, the structure of the different objects involved, and what the different commands do. There are some optional questions - it is probably best if you try and answer them in the end, if you still have time.   
**We are available to help if you have any question or if anything is unclear.**  
***

## 1. Preparation

### Install the software

DESeq2 runs in R and is available as a package via [Bioconductor](https://www.bioconductor.org/), which is a large-scale project to develop, support, and disseminate open source software for bioinformatic data analysis. Many tools used by computational biologists are available there. 

You should have [RStudio](https://posit.co/downloads/) installed on your computer. Google 'DESeq2 Bioconductor' to find out how to install DESeq2. Once DESeq2 is installed, you need R to load the DESeq2 package by typing `library("DESeq2")`. *Remember, code is ALWAYS case-sensitive*. Once this is done, move to the next step.

### Get the data

Bioconductor also contains some 'pre-packaged' datasets that can be easily downloaded and used as examples. Today, we will work on the 'pasilla' dataset ([Brooks et al., Genome Research 2011](https://pubmed.ncbi.nlm.nih.gov/20921232/)) which explores the effects of RNAi knockdown of *pasilla*, a nuclear RNA binding protein implicated in splicing and ortholog of mammalian NOVA1 and NOVA2, on the transcriptome of cultured *Drosophila melanogaster* cells.

The data from this experiment is available as an R package. As above, use Google to find out how to install the `pasilla` package from Bioconductor, and once this is done, don't forget to type `library("pasilla")` to load the package and the data. 

### Additional packages to install

We will use [tidyverse](https://www.tidyverse.org/) to manipulate data. Tidyverse is available via CRAN, which means that you can install it from the 'Packages' tab on the top right of your screen in RStudio. Then, load the package. 

[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) might be an old friend of yours. If it's already installed on your computer you can just load the package. Otherwise, you can also install `ggplot2` from CRAN. Then, load the package.  

[apeglm](https://bioconductor.org/packages/release/bioc/html/apeglm.html) is an additional package needed for LFC shrinkage (more details below). Install the package from Bioconductor, and load it.  

We will install additional packages for gene set enrichment analysis later. 

<br/>

## 2. Generate a *DESeqDataSet*

### Extract the data
A *DESeqDataSet* is an object class in R that stores RNA-seq read counts and other metadata, as well as the results of intermediate calculations, in a single R object. We need to create a new *DESeqDataSet* from the 'pasilla' data. We first create an object containing the **counts data as a matrix** from the 'pasilla' package that we have previously installed and loaded:
```
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE) # generate a path to where the data is located on your computer
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id")) # make a matrix with the counts data
```
And we also create a separate object with the **sample information**:
```
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE) # as above, we generate a path to where the sample annotation is located on your computer
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")] # we only keep the relevant column from this file
```
We can now visualise the **count matrix**, to make sure it actually contains our data.

>The count matrix likely contains several thousands of rows (one per gene). What command would you use to display only the first few rows?  
>Can you tell how many rows and columns this table has? If you don't know how to, use Google  

We can also visualise the **sample information** (it is a small table so we can just display it entirely). 

It is essential that the column headers of the **counts matrix** correspond to the sample names in the **sample information** table, and are arranged in the same order. We can easily see that this is not the case: the sample names in the **sample information** have an extra 'fb' that needs to be removed, and the columns in the **counts matrix** needs to be rearranged. So we need to edit the sample names
```
rownames(coldata) <- sub("fb", "", rownames(coldata))
```
and arrange the order of the columns in the **counts matrix**
```
cts <- cts[, rownames(coldata)]
```

> What do these commands do?  
> Visualise the counts matrix and sample info again. What has changed?

### Build the *DESeqDataSet* 
We can now load the **count matrix** and **sample info** into a new *DESeqDataSet*, which we call `dds`. For this, we use the *DESeqDataSetFromMatrix* function, which is part of the DESeq2 package:
```
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = coldata,
                              design = ~ condition) # the design element is essential to tell the software which samples need to be compared against each other. In our case, we want to use the 'condition' column from the sample info table. 
```

>Here, we created a new object that contains count data and sample info. Type `dds` to display some information about the new object that you created. What kind of information do you get?  
>DESeq2 allows easy access to the key elements of a *DESeqDataSet*. For example, you can type `counts(dds)` to quickly view the **count matrix**.

<br/>

## 3. Pre-filtering and changing reference levels
The genes with the lowest count numbers are likely to be almost absent from the cells that we have sequenced. Removing these genes is not crucial, however having less genes might speed up the analysis. Here, we only do minimal filtering by removing genes for which we have less than 10 counts across all samples.
```
keep <- rowSums(counts(dds)) >= 10 #gives a TRUE or FALSE value for each gene (row): TRUE if there are more than 10 reads, FALSE otherwise.
dds <- dds[keep,] #keep only genes with a TRUE value
```

>How many genes were removed?  

**Important** By default, conditions are considered in alphabetical order, with the first one being assigned as the 'reference'. In our case, conditions are labelled 'untreated' and 'treated', so by default the software will assign the 'treated' group as the control, which would be wrong. Conditions are saved as *factors* in R and if you type `dds$condition`, you will notice that indeed 'treated' comes first in the list of 'levels'. Use the following command to make sure that 'untreated' is the first level
```
dds$condition <- relevel(dds$condition, ref = "untreated") # the 'relevel' function can be used to assign a new 'reference' level 
```
>Type `dds$condition` again. What has happened?

<br/>

## 4. Running the differential expression analysis

The actual command to run DESeq2 is pretty simple:
```
dds <- DESeq(dds) # more options can be added, but for this exercise we will only use default parameters
```
This can run for up to a minute, depending on your computer. The results are saved in the `dds` object and you can visualise a preview of the results by typing 
```
res <- results(dds)
res
```

>Take a moment to look at the table in front of you. What information is in the header? Can you make sense of what each column represents?  
>If you are struggling, the command `mcols(res)$description` conveniently provides more information about each column.  
>Can you find the gene with the highest/lowest *log2FoldChange*? And the gene with the lowest adjusted p-value? Hint: you can sort the table according to the values of selected columns.  
>Why are some *log2FoldChange* values negative?  
>OPTIONAL - to avoid false positives (type I errors), the p-value is corrected with the Benjamini-Hochberg (BH) method. Use Google to find more about the BH method.  

We can add options to the `results` command. For example, `contrast` can be used to use a different reference sample. 

>Try the two following commands. What effect does this have on the results?  
>`results(dds, contrast=c("condition","treated","untreated"))`  
>`results(dds, contrast=c("condition","untreated","treated"))`

`summary(res, alpha=0.05)` provides some interesting statistics about the data

### Size factors

You probably remember from the transcriptomics lecture that DESeq2 calculates 'size factors' to normalise the data. Use Google to find out how to display the size factors for each sample. 

You can access the raw counts and normalised counts by typing `counts(dds)` and `counts(dds, normalized=TRUE)`, respectively. 

> Can you verify that the normalised counts are indeed the raw counts divided by the size factors? 
> OPTIONAL: write your own R script to calculate size factors

<br/>

## 5. Visualising the data

A large table can be very difficult to read for the human eye. Therefore, we have to use graphics to visualise the data. Below are just a few examples. 

### MA-plot
DESeq2 offers simple commands to visualise the data. For example, the `plotMA` function displays *log2FoldChange* over the mean of normalised counts, while highlighting genes with a significant adjusted p-value.
```
plotMA(res, ylim=c(-4,4), alpha=0.05)
```
>What does this plot tell you? What happens if you use `alpha=0.01`?  
>OPTIONAL - Can you recreate this plot using ggplots? 

### Volcano plot
A volcano plot is a scatter plot representing the negative log of the *padj* over the *log2FoldChange*. As for the MA-plot, each gene is represented by a point. DESeq2 does not include a function to make volcano plots, but we can easily make one using `ggplot`. 
```
#We start with a very basic representation
p <- ggplot(as.data.frame(res), aes(x=log2FoldChange, y=-log(padj))) # 'res' needs to be transformed into a data.frame. Then, we tell ggplot what information from this table needs to be on each axis. 
p + geom_point() # this command represents the data as points.

# And make it look a bit nicer
p + geom_point(color=ifelse(res$padj<0.05 & abs(res$log2FoldChange)>0.5, "red", ifelse(abs(res$log2FoldChange)>0.5, "steelblue1", ifelse(res$padj<0.05, "steelblue4", "grey")))) + # adding some colour: if padj<0.05 and |log2FoldChange|>0.5, the points are red, if only one of the conditions is true, they are blue, otherwise they are grey.
  lims(y=c(0,40), x=c(-2.5,2.5)) # cropping X and Y axes

# OPTIONAL: we can also add gene labels, using the `ggrepel` package
install.packages("ggrepel")
library("ggrepel")
p + geom_point(color=ifelse(res$padj<0.05 & abs(res$log2FoldChange)>0.5, "red", ifelse(abs(res$log2FoldChange)>0.5, "steelblue1", ifelse(res$padj<0.05, "steelblue4", "grey")))) + # add the colour
    lims(y=c(0,40), x=c(-2.5,2.5)) + # cropping X and Y axes
    geom_text_repel(aes(x=log2FoldChange, y=-log(padj), label=ifelse(rownames(res)=="FBgn0038198", rownames(res), ""))) # add text: if the gene is called "FBgn0038198", then print its name, otherwise print nothing.
```

### P-value histogram
Now let's plot a histogram of the adjusted p-values (in base R, but this can also be done with ggplot).
```
hist(res$padj,breaks = 100); abline(v=0.05,col="red")
```
>What does this plot tell us?  

### Plot normalised counts for specific genes of interest
Using the following command, you can plot the normalised counts between conditions for any gene you like (just replace *XXX* with the name of your gene of choice). 
```
plotCounts(dds, gene=XXX, intgroup="condition")
```
To plot the gene with the lowest *padj*, you can replace *XXX* with `which.min(res$padj)`. 

>What command would you use to plot the gene with highest *log2FoldChange*?  
>Can you confirm that the knock down of the *pasilla* gene was indeed successful? Note: you need to search [FlyBase](flybase.org) to find the FlyBase ID (FBgn) for the *pasilla* gene.  
>OPTIONAL - how would you plot this with ggplot?

### Heatmap
Often, a heatmap is an ideal way of representing variation in expression across several genes in a single figure panel. Here is an example (using `ggplot`) for the top 10 genes with the highest *log2FoldChange* and significant adjusted p-value. 
```
res.sig <- res[!is.na(res$padj) & res$padj<0.05,] # Filtering only significant results
top10genes <- rownames(res.sig[rev(order(res.sig$log2FoldChange)),][1:10,]) # Selecting the genes with the highest log2FoldChange
top10counts <- counts(dds,normalized=TRUE)[top10genes,] # extract normalised counts for each gene
top10counts <- as.data.frame(t(scale(t(top10counts)))) # calculate z-scores **note: to use the 'scale' function we need to transpose the table (invert rows and columns), and transpose it back again. 
top10counts <- cbind(top10counts, Gene=rownames(top10counts)) # add extra column with gene name
top10counts <- gather(top10counts, Treatment, Z_score, -Gene) # reshape data to use with ggplot

ggplot(top10counts, aes(x=Treatment, y=Gene, fill=Z_score)) + # create a plot and provide information about axes. 'fill' indicates what information will be used to fill each tile of the heatmap
    geom_tile() + # add tiles
    scale_fill_gradient2(low="navy", mid="linen", high="darkred", na.value="transparent") # add colour scheme
```
> Can you recreate this heatmap for the 10 most downregulated genes? 
> OPTIONAL - can you understand every part of the code above? 

### PCA
Principal component analysis can be used to establish how different samples are from each other. Conveniently, the `plotPCA` function is included in DESeq2, but first, we need to transform the raw count data using *variance stabilising transformations* (*VST*), which produces normalised, log2 scale values. 
```
vsd <- vst(dds, blind=FALSE) # perform the variance stabilising transformations. We use blind=FALSE to take the 'within group' variabliity into account.
head(assay(vsd), 3) # print the calculated values for the first three genes
plotPCA(vsd, intgroup=c("condition", "type"))
```
> What can you conclude from this plot?

<br/>

## 6. Log fold change shrinkage
Genes with low counts are more likely to have high *log2FoldChange* values because the natural variation between samples may create artificial differences. For example, for the data in the table below, the *mean* count value for the treated samples is more than 4x lower than for the untreated sample for **Gene A**, but almost equal for **Gene B**, even though the *difference* in read counts between both samples is the same. This leads to a bias towards low-count genes having a high *log2FoldChange*, with **Gene A** having a much higher *log2FoldChange* than **Gene B**. 

|Gene  |untreated 1|untreated 2|untreated 3|untreated 4|treated 1|treated 2|treated 3|
|------|-----------|-----------|-----------|-----------|---------|---------|---------|
|Gene A|10         |8          |4          |5          |1        |3        |1        |
|Gene B|1010       |1008       |1004       |1005       |1001     |1003     |1001     |

To correct for this, we apply a DESeq2 function called `lfcShrink`. This function requires the user to enter a `coef`, which is an argument specifying the comparison to extract from the data. The `resultsNames(dds)` command displays available `coef` options. We will select `"condition_treated_vs_untreated"`. 

We can now run the `lfcShrink` function. The 'type' argument specifies the LFC shrinkage method. For our analysis, the more recent `apeglm` method is recommended.
```
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
plotMA(resLFC, ylim=c(-4,4), alpha=0.05)
```
>What has changed between this plot and the previous one? 
>Repeat some of the plots above using the LFC-shrunk data. Can you see any differences? 

<br/>

## 7. Gene-set enrichment analysis (GSEA)

We can now tell which genes are differently regulated in *pasilla* knock-down cells, however this in itself does not provide much biological information. What would be interesting at this stage is to switch the focus of our analysis from individual genes onto biological pathways. For this, we can use gene-set enrichment analysis (GSEA), also called pathway enrichment analysis.  

GSEA uses existing databases that contain genome-wide information about the characteristics of each gene, and therefore categorises genes into specific groups of similar function. The best known example of such database is the Gene Ontology (GO), which categorises (almost) all known genes according to their **cellular component, molecular function, or biological process**.  

Here, we will use the `clusterProfiler` package to perform GSEA, the `org.Dm.eg.db` database of genomic information, which contains GO information for all *Drosophila* genes, and the `enrichplot` package, which will allow us to plot gene networks.

First, install (from Bioconductor) and load these three packages.  

### Prepare the data
GO enrichment analysis requires a list of genes considered as significantly differently expressed. Here, we will select the genes with *log2FoldChange*>0.5 (up-regulated) and *padj*<0.05 from our DESeq2 analysis.  
```
res_filter <- filter(as.data.frame(res), log2FoldChange>0.6 & padj<0.05) #filtering the data into a new object
gene <- na.omit(rownames(res_filter)) # we extract the list of differently expressed genes. 
```
> How many genes are in this list?
> How would you prepare a similar list of downregulated genes?

### Run the enrichment analysis

We use the `enrichGO` function to search for enriched categories. Using the `ont=BP` option, we specifically search trough **biological process** categories.  
```
ego <- clusterProfiler::enrichGO(gene          = gene, # the list of genes we generated above
                                 OrgDb         = org.Dm.eg.db, # the database with GO information
                                 ont           = "BP", # we are searching the 'biological process' category
                                 keyType       = "ENSEMBL") # our genes are in ENSEMBL (= FlyBase ID) format
head(ego@result) # print the first few lines of the results
```
> What information is contained in each column of this table?
> Repeat the analysis for the cellular compartment (CC) and molecular function (MF) categories. What can you see?

### Plot the results

We can visualise these results using a variety of plots. Try the two options below:
```
barplot(ego, showCategory=25)
```
```
lfc_vector <- res$log2FoldChange # creating a vector with log2FoldChange values from the DESeq2 analysis
names(lfc_vector) <- rownames(res) # naming elements of this vector with the names of the corresponding genes
cnetplot(ego, showCategory=10, foldChange=lfc_vector)
```
> What information is displayed in these two plots? 
> Can you produce similar plots for the cellular compartment (CC) and molecular function (MF) categories?

<br/>

## 8. OPTIONAL - Differential isoform expression

The *pasilla* gene encodes a nuclear RNA binding protein implicated in mRNA splicing. Therefore, we expect that knocking down this gene will have an effect on splicing, and that differential splicing events will occur in the treated sample. To investigate this, we can use the DEXSeq package.  
Download DEXSeq from Bioconductor, and follow the steps in the [DEXSeq vignette](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html) (you can start from step 3.1.). 
> How many genes can you find with significantly different exon usage?  
> Read the [original *pasilla* paper](https://genome.cshlp.org/content/21/2/193.long). Does your analysis help you understand the findings from this paper?

<br/>
<br/>

<h1 align="center">The end</h1>
<p align="center">(you are now an expert!)</p>
