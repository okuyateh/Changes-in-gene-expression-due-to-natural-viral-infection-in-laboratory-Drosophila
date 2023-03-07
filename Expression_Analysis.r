
#Title: Gene Expression Analysis

#Author: Oumie Kuyateh

#Date written: 2023

#Aim: The aim of this analysis is to do Differential Expression (DE) analyses using dream. dream is built on limma/voom but allows for more than one random effect to be put in the model.
Differential expression for repeated measures (dream) uses a linear model model to increase power and decrease false positives for RNA-seq datasets with multiple measurements per individuals. Dream achieves this by combining multiple statistical concepts into a single statistical model. The model includes:

#- flexible modeling of repeated measures gene expression data
#- precision weights to model measurement error in RNA-seq counts
#- ability to model multiple random effects
#- random effects estimated separately for each gene
#- hypothesis testing for fixed effects in linear mixed models
#- small sample size hypothesis test


```{r}
#Firstly, Remove all files from R's environment

rm(list = ls())
```


```{r}
#Now load all the libraries I will need
library(edgeR)#edgeR package loads limma as a dependency)
library(ggplot2)
library(dplyr)
library(tibble) 
library(variancePartition)
library(BiocParallel)
library(EnhancedVolcano) #To plot the Volcano plot
library(ggfortify) #To convert the FBgn numbers to gene symbols
```

```{r}
#Now it's time to read in the Count_Matrix. We need to  transform the count matrix since we want the gene names as rownames and sample names as column names
Count_Matrix <- data.frame (t (read.table( "Processed_Count_Matrix_PRJNA258012_Yanzhu_11_11_2020.txt", sep = "\t", row.names = 1, header= TRUE)))


```



```{r}

#Our Count_Matrix data has FlyBase identifier (FBgn#) as the rownames. To get the GeneSymbol, we use data from the the "org.Dm.eg.db" AnnotationDb object which was loaded in the beginning of our analysis. To do this, we use the "mapIds()" function, to look in the rownames of the Count_Matrix  which contain the "FLYBASE" keytype. For each"FLYBASE"keytype return the Gene SYmbol. If multiple Symbols are available, give me the first one. Save this in a data.frame called "FBgn_to_Symbol".


FBgn_to_Symbol <- data.frame( mapIds(org.Dm.eg.db, row.names (Count_Matrix),  column="SYMBOL" , keytype="FLYBASE", multiVals="filter"))

#First order the FBgn_to_Symbol and the Count_Matrix by rownames 


FBgn_to_Symbol<-  FBgn_to_Symbol[ order(row.names( FBgn_to_Symbol)),, drop =FALSE]

Count_Matrix<-  Count_Matrix[ order(row.names( Count_Matrix)),, drop =FALSE]

#Check that they are in the same order

print (paste (" The row.names of FBgn_to_Symbol and Count_Matrix are in the same order:", identical(row.names(Count_Matrix), row.names (FBgn_to_Symbol) )))

#Now attach a column in Count_Matrix containing the gene_Symbols

Count_Matrix$Gene_Symbols <- FBgn_to_Symbol[,1]

#Some FBgn# do not have corresponding gene Symbols and will be put as NA. However, we cannot have NA as ronames so remove them  

Count_Matrix<- na.omit(Count_Matrix)

#Now make the gene names as row.names 

#First remove the current row.names which is the FBgn numbers

row.names(Count_Matrix) <- NULL

#Now make the gene symbols the rownames

Count_Matrix <- column_to_rownames(Count_Matrix, var = "Gene_Symbols")


#As we can see after removing the FBgn without corresponding gene symbols, the number of genes in the count matrix decreased from 16613 t0 16595
```


```{r}

Meta_data <- read.table( "Processed_Metadata_PRJNA258012_Yanzhu_11_11_2020.txt", sep = "\t", row.names = 1, header = TRUE)

#Now remove samples from the Meta_data where the Author_Sex is different from the Feature_Count_Sex and call this Meta_data. The Author_Sex is the sex that was given to the samples by the authors whilst the Feature_Count_Sex was the sex given based on total number of male and female features obtained from featureCount. 

Meta_data <- Meta_data [  Meta_data$Author_Sex == Meta_data$featureCounts_Sex,] 
  
#Select the columns we require from the Meta_data data frame for our analysis which are:"DGRP_Line", "Environment" ,  "Author_Sex" ,  "Galbut_150_Presence", "Motts_Mill_150_Presence". The 150 in the "Galbut_150_Presence" & "Motts_Mill_150_Presence" represent the treshold that was used to estimate barcode switching.

Meta_data<- dplyr::select(Meta_data, c( "DGRP_Line", "Environment", "Author_Sex" ,  "Motts_Mill_Presence", "Galbut_Presence", "Galbut_Presence"))

#Now look at the number of samples with Galbut, MottsMill and both in the Meta_data data frame. #The Meta_data  will be used in our dream analysis. 

print(paste("Out of", nrow (subset(Meta_data)), "samples", nrow (subset(Meta_data, Galbut_Presence == "Yes" &Authr)), "are infected with Galbut only"))

print(paste("Out of", nrow (subset(Meta_data)), "samples", nrow (subset(Meta_data, Motts_Mill_Presence == "Yes")), "are infected with Motts_Mill only"))

print(paste("Out of",nrow (subset(Meta_data)), "samples", nrow (subset(Meta_data, Galbut_Presence == "Yes" & Motts_Mill_Presence == "Yes")), "are infected with both Galbut & Motts_Mill")) 

```

```{r}

Meta_data <- read.table( "Processed_Metadata_PRJNA258012_Yanzhu_11_11_2020.txt", sep = "\t", row.names = 1, header = TRUE)

```



```{r}
#Now keep only samples in the Count_Matrix that are in the Meta_data. The t() function is used since the subset works on rows and the Count_Matrix and Meta_data are in different orientation.

Count_Matrix <- subset(Count_Matrix, row.names(Count_Matrix) %in% row.names(Meta_data) )

#The countmatrix has the gene IDs as columnanmes and sample names as rownames. However, for dream analysis we want the gene IDs as rowanems and sample names as column names. To achieve this we need to use the transpose function t().

Count_Matrix <- data.frame(t(Count_Matrix))

#Now the columnames of the Count_Matrix  is the same as the rownames of the Meta_data. However, we need to test that they are identical as a QC step. To do this, we first order the Meta_data and the Count_Matrix by rownames and columnnames respectively.


Meta_data <- Meta_data [order(row.names(Meta_data)),]


Count_Matrix <- Count_Matrix[ , order(colnames(Count_Matrix))]


#Now check that the rownames of the Meta_data and colnames of Count_Matrix are in the same order


print (paste ("row.names of Meta_data and colnames of Count_Matrix are identical:", identical(row.names( Meta_data), colnames (Count_Matrix) )))
```


```{r}

# Filter genes by number of counts in the Count_Matrix data frame  by:

#Summing the counts per million (cpm) that are above 0.1 for each gene and then select the genes whose sum is greater than or equal to 5 cpm. Save these genes in isexpr


isexpr <- rowSums(cpm(Count_Matrix) > 0.1) >= 5


#Now create a DGEList object from the Count_Matrix of all the genes in isexpr that are in Count_Matrix and save this as "geneExpr. 


geneExpr <- DGEList( Count_Matrix[isexpr,] )

#Now calculate normalization facotors from the "geneExpr" DGEList object using calcNormFactors. NB calcNormFactors doesn't normalize the data, it just calculates normalization factors for use downstream.


geneExpr <-calcNormFactors( geneExpr )

```

#The dream method replaces two core functions of limma with a linear mixed model.

#1. voomWithDreamWeights() replaces voom() to estimate precision weights
#2. dream() replaces  lmFit() to estimate regression coefficients
#Otherwise dream uses the same workflow as limma with topTable(), since any statistical differences are handled behind the scenes.
#In limma/voom, the voom() function prepares our RNAseq data for a format suitable for linear modelling using limma. Voom is an acronym for "mean-variance modelling at the observational level". It is meant to calculate the mean-variance relationship in the data and use this to assign weights to each observation based on this. Count-data have a significant (non-trivial) mean-variance relationship. Raw counts obtained from featureCounts show increasing variance with increasing count-size increases. On the other hand, log-counts typically show a decreasing mean-variance trend i.e. the higher the logCount, the lower the variance. The voom finction estimates the mean-variance trend for log-counts, then assigns a weight to each observation based on its predicted variance. The weights are then used in the linear modeglling process to adjust for heteroscedasticity. Heteroscedasticity (also spelled heteroskedasticity) refers to the circumstance in which the variability of a variable is unequal across the range of values of a second variable that predicts it.

#The voom function works in the following ways:
#1.Transforms the counts  to log2 counts per million reads (CPM), where "per million reads" is defined based on the normalization factors we calculated earlier
#2. A linear model is then  fitted to the log2 CPM for each gene, and the residuals (erros) calculated
#3.A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot below)
#4. The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

#As mentioned earlier, voomWithDreamWeights  is the same as what described above for voom(), except that it allows random effects in the formula.

```{r}
# Specify parallel processing parameters. This is used implicitly by dream() to run in parallel. It allows three types of clusters to be used. Possible values are SOCK (default) and MPI. Instead of type=FORK use MulticoreParam.type=c("SOCK", "MPI", "FORK"),

param <- SnowParam(7, "SOCK", progressbar=TRUE)

register(param)

# Specify the model to be tested and save this as form. NB the DGRP line is a random effect

form <- ~ Author_Sex  +  Motts_Mill_Presence + Galbut_Presence +  Motts_Mill_Presence:Author_Sex  + Galbut_Presence:Author_Sex + (1|DGRP_Line) + (1|Environment)


# estimate weights using linear mixed model of dream

vobjDream <- voomWithDreamWeights( geneExpr, form, Meta_data )

```


```{r}

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test


fitmm <- dream( vobjDream, form, Meta_data)


```


```{r}
#Now look at the coefficients of the model

head(fitmm$coefficients)

```


```{r}
#In our model Females and No Infection are the baselines
#Extract results for each coefficient from fitmm using the topTable() function and save them in data frames. Put all the for each coefficient in a data frame, n=INF means everything

#look at the sex the Sex Differences in Gene Expression in The Absence of Virus Infection and save this in a data frame called "Res_Sex". n=INF means infinity i.e. results for all the genes

Res_Sex <- topTable( fitmm , coef="Author_SexM", n = Inf )

#Now collect data for the individual sexes. First, look at the changes in gene expression in response to Galbut Infection In Females and save this as "Res_Galbut_Females"

Res_Galbut_Females <- topTable( fitmm , coef="Galbut_PresenceYes", n = Inf )

#Res_Galbut_and_Chaq_Females <- topTable( fitmm , coef="Galbut_and_Chaq_PresenceYes", n = Inf )


#Now I want to look at the changes in gene expression in response to Motts Mill Infection In Females and save this in "Res_Motts_Mill_Females"

Res_Motts_Mill_Females <- topTable( fitmm , coef="Motts_Mill_PresenceYes", n = Inf )


#Now I want to look at Sex Differences in Response to Gene Expression for Galbut. 

Res_Sex_Differences_Galbut <-topTable( fitmm , coef="Author_SexM:Galbut_PresenceYes", n = Inf )

#Res_Sex_Differences_Galbut_and_Chaq <-topTable( fitmm , coef="Author_SexM:Galbut_and_Chaq_PresenceYes", n = Inf )


#Now I want to look at Sex Differences in Response to Gene Expression for Motts Mill. 

Res_Sex_Differences_Motts_Mill <-topTable( fitmm , coef="Author_SexM:Motts_Mill_PresenceYes", n = Inf )

```


```{r}
#Now look at the plot for the change in gene expression in response to Motts_mill infectionin females
 EnhancedVolcano(Res_Galbut,
    lab = rownames(Res_Galbut),
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Female Gene Expression in Response to Motts Mill Infection',
   subtitle = "(PRJNA258012_Yanzhu)",
    titleLabSize = 10,
    subtitleLabSize = 10,
  captionLabSize = 10,
    pCutoff = 0.05,
     FCcutoff = 0.4,
    axisLabSize = 10,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1, 
   cutoffLineType = 'blank',
     xlim=c(-4, 4),
#ylim=c(0, 500),
pointSize = 1.0,
    labSize = 3.0)

```



#Now that the dream analysis is done. We need to plot volcano plots to look at our results. 

```{r}
#Firstly, we will look at the plot for the overall Sex differences in immune response in the absence of Virus Infection. Here, we are going to highlight the sex-specific genes that were used in my barcode switching estimation. We need to know the gene symbols that corresponds to theose FBgn in FBgn_to_Symbol data.frame.

subset (FBgn_to_Symbol, grepl("FBgn0004181|FBgn0000356|FBgn0000358|FBgn0000359|FBgn0000360|FBgn0011669|FBgn0011694|FBgn0046294|FBgn0053340|FBgn0250832|FBgn0259795|FBgn0259971|FBgn0259975|FBgn0262099|FBgn0262623|FBgn0270925", row.names(FBgn_to_Symbol)))
```


```{r}
#Now that we know them we can draw the volcano plot and highlight these genes.  For this, I am making  xlim=c(-11, 11) so as to capture them in the graph.
EnhancedVolcano(Res_Sex,
    lab = rownames(Res_Sex),
    selectLab = c('Cp16','Cp19', 'Cp36', 'Cp38', 'Ebp', 'Mst57Db', '	EbpII', 'CG12699', 'CG33340', 'Dup99B', 'loopin-1', 'CG42481', 'Sfp87B', 'CG42852', 'CG43147', 'CG4836'),
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Sex Differences in Gene Expression in the Absence of Virus Infection',
    subtitle = "(PRJNA258012_Yanzhu)",
    titleLabSize = 10,
    subtitleLabSize = 10,
  captionLabSize = 10,
    pCutoff = 0.05,
     FCcutoff = 2,
    axisLabSize = 10,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1, 
   cutoffLineType = 'blank',
drawConnectors = TRUE,
     xlim=c(-15, 15),
#ylim=c(0, 500),
pointSize = 1.0,
    labSize = 2.0,
    colConnectors = 'black', widthConnectors = 0.5, boxedLabels = TRUE)
```


```{r}
#Now look at the plot for the change in gene expression in response to Galbut  only  infection in females
 EnhancedVolcano(Res_Galbut_Females,
    lab = rownames(Res_Galbut_Females),
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Female Gene Expression in Response to Galbut only Infection',
    subtitle = "(PRJNA258012_Yanzhu)",
    titleLabSize = 10,
    subtitleLabSize = 10,
  captionLabSize = 10,
    pCutoff = 0.05,
     FCcutoff = 0.4,
    axisLabSize = 10,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1, 
   cutoffLineType = 'blank',
     xlim=c(-6, 6),
#ylim=c(0, 500),
pointSize = 1.0,
    labSize = 3.0)

```


```{r}
#Now look at the plot for the change in gene expression in response to Motts_mill infectionin females
 EnhancedVolcano(Res_Motts_Mill_Females,
    lab = rownames(Res_Motts_Mill_Females),
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Female Gene Expression in Response to Motts Mill Infection',
   subtitle = "(PRJNA258012_Yanzhu)",
    titleLabSize = 10,
    subtitleLabSize = 10,
  captionLabSize = 10,
    pCutoff = 0.05,
     FCcutoff = 0.4,
    axisLabSize = 10,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1, 
   cutoffLineType = 'blank',
     xlim=c(-4, 4),
#ylim=c(0, 500),
pointSize = 1.0,
    labSize = 3.0)


```




```{r}
#Now look at the plot for the sex differences in gene expression in response to Galbut only infection 
EnhancedVolcano(Res_Sex_Differences_Galbut ,
    lab = rownames(Res_Sex_Differences_Galbut ),
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Sex Differences in Response to Galbut  Infection',
    subtitle = "(PRJNA258012_Yanzhu)",
    titleLabSize = 10,
    subtitleLabSize = 10,
  captionLabSize = 10,
    pCutoff = 0.05,
     FCcutoff = 0.4,
    axisLabSize = 10,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1, 
   cutoffLineType = 'blank',
     xlim=c(-3, 3),
#ylim=c(0, 500),
pointSize = 1.0,
    labSize = 3.0)


```




```{r}
#Now look at the plot for the sex differences in gene expression in response to Motts Mill infection 
 EnhancedVolcano(Res_Sex_Differences_Motts_Mill,
    lab = rownames(Res_Sex_Differences_Motts_Mill),
    x = 'logFC',
    y = 'adj.P.Val',
    title = 'Sex Differences in Response to Motts Mill Infection',
     subtitle = "(PRJNA258012_Yanzhu)",
    titleLabSize = 10,
    subtitleLabSize = 10,
  captionLabSize = 10,
    pCutoff = 0.05,
     FCcutoff = 0.4,
    axisLabSize = 10,
    legendPosition = 'right',
    legendLabSize = 10,
    legendIconSize = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1, 
   cutoffLineType = 'blank',
     xlim=c(-6, 6),
pointSize = 1.0,
    labSize = 3.0)
```





























