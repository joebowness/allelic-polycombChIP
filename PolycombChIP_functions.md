## Useful R functions for allelic Polycomb modification ChIP-seq data analysis
**CentreColandNames**

This function is designed to run directly on a nested read_delim() command.
The read_delim() loads a .txt file produced from runnning the custom Python script *ExtractInfoFrombedGraph_AtBed.py* (https://github.com/guifengwei) on a bedgraph of the ChIP library with a .bed file of windows of a given length (eg. 250kb) of Chromosome X.

The function names columns and creates an additional column of the centre of the window (for plotting).

```{CentreColandNames}
CentreColandNames <- function(table){
  table <- cbind(table[1:3],0.5*(table[2] + table[3]), table[4])
  colnames(table)<- c("Chr","Start","End","Centre","Value")
  return(table)
}

```

Here is an example of this command being used:

```
 MyFile <- CentreColandNames(read_delim("~/PATH/TO/MyFile.sort.Norm.chrX.bedGraph.250KB.txt",
                            "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_skip(), X6 = col_skip()), trim_ws = TRUE))
```

**IPoverInput**

This normalises IP files to an appropriate input to calculate enrichment for each window over the chromosome.

```{IPoverInput}
IPoverInput <- function(IP, input){
  Enrichment <- as.numeric(IP[,5] / input[,5])
  Enrichment[is.nan(Enrichment)] <- NA
  IPoverinput_table <- cbind(IP[1:4], Enrichment)
  return(IPoverinput_table)
}
```

**DoxMinusNoDox**

For non-allelic analysis, this function subtracts a NoDox sample from a Dox sample. This enables plots of the distribution of Xist-specific gain of Polycomb enrichment over the chromosome.

```{DoxMinusNoDox}
DoxMinusNoDox <- function(Dox,NoDox){
  DoxMinusNoDox <- cbind(Dox[,1:4], Dox[,5] - NoDox[,5])
  colnames(DoxMinusNoDox)[5] <- "DoxMinusNoDox"
  return(DoxMinusNoDox)
}

```

**XiMinusXa**

For allelic analysis, this function subtracts the distribution pattern of enrichment over the active X (Xa) from the inactive X (Xi) chromosome. Thus this enables plotting of Xi-specific Polycomb enrichment 'internally normalized' within one sample, so is more robust to technical variability (eg. in ChIP efficiency) between samples than the non-allelic DoxMinusNoDox approach (see above). 

In iXist-ChrX-Dom lines, genome1 = *Castaneous* = **Xa** and genome2 = *Domesticus*/129 = **Xi**. 

```{XiMinusXa}
XiMinusXa <- function(genome2,genome1){
  XiMinusXa <- cbind(genome2[,1:4], genome2[,5] - genome1[,5])
  colnames(XiMinusXa)[5] <- "XiMinusXa"
  return(XiMinusXa)
}

```

**poormappability.blacklist**

This function defines windows with outlier signal in an non-allelic 'input' sample, which are often the result of poor mappability in repetitive regions of the genome. We suggest to define ‘poor mappability’ regions as windows with +/- 2.5 median absolute deviation (from visual inspection of plots), although this mad threshold is adjustable.


```{poormappability.blacklist}
poormappability.blacklist <- function(input,madfactor){
  values <- unlist(input$Value) 
  values[values == 0] <- NA
  values.mad <- mad(values, na.rm=TRUE)
  values.median <- median(values, na.rm=TRUE)
  upper <<- values.median + madfactor*values.mad
  lower <<- values.median - madfactor*values.mad
  blacklist <- input[which(
    input$Value > upper | input$Value < lower),1:3]
  return(blacklist)
}
```

**lowallelic.blacklist**

This function defines windows for which there are not sufficient allelic reads (eg. poorly mappable, repetitive, or very few strain-specific SNPs) to confidently assess allele-specific enrichment. It takes as allelic input files (in iXist-ChrX-Dom lines, g1 = *Castaneous* and g2 = *Domesticus*/129), and marks windows ranking in the lowest X% of signal, where the threshold X is adjustable. We recommend 5% (from visual inspection of plots). 

```{lowallelic.blacklist}
lowallelic.blacklist <- function(input_g1,input_g2,threshold){
  #g1
  values <- unlist(input_g1$Value) 
  values[values == 0] <- NA
  lower.g1 <<- quantile(values,threshold, na.rm=TRUE)
  #g2
  values <- unlist(input_g2$Value) 
  values[values == 0] <- NA
  lower.g2 <<- quantile(values,0.1, na.rm=TRUE)
  
  blacklist_g1 <- input_g1[which(input_g1$Value < lower.g1),]
  blacklist_g2 <- input_g2[which(input_g2$Value < lower.g2),]
  blacklist <- rbind(blacklist_g1,
                             blacklist_g2)
  blacklist <- unique(blacklist[,1:3])
  return(blacklist)  
}
```

**shade.blacklist.regions**

This function overlays semi-transparent bars at coordinates of blacklisted regions over a line graph of ChIP enrichment over a chromosome. The first parameter 'input' is arbitrary - it is just to collect coordinates so can be any file of the same structure (same windows as rows) as the graph to overlay. The 'blacklist' parameter can be a 'low mappability' blacklsit (for non-allelic graphs), a 'low allelic' blacklist (for allelic graphs) or a combined blacklsit by both criteria. ylim should be fit to the y axis of the graph the blacklist is to be overlain on top of. Unfortunately it was tricky to extend shahed blacklists to negative 'y' values. If necessary, these can be extended manually in a figure editing software (eg. Abode Illustrator, Affinity Designer) after generating pdfs of the graphs. 


```{shade.blacklist.regions}
shade.blacklist.regions <- function(input,blacklist,ylim){
  blacklist_binary <- input
  blacklist_binary[which(
    blacklist_binary$End %in% blacklist$End),]$Value <- 1
  blacklist_binary[which(
    blacklist_binary$End %notin% blacklist$End),]$Value <- 0
  y <- (blacklist_binary$Value)*ylim
  x <- (blacklist_binary$Start)/1000000
  y2 <- rep(y, each=2)
  y2 <- y2[-length(y2)]
  x2 <- rep(x, each=2)[-1]
  x3 <- c(min(x2), x2, max(x2))
  y3 <- c(0, y2, 0)
  polygon(x3, y3, border=NA, col=rgb(240, 240, 240, max = 255, alpha = 180)) #to add shaded regions
```

**calculateCorr**

This function calculates the correlation coefficent (R) between the distribution patterns of two samples. It also takes a 'blacklist' parameter to discount blacklisted regions from the correlation calculation.

```{calculateCorr}
calculateCorr <- function(table1,table2,blacklist){
  table1_scale <- table1
  table1_scale[which(table1_scale$End%in% blacklist$End),][,5] <- NA
  table1_scale <- table1_scale[,5]
  table2_scale <- table2
  table2_scale[which(table2_scale$End%in% blacklist$End),][,5] <- NA
  table2_scale <- table2_scale[,5]
  mat <- as.matrix(cbind(table1_scale,
                         table2_scale))
  correlation <- cor(mat, use = "complete.obs", method="pearson")
  return(correlation)  
}
```
