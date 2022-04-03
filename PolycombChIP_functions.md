## Useful R functions for allelic Polycomb modification ChIP-seq data analysis
**CentreColandNames**

This function is designed to run directly on a nested read_delim() command. 
The read_delim() loads a .txt file produced from runnning the custom Python script ExtractInfoFrombedGraph_AtBed.py (https://github.com/guifengwei) on a bedgraph of the ChIP library with a .bed file of windows of a given length (eg. 250kb) of Chromosome X. 
The function names columns and creates an additional column of the centre of the window (for plotting).

```{CentreColandNames}
CentreColandNames <- function(table){
  table <- cbind(table[1:3],0.5*(table[2] + table[3]), table[4])
  colnames(table)<- c("Chr","Start","End","Centre","Value")
  return(table)
}

```
