## Useful R functions for allelic Polycomb modification ChIP-seq data analysis
**LoadingCountFile**

This function is designed to run directly on a nested read_delim() command. 
The read_delim() loads a .txt file produced from runnning the custom Python script ExtractInfoFrombedGraph_AtBed.py (https://github.com/guifengwei) ![image](https://user-images.githubusercontent.com/32739123/161447460-4976325a-5644-4180-ae67-d4f965b0fbf5.png) on a bedgraph of the ChIP library and a .bed files of windows of a given length (eg. 250kb) of Chromosome X. 

```
