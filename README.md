# allelic-polycombChIP
R functions and scripts for allelic analysis of ChIP-seq of Polycomb modifications, H3K27me3 and H2AK119ub1. Input files are of the form of tables of enrichment over 250KB windows across ChrX (as .txt files), produced by the custom python script *ExtractInfoFrombedGraph_AtBed.py* (https://github.com/guifengwei).



These analysis methods are broadly applicable to ChIP-seq data of other histone modifications or factors allelically-enriched over Xi which have a broad/blanket pattern of distribution, and less applicable for data types with defined enrichment peaks.
