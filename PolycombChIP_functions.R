#### Packages ####
library(readr)
library(dplyr)
library(mefa4)
library(Hmisc)
library(stringr)
#### Functions ####
CentreColandNames <- function(table){
  table <- cbind(table[1:3],0.5*(table[2] + table[3]), table[4])
  colnames(table)<- c("Chr","Start","End","Centre","Value")
  return(table)
}
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
  #to add lines to borders of blacklist regions
  
}
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
IPoverInput <- function(IP, input){
  Enrichment <- as.numeric(IP[,5] / input[,5])
  Enrichment[is.nan(Enrichment)] <- NA
  IPoverinput_table <- cbind(IP[1:4], Enrichment)
  return(IPoverinput_table)
}
DoxMinusNoDox <- function(Dox,NoDox){
  DoxMinusNoDox <- cbind(Dox[,1:4], Dox[,5] - NoDox[,5])
  colnames(DoxMinusNoDox)[5] <- "DoxMinusNoDox"
  return(DoxMinusNoDox)
}
XiMinusXa <- function(genome2,genome1){
  XiMinusXa <- cbind(genome2[,1:4], genome2[,5] - genome1[,5])
  colnames(XiMinusXa)[5] <- "XiMinusXa"
  return(XiMinusXa)
}
InputNorm <- function(DoxMinusNoDox,input){
  IP_vals <- unlist(DoxMinusNoDox[,5])
  IP_vals[IP_vals == 0] <- NA
  input_vals <- unlist(input[,5])
  input_vals[input_vals == 0] <- NA
  DoxMinusNoDox_Norm <- cbind(DoxMinusNoDox[1:4],
                                            (IP_vals)*(mean(input_vals,na.rm=TRUE)/input_vals))
  DoxMinusNoDox_Norm[,5] <- as.numeric(DoxMinusNoDox_Norm[,5])
  colnames(DoxMinusNoDox_Norm)[5] <- "DoxMinusNoDox_Norm"
  
  return(DoxMinusNoDox_Norm)
}
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
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
