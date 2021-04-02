#!/usr/bin/env Rscript

library("optparse")
library("stringr")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name of format .fa/.fasta/.csv with just id, size and sequence. If .csv, only ',' accepted", metavar="character"),
  make_option(c("-x", "--first_nt_position"), type="numeric", default=3, 
              help="first nucleotide position, 3 by default", metavar="numeric"),
  make_option(c("-o", "--output_name"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



if (is.null(opt$input)){
  print_help(opt_parser)
  stop("input_file needed", call.=FALSE)
}
if (is.null(opt$output_name)){
  print_help(opt_parser)
  stop("output_file needed", call.=FALSE)
}

if (grepl('.fasta$',opt$input) | grepl('.fa$',opt$input)) {
  
  inicial <- read.csv(opt$input, sep=';', header = F, stringsAsFactors = F)
  # inicial <- read.csv('~/Nextcloud/1_tesi_Adrià/test_DnoisE/PHY1subset_final.fa',sep=';', header = F,stringsAsFactors = F)
  info <- inicial[seq(1,dim(inicial)[1],2),]
  seqs <- inicial[seq(2,dim(inicial)[1],2),]
  inicial <- cbind(info[,1:2],seqs[,1])
  inicial[,2] <- as.numeric(str_extract(inicial[,2],regex('[0-9].*')))
  rm('info','seqs')
  names(inicial)<-c("id","nreads","seq")
  inicial$seq <- as.character(inicial$seq)
  inicial <- dplyr::arrange(inicial,inicial$id)
  
} else if (grepl('.csv$',opt$input)) {
  inicial <- read.csv(opt$input,stringsAsFactors = F)
  names(inicial)<-c("id","nreads","seq")
  inicial <- dplyr::arrange(inicial,inicial$id)
} else {
  stop("unknown input file format!" )
}





library(entropy)

entrpy<-function(seqsbones)
{
  # message("computing entropy at alpha= ",alpha[i])
  
  seqsbones<-seqsbones[nchar(seqsbones$seq)==seq_length,]
  
  #el primer nucleotid es posicio 2 del codo
  pos1<-seq(2,seq_length,by=3)
  pos2<-seq(3,seq_length,by=3)
  pos3<-seq(1,seq_length,by=3)
  
  
  long<-seq_length
  n_seqs<-vector("integer",dim(seqsbones)[1])
  total<-matrix(as.numeric("0"),4,313)
  rownames(total)<-c("A","C","G","T")
  
  
  seqs<-seqsbones$seq
  seqqs<-matrix("character",length(seqs),long)
  for (j in 1:length(seqs)) seqqs[j,]<-substring(seqs[j], seq(1, nchar(seqs[j]), 1), seq(1, nchar(seqs[j]), 1))
  
  rebost<-matrix(as.numeric("0"),4,long)
  rownames(rebost)<-c("A","C","G","T")
  for (j in 1:long)
    for (h in 1:length(tapply(seqsbones$nreads,seqqs[,j],sum)))
      rebost[rownames(rebost)==names(tapply(seqsbones$nreads,seqqs[,j],sum)[h]),j]<-rebost[rownames(rebost)==names(tapply(seqsbones$nreads,seqqs[,j],sum)[h]),j]+tapply(seqsbones$nreads,seqqs[,j],sum)[h]
  total<-total+rebost
  
  
  
  #generem fitxers amb resultats
  
  
  entropy1<-vector("numeric",length(pos1))
  entropy2<-vector("numeric",length(pos2))
  entropy3<-vector("numeric",length(pos3))
  
  
  for (q in 1:length(pos1)) entropy1[q]<-entropy(total[,pos1[q]])
  for (q in 1:length(pos2)) entropy2[q]<-entropy(total[,pos2[q]])
  for (q in 1:length(pos3)) entropy3[q]<-entropy(total[,pos3[q]])
  
  
  meanent1<-mean(entropy1)
  meanent2<-mean(entropy2)
  meanent3<-mean(entropy3)
  out<-data.frame('E1'=meanent1,'E2'=meanent2,'E3'=meanent3)
  return(out)
  
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#main program

seq_length=getmode(str_count(inicial$seq))

inicial <- inicial[str_count(inicial$seq)==seq_length,]

entropia_inicial <- entrpy(inicial)
write.csv(entropia_inicial,opt$output_name, row.names = F)

