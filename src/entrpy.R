#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-x", "--first_nt_position"), type="numeric", default=3, 
              help="first nucleotide position, 3 by default", metavar="numeric"),
  make_option(c("-f", "--fasta_input"), type="logical", default=FALSE, 
              help="input file is a .fasta file (TRUE) or .csv (FALSE, default)", metavar="logical"),
  make_option(c("-o", "--output_name"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("input_file needed", call.=FALSE)
}

library(entropy)

entrpy<-function(seqsbones)
{
  # message("computing entropy at alpha= ",alpha[i])
  #per si de cas comprovo que nomes em quedo  amb  seqs de 313
  message(i)
  
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
  out<-c(meanent1,meanent2,meanent3)
  return(out)
  
}

#main program

seq_length=313

inicial <- read.csv('PHY1bis_final.csv',stringsAsFactors = F)
names(inicial)<-c("id","nreads","seq")
inicial <- dplyr::arrange(inicial,inicial$id)
inicial2 <- read.csv('PHY1bis_final_subset.csv',stringsAsFactors = F)
names(inicial2)<-c("id","nreads","seq")
inicial2 <- dplyr::arrange(inicial2,inicial2$id)

# input_files <- c('PHY1bis_final.csv_no_Adcorr_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_no_Adcorr_Adcorr_denoised_ratio_d.csv',
#                  'PHY1bis_final.csv_Adcorr_nou_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adcorr_nou_Adcorr_denoised_ratio_d.csv',
#                  'PHY1bis_final.csv_Adcorr_antic_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adcorr_antic_Adcorr_denoised_ratio_d.csv',
#                  'PHY1bis_final.csv_seqs_juntes_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_seqs_juntes_Adcorr_denoised_ratio_d.csv',
#                  'PHY1bis_final.csv_inexplicades_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_inexplicades_Adcorr_denoised_ratio_d.csv',
#                  'PHY1bis_final.csv_quarta_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_quarta_Adcorr_denoised_ratio_d.csv',
#                  'PHY1bis_final.csv_Adsum_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adsum_Adcorr_denoised_ratio_d.csv')
# 
# pos_daugthers <-c()
# for (i in c(1:(length(input_files)/2))) {
#   file_info <- read.csv(input_files[i*2-1])
#   file_info <-  dplyr::arrange(file_info,file_info$daughter)
#   file2 <- read.csv(input_files[i*2],stringsAsFactors = F)
#   names(file2)<-c("id","nreads","seq")
#   if (i==1) {
#     entropia_ini <- c()
#   }
#   # if (i==7) {
#   #   entropia_ini <- c()
#   # }
#   pos_daugthers_1 <- c()
#   for (a in seq(0,100,20)) {
#     if (i==1) {
#       entropia_inicial <- entrpy(inicial[inicial$nreads>=a,])
#       names(entropia_inicial) <- c("meanent1_ini","meanent2_ini","meanent3_inic",
#                                "entropy_ratio_ini",
#                                "seent1_ini","seent2_ini","seent3_ini")
#       entropia_ini <- rbind(entropia_ini, entropia_inicial)
#     }
#     # if (i==7) {
#     #   entropia_inicial <- entrpy(inicial2[inicial2$nreads>=a,])
#     #   names(entropia_inicial) <- c("meanent1_ini","meanent2_ini","meanent3_inic",
#     #                                "entropy_ratio_ini",
#     #                                "seent1_ini","seent2_ini","seent3_ini")
#     #   entropia_ini <- rbind(entropia_ini, entropia_inicial)
#     # }
#     entropia_final <- entrpy(file2[file2$nreads>=a,])
#     
#     names(entropia_final) <- c("meanent1","meanent2","meanent3",
#                                "entropy_ratio",
#                                "seent1","seent2","seent3")
#     
#     pos <- c('file'=input_files[i*2-1], 'num' = a, 
#              'pos1' = sum(file_info$difpos1[inicial$nreads>=a], na.rm = T), 
#              'pos2' = sum(file_info$difpos2[inicial$nreads>=a], na.rm = T), 
#              'pos3' = sum(file_info$difpos3[inicial$nreads>=a], na.rm = T), 
#              'ASV' = dim(file2[file2$nreads>=a,])[1], 
#              'daughters' = sum(file_info$d[inicial$nreads>=a]>0, na.rm = T), 
#              entropia_final)
#     pos_daugthers_1 <- rbind(pos_daugthers_1,pos)
#   }
#   pos_daugthers_1 <- cbind(pos_daugthers_1, entropia_ini)
#   pos_daugthers <- rbind(pos_daugthers,pos_daugthers_1)
#   
# }


input_files <- c('PHY1bis_final.csv_no_Adcorr_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_no_Adcorr_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_Adcorr_nou_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_Adcorr_nou_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_Adcorr_antic_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_Adcorr_antic_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_seqs_juntes_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_seqs_juntes_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_inexplicades_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_inexplicades_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_quarta_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_quarta_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_Adsum_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_Adsum_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_Ad334_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_Ad334_Adcorr_denoised_ratio_d.csv',
                 'PHY1bis_final.csv_Adcorr1001_Adcorr_denoising_info.csv',
                 'PHY1bis_final.csv_Adcorr1001_Adcorr_denoised_ratio_d.csv')

pos_daugthers_3 <-c()
for (i in c(1:(length(input_files)/2))) {
  file_info <- read.csv(input_files[i*2-1])
  file_info <-  dplyr::arrange(file_info,file_info$daughter)
  file2 <- read.csv(input_files[i*2],stringsAsFactors = F)
  names(file2)<-c("id","nreads","seq")
  pos_daugthers_1 <- c()
  for (a in c(seq(2,5,1),10,15,20)) {
    
    pos <- data.frame('file'=stringr::str_remove(stringr::str_remove(input_files[i*2-1],'PHY1bis_final.csv_'),'_Adcorr_denoising_info.csv'), 
             'reads' = a, 
             'pos1' = sum(file_info$difpos1[inicial$nreads>=a], na.rm = T), 
             'pos2' = sum(file_info$difpos2[inicial$nreads>=a], na.rm = T), 
             'pos3' = sum(file_info$difpos3[inicial$nreads>=a], na.rm = T),
             'pos1 d=1' = sum(file_info$difpos1[inicial$nreads>=a & file_info$difpos1 == 1 & file_info$difpos2 == 0 & file_info$difpos3 == 0], na.rm = T), 
             'pos2 d=1' = sum(file_info$difpos2[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 1 & file_info$difpos3 == 0], na.rm = T), 
             'pos3 d=1' = sum(file_info$difpos3[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 0 & file_info$difpos3 == 1], na.rm = T), 
             'pos1 per reads' = sum(file_info$difpos1[inicial$nreads>=a]*inicial$nreads[inicial$nreads>=a], na.rm = T), 
             'pos2 per reads' = sum(file_info$difpos2[inicial$nreads>=a]*inicial$nreads[inicial$nreads>=a], na.rm = T), 
             'pos3 per reads' = sum(file_info$difpos3[inicial$nreads>=a]*inicial$nreads[inicial$nreads>=a], na.rm = T),
             'pos1 d=1 per reads' = sum(file_info$difpos1[inicial$nreads>=a & file_info$difpos1 == 1 & file_info$difpos2 == 0 & file_info$difpos3 == 0] * inicial$nreads[inicial$nreads>=a & file_info$difpos1 == 1 & file_info$difpos2 == 0 & file_info$difpos3 == 0], na.rm = T), 
             'pos2 d=1 per reads' = sum(file_info$difpos2[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 1 & file_info$difpos3 == 0] * inicial$nreads[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 1 & file_info$difpos3 == 0], na.rm = T), 
             'pos3 d=1 per reads' = sum(file_info$difpos3[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 0 & file_info$difpos3 == 1] * inicial$nreads[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 0 & file_info$difpos3 == 1], na.rm = T), 
             'ASV' = dim(file2[file2$nreads>=a,])[1], 
             'daughters' = sum(file_info$d[inicial$nreads>=a]>0, na.rm = T))
    pos_daugthers_1 <- rbind(pos_daugthers_1,pos)
  }
  pos_daugthers_3 <- rbind(pos_daugthers_3,pos_daugthers_1)
  
}

# input_files <- c('PHY1bis_final.csv_no_Adcorr_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adcorr_nou_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adcorr_antic_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_seqs_juntes_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_inexplicades_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_quarta_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adsum_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Ad334_Adcorr_denoising_info.csv',
#                  'PHY1bis_final.csv_Adcorr1001_Adcorr_denoising_info.csv')
# 
# pos_daugthers2 <-c()
# for (i in c(1:(length(input_files)))) {
#   file_info <- read.csv(input_files[i])
#   file_info <-  dplyr::arrange(file_info,file_info$daughter)
#   # file2 <- read.csv(input_files[i*2],stringsAsFactors = F)
#   # names(file2)<-c("id","nreads","seq")
#   for (a in seq(0,100,20)) {
#     # if (i == 1 |i == 2 |i == 3 |i == 4 |i == 5 |i == 6 ) {
#       pos <- c('file'=input_files[i], 'num' = a, 
#                'pos1' = sum(file_info$difpos1[inicial$nreads>=a & file_info$difpos1 == 1 & file_info$difpos2 == 0 & file_info$difpos3 == 0], na.rm = T), 
#                'pos2' = sum(file_info$difpos2[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 1 & file_info$difpos3 == 0], na.rm = T), 
#                'pos3' = sum(file_info$difpos3[inicial$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 0 & file_info$difpos3 == 1], na.rm = T))
#       
#     # } else {
#     #   pos <- c('file'=input_files[i], 'num' = a, 
#     #          'pos1' = sum(file_info$difpos1[inicial2$nreads>=a & file_info$difpos1 == 1 & file_info$difpos2 == 0 & file_info$difpos3 == 0], na.rm = T), 
#     #          'pos2' = sum(file_info$difpos2[inicial2$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 1 & file_info$difpos3 == 0], na.rm = T), 
#     #          'pos3' = sum(file_info$difpos3[inicial2$nreads>=a & file_info$difpos1 == 0 & file_info$difpos2 == 0 & file_info$difpos3 == 1], na.rm = T))
#     # 
#     # }
#     
#     
#     pos_daugthers2 <- rbind(pos_daugthers2,pos)
#   }
# }
# write.csv(pos_daugthers,'pos_daugthers.csv')
write.csv(pos_daugthers_3,'pos_daugthers3.csv')
# write.csv(pos_daugthers2,'pos_daugthers2.csv')
