arg <- commandArgs(T)

if (length(arg)<3) {
  stop("WAAAAIIITTT!!!! something is missing")
} else{
  parent_motu <- arg[1]
  data_lulu_input <- arg[2]
  output_name <- arg[3]
}


data_lulu <- read.csv2(data_lulu_input)
#parent_motu <- "PHY1_000003276"
#output_name <- "PHY1_000003276_lulu"

daughters_motu <- data_lulu$id_removed[grep(parent_motu,data_lulu$parent_id)]

write.table(daughters_motu,file = output_name , quote = F,col.names = F,row.names = F)

  
  
  