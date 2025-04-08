
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')


args<- commandArgs(TRUE)
args_v<-NA
if(length(args)>1){
   args_v <- args[c(FALSE, TRUE)]
   names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
   print(args_v)
}
args<-NULL

csv_dir <- ifelse(is.na(args_v["csv"]), "_csv_", args_v["csv"])
out_dir <- ifelse(is.na(args_v["out"]), "_out_", args_v["out"])


filenames <- list.files(csv_dir, pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)

cluster.all.reads <- Reduce(function(x, y) merge(x, y, by="PasID", all=T), ldf)
cluster.all.reads[is.na(cluster.all.reads)] <- 0

cluster.all.reads <- cluster.all.reads %>% separate(PasID, sep = ":", into = c("chromosome","strand","position"), remove = T)

write.csv(cluster.all.reads, out_dir, col.names = T, row.names = F, quote = F)





