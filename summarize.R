#! /usr/bin/Rscript-4.1

# making associated loci given significant variants

args = commandArgs(TRUE)
argsLen <- length(args)

if(argsLen != 6) stop('error: wrong number of arguments (5e-08_list.txt distance(kb) flanking(kb) outdir)')

library(dplyr)

# arguments #####
var_list = args[1] #All.results.5e-8.txt
dist = as.numeric(args[2])*1000
flank = as.numeric(args[3])*1000
outdir=args[4]
file_n=args[5]
thresh=args[6]


# reading significant variants ####
all_sig_variants = read.table(var_list, header=F, sep="\t")
colnames(all_sig_variants)= c("SNP","CHR","POS","P")

# making associated loci ####
loci_results = matrix(NA,nrow=0, ncol=8)
colnames(loci_results)= c("LocusNo","START","END","START_pos","END_pos","START_pos2","END_pos2","CHR")

k=1
for (chr in c(1:22,"X")) {
#for (chr in 1:22) {
  temp = all_sig_variants[which(all_sig_variants[,2]==chr),] # extracting i-th chromosome
  temp = temp[order(temp[,3],decreasing = F),] # ordering by position
  Norow= nrow(temp) # the number of variants in i-th chromosome
  t=1
  while(t <= Norow ) { # scanning variants
    hold = c(paste("Chr",chr,"_Locus",k,sep=""), as.character(temp[1,1]), as.character(temp[1,1]), #LocusNo, start, end
             temp[1,3], temp[1,3], # start, end
             temp[1,3]-flank, temp[1,3]+flank, chr) # start_with_flaking 2, end_with_flaking 2, chromosome
    m=1
    i=2 # the next variant
    t=t+1 #
    
    if(t > Norow) {
      m=2
    }
    
    while (m == 1) {
      if(temp[(i-1),3] + dist < temp[i,3]) { # check the distance between this and the next variant
        m=2
      } else {
        hold[3] = as.character(temp[i,1]); hold[5] = temp[i,3]; hold[7] = temp[i,3] + flank # update end information with the next variant
        i=i+1
        if(t == Norow) {
          m=2
        }
        t=t+1
      }
    }
    
    # print(hold) # checking output
    loci_results = rbind(loci_results, hold) # merging with header information and associated loci
    colnames(loci_results)= c("LocusNo","START","END","START_pos","END_pos","START_pos2","END_pos2","CHR") # attaching header
  # RA_results
    temp = temp[-c(1:(i-1)),]
    k=k+1
  }
}

leads = c()
for(i in 1:nrow(loci_results)){
  t_chr = loci_results[i,8]; t_start = loci_results[i,4]; t_end = loci_results[i,5]
  temp_1 = all_sig_variants %>% filter(CHR == t_chr, POS >= t_start, POS <= t_end) %>% arrange(P) %>% .[1,]
  leads <- rbind(leads, temp_1)
}

colnames(leads)= c("Lead_SNP","Lead_CHR","Lead_POS","P")
loci_results =  cbind(loci_results,leads)


# saving the result
write.table(loci_results, file=paste(outdir,"/associated_region.summary_",file_n,"_",as.numeric(args[2]),"-",thresh,sep=""), 
            quote=F, col.names = T, row.names = F, sep=" ")
