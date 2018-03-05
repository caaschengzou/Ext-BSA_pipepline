#!/usr/bin/env bash
# The pipeline script for Ext-BSA

# Predefined parameters
FREBAYES="/programs/freebayes/bin"
REF="PacBio_cinerea_canu1.contigs.fasta.upper.fa"  # genome reference for SNP calling
DPlow=20    # the lowest sequencing depth threshold
DPhigh=60   # the highest sequencing depth threhold

# Step1: Variant calling
#-----------------------------------------------------------------------------------
echo -e "-----------------------STEP 1------------------------------------\n"
echo -e "Variant calling for BAM files prepared...\n"
echo -e "-----------------------------------------------------------------\n\n"

# Index reference fasta file
samtools faidx ${REF}

# Call variants using bam
~/bin/bamaddrg/bamaddrg -b High.bulk.sorted.bam -s HighBulk -r HighBulk.s1.1 \
     -b Low.bulk.sorted.bam -s LowBulk -r LowBulk.s1.1 |\
     ${FREBAYES}/freebayes -f ${REF} -K -m 30 -q 20 -C 2 -3 10 -G 10 --stdin > F00.freebayes.BSA.vcf

egrep '#|TYPE=snp' F00.freebayes.BSA.vcf | egrep -v "scaffold|Mt|Pt" | awk '$10 != "." && $11 != "."'> F01.ensure.gt.vcf

# Generate text format for next step
bcftools  query\
     -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%DP[\t%SAMPLE;%GT;%DP;%AD]\n' F01.ensure.gt.vcf\
     > F02.bcftools.format.vcf.txt

# Define the sequencing depth thresholds
awk -v dplow=${DPlow} -v dphigh=${DPhigh} '$5 >= dplow && $5 <= dphigh ' F02.bcftools.format.vcf.txt |\
     sed 's/[,;]/\t/g' |\
     awk 'length($3) == 1 && length($4) ==1 && NF==15'|\
     awk '{printf"%s\t%s\t%s\t%s\t%s:%s\t%.4f\t%s:%s\t%.4f\n",$1,$2,$3,$4,$6,$7,$10/$8,$11,$12,$15/$13}'\
     > F03.allele.frequency.txt

# Format allele frequency text for Rscript testing
cut -f1,2,6,8 F03.allele.frequency.txt | sort -n -k 1 -k 2 > F04.ratio.for.test.sort.txt

echo -e "Step 1: Variant calling finished!\n\n"

# Step2: Ext-BSA testing
#-----------------------------------------------------------------------------------
echo -e "-----------------------STEP 2------------------------------------\n"
echo -e "Starting Ext-BSA...\n"
echo -e "-----------------------------------------------------------------\n\n"

cat << EOF > Ext-BSA.script.R
## Ext-BSA -- step 1
## Import allele frequency file from stdin

Args <- commandArgs()
file <- Args[6]
t <- read.table(file, header = F)
ouf_data1 <- "Ext_BSA_extension_result.txt"
ouf_data2 <- "Ext_BSA_T_test_results.txt"

## Ext-BSA -- step 2
## Define data format and significant level

diffr <- abs(t\$V3 - t\$V4)
s.diffr <- sort(diffr)
r1 <- rank(diffr)
y <- -log10(1 - (r1-0.01)/length(r1))
index.y <- which(r1 > length(r1)*0.99)
threshold0 <-  -log10(1 -  length(r1)*0.99/length(r1))
print(threshold0)

## Ext-BSA -- step 3
## Construct test function and run it

walking <- function(start, step = 5, threshold = threshold0) {
  score = y[start]
  left = list()
  left[[1]] = c(start,score)
  i = 1
  while(score > threshold) {
    tinterval = y[(start - step*i):(start - step*(i - 1) - 1)]
    score = score + sum(tail(sort(tinterval), 2)) - step
    i = i+1
    inp = c(start - step*(i - 1), score)
    left[[i]] = inp
  }
  
  start_score = y[start]
  right = list()
  right[[1]] = c(start,score)
  i = 1
  while(score > threshold) {
    tinterval = y[(start + step*(i - 1) + 1):(start + step*i)]
    score = score + sum(tail(sort(tinterval), 2)) - step
    i = i + 1
    inp = c(start-step*(i - 1), score)
    right[[i]] = inp
    
  }
  return(c(left, right))
}

index.y <- subset(index.y, index.y > 5)
lout <- lapply(index.y, walking)
outm <- matrix(unlist(lout), ncol = 2, byrow = T)
colnames(outm) <- c("ind", "value")
df <- data.frame(outm)
df2 <- data.frame( t[df\$ind, 2], t[df\$ind, 1], t[df\$ind, 2], df\$value )
colnames(df2) <- c("SNP", "CHR", "BP", "P")


## Ext-BSA -- step 4 
## Remove replicate and save results
df2[] <- lapply(df2, function(x) type.convert(as.character(x)))
final_result <- aggregate(P~CHR+BP, df2, max)
write.table(final_result, file = ouf_data1, col.names = T, row.names = F, quote = F, sep = "\t")


## T-test -- step5
## Bulked sample number
n1 <- 30
n2 <- 30

raw.data <- t
t.stat  <- matrix(NA, ncol=1, nrow=nrow(raw.data))
p.value <- matrix(NA, ncol=1, nrow=nrow(raw.data))

## Calculate the t.stat
raw.data[,5] <- (raw.data[,3] * n1 + raw.data[,4] * n2)/(n1 + n2)
t.stat[, 1] <- (raw.data[,3] - raw.data[,4])/sqrt((1.0/(2.0 * n1) + 1.0/(2.0 * n2)) * raw.data[,5] * (1.0 - raw.data[,5]))

## Calculate the pvalues and define the DF of Test

p.value <- 2 * pt(abs(t.stat),  n1 - 1, lower.tail = F)
t_result <- data.frame(raw.data[,c(2,1:2)], p.value)
colnames(t_result) <- c("SNP", "CHR", "BP", "P")

## Save T-test results 
write.table(t_result, file = ouf_data2, quote=F, row.names=F, col.names = T, sep = "\t")

EOF

# Run R script

Rscript  Ext-BSA.script.R  F04.ratio.for.test.sort.txt  

echo -e "Step 2: Ext-BSA finished!\n\n"
#-----------------------------------------------------------------------------------
