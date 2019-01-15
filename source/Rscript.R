########################################################
# RNA-seq and Methyl-seq crosstalk to find driver gene #
# Author: Ran Yin                                      #
# Date: 12-10-2018                                     #
########################################################

require(data.table)
library(readr)
# Load data and find gene significantly changes in promoter and intergenic region#
dtm2 <- read.csv(file = "./data/wk2_pro_inter.csv", 
            header = TRUE, sep ="\t", col.names = )
colnames(dtm2)

#Make hyper and hypo methylated list for week 2
lsm2_hyper <- dtm2$geneId[dtm2$X2w_UVB..2w_UA.diff >= 0.3]
lsm2_hypo <- dtm2$geneId[dtm2$X2w_UVB..2w_UA.diff <= -0.3]

summary(lsm2_hyper)
summary(lsm2_hypo)
# find unique genes from both lists
lsm2_hyper <- unique(lsm2_hyper)
lsm2_hypo <- unique(lsm2_hypo)

summary(lsm2_hyper)
summary(lsm2_hypo)

View(lsm2_hyper)

# load RNA-seq data
dtr2 <- read_csv("data/RNA/Wk2 UA Ran version.csv", 
                               col_types = cols(log2FoldChange = col_number(), 
                                                padj = col_number(), pvalue = col_number()))

# find unique genes from up- and down-regulated genes in RNA
lsr2_up <- dtr2$gene[dtr2$log2FoldChange >= 2.0]
lsr2_down <- dtr2$gene[dtr2$log2FoldChange <= -2.0]


summary(lsr2_up)
View(lsr2_up)
summary(lsr2_down)
#Paired hyper-methylation with down-regulation gene
hyper_down <- unique(lsr2_down[lsr2_down %in% lsm2_hyper])
hypo_up <- unique(lsr2_up[lsr2_up %in% lsm2_hypo])

write.csv(hyper_down, file = "./tmp/hyper_down_wk2_UA.csv")
write.csv(hypo_up, file = "./tmp/hypo_up_wk2_UA.csv")
#From the summary, 305 genes find in hyper_down, and 53 genes find in hypo_up

#------------------------------------------------------------------------------#
# Load week 25 tumor data
dtm25t <- read.csv(file = "./data/wk25t_pro_inter.csv", 
                 header = TRUE, sep ="\t", col.names = )

#Make hyper and hypo methylated list for week 25t
lsm25t_hyper <- dtm25t$geneId[dtm25t$X25t_UVB..25t_UA.diff >= 0.3]
lsm25t_hypo <- dtm25t$geneId[dtm25t$X25t_UVB..25t_UA.diff <= -0.3]

summary(lsm25t_hyper)
summary(lsm25t_hypo)
# find unique genes from both lists
lsm25t_hyper<- unique(lsm25t_hyper)
lsm25t_hypo<- unique(lsm25t_hypo)

# load RNA-seq data
dtr25t <- read_csv("data/RNA/Tumor UA Ran version.csv", 
                 col_types = cols(log2FoldChange = col_number(), 
                                  padj = col_number(), pvalue = col_number()))

# find unique genes from up- and down-regulated genes in RNA
lsr25t_up <- dtr25t$gene[dtr25t$log2FoldChange >= 2.0]
lsr25t_down <- dtr25t$gene[dtr25t$log2FoldChange <= -2.0]

#Paired hyper-methylation with down-regulation gene
hyper_down_25t <- unique(lsr25t_down[lsr25t_down %in% lsm25t_hyper])
hypo_up_25t <- unique(lsr25t_up[lsr25t_up%in% lsm25t_hypo])

summary(hyper_down_25t)
summary(hypo_up_25t)

write.csv(hyper_down_25t, file = "./tmp/hyper_down_wk25t_UA.csv")
write.csv(hypo_up_25t, file = "./tmp/hypo_up_wk25t_UA.csv")

#section off




