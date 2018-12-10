########################################################
# RNA-seq and Methyl-seq crosstalk to find driver gene #
# Author: Ran Yin                                      #
# Date: 12-10-2018                                     #
########################################################

require(data.table)



# find gene significantly changes in promoter region#
dtm2 <- read.csv(file = "./data/wk2_pro_inter.csv", 
            header = TRUE, sep ="\t", col.names = )

