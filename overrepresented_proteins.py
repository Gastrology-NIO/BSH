'''
# counts.csv
  protein,sample01,sample02,sample03,...,sample40
  P001,120,98,134,...,210
  P002,45,51,39,...,20
  P003,2300,1900,2500,...,3100

# metadata.csv
  sample,group
  sample01,healthy
  sample02,healthy
  ...
  sample21,diseased
  ...
  sample40,diseased
library(DESeq2)

counts <- read.csv("counts.csv", row.names = 1)
coldata <- read.csv("metadata.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treatment", "control"))

res <- res[order(res$padj), ]

write.csv(as.data.frame(res), "DESeq2_results.csv")



'''
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import gammaln
import pandas as pd

counts = pd.read_csv("counts.csv", index_col=0)
metadata = pd.read_csv("metadata.csv", index_col=0)

print(counts.shape)
print(metadata)

#count ~ disease_status
#zrobic tu na zewnątrz deseq2

results = pd.read_csv("deseq2_results.csv")

significant = results[
    (results["padj"] < 0.05) &
    (results["log2FoldChange"] > 1)
]

significant = significant.sort_values("padj")
